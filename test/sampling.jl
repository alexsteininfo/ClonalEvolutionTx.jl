@isdefined(grown) || include("fixtures.jl")

# id => value maps over the living cells
byid(pop, v) = Dict(zip(cell_ids(pop), v))
treegrown(; N = 300, seed = 1) = grown(:tree; N, seed)
# The reference draw, re-derived independently of the package code.
function reference_draw(ids, n, seed)
    rng, perm = StableRNG(seed), collect(1:length(ids))
    for i in 1:n
        j = rand(rng, i:length(ids))
        perm[i], perm[j] = perm[j], perm[i]
    end
    return ids[perm[1:n]]
end
tree_snapshot(pop) = (copy(pop.tree.id), copy(pop.tree.parent), copy(pop.alive), copy(pop.clones.size),
                      mutations_per_cell(pop), site_frequency_spectrum(pop))

@testset "GOLDEN: the draw recipe (StableRNG + partial Fisher–Yates over the living order)" begin
    pop = treegrown()
    for (seed, n) in ((1, 1), (1, 5), (42, 30), (7, popsize(pop)))
        @test cell_ids(sample_population(pop, n; seed)) == reference_draw(cell_ids(pop), n, seed)
    end
    cp = grown(:cells; N = 300)
    @test cell_ids(sample_population(cp, 20; seed = 9)) == reference_draw(cell_ids(cp), 20, 9)
    @test ClonalEvolutionTx._draw_indices(StableRNG(1), 10, 4) == reference_draw(collect(1:10), 4, 1)
end

@testset "n = N_full reproduces the source tree" begin
    pop = treegrown()
    s = sample_population(pop, popsize(pop); seed = 1)
    @test sort(cell_ids(s)) == sort(cell_ids(pop))
    @test byid(s, mutations_per_cell(s)) == byid(pop, mutations_per_cell(pop))
    @test byid(s, divisions_per_cell(s)) == byid(pop, divisions_per_cell(pop))
    @test site_frequency_spectrum(s) == site_frequency_spectrum(pop)
    @test issubset(s.tree.id, pop.tree.id)   # the source also holds ancestors of dead cells
    @test s.clones.size == pop.clones.size
end

@testset "sampling is deterministic in the seed" begin
    pop = treegrown()
    a, b = sample_population(pop, 40; seed = 42), sample_population(pop, 40; seed = 42)
    @test cell_ids(a) == cell_ids(b) && a.tree.id == b.tree.id && a.tree.parent == b.tree.parent
    @test sample_population(pop, 40; seed = UInt64(42)).alive == a.alive
end

@testset "sample_population errors on out-of-range n" begin
    for pop in (treegrown(N = 10), grown(:cells; N = 10))
        @test_throws ArgumentError sample_population(pop, 0; seed = 1)
        @test_throws ArgumentError sample_population(pop, popsize(pop) + 1; seed = 1)
    end
    @test_throws ArgumentError sample_population(grown(:tree; N = 10), 11; seed = 1)
end

@testset "a single-cell population samples to itself" begin
    pop = initialize_population(1)
    s = sample_population(pop, 1; seed = 1)
    @test cell_ids(s) == [1] && length(s.tree) == 1 && s.alive == Int32[1]
    @test_throws ArgumentError sample_population(pop, 2; seed = 1)
end

@testset "INVARIANCE: burden and depth are the cell's full-tree values" begin
    pop = treegrown(N = 400)
    burden, depth = byid(pop, mutations_per_cell(pop)), byid(pop, divisions_per_cell(pop))
    for n in (1, 5, 40, popsize(pop))
        s = sample_population(pop, n; seed = n)
        ids = cell_ids(s)
        @test mutations_per_cell(s) == [burden[i] for i in ids]
        @test divisions_per_cell(s) == [depth[i] for i in ids]
        @test issubset(ids, cell_ids(pop)) && allunique(ids) && length(ids) == n
    end
end

@testset "every retained node has a sampled leaf below it, and vice versa" begin
    pop = treegrown()
    s = sample_population(pop, 25; seed = 5)
    t = s.tree
    drawn = Set(cell_ids(s))
    @test Set(t.id[i] for i in eachindex(t.id) if ClonalEvolutionTx._isleaf(t, i)) == drawn
    for r in roots(s), node in PreOrderDFS(r)
        @test !isempty(intersect(Set(l.id for l in Leaves(node)), drawn))
    end
    @test length(t) == length(Set(n.id for r in roots(s) for n in PreOrderDFS(r)))
end

@testset "n = 1 gives a path, and the founder stays the root" begin
    pop = treegrown(N = 100)
    seen = Set{Int64}()
    for seed in 1:30
        s = sample_population(pop, 1; seed)
        push!(seen, cell_ids(s)[1])
        @test length(roots(s)) == 1 && roots(s)[1].id == 1
        @test s.tree.parent[1] == ClonalEvolutionTx.NOPARENT
        @test count(i -> ClonalEvolutionTx._isleaf(s.tree, i), 1:length(s.tree)) == 1
        @test all(i -> (s.tree.left[i] != 0) + (s.tree.right[i] != 0) <= 1, 1:length(s.tree))
        @test length(s.tree) == divisions_per_cell(s)[1] + 1
    end
    @test length(seen) > 10
end

@testset "left/right slots are preserved and links are consistent" begin
    pop = treegrown(N = 200)
    s = sample_population(pop, 30; seed = 3)
    t = s.tree
    for i in 1:length(t)
        for c in (t.left[i], t.right[i])
            c == 0 || @test t.parent[c] == i
        end
        p = t.parent[i]
        p == 0 || @test t.left[p] == i || t.right[p] == i
        j = searchsortedfirst(pop.tree.id, t.id[i])
        p == 0 || @test (t.left[p] == i) == (pop.tree.left[pop.tree.parent[j]] == j)
    end
    @test issorted(t.id)
end

@testset "a forest keeps each sampled lineage's root" begin
    pop = simulate!(initialize_population(3), BirthDeathBlock(birth = s -> 1.0, death = s -> 0.2, μ = 1.0, stop_above = 150,
                    restart_on_extinction = true), StableRNG(5))
    @test length(roots(pop)) == 3
    s = sample_population(pop, 60; seed = 2)
    @test 1 <= length(roots(s)) <= 3 && issubset(Set(r.id for r in roots(s)), Set(r.id for r in roots(pop)))
    full = sample_population(pop, popsize(pop); seed = 2)
    @test length(roots(full)) == 3
    @test sort(mutations_per_cell(s)) == sort(Int[byid(pop, mutations_per_cell(pop))[i] for i in cell_ids(s)])
end

@testset "the source population is never mutated" begin
    pop = treegrown()
    before = tree_snapshot(pop)
    rngstate = copy(pop.neutral_rng)
    for n in (1, 10, popsize(pop))
        sample_population(pop, n; seed = n)
    end
    @test tree_snapshot(pop) == before && pop.neutral_rng == rngstate
end

@testset "node data equals the source's, in a separate tree" begin
    pop = treegrown()
    s = sample_population(pop, 30; seed = 1)
    @test s.tree !== pop.tree && s.clones !== pop.clones && s.log !== pop.log && s.neutral_rng !== pop.neutral_rng
    for i in 1:length(s.tree)
        j = searchsortedfirst(pop.tree.id, s.tree.id[i])
        @test pop.tree.id[j] == s.tree.id[i]
        @test (s.tree.birthtime[i], s.tree.clone[i], s.tree.neutral[i], s.tree.total[i]) ==
              (pop.tree.birthtime[j], pop.tree.clone[j], pop.tree.neutral[j], pop.tree.total[j])
    end
    @test (s.t, s.m, s.clonal_mutations, s.next_id, s.neutral_rng, s.log) ==
          (pop.t, pop.m, pop.clonal_mutations, pop.next_id, pop.neutral_rng, pop.log)
end

@testset "different seeds differ; one seed gives nested draws across n" begin
    pop = treegrown(N = 200)
    a = cell_ids(sample_population(pop, 10; seed = 42))
    @test a != cell_ids(sample_population(pop, 10; seed = 2))
    @test cell_ids(sample_population(pop, 4; seed = 42)) == a[1:4]
    cp = grown(:cells; N = 200)
    @test cell_ids(sample_population(cp, 4; seed = 42)) == cell_ids(sample_population(cp, 10; seed = 42))[1:4]
end

@testset "SFS of a sample sums to the sample's total burden" begin
    pop = treegrown(N = 150)
    for n in (1, 10, 150)
        s = sample_population(pop, n; seed = 100 + n)
        sfs = site_frequency_spectrum(s)
        @test length(sfs) == n
        @test sum(k * sfs[k] for k in 1:n) == sum(mutations_per_cell(s) .- s.clonal_mutations)
    end
end

@testset "clone sizes are recounted over the sample" begin
    pop = grown(:tree; N = 2000)
    s = sample_population(pop, 100; seed = 1)
    @test popsize(s) == 100 && sum(clone_sizes(s)) == 100 && n_clones(s.clones) == n_clones(pop.clones)
end

@testset "the sampled clone table is internally consistent" begin
    for mode in (:tree, :cells)
        pop = grown(mode; N = 2000, ν = 0.05)
        @test n_clones(pop.clones) > 3
        for n in (1, 7, 100, popsize(pop))
            ct = sample_population(pop, n; seed = n).clones
            counts = zeros(Int, n_clones(ct))
            for k in cell_clones(sample_population(pop, n; seed = n)); counts[k] += 1; end
            @test ct.size == counts
            @test sort(ct.alive) == findall(>(0), counts)
            @test all(ct.alive_pos[k] == findfirst(==(k), ct.alive) for k in ct.alive)
            @test all(ct.alive_pos[k] == 0 for k in 1:n_clones(ct) if counts[k] == 0)
            @test ct.sumNβ ≈ sum(counts .* ct.β) && ct.sumNα ≈ sum(counts .* ct.α)
            @test ct.maxβ == maximum(ct.β[ct.alive]) && ct.β[ct.argmaxβ] == ct.maxβ
            @test ct.maxα == maximum(ct.α[ct.alive]) && ct.α[ct.argmaxα] == ct.maxα
            @test ct.ever_alive == pop.clones.ever_alive
            @test ct.parent == pop.clones.parent && ct.origin_time == pop.clones.origin_time
        end
        # the sampled table keeps working: a clone alive before may go extinct and be resized
        ct = sample_population(pop, 10; seed = 1).clones
        k = ct.alive[1]
        ClonalEvolutionTx._change_size!(ct, k, -ct.size[k])
        @test ct.size[k] == 0 && !(k in ct.alive)
    end
end

@testset "cell mode samples carry the cells' burdens" begin
    pop = grown(:cells; N = 2000)
    s = sample_population(pop, 100; seed = 2)
    @test s isa CellPopulation
    @test issubset(cell_ids(s), cell_ids(pop))
    byid = Dict(zip(cell_ids(pop), mutations_per_cell(pop)))
    @test mutations_per_cell(s) == [byid[i] for i in cell_ids(s)]
    @test s.cells == [pop.cells[findfirst(c -> c.id == i, pop.cells)] for i in cell_ids(s)]
end

@testset "tree and cell samples of the same seed draw the same cells" begin
    t, c = grown(:tree; N = 300), grown(:cells; N = 300)
    @test cell_ids(t) == cell_ids(c)
    @test cell_ids(sample_population(t, 50; seed = 8)) == cell_ids(sample_population(c, 50; seed = 8))
end

@testset "snapshots evaluate statistics on samples" begin
    acc = MeasurementAccumulator(MeasurementSpec(triggers = [AtBlockEnd()], stats = [popsize, SFS()], samples = [50, 200], sample_seed = 3))
    simulate!(initialize_population(1), BirthDeathBlock(birth = s -> 1.0, death = s -> 0.0, μ = 1.0, stop_above = 1000), StableRNG(4); accumulator = acc)
    snap = finalize_measurements(acc).snapshots[end]
    @test [x.n for x in snap.samples] == [50, 200] && snap.samples[1].stats.popsize == 50
    @test length(snap.samples[2].stats.sfs) == 200
    @test_throws ArgumentError sample_population(grown(:tree; N = 10), 11; seed = 1)
end

# Check every snapshot's samples against samples drawn independently from the population
# the snapshot captured (`:pop => copy`), with i the snapshot's position in the final vector.
function check_snapshot_samples(snaps, ns, base)
    for (i, snap) in enumerate(snaps)
        @test [x.n for x in snap.samples] == ns
        for x in snap.samples
            expected = sample_population(snap.stats.pop, x.n; seed = derive_seed(base, i, x.n))
            @test cell_ids(x.stats.pop) == cell_ids(expected)
            @test x.stats.popsize == x.n
        end
    end
end

@testset "snapshot samples use derive_seed(sample_seed, snapshot_index, n)" begin
    for mode in (:tree, :cells)
        spec = MeasurementSpec(triggers = [AtSize(100), AtSize(200), AtBlockEnd()],
                               stats = [:pop => copy, popsize], samples = [10, 30], sample_seed = 11)
        acc = MeasurementAccumulator(spec)
        simulate!(initialize_population(1; mode), BirthDeathBlock(birth = s -> 1.0, death = s -> 0.0, stop_above = 200),
                  StableRNG(1); accumulator = acc)
        snaps = finalize_measurements(acc).snapshots
        @test length(snaps) == 3
        check_snapshot_samples(snaps, [10, 30], 11)
        # the index and the sizes both matter: other seeds draw other cells
        @test cell_ids(snaps[1].samples[2].stats.pop) !=
              cell_ids(sample_population(snaps[1].stats.pop, 30; seed = derive_seed(11, 2, 30)))
    end
    @test all(isempty(s.samples) for s in
              finalize_measurements(MeasurementAccumulator(MeasurementSpec())).snapshots)
end

@testset "snapshot sample indices survive extinction restarts" begin
    spec = MeasurementSpec(triggers = [AtSize(3), AtBlockEnd()], stats = [:pop => copy, popsize],
                           samples = [2, 3], sample_seed = 5)
    restarted = 0
    for seed in 1:40
        acc = MeasurementAccumulator(spec)
        pop = simulate!(initialize_population(1), BirthDeathBlock(birth = s -> 1.0, death = s -> 0.9, stop_above = 12,
                        restart_on_extinction = true), StableRNG(seed); accumulator = acc)
        r = sum(b.restarts for b in pop.log)
        r > 0 || continue
        restarted += 1
        snaps = finalize_measurements(acc).snapshots
        @test length(snaps) == 2
        check_snapshot_samples(snaps, [2, 3], 5)
    end
    @test restarted > 0
end

@testset "a snapshot sample larger than the population is nothing" begin   # final review F1
    acc = MeasurementAccumulator(MeasurementSpec(triggers = [AtBlockEnd()], stats = [popsize], samples = [500, 20]))
    pop = simulate!(initialize_population(1), BirthDeathBlock(birth = s -> 1.0, death = s -> 0.0, stop_above = 50), StableRNG(1); accumulator = acc)
    s = only(finalize_measurements(acc).snapshots).samples
    @test s[1] == (n = 500, stats = nothing) && s[2] == (n = 20, stats = (popsize = 20,))
    @test popsize(pop) == 50 && pop.log[end].end_reason == :above
end
