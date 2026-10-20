@testset "single founder" for mode in (:tree, :cells)
    pop = initialize_population(1; mode)
    @test popsize(pop) == 1 && pop.t == 0.0 && pop.m == 0.0
    @test living_clones(pop.clones) == [1] && traits(pop.clones, 1) == Traits(0.0, 1.0)
    @test mutations_per_cell(pop) == [0] && divisions_per_cell(pop) == [0]
end

@testset "several founders and clonal mutations" begin
    pc = initialize_population(5; mode = :cells, s = 0.1, α = 0.9, clonal_mutations = 7)
    @test pc isa CellPopulation && mutations_per_cell(pc) == fill(7, 5) && cell_ids(pc) == 1:5
    pt = initialize_population(5; clonal_mutations = 7)
    @test pt isa TreePopulation && length(roots(pt)) == 5 && mutations_per_cell(pt) == fill(7, 5)
    @test pt.next_id == 6 && pc.next_id == 6 && cell_ids(pt) == 1:5
end

@testset "clone specification" begin
    pop = initialize_population([(n = 900, s = 0.0, α = 1.0), (n = 100, s = 0.0, α = 0.8)])
    @test pop.clones.size == [900, 100] && pop.clones.parent == [0, 0]
    @test pop.clones.kind == [:founder, :founder] && length(roots(pop)) == 1000
    @test count(==(2), cell_clones(pop)) == 100
end

@testset "validation" begin
    @test_throws ArgumentError initialize_population(0)
    @test_throws ArgumentError initialize_population(1; mode = :foo)
    @test_throws ArgumentError initialize_population(1; α = 1.5)
end

@testset "copy is independent" begin   # Review Focus 3
    pop = initialize_population(3)
    c = copy(pop)
    @test c.clones !== pop.clones && c.tree !== pop.tree && c.alive !== pop.alive && c.log !== pop.log
    push!(c.alive, c.alive[1]); c.clones.size[1] += 1
    @test length(pop.alive) == 3 && pop.clones.size[1] == 3
    pc = initialize_population(3; mode = :cells)
    cc = copy(pc)
    @test cc.cells !== pc.cells && cc.neutral_rng !== pc.neutral_rng
    @test rand(cc.neutral_rng) == rand(pc.neutral_rng)
end

@testset "simulating a copy leaves the original untouched" for mode in (:tree, :cells)   # Review Focus 3
    orig = grown(mode; N = 200)
    state(p) = (mode === :tree ? (ntuple(i -> copy(getfield(p.tree, i)), fieldcount(typeof(p.tree))), copy(p.alive)) : copy(p.cells),
                ntuple(i -> (x = getfield(p.clones, i); x isa AbstractVector ? copy(x) : x), fieldcount(CloneTable)),
                copy(p.log), p.t, p.m, p.next_id, p.clonal_mutations, rand(copy(p.neutral_rng), UInt64))
    before = state(orig)
    c = copy(orig)
    churn = BirthDeathBlock(birth = s -> 1.0, death = s -> 0.9, μ = 2.0, stop_above = 400, tmax = orig.t + 50.0,
                            mutations = (Resistance(ν = 0.05, α = 0.3), Drivers(ν = 0.05, Δs = 0.1)))
    simulate!(c, churn, StableRNG(5))
    @test c.log[end].accepted >= 300 && n_clones(c.clones) > n_clones(orig.clones)   # precondition: real work
    @test state(orig) == before
end
