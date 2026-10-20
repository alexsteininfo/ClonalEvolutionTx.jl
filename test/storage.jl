@testset "round trip in both modes" for mode in (:tree, :cells)
    pop = grown(mode; N = 1000)
    path = tempname()
    save_population(path, pop)
    q = load_population(path)
    @test typeof(q) == typeof(pop) && q.t == pop.t && q.clonal_mutations == pop.clonal_mutations
    @test cell_ids(q) == cell_ids(pop) && mutations_per_cell(q) == mutations_per_cell(pop)
    @test q.clones.size == pop.clones.size && q.clones.α == pop.clones.α && q.clones.kind == pop.clones.kind
    simulate!(q, BirthDeathBlock(birth = s -> 1.0, death = s -> 0.3, tmax = q.t + 0.1), StableRNG(1))
    @test q.t ≈ pop.t + 0.1
end

@testset "bad files are rejected" begin
    path = tempname(); write(path, "not a population")
    @test_throws ArgumentError load_population(path)
end

@testset "newick of the hand-built tree" begin
    io = IOBuffer(); write_newick(io, small_tree())
    @test String(take!(io)) == "(2:1,3:1):0;\n"
end

@testset "save_tree writes a subtree" begin
    pop = grown(:tree; N = 200)
    node = first(children(roots(pop)[1]))
    path = tempname(); save_tree(path, pop, node)
    @test popsize(load_population(path)) == length(collect(Leaves(node)))
end

@testset "truncated, versioned and wrong-mode files" begin
    pop = grown(:tree; N = 100)
    path = tempname(); save_population(path, pop)
    bytes = read(path)
    for cut in (0, 5, 20, length(bytes) ÷ 2, length(bytes) - 1)
        p = tempname(); write(p, bytes[1:cut])
        @test_throws ArgumentError load_population(p)
    end
    bad = copy(bytes); bad[9] = 0x07                  # version
    p = tempname(); write(p, bad); @test_throws ArgumentError load_population(p)
    bad = copy(bytes); bad[13] = 0x09                 # mode
    p = tempname(); write(p, bad); @test_throws ArgumentError load_population(p)
end

@testset "saving never changes the population; removed rows are left out" begin
    pop = grown(:tree; N = 300)
    simulate!(pop, BirthDeathBlock(birth = s -> 1.0, death = s -> 1.2, tmax = pop.t + 1.0,
                                   stop_below = 50), StableRNG(2))
    before = (copy(pop.tree.parent), copy(pop.alive), pop.tree.ndead)
    path = tempname(); save_population(path, pop)
    @test (pop.tree.parent, pop.alive, pop.tree.ndead) == before
    q = load_population(path)
    @test q.tree.ndead == 0 && cell_ids(q) == cell_ids(pop)
    @test divisions_per_cell(q) == divisions_per_cell(pop)
    @test length(q.tree) <= length(pop.tree)
end

@testset "newick options and shapes" begin
    pop = small_tree()
    pop.t = 3.0                        # the fixture's cells are born at 1.0
    io = IOBuffer(); write_newick(io, pop; branch_lengths = :time)
    @test String(take!(io)) == "(2:2.0,3:2.0):1.0;\n"
    @test_throws ArgumentError write_newick(IOBuffer(), pop; branch_lengths = :nope)
    one = initialize_population(1); io = IOBuffer(); write_newick(io, one)
    @test String(take!(io)) == "(1:0):0;\n"
    many = grown(:tree; N = 50); many2 = initialize_population(3)
    io = IOBuffer(); write_newick(io, many2)
    @test String(take!(io)) == "(1:0):0;\n(2:0):0;\n(3:0):0;\n"
    io = IOBuffer(); write_newick(io, many)
    out = String(take!(io))
    @test count(==(';'), out) == length(roots(many)) && count(==('\n'), out) == length(roots(many))
    path = tempname(); write_newick(path, pop)
    @test read(path, String) == "(2:1,3:1):0;\n"
end

@testset "save_tree rejects foreign nodes" begin
    a, b = grown(:tree; N = 50), grown(:tree; N = 50, seed = 2)
    @test_throws ArgumentError save_tree(tempname(), a, roots(b)[1])
end

@testset "newick of a grown tree: every living cell once, deep chain safe" begin
    pop = grown(:tree; N = 200)
    io = IOBuffer(); write_newick(io, pop)
    out = String(take!(io))
    ids = [parse(Int, m.captures[1]) for m in eachmatch(r"[(,](\d+):", out)]
    @test sort(ids) == sort(cell_ids(pop))
    # a long chain of single-child nodes collapses to one stem
    deep = initialize_population(1)
    for k in 2:20000
        i = ClonalEvolutionTx._push_node!(deep.tree, deep.alive[1], k, Float64(k), 1, 0, 0)
        deep.tree.left[deep.alive[1]] = i
        deep.alive[1] = i
    end
    io = IOBuffer(); write_newick(io, deep)
    @test String(take!(io)) == "(20000:0):19999;\n"
end

# ── Final review fixes ────────────────────────────────────────────────────────

# A hand-built tree with dyadic times, at t = 4: founder 1 (born 0) has one surviving
# daughter 2 (born 0.5); 2 divides at 1 into 4 and 5; 4 has one surviving daughter 6
# (born 2); 5 divides at 3 into 8 and 9. Living: 6, 8, 9.
function timed_tree()
    pop = initialize_population(1)
    tree = ClonalEvolutionTx.LineageTree()
    push = ClonalEvolutionTx._push_node!
    n1 = push(tree, Int32(0), 1, 0.0, 1, 0, 0)
    n2 = push(tree, n1, 2, 0.5, 1, 0, 0); tree.left[n1] = n2
    n4 = push(tree, n2, 4, 1.0, 1, 0, 0)
    n5 = push(tree, n2, 5, 1.0, 1, 0, 0); tree.left[n2], tree.right[n2] = n4, n5
    n6 = push(tree, n4, 6, 2.0, 1, 0, 0); tree.left[n4] = n6
    n8 = push(tree, n5, 8, 3.0, 1, 0, 0)
    n9 = push(tree, n5, 9, 3.0, 1, 0, 0); tree.left[n5], tree.right[n5] = n8, n9
    pop.tree, pop.alive, pop.t, pop.next_id = tree, Int32[n6, n8, n9], 4.0, 10
    ClonalEvolutionTx._change_size!(pop.clones, 1, 2)
    return pop
end

# Parse one Newick tree as written by write_newick: each leaf's own branch length and the
# groups (internal nodes; the outermost is the root, whose length is the stem) above it.
function parse_newick(line::AbstractString)
    leaf_len, above, glen, open = Dict{Int,Float64}(), Dict{Int,Vector{Int}}(), Float64[], Int[]
    i = 1
    token() = (j = i; while i <= lastindex(line) && !(line[i] in ",();\n:"); i += 1; end; line[j:i-1])
    while i <= lastindex(line)
        c = line[i]
        if c == '('
            push!(glen, NaN); push!(open, length(glen)); i += 1
        elseif c == ')'
            g = pop!(open); i += 2                       # ")" and ":"
            glen[g] = parse(Float64, token())
        elseif c in ",;\n"
            i += 1
        else
            id = parse(Int, token()); i += 1             # ":"
            leaf_len[id] = parse(Float64, token()); above[id] = copy(open)
        end
    end
    return leaf_len, above, glen
end
root_to_tip(nw, id) = nw[1][id] + sum(nw[3][g] for g in nw[2][id])
to_ancestor(nw, a, b) = nw[1][a] + sum((nw[3][g] for g in setdiff(nw[2][a], nw[2][b])); init = 0.0)
newick_time(pop) = (io = IOBuffer(); write_newick(io, pop; branch_lengths = :time); String(take!(io)))

@testset "newick time lengths: splits at divisions, leaves at pop.t" begin   # final review F3
    pop = timed_tree()
    @test newick_time(pop) == "(6:3.0,(8:1.0,9:1.0):2.0):1.0;\n"
    io = IOBuffer(); write_newick(io, pop)
    @test String(take!(io)) == "(6:2,(8:1,9:1):1):1;\n"            # divisions unchanged
    one = initialize_population(1); one.t = 2.5
    @test newick_time(one) == "(1:0):2.5;\n"                        # stem = pop.t − founder birth
end

@testset "newick time trees are ultrametric and agree with coalescence_times" begin   # F3
    pop = grown(:tree; N = 120)
    simulate!(pop, BirthDeathBlock(birth = s -> 1.0, death = s -> 1.0, tmax = pop.t + 0.5), StableRNG(3))
    nw = parse_newick(newick_time(pop))
    ids = cell_ids(pop)
    @test sort(collect(keys(nw[1]))) == sort(ids)
    @test all(id -> isapprox(root_to_tip(nw, id), pop.t; atol = 1e-9), ids)   # founder born at 0
    ct = coalescence_times(pop)
    pairs = [(ids[x], ids[y]) for x in eachindex(ids) for y in (x + 1):lastindex(ids)]
    @test all(isapprox(to_ancestor(nw, a, b), ct[k]; atol = 1e-9) &&
              isapprox(to_ancestor(nw, b, a), ct[k]; atol = 1e-9) for (k, (a, b)) in enumerate(pairs))
    forest = simulate!(initialize_population(3), BirthDeathBlock(birth = s -> 1.0, death = s -> 0.3, tmax = 3.0), StableRNG(4))
    lines = split(newick_time(forest), '\n'; keepempty = false)
    @test length(lines) == length(roots(forest))
    @test all(all(id -> isapprox(root_to_tip(nw, id), forest.t; atol = 1e-9), keys(nw[1]))
              for nw in parse_newick.(lines))
    # a saved subtree: its root was born after 0, so the stem starts at that birth time
    path = tempname(); save_tree(path, pop, first(children(first(children(roots(pop)[1])))))
    q = load_population(path)
    nwq = parse_newick(newick_time(q))
    @test q.tree.birthtime[1] > 0 && all(id -> isapprox(root_to_tip(nwq, id), q.t - q.tree.birthtime[1]; atol = 1e-9), cell_ids(q))
end

@testset "save_tree: the subtree root carries the inherited mutations" begin   # final review F4
    pop = grown(:tree; N = 400)
    pop.clonal_mutations = 2
    node = first(n for n in PreOrderDFS(roots(pop)[1]) if n.total > n.neutral && length(collect(Leaves(n))) >= 20)
    path = tempname(); save_tree(path, pop, node)
    q = load_population(path)
    sfs = site_frequency_spectrum(q)
    @test sum(k * x for (k, x) in enumerate(sfs)) == sum(mutations_per_cell(q))
    @test sfs[end] == clonal_mutations(q)
    burden = Dict(zip(cell_ids(pop), mutations_per_cell(pop)))
    @test mutations_per_cell(q) == [burden[id] for id in cell_ids(q)]
    @test q.tree.neutral[1] == q.tree.total[1] == node.total
    depth = let j = ClonalEvolutionTx._index(node), d = 0
        while pop.tree.parent[j] != 0; j = pop.tree.parent[j]; d += 1; end
        d
    end
    divs = Dict(zip(cell_ids(pop), divisions_per_cell(pop)))
    @test divisions_per_cell(q) == [divs[id] - depth for id in cell_ids(q)]   # from the subtree root
    repaint_neutral!(q, 0.0, StableRNG(1))       # roots keep theirs: the inherited count stays
    @test all(==(node.total + q.clonal_mutations), mutations_per_cell(q))
end
