# Shared test fixtures.

# A population grown from one cell to `N` cells, with neutral and resistance mutations.
grown(mode; N = 500, seed = 1, μ = 1.0, ν = 0.01) =
    simulate!(initialize_population(1; mode),
              BirthDeathBlock(birth = s -> 1.0, death = s -> 0.3, μ = μ,
                              mutations = (Resistance(ν = ν, α = 0.5),), stop_above = N,
                              restart_on_extinction = true),
              StableRNG(seed))

# A hand-built tree population: root 1 (neutral 2, total 2) with two living children,
# `a` = id 2 (neutral 1, total 3) and `b` = id 3 (neutral 3, total 5); one clonal mutation.
function small_tree()
    pop = initialize_population(1)
    tree = ClonalEvolutionTx.LineageTree()
    push = ClonalEvolutionTx._push_node!
    r = push(tree, Int32(0), 1, 0.0, 1, 2, 2)
    a = push(tree, r, 2, 1.0, 1, 1, 3)
    b = push(tree, r, 3, 1.0, 1, 3, 5)
    tree.left[r], tree.right[r] = a, b
    pop.tree = tree
    pop.alive = Int32[a, b]
    pop.next_id = 4
    pop.clonal_mutations = 1
    ClonalEvolutionTx._change_size!(pop.clones, 1, 1)
    return pop
end
