# Uniform sampling of living cells (SPEC §8). The draw is a partial Fisher–Yates shuffle of
# the living-cell order driven by `StableRNG(seed)`, so one seed gives nested draws across
# `n`. Tree mode returns the induced subtree: the sampled leaves plus every ancestor, unary
# nodes kept, so a sampled cell's burden and depth are its full-tree values.

# The first `n` entries of a uniformly random permutation of `1:N` (partial Fisher–Yates).
# Uses only `rand(rng, a:b)`, whose stream `StableRNG` guarantees across Julia versions
# (unlike `randperm`). Changing this function invalidates every stored draw.
function _draw_indices(rng::StableRNG, N::Int, n::Int)
    perm = collect(1:N)
    for i in 1:n
        j = rand(rng, i:N)
        perm[i], perm[j] = perm[j], perm[i]
    end
    return perm[1:n]
end

function _checked_draw(pop::AbstractPopulation, n::Integer, seed::Integer)
    N = popsize(pop)
    1 <= n <= N || throw(ArgumentError("cannot draw n = $n cells from a population of $N cells"))
    return _draw_indices(StableRNG(seed % UInt64), N, Int(n))
end

# Copy the rows `kept` (sorted indices) into a new tree, in the same order. Left/right slots
# are preserved, so a node whose left lineage was dropped keeps no left child; a node whose
# parent was not kept (a root, or the top of a saved subtree) becomes a root.
function _copy_rows(tree::LineageTree, kept::Vector{Int32})
    new = LineageTree()
    sizehint!(new, length(kept))
    for (k, i) in enumerate(kept)
        p0 = tree.parent[i]
        j = p0 == NOPARENT ? 0 : searchsortedfirst(kept, p0)
        p = (j != 0 && j <= length(kept) && kept[j] == p0) ? Int32(j) : NOPARENT
        _push_node!(new, p, tree.id[i], tree.birthtime[i], tree.clone[i], tree.neutral[i],
                    tree.total[i])
        if p != NOPARENT
            tree.left[p0] == i ? (new.left[p] = Int32(k)) : (new.right[p] = Int32(k))
        end
    end
    return new
end

# A copy of `ct` whose sizes are the counts of the sampled cells' clones.
function _recounted_clones(ct::CloneTable, clones::AbstractVector{<:Integer})
    new = copy(ct)
    fill!(new.size, 0)
    for k in clones
        new.size[k] += 1
    end
    return _rebuild_derived!(new)
end

"""
    sample_population(pop, n; seed) -> population

Draw `n` of the living cells uniformly without replacement and return them as a new
population of the same kind. `n` must satisfy `1 ≤ n ≤ popsize(pop)` (`ArgumentError`
otherwise); `pop` is left untouched.

The draw is a partial Fisher–Yates shuffle of the living-cell order (`pop.alive` or
`pop.cells`) driven by `StableRNG(seed)`; it is a pure function of `(pop, n, seed)`, stable
across Julia versions. For one `seed`, draws at different `n` are nested: the `n`-draw is the
first `n` cells of any larger draw. Distinct draws need distinct seeds
([`derive_seed`](@ref) makes them). The sampled cells are in draw order.

A [`TreePopulation`](@ref) yields the induced tree: the sampled leaves plus every ancestor of
a sampled leaf, unary nodes kept (not collapsed), and the root of each sampled lineage stays
a root (so the founder is the root, not the sample's MRCA, and a forest stays a forest).
Ids, birth times, clones and mutation counts are those of the source cells, hence burdens and
divisions of the sampled cells equal their full-tree values and every statistic applies
unchanged. A [`CellPopulation`](@ref) yields the drawn cells, copied as values.

Either way the clone table keeps every clone (also those with no cell in the sample) with
sizes recounted over the sample, and `t`, `m`, `clonal_mutations`, `next_id`, the neutral
stream and the log are copied.
"""
function sample_population(pop::TreePopulation, n::Integer; seed::Integer)
    idx = _checked_draw(pop, n, seed)
    tree = pop.tree
    marked = falses(length(tree))
    kept = Int32[]
    for k in idx
        i = pop.alive[k]
        while !marked[i]
            marked[i] = true
            push!(kept, i)
            i = tree.parent[i]
            i == NOPARENT && break
        end
    end
    sort!(kept)
    newtree = _copy_rows(tree, kept)
    alive = Int32[searchsortedfirst(kept, pop.alive[k]) for k in idx]
    clones = _recounted_clones(pop.clones, Int32[newtree.clone[i] for i in alive])
    return TreePopulation(newtree, alive, clones, pop.t, pop.m, pop.clonal_mutations,
                          pop.next_id, copy(pop.neutral_rng), copy(pop.log))
end

function sample_population(pop::CellPopulation, n::Integer; seed::Integer)
    idx = _checked_draw(pop, n, seed)
    cells = Cell[pop.cells[k] for k in idx]
    clones = _recounted_clones(pop.clones, Int32[c.clone for c in cells])
    return CellPopulation(cells, clones, pop.t, pop.m, pop.clonal_mutations, pop.next_id,
                          copy(pop.neutral_rng), copy(pop.log))
end
