# The per-event primitives of the thinning loop: starting a block, dividing a cell and
# killing one. All three are allocation-free in the steady state, and a division draws
# the same random numbers in tree and cell mode.

# Treatment level of the population's current state under the block's dose (at block
# start, and after every event for the state-dependent doses). B(N), C(N) are checked.
_dose_level(d::FixedDose, block::BirthDeathBlock, pop) = d.m
_dose_level(d::SizeDose, block::BirthDeathBlock, pop) = _size_level(d, popsize(pop))
_dose_level(d::TimeDose, block::BirthDeathBlock, pop) = _time_level(d, 0.0)
function _dose_level(d::StabilisingDose, block::BirthDeathBlock, pop)
    N = popsize(pop)
    N == 0 && return 0.0                     # nothing to stabilise; avoids calling B(0), C(0)
    ct = pop.clones
    return _stabilising_level(block.response, d.m_max;
                              B = _density(block.birth_density, N, "birth_density B(N)"),
                              C = _density(block.death_density, N, "death_density C(N)"),
                              N, sumNβ = ct.sumNβ, sumNδ = ct.sumNδ, sumNα = ct.maxα == 0.0 ? 0.0 : ct.sumNα)   # maxα is exact; the running sum drifts
end

function _start_block!(pop::AbstractPopulation, block::BirthDeathBlock, rng::AbstractRNG)
    _bind_rates!(pop.clones, block.birth, block.death)
    pop.neutral_rng = Xoshiro(rand(rng, UInt64))
    pop.m = _dose_level(block.dose, block, pop)
    return pop
end

# Apply each mutation kind in tuple order to a daughter; every new clone nests under the
# daughter's current clone. No random number is drawn once `max_clones` is reached.
@inline _mutate_daughter!(ct::CloneTable, ::Tuple{}, rng, max_clones, t, id, cur) = cur
@inline function _mutate_daughter!(ct::CloneTable, muts::Tuple, rng, max_clones, t, id, cur)
    kind = first(muts)
    if n_clones(ct) < max_clones && rand(rng) < kind.ν
        cur = _add_clone!(ct, cur, mutate(kind, traits(ct, cur), rng), mutation_kind(kind),
                          t, id)
    end
    return _mutate_daughter!(ct, Base.tail(muts), rng, max_clones, t, id, cur)
end

@inline _neutral_count(block::BirthDeathBlock, rng) =
    block.μ == 0.0 ? 0 : Int(rand(rng, Poisson(block.μ)))

# Draw one daughter: neutral count, then mutations. Returns (final clone, neutral count).
@inline function _draw_daughter!(pop, block, t, rng, id, parent_clone)
    count = _neutral_count(block, pop.neutral_rng)
    clone = _mutate_daughter!(pop.clones, block.mutations, rng, block.max_clones, t, id,
                              parent_clone)
    return clone, count
end

# Parent clone −1 plus the daughters that stayed in it, in one call (so a size-1 clone
# whose daughter stays is never seen at size 0); then +1 for each daughter's new clone.
@inline function _update_sizes!(ct::CloneTable, pc, c1, c2)
    net = -1 + (c1 == pc) + (c2 == pc)
    net == 0 || _change_size!(ct, pc, net)
    c1 == pc || _change_size!(ct, c1, 1)
    c2 == pc || _change_size!(ct, c2, 1)
    return nothing
end

function _divide!(pop::TreePopulation, k::Int, t::Float64, block::BirthDeathBlock, rng::AbstractRNG)
    tree = pop.tree
    node = pop.alive[k]
    pc = tree.clone[node]
    id1 = pop.next_id
    id2 = id1 + 1
    c1, n1 = _draw_daughter!(pop, block, t, rng, id1, pc)
    c2, n2 = _draw_daughter!(pop, block, t, rng, id2, pc)
    pop.next_id = id1 + 2
    _update_sizes!(pop.clones, pc, c1, c2)
    tot = Int(tree.total[node])
    child1 = _push_node!(tree, node, id1, t, c1, n1, tot + n1)
    child2 = _push_node!(tree, node, id2, t, c2, n2, tot + n2)
    tree.left[node] = child1
    tree.right[node] = child2
    pop.alive[k] = child1
    push!(pop.alive, child2)
    return nothing
end

function _divide!(pop::CellPopulation, k::Int, t::Float64, block::BirthDeathBlock, rng::AbstractRNG)
    parent = pop.cells[k]
    pc = parent.clone
    id1 = pop.next_id
    id2 = id1 + 1
    c1, n1 = _draw_daughter!(pop, block, t, rng, id1, pc)
    c2, n2 = _draw_daughter!(pop, block, t, rng, id2, pc)
    pop.next_id = id1 + 2
    _update_sizes!(pop.clones, pc, c1, c2)
    pop.cells[k] = Cell(id1, t, c1, parent.divisions + Int32(1), parent.neutral + n1)
    push!(pop.cells, Cell(id2, t, c2, parent.divisions + Int32(1), parent.neutral + n2))
    return nothing
end

# Compaction thresholds: compact once the dead rows are at least this many and at least
# `_COMPACT_FRACTION` of the tree.
const _COMPACT_MIN = 4096
const _COMPACT_FRACTION = 0.5

function _die!(pop::TreePopulation, k::Int)
    tree = pop.tree
    alive = pop.alive
    node = alive[k]
    _change_size!(pop.clones, tree.clone[node], -1)
    alive[k] = alive[end]
    pop!(alive)
    _prune!(tree, node)
    if tree.ndead >= _COMPACT_MIN && tree.ndead >= _COMPACT_FRACTION * length(tree)
        _compact_population!(pop)
    end
    return nothing
end

# Rare: drop the removed rows and translate the living-cell indices.
@noinline function _compact_population!(pop::TreePopulation)
    remap = _compact!(pop.tree, Int32(1))
    alive = pop.alive
    @inbounds for i in eachindex(alive)
        alive[i] = remap[alive[i]]
    end
    return nothing
end

function _die!(pop::CellPopulation, k::Int)
    cells = pop.cells
    _change_size!(pop.clones, cells[k].clone, -1)
    cells[k] = cells[end]
    pop!(cells)
    return nothing
end
