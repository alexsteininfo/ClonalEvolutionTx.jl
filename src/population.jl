# Populations: the state `simulate!` evolves and every statistic reads. A tree population
# keeps the lineage tree; a cell population keeps only a flat vector of cells.

"""
    Cell(id, birthtime, clone, divisions, neutral)

One living cell of a [`CellPopulation`](@ref): its unique `id`, `birthtime`, `clone`
index, number of `divisions` in its ancestry and running total of `neutral` mutations.
"""
struct Cell
    id::Int64
    birthtime::Float64
    clone::Int32
    divisions::Int32
    neutral::Int32
end

"""
    BlockSummary

Log entry of one simulated block: `label`, start and end time (`t_start`, `t_end`),
`end_reason` (`:stopped`, `:above`, `:below`, `:extinct` or `:tmax`), the number of
`ticks`, `accepted` events and `restarts`, the `dose` description and the drug level
`m_end` at the end (for a [`TimeDose`](@ref), its level at `t_end`).
"""
struct BlockSummary
    label::String
    t_start::Float64
    t_end::Float64
    end_reason::Symbol
    ticks::Int
    accepted::Int
    restarts::Int
    dose::String
    m_end::Float64
end

"""
    AbstractPopulation

Supertype of [`TreePopulation`](@ref) and [`CellPopulation`](@ref).
"""
abstract type AbstractPopulation end

"""
    TreePopulation

A population that records the full lineage tree. Fields: `tree`, `alive` (node index of
each living cell, in sampling order), `clones`, time `t`, drug level `m`,
`clonal_mutations` shared by all cells, `next_id`, the `neutral_rng` stream and the
block `log`. Create with [`initialize_population`](@ref).
"""
mutable struct TreePopulation <: AbstractPopulation
    tree::LineageTree
    alive::Vector{Int32}
    clones::CloneTable
    t::Float64
    m::Float64
    clonal_mutations::Int
    next_id::Int64
    neutral_rng::Xoshiro
    log::Vector{BlockSummary}
end

"""
    CellPopulation

A population that keeps only its living cells (no genealogy). Fields: `cells`, `clones`,
time `t`, drug level `m`, `clonal_mutations`, `next_id`, the `neutral_rng` stream and the
block `log`. Create with [`initialize_population`](@ref).
"""
mutable struct CellPopulation <: AbstractPopulation
    cells::Vector{Cell}
    clones::CloneTable
    t::Float64
    m::Float64
    clonal_mutations::Int
    next_id::Int64
    neutral_rng::Xoshiro
    log::Vector{BlockSummary}
end

"""
    initialize_population(n = 1; mode = :tree, s = 0.0, α = 1.0, clonal_mutations = 0)
    initialize_population(spec; mode = :tree, clonal_mutations = 0)

A population at time 0 of `n` founder cells with traits `(s, α)` (one founder clone), or
one founder clone per entry `(n, s, α)` of the vector `spec`. `mode` is `:tree` (a
[`TreePopulation`](@ref)) or `:cells` (a [`CellPopulation`](@ref)). Every initial cell is
its own root with no neutral mutations of its own; `clonal_mutations` are carried by all.
"""
function initialize_population(n::Integer = 1; mode::Symbol = :tree, s::Real = 0.0,
                               α::Real = 1.0, clonal_mutations::Integer = 0)
    return initialize_population([(n = n, s = s, α = α)]; mode, clonal_mutations)
end

function initialize_population(spec::AbstractVector{<:NamedTuple}; mode::Symbol = :tree,
                               clonal_mutations::Integer = 0)
    mode in (:tree, :cells) || throw(ArgumentError("mode must be :tree or :cells, got :$mode"))
    clonal_mutations >= 0 || throw(ArgumentError("clonal_mutations must be ≥ 0, got $clonal_mutations"))
    isempty(spec) && throw(ArgumentError("the clone specification is empty"))
    for e in spec
        e.n >= 1 || throw(ArgumentError("each clone needs n ≥ 1 cells, got n = $(e.n)"))
    end
    clones = CloneTable()
    ntotal = sum(e -> Int(e.n), spec)
    for e in spec
        k = _add_clone!(clones, 0, Traits(e.s, e.α), :founder, 0.0, Int64(0))
        _change_size!(clones, k, Int(e.n))
    end
    cloneof = Vector{Int32}(undef, ntotal)
    pos = 0
    for (k, e) in enumerate(spec), _ in 1:Int(e.n)
        cloneof[pos += 1] = Int32(k)
    end
    rng = Xoshiro(0)
    log = BlockSummary[]
    if mode === :tree
        tree = LineageTree()
        sizehint!(tree, ntotal)
        alive = Vector{Int32}(undef, ntotal)
        for i in 1:ntotal
            alive[i] = _push_node!(tree, NOPARENT, i, 0.0, cloneof[i], 0, 0)
        end
        return TreePopulation(tree, alive, clones, 0.0, 0.0, Int(clonal_mutations),
                              Int64(ntotal + 1), rng, log)
    else
        cells = [Cell(i, 0.0, cloneof[i], 0, 0) for i in 1:ntotal]
        return CellPopulation(cells, clones, 0.0, 0.0, Int(clonal_mutations),
                              Int64(ntotal + 1), rng, log)
    end
end

"""
    popsize(pop) -> Int

Number of living cells.
"""
popsize(pop::TreePopulation) = length(pop.alive)
popsize(pop::CellPopulation) = length(pop.cells)

"""
    cell_ids(pop) -> Vector{Int64}

Ids of the living cells, in living-cell order.
"""
cell_ids(pop::TreePopulation) = Int64[pop.tree.id[i] for i in pop.alive]
cell_ids(pop::CellPopulation) = Int64[c.id for c in pop.cells]

"""
    cell_clones(pop) -> Vector{Int32}

Clone index of each living cell, in living-cell order.
"""
cell_clones(pop::TreePopulation) = Int32[pop.tree.clone[i] for i in pop.alive]
cell_clones(pop::CellPopulation) = Int32[c.clone for c in pop.cells]

"""
    mutations_per_cell(pop) -> Vector{Int}

Neutral mutations carried by each living cell, including the `clonal_mutations`.
"""
mutations_per_cell(pop::TreePopulation) =
    Int[Int(pop.tree.total[i]) + pop.clonal_mutations for i in pop.alive]
mutations_per_cell(pop::CellPopulation) =
    Int[Int(c.neutral) + pop.clonal_mutations for c in pop.cells]

"""
    divisions_per_cell(pop) -> Vector{Int}

Number of divisions in the ancestry of each living cell.
"""
function divisions_per_cell(pop::TreePopulation)
    parent = pop.tree.parent
    return Int[let d = 0, j = i
                   while parent[j] != NOPARENT
                       j = parent[j]; d += 1
                   end
                   d
               end for i in pop.alive]
end
divisions_per_cell(pop::CellPopulation) = Int[Int(c.divisions) for c in pop.cells]

"""
    roots(pop::TreePopulation) -> Vector{CellNode}

The roots of the lineage trees of `pop` that have not been removed.
"""
roots(pop::TreePopulation) = CellNode[CellNode(pop.tree, i) for i in _rootidxs(pop.tree)]

"""
    copy(pop)

Deep copy of a population (tree, living cells, clone table, log and neutral stream);
simulating the copy leaves the original untouched.
"""
Base.copy(pop::TreePopulation) =
    TreePopulation(copy(pop.tree), copy(pop.alive), copy(pop.clones), pop.t, pop.m,
                   pop.clonal_mutations, pop.next_id, copy(pop.neutral_rng), copy(pop.log))
Base.copy(pop::CellPopulation) =
    CellPopulation(copy(pop.cells), copy(pop.clones), pop.t, pop.m, pop.clonal_mutations,
                   pop.next_id, copy(pop.neutral_rng), copy(pop.log))
