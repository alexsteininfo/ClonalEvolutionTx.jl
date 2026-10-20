# Concepts

## Cells, clones and traits

A **cell** carries only a clone index plus bookkeeping: a unique id, its birth time, the
number of divisions in its ancestry and its neutral mutations. Cells within a clone are
exchangeable.

A **clone** carries the heritable traits ([`Traits`](@ref)):

- the **fitness** `s`, a finite number of either sign that the user's `birth` and `death`
  functions turn into rates, and
- the **drug resistance** `α ∈ [0, 1]`. `α = 1` is fully drug-sensitive, `α < 1` is
  "resistant", `α = 0` is fully resistant.

A new clone is created when a daughter cell is hit by a driver, resistance or custom
mutation (see [Mutations](@ref)). Clones are numbered in order of creation and never
removed: an extinct clone stays in the [`CloneTable`](@ref) with size 0, so the history of
the clone structure is kept, but it leaves the list of living clones. Clone-level
statistics ([`clone_sizes`](@ref), [`hill_number`](@ref) of the clone sizes) read the living
clones.

The **environment** is the population size `N`, the time `t` and the treatment level
`m ∈ [0, 1]`. The treatment level belongs to the block or schedule, never to a cell.

## Rates

For a cell of clone `k`, in a population of `N` cells at treatment level `m`:

```
birth rate   b_k(N)    = β(s_k) · B(N)
death rate   d_k(N, m) = δ(s_k) + C(N) + α_k · D(m)
```

`β(s)` and `δ(s)` are user functions of the fitness, evaluated once per clone when it is
created. `B(N)` and `C(N)` are user functions of the population size shared by all clones
(defaults `B ≡ 1`, `C ≡ 0`), and `D(m)` is the drug response. A negative rate is an error.
See [The simulation block](@ref) and [Treatment and schedules](@ref).

## Tree mode and cell mode

A population is created by [`initialize_population`](@ref) in one of two storage modes.

- `mode = :tree` (the default) gives a [`TreePopulation`](@ref): the full lineage tree of
  the survivors in a compact struct-of-arrays [`LineageTree`](@ref). It supports every
  statistic, including the site frequency spectrum, pairwise distances and Newick export.
  Cells whose lineage has died out are pruned from the tree.
- `mode = :cells` gives a [`CellPopulation`](@ref): only a flat vector of living cells.
  It is faster and uses less memory, and supports burdens, divisions per cell and all clone
  statistics, but no genealogy.

Both modes keep the living cells in the same order and consume the same random numbers, so
one seed gives identical living cells in both. Every initial cell is its own root (the
tree is a forest of founders).

```@example concepts
using ClonalEvolutionTx, Random

block = BirthDeathBlock(birth = s -> 1.0, death = s -> 0.5, μ = 1.0,
                        stop_above = 300, restart_on_extinction = true)

tree_pop = initialize_population(1; mode = :tree)
cell_pop = initialize_population(1; mode = :cells)
simulate!(tree_pop, block, Xoshiro(3))
simulate!(cell_pop, block, Xoshiro(3))

(cell_ids(tree_pop) == cell_ids(cell_pop),
 mutations_per_cell(tree_pop) == mutations_per_cell(cell_pop))
```

The per-cell burdens agree as well: sisters share their parent's burden, so burdens are
correct jointly and not just marginally.

## Random streams

A population owns two random streams.

1. The **dynamics stream** is the `rng` that you pass to [`simulate!`](@ref). It drives the
   waiting times, the choice of cell, the event type and the mutation effects.
2. The **neutral stream** (`pop.neutral_rng`) draws the Poisson numbers of neutral
   mutations. It is seeded from `rng` at the start of each block.

Because they are separate, the same seed gives the same living cells for any neutral
mutation rate `μ` (including `μ = 0`) and in both storage modes, and measurements at
different times of one run see the same mutations on shared branches.

```@example concepts
a = initialize_population(1; mode = :cells)
b = initialize_population(1; mode = :cells)
simulate!(a, remake(block; μ = 0.0), Xoshiro(5))
simulate!(b, remake(block; μ = 2.5), Xoshiro(5))
(cell_clones(a) == cell_clones(b), a.t == b.t)
```

There is no global RNG and no global mutable state: one population belongs to one task, so
runs are thread-safe. For ensembles, [`derive_seed`](@ref)`(base, i)` turns a base seed and
an index into independent-looking seeds, stable across Julia versions.

## Initial conditions

[`initialize_population`](@ref) accepts

- `n` founders with the same traits, `initialize_population(n; s, α)`;
- a clone specification, one founder clone per entry, e.g.
  `[(n = 900, s = 0.0, α = 1.0), (n = 100, s = 0.0, α = 0.8)]`;
- any population produced by an earlier `simulate!` (chaining, see
  [The simulation block](@ref)) or loaded from a file (see [Sampling](@ref)).

The clones of a specification are siblings of each other in the clone table (no nesting).
`clonal_mutations` is a population-wide count of neutral mutations carried by every cell,
for instance mutations acquired before the simulated tumour began.

```@example concepts
pop = initialize_population([(n = 90, s = 0.0, α = 1.0), (n = 10, s = 0.0, α = 0.8)];
                            mode = :cells, clonal_mutations = 50)
(popsize(pop), clone_sizes(pop), n_sensitive(pop), n_resistant(pop), mean_burden(pop))
```
