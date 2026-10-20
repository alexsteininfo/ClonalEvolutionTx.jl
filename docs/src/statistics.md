# Tree statistics

Statistics from neutral mutations need **tree mode** (`mode = :tree`); in cell mode they throw
an `ArgumentError`. Clone-level statistics work in both modes.

```@example stat
using ClonalEvolutionTx, Random

block = BirthDeathBlock(birth = s -> 1.0, death = s -> 0.5, μ = 1.0,
                        mutations = [Resistance(ν = 0.02, α = 0.8)],
                        response = LinearResponse(4.0),
                        stop_above = 400, restart_on_extinction = true)
pop = initialize_population(1; clonal_mutations = 5)
simulate!(pop, block, Xoshiro(1))
popsize(pop)
```

## Site frequency spectrum and burdens

[`site_frequency_spectrum`](@ref)`(pop)` returns `sfs`, where `sfs[k]` is the number of neutral
mutations carried by exactly `k` living cells; the vector has length `N` and the
`clonal_mutations` are counted at `k = N`. [`total_mutational_burden`](@ref) is the number of
distinct mutations, `sum(sfs)`. [`mutations_per_cell`](@ref) (both modes) is the burden of each
living cell, including the clonal mutations, and [`mean_burden`](@ref) its mean.
[`divisions_per_cell`](@ref) is the number of divisions in each cell's ancestry (the leaf
depth in tree mode).

```@example stat
sfs = site_frequency_spectrum(pop)
(sfs[1], sfs[end], total_mutational_burden(pop), round(mean_burden(pop), digits = 1),
 clonal_mutations(pop))
```

[`clonal_mutations`](@ref)`(pop)` is the number of mutations carried by every living cell: the
population's `clonal_mutations` plus the path total of the most recent common ancestor of all
living cells.

## Distances and coalescence times

[`pairwise_distances`](@ref) returns, for every pair of living cells `i < j`, the number of
differing neutral mutations, and [`coalescence_times`](@ref) the time from `pop.t` back to the
division of their most recent common ancestor. Both are quadratic in the population size; use a
sample for large populations (see [Sampling](@ref)).

```@example stat
small = sample_population(pop, 50; seed = 1)
d = pairwise_distances(small)
c = coalescence_times(small)
(length(d), round(sum(d) / length(d), digits = 1), round(maximum(c), digits = 2))
```

## Clone statistics

- [`clone_sizes`](@ref): sizes of the living clones (non-overlapping), in
  [`living_clones`](@ref) order.
- [`n_sensitive`](@ref), [`n_resistant`](@ref): cells in clones with `α == 1` and with
  `α < 1`.
- [`mutation_frequencies`](@ref): for every non-founder clone, in id order, its size plus the
  sizes of all descendant clones: the number of cells that carry that driver or resistance
  mutation.

```@example stat
(clone_sizes(pop), n_sensitive(pop), n_resistant(pop), mutation_frequencies(pop))
```

## Hill numbers and Rényi entropies

[`hill_number`](@ref)`(counts, q)` is the effective number of types of the distribution
`p = counts / sum(counts)`: `q = 0` the number of nonzero entries, `q = 1` `exp(Shannon)`,
`q = 2` the inverse Simpson index, `q = Inf` `1 / max p` and `q = -Inf` `1 / min p` (the
limits of `q → ±∞`, which large finite `|q|` approach). The Hill number of an **empty
distribution** (all counts zero) is `0` for every `q`. [`renyi_entropy`](@ref) is its
logarithm. [`HillNumbers`](@ref)`(q = [0, 1, 2, Inf])` is the statistic of the living clone
sizes.

```@example stat
(hill_number([50, 30, 20], 0), hill_number([50, 30, 20], 2), hill_number([0, 0], 1),
 round(renyi_entropy([50, 30, 20], 1), digits = 3))
```

### Hill numbers of an SFS: two definitions

[`sfs_hill_number`](@ref)`(sfs, q; frequencies = :normalised)` treats the mutations of a site
frequency spectrum as the types, under one of two definitions chosen by an explicit keyword.

- `frequencies = :normalised` (the default): `p = k / A` with `A = Σₖ k·Sₖ` (all mutation
  copies). A proper probability distribution, so the family is consistent for every `q`
  (`q → 1` gives Shannon), `∞D = A / k_max` and `₋∞D = A / k_min`.
- `frequencies = :cell_fraction`: `p = k / N`, as in the v1 papers and in
  BirthDeathMutation_AT's `genetic_hill`. For `q = 1` this is `exp(−Σ Sₖ (k/N) ln(k/N))`, and
  `∞D = N / k_max`, `₋∞D = N / k_min`. Not a proper Hill family: the `p` sum to the mean burden.

`q = 0` is the number of mutations in both. With the mean burden `m̄ = A/N`, the conversions are
`²D_norm = m̄²·²D_cell` and `H₁,norm = H₁,cell/m̄ + ln m̄`. [`sfs_renyi_entropy`](@ref) is the
logarithm, and [`SFSHillNumbers`](@ref)`(; q, frequencies)` the statistic (tree mode).

```@example stat
D_norm = sfs_hill_number(sfs, 2)
D_cell = sfs_hill_number(sfs, 2; frequencies = :cell_fraction)
A = sum(k * sfs[k] for k in eachindex(sfs))      # all mutation copies
m̄ = A / popsize(pop)                              # mean burden, clonal mutations included
(round(D_norm, digits = 3), round(D_cell, digits = 3), D_norm ≈ m̄^2 * D_cell)
```

## Re-painting

[`repaint_neutral!`](@ref)`(pop, μ′, rng)` redraws the neutral mutations of a tree-mode population
with a new rate `μ′`, leaving the dynamics, clones and `clonal_mutations` untouched.
