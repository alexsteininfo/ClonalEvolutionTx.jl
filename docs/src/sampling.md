# Sampling

[`sample_population`](@ref)`(pop, n; seed)` draws `n` of the living cells uniformly without
replacement and returns them as a new population of the same kind; `pop` is left untouched.
All statistics apply to the sample unchanged.

```@example samp
using ClonalEvolutionTx, Random

block = BirthDeathBlock(birth = s -> 1.0, death = s -> 0.5, μ = 1.0,
                        stop_above = 1_000, restart_on_extinction = true)
pop = initialize_population(1)
simulate!(pop, block, Xoshiro(1))

s100 = sample_population(pop, 100; seed = 1)
(popsize(s100), round(mean_burden(s100), digits = 1), round(mean_burden(pop), digits = 1))
```

## Reproducibility

The draw is a partial Fisher–Yates shuffle of the living-cell order driven by
`StableRNG(seed)`, so it is a pure function of `(pop, n, seed)`, stable across Julia versions.
For one seed the draws at different `n` are **nested**: the `n`-draw is the first `n` cells of
any larger draw. Distinct draws need distinct seeds, which [`derive_seed`](@ref)`(base, i, ...)`
provides.

```@example samp
a = sample_population(pop, 20; seed = 7)
b = sample_population(pop, 50; seed = 7)
Set(cell_ids(a)) ⊆ Set(cell_ids(b))
```

## What a sample contains

- **Tree mode** gives the *induced tree*: the sampled leaves plus every ancestor of a sampled
  leaf, unary nodes kept (not collapsed), and the root of each sampled lineage stays a root
  (the founder is the root, not the sample's most recent common ancestor). Ids, birth times,
  clones and mutation counts are those of the source cells, so a sampled cell's burden and
  division count equal its full-tree values.
- **Cell mode** gives the drawn cells, copied as values.
- Either way the clone table keeps every clone (also those with no cell in the sample), with
  sizes recounted over the sample, and `t`, `m`, `clonal_mutations` and the log are copied.

`n` must satisfy `1 ≤ n ≤ popsize(pop)`; otherwise an `ArgumentError` is thrown. A snapshot
that asks for a sample larger than the population does not throw: it records
`(n = n, stats = nothing)` for that size (see below).

## Samples inside measurements

A [`MeasurementSpec`](@ref) with `samples = [100, 1000]` evaluates every statistic of every
snapshot on samples of those sizes as well. For the `i`-th snapshot of the accumulator the
sample of size `n` is `sample_population(pop, n; seed = derive_seed(sample_seed, i, n))`.
When the population has fewer than `n` cells, for example at extinction, that entry is
`(n = n, stats = nothing)` and no sample is drawn; the other sizes keep their seeds.

```@example samp
spec = MeasurementSpec(triggers = [AtBlockEnd()], stats = [SFS(), MeanBurden()],
                       samples = [50, 200], sample_seed = 1)
acc = MeasurementAccumulator(spec)
grow = initialize_population(1)
simulate!(grow, block, Xoshiro(2); accumulator = acc)
sn = finalize_measurements(acc).snapshots[1]
[(s.n, length(s.stats.sfs)) for s in sn.samples]
```

## Saving and loading

[`save_population`](@ref)`(path, pop)` writes a whole population (tree or cell mode, with the clone
table) in a flat, versioned, little-endian binary format that does not depend on Julia's
`Serialization`; [`load_population`](@ref)`(path)` reads it back. The neutral stream and the log are
not stored. A population loaded from a file lists its living clones in clone-id order, so
`clone_sizes` may come in a different order than before saving (the multiset is the same).

A loaded population can be analysed, sampled or continued with `simulate!` with a new `rng`. A
continued run is **not** bit-identical to an uninterrupted one (see
[Limitations and open questions](@ref)).

[`save_tree`](@ref)`(path, pop, node)` takes the **population** and a [`CellNode`](@ref) of its tree
and writes the subtree under `node` as a single-root tree-mode file. The subtree root carries
all neutral mutations inherited from above, so in the loaded file they count at `k = n` in the
site frequency spectrum and in `clonal_mutations`, and every cell keeps its burden;
`divisions_per_cell` counts divisions from the subtree root.

```@example samp
path = joinpath(mktempdir(), "pop.cepop")
save_population(path, pop)
back = load_population(path)
(cell_ids(back) == cell_ids(pop), site_frequency_spectrum(back) == site_frequency_spectrum(pop),
 sort(clone_sizes(back)) == sort(clone_sizes(pop)))
```

```@example samp
root = first(roots(pop))
sub = joinpath(dirname(path), "subtree.cepop")
save_tree(sub, pop, root)
popsize(load_population(sub)) <= popsize(pop)
```

## Newick export

[`write_newick`](@ref)`(io_or_path, pop; branch_lengths = :divisions)` writes the genealogy of
the living cells as Newick (one tree per founder, one line each), for R/ape or ete3. Leaves are
labelled by cell id, chains of single-child nodes are collapsed, and each tree is rooted at the
most recent common ancestor of its cells, with a stem from the founder. A branch length is the
number of divisions along the chain (`:divisions`) or a time (`:time`). With `:time` every node
sits on the time axis: a branching node at its division, a leaf at `pop.t`, and the stem starts
at the founder's birth. The tree is then ultrametric, and the length from two leaves up to their
common ancestor is their [`coalescence_times`](@ref) entry.

```@example samp
io = IOBuffer()
write_newick(io, sample_population(pop, 5; seed = 3))
String(take!(io))
```
