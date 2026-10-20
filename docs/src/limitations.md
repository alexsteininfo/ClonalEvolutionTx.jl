# Limitations and open questions

## Out of scope

- Non-exponential cell-cycle timing: use
  [NonMarkovEvolution.jl](https://github.com/alexsteininfo/NonMarkovEvolution.jl).
- Spatial structure, inference of any kind (this package simulates; it does not fit), and
  sequencing-noise models.
- Ensemble runners, result-file formats and scenario files. The package only makes them easy:
  all randomness comes from the passed RNGs, [`derive_seed`](@ref) derives per-run seeds,
  populations are thread-safe, and measurement output is small and table-friendly.
- Cluster or SLURM tooling.
- Paper-specific parameter sets, sweeps, results, figures and power-law fits.
- Planned for a later milestone and not precluded by the design: cytostatic treatment (a drug
  acting on birth), clone-specific competition coefficients, a Moran block, a time-based
  neutral mutation clock and a rejection-free algorithm.

## Known constraints

- **Common density terms only.** Density dependence enters through `B(N)` and `C(N)`, which are
  shared by all clones. Clone-specific competition coefficients (`l_S ≠ l_R`) are not covered.
- **Cytotoxic treatment only.** The drug acts on the death rate, `α·D(m)`.
- **No shared founder mutations.** A clone specification sets sizes and traits only; every
  initial cell is its own root, and the founder clones are siblings in the clone table. Shared
  neutral history comes only from `clonal_mutations`. To start from a realistic genealogy, grow
  a population with a block and chain.
- **Continued runs after loading are not bit-identical.** A population read by
  [`load_population`](@ref) gets a fresh neutral stream and an empty log, and a continued run
  uses a new `rng`, so it is statistically equivalent to, but not identical with, the
  uninterrupted run.
- A block that can never stop (`restart_on_extinction` for a population certain to die out)
  restarts forever; a block with no possible event and `tmax = Inf` throws. A schedule with
  the default `t_max = Inf` (`MTD`, `Metronomic`, `CustomDose`) whose population can neither
  reach `N_prog` nor die out runs forever: give it a finite `t_max`.
- Pairwise distances and coalescence times are quadratic in the population size; use a sample.

## Conventions

- **`tmax` versus `t_max`.** A block's `tmax` is an absolute time. The `t_max` of every schedule
  is measured from the schedule start, `t0 + t_max`, **including `Stabilising`**. This differs
  from BirthDeathMutation_AT, which measured the maintenance `t_max` of its stabilising ("maintenance")
  schedule from the start of the maintenance phase. Metronomic intervals are relative to the
  schedule start as well.
- **Stops are inclusive:** a block stops at `N ≥ stop_above` and `N ≤ stop_below`, and the
  built-in schedules progress at `N ≥ N_prog`.
- **Two Hill-number definitions for an SFS.** `:normalised` (default, a proper distribution)
  and `:cell_fraction` (the v1 papers' definition, not a proper Hill family); see
  [Tree statistics](@ref). Which one a paper uses is decided there, not here.
- **`n_sensitive` means `α == 1`;** every cell in a clone with `α < 1` counts as resistant.
- **The Hill number of an empty distribution is 0**, and the mean burden of an empty population
  is `NaN`.
- **A sample larger than the population is `nothing`.** A snapshot that asks for a sample of
  `n` cells from a population with fewer cells (for example at extinction) records
  `(n = n, stats = nothing)` instead of throwing.
- **Newick time lengths** place every branching node at its division time and every leaf at
  `pop.t`, so a `:time` tree is ultrametric and agrees with [`coalescence_times`](@ref).
- **Ordering.** Loaded or sampled populations list their living clones in clone-id order, and
  a sample lists its cells in draw order.
- **Resistance only increases.** The child's `α` is `min(α_drawn, α_parent)`.

## Open questions

- Whether a rejection-free algorithm is needed depends on the acceptance fractions measured on
  the [Performance](@ref) page.
- The performance targets (cell mode within about 0.5 µs per created cell at 10⁶ cells; tree
  mode within about 1.5× of NonMarkovEvolution.jl) are to be verified by the benchmark.
