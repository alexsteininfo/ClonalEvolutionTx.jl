# ClonalEvolutionTx.jl

Stochastic, **Markovian** (exponential waiting time) birth–death simulation of cell
populations in which **clones** are tracked. Every clone carries two traits, a fitness `s`
and a drug resistance `α`, and **treatment** is a built-in, composable part of the dynamics:
a drug level `m(t)` raises the death rate of every cell in proportion to its sensitivity
`α`, and the level can be fixed, follow a schedule, depend on the population size or on
time, or be set after every event so that the tumour just stops growing. Neutral mutations
are recorded either on the full lineage tree of the survivors (*tree mode*) or as per-cell
counts on the living cells only (*cell mode*).

The simulation is exact: waiting times are drawn by thinning (uniformisation) against a
rate bound that is cheap to maintain and never too small, stops never overshoot, and any
rate that would be negative is an error rather than a silent clamp. Statistics (site
frequency spectra, clone sizes, Hill numbers, ...) can be computed while the population is
in memory, on the whole population or on small samples, so that only small results need to
be stored.

This package is a **simulator**. It does not infer anything from data, fit parameters, or
model sequencing noise. See [Limitations and open questions](@ref) for what is out of scope.

## Why it exists

ClonalEvolutionTx.jl is the clone-tracking, Markovian, low-mutation-rate counterpart of
[NonMarkovEvolution.jl](https://github.com/alexsteininfo/NonMarkovEvolution.jl) (same
author), which covers non-exponential cell-cycle times, per-cell fitness and high mutation
rates without clone labels. The two read alike: blocks, `initialize_population`,
`simulate!`, measurement specs, `restart_on_extinction` and sampling follow the same
conventions.

The code descends from
[SomaticEvolution.jl](https://github.com/jessierenton/SomaticEvolution.jl) (Jessie
Renton, MIT) via BirthDeathMutation.jl. It is the modelling package for the manuscript
*The dynamics of genetic diversity during evolutionary therapies* (Stein, Staňková,
Werner). Paper-specific parameter sets, sweeps and figures are not part of the package.

## Installation

Requires Julia 1.10 or newer.

```julia
using Pkg
Pkg.add(url = "https://github.com/alexsteininfo/ClonalEvolutionTx.jl")
```

## Quickstart

Grow a tumour from one cell to 2000 cells, then treat it with adaptive therapy (full dose
until the population has fallen to 1000 cells, no drug until it is back at 2000, and so
on), recording a snapshot at every treatment switch.

```@example quickstart
using ClonalEvolutionTx, Random

base = BirthDeathBlock(
    birth     = s -> 1.0,                              # β(s): birth rate 1
    death     = s -> 0.5,                              # δ(s): death rate 0.5
    response  = LinearResponse(4.0),                   # drug death rate D(m) = 4m (times α)
    mutations = [Resistance(ν = 0.002, α = 0.8)],       # 0.2% of daughters become resistant
    μ         = 1.0,                                   # neutral mutations per daughter
    label     = "growth",
)

rng = Xoshiro(2)
pop = initialize_population(1)                         # one founder, tree mode
simulate!(pop, remake(base; stop_above = 2_000, restart_on_extinction = true), rng)
(popsize(pop), round(pop.t, digits = 2))
```

```@example quickstart
schedule = AdaptiveOnOff(base; N_low = 1_000, N_up = 2_000, N_prog = 4_000, t_max = 5.0)
acc = MeasurementAccumulator(MeasurementSpec(
    triggers = [AtBlockEnd()],
    stats    = [popsize, n_resistant, MeanBurden()],
))
simulate!(pop, schedule, rng; accumulator = acc)
m = finalize_measurements(acc)
schedule_outcome(pop)
```

The snapshots, one per block:

```@example quickstart
[(s.label, round(s.t, digits = 2), s.stats.popsize, s.stats.n_resistant) for s in m.snapshots]
```

The treatment log of the population has one entry per block:

```@example quickstart
[(b.label, round(b.t_start, digits = 2), round(b.t_end, digits = 2), b.end_reason) for b in pop.log]
```

## Where to go next

- [Concepts](@ref): cells, clones, traits, the two storage modes and the random streams.
- [The simulation block](@ref): rates, density dependence, stops, restarts and chaining.
- [Mutations](@ref) and [Treatment and schedules](@ref): the model's moving parts.
- [Measurements and output](@ref), [Tree statistics](@ref) and [Sampling](@ref): getting
  results out.
- The `examples/` folder holds eight runnable scripts, from a single-cell expansion to
  cyclic and stabilising therapy.
