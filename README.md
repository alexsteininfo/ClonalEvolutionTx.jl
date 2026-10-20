# ClonalEvolutionTx.jl

Stochastic, **Markovian** (exponential waiting time) birth–death simulation of cell
populations in which **clones** are tracked. Every clone carries two traits, a fitness `s`
and a drug resistance `α`, and **treatment** is a built-in, composable part of the dynamics.
Neutral mutations are recorded either on the full lineage tree of the survivors (tree mode)
or as per-cell counts on the living cells only (cell mode).

The simulation is exact (thinning against a maintained rate bound, exact stops, no silent
clamping of rates). Statistics such as site frequency spectra, clone sizes and Hill numbers
can be computed during a run, on the whole population or on samples, so that only small
results need to be stored.

This package is a **simulator**. It does not infer anything from data, fit parameters, or
model sequencing noise.

## Installation

Requires Julia 1.10 or newer.

```julia
using Pkg
Pkg.add(url = "https://github.com/alexsteininfo/ClonalEvolutionTx.jl")
```

## Quickstart

Grow a tumour from one cell to 2000 cells, then treat it with adaptive therapy and record a
snapshot at every treatment switch.

```julia
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

schedule = AdaptiveOnOff(base; N_low = 1_000, N_up = 2_000, N_prog = 4_000, t_max = 5.0)
acc = MeasurementAccumulator(MeasurementSpec(
    triggers = [AtBlockEnd()],
    stats    = [popsize, n_resistant, MeanBurden()],
))
simulate!(pop, schedule, rng; accumulator = acc)
m = finalize_measurements(acc)
schedule_outcome(pop)

[(s.label, round(s.t, digits = 2), s.stats.popsize, s.stats.n_resistant) for s in m.snapshots]

[(b.label, round(b.t_start, digits = 2), round(b.t_end, digits = 2), b.end_reason) for b in pop.log]
```

## Documentation

The documentation (concepts, the simulation block, mutations, treatment and schedules,
measurements, statistics, sampling, performance, limitations and the API reference) is built
with Documenter.jl from `docs/`:

```bash
julia --project=docs -e 'using Pkg; Pkg.instantiate()'
julia --project=docs docs/make.jl     # output in docs/build/
```

Eight runnable examples are in `examples/` (`julia --project=. examples/01_SingleCellExpansion.jl`).

## Development

```bash
# one-time test environment setup
julia --project=test -e 'using Pkg; Pkg.instantiate()'

# whole test suite
julia --project=. -e 'using Pkg; Pkg.test()'

# selected test files (names without .jl)
julia --project=test test/runtests.jl package traits

# benchmarks (local only; see the Performance page)
julia --project=. benchmark/benchmark.jl 1000 10000 100000
bash benchmark/run_benchmarks.sh
```

## Related packages and provenance

- [NonMarkovEvolution.jl](https://github.com/alexsteininfo/NonMarkovEvolution.jl) (same
  author) covers the complementary regime: non-exponential cell-cycle times, per-cell
  fitness, high mutation rates, no clone labels. The two packages read alike.
- The code descends from [SomaticEvolution.jl](https://github.com/jessierenton/SomaticEvolution.jl)
  (Jessie Renton, MIT) via BirthDeathMutation.jl.

## License

MIT, see `LICENSE`. Copyright (c) 2021 Jessie Renton and contributors; copyright (c) 2026
Alexander Stein.
