# Measurements and output

Measurements are taken **during** a run, while the population is in memory, because the tree
later loses extinct lineages. A [`MeasurementSpec`](@ref) says what to record, a
[`MeasurementAccumulator`](@ref) records it while [`simulate!`](@ref) runs (of a block or a
schedule), and [`finalize_measurements`](@ref) returns a [`Measurements`](@ref) with two
plain vectors, `trajectory` and `snapshots`. Plain structs and `NamedTuple`s: a
`Vector{NamedTuple}` is already a Tables.jl table, and no file is written.

```@example meas
using ClonalEvolutionTx, Random

base = BirthDeathBlock(birth = s -> 1.0, death = s -> 0.5, μ = 1.0,
                       response = LinearResponse(4.0),
                       mutations = [Resistance(ν = 0.01, α = 0.8)])

spec = MeasurementSpec(
    trajectory_dt = 1.0,
    triggers      = [AtTime(4.0), AtSize(400), AtBlockEnd()],
    stats         = [popsize, MeanBurden(), n_resistant],
)
acc = MeasurementAccumulator(spec)
pop = initialize_population(1)
simulate!(pop, remake(base; stop_above = 800, restart_on_extinction = true), Xoshiro(2);
          accumulator = acc)
meas = finalize_measurements(acc)

[(sn.trigger, round(sn.t, digits = 2), sn.N) for sn in meas.snapshots]
```

## Triggers (when)

| trigger | semantics |
|:--|:--|
| [`AtTime`](@ref)`(t)` | the state exactly at the absolute time `t`: after every event at or before `t`, before any later one. Fires once; `t == tmax` fires at the end of that block. An `AtTime` before the population's current time throws, since that state is gone |
| [`AtSize`](@ref)`(N; direction = :up, every = false)` | first (or every) arrival at size `N` from below (`:up`), above (`:down`) or either side |
| [`Every`](@ref)`(Δt)` | full snapshots on the grid `t₀ + kΔt` |
| [`When`](@ref)`(f; every = false)` | when `f(pop)` returns `true`, checked at block start and after each event; keep it cheap |
| [`AtBlockEnd`](@ref)`()` | at the end of every block, whatever ended it. Treatment switches are block ends |

`trajectory_dt = Δt` in the spec adds cheap [`TrajectoryPoint`](@ref)s on the time grid
`t₀ + kΔt` (`t₀` is the time at which the accumulator first starts a block). With
`trajectory_clone_sizes = true` each point also carries the sizes of all living clones.

Measurements record the state **exactly** at the measurement time: the clock `pop.t`
that statistics see equals the item's time, and the recorded drug level `m` is the dose's
level at that time (for a [`TimeDose`](@ref), the level `m(t)` itself). At one instant, time
items come first, then `AtSize`/`When`, then `AtBlockEnd`.

One accumulator can be carried across blocks and schedules: its grids continue, each
one-shot trigger fires at most once, and trigger firing does not repeat across chained
blocks. Grid times that pass while the population is simulated *without* the accumulator
are skipped.

## What is recorded

A **trajectory point** has `t`, `N`, the drug level `m`, `n_sensitive` (cells in clones with
`α == 1`), `n_resistant` (all other cells) and optionally `clone_sizes`.

```@example meas
p = meas.trajectory[end]
(p.t, p.N, p.n_sensitive, p.n_resistant, p.m)
```

A **snapshot** has `t`, `N`, `m`, the `trigger`, the index `block` of the block in the
population's log, its `label`, the requested statistics `stats` as a `NamedTuple`, and
`samples` (see [Sampling](@ref)), with `stats = nothing` for a sample size larger than the
population.

```@example meas
sn = meas.snapshots[1]
(sn.label, sn.stats.popsize, round(sn.stats.mean_burden, digits = 1), sn.samples)
```

## Statistics (what)

An entry of `stats` is

- a statistic type: [`SFS`](@ref), [`MutationsPerCell`](@ref), [`MeanBurden`](@ref),
  [`CloneSizes`](@ref), [`MutationFrequencies`](@ref), [`HillNumbers`](@ref),
  [`SFSHillNumbers`](@ref), [`Divisions`](@ref), [`SensitiveResistant`](@ref);
- a **named function** `pop -> value`, such as `popsize` or `n_resistant`, keyed by its name;
- or a pair `:name => f` for anything else.

**Anonymous functions are rejected** (they have no usable key): write `:name => (pop -> ...)`
instead. Two statistics with the same name are an error. A new statistic type subtypes
[`AbstractStatistic`](@ref) and defines [`measure`](@ref)`(stat, pop)` and, optionally,
[`statistic_name`](@ref)`(stat)`.

```@example meas
n_big_clones(pop) = count(>(50), clone_sizes(pop))
spec2 = MeasurementSpec(triggers = [AtBlockEnd()],
                        stats = [SFS(), :resistant_fraction => (pop -> n_resistant(pop) / popsize(pop)),
                                 n_big_clones])
keys(measure(pop, spec2.stats))
```

## Direct measurement

[`measure`](@ref)`(pop, stats)` evaluates a vector of statistics between `simulate!` calls and
returns a `NamedTuple`. It never changes the population.

```@example meas
measure(pop, [popsize, MeanBurden(), HillNumbers(q = [0, 1, 2])])
```

## The treatment log

Every population carries a log with one [`BlockSummary`](@ref) per simulated block: label,
`t_start`, `t_end`, `end_reason`, ticks, accepted events, restarts, the dose description and
the final drug level. Figures that shade treatment periods read it.

```@example meas
[(b.label, b.dose, b.end_reason) for b in pop.log]
```

## Measuring across a schedule

```@example meas
sched = AdaptiveOnOff(base; N_low = 400, N_up = 800, N_prog = 1_600, t_max = 3.0)
acc2 = MeasurementAccumulator(MeasurementSpec(trajectory_dt = 0.5, triggers = [AtBlockEnd()],
                                              stats = [popsize, n_resistant]))
simulate!(pop, sched, Xoshiro(2); accumulator = acc2)
m2 = finalize_measurements(acc2)
(schedule_outcome(pop), length(m2.trajectory), [(s.label, s.stats.popsize) for s in m2.snapshots])
```

## Samples larger than the population, and errors

A snapshot of a population smaller than a requested sample size `n` (for example the
[`AtBlockEnd`](@ref) snapshot of a block that ended in extinction) records
`(n = n, stats = nothing)` for that size: no sample is drawn, and the other sizes are drawn
with their usual seeds. Invalid triggers, statistics or spec options throw at construction.
