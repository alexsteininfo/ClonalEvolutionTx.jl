# Treatment and schedules

Treatment is cytotoxic: the drug raises the death rate of a cell of clone `k` by
`α_k · D(m)`, where `m ∈ [0, 1]` is the treatment level and `D` the drug response. A fully
resistant clone (`α = 0`) is unaffected; a fully sensitive one (`α = 1`) feels all of it.

## Responses

A response `D(m)` must be monotone non-decreasing with `D(0) = 0`.

- [`LinearResponse`](@ref)`(d_S)`: `D(m) = d_S·m`.
- [`HillResponse`](@ref)`(d_max, n, EC50)`: `D(m) = d_max·mⁿ / (EC50ⁿ + mⁿ)`.

[`drug_death`](@ref)`(response, m)` evaluates `D`, and [`inverse_drug_death`](@ref)`(response, x)`
its inverse (closed form for both built-ins). A custom response subtype of
[`AbstractResponse`](@ref) defines `drug_death` and, optionally, `inverse_drug_death`;
without the latter a bisection fallback (to `1e-12`) is used.

```@example treatment
using ClonalEvolutionTx, Random

hill = HillResponse(4.0, 2.0, 0.5)
x = drug_death(hill, 0.7)
(x, inverse_drug_death(hill, x))
```

## Doses

A block's `dose` says how `m` is set. All four forms work in any block, inside the built-in
schedules and in user schedules, and the level is recorded in trajectories and snapshots.

| dose | level |
|:--|:--|
| [`FixedDose`](@ref)`(m)` | constant `m` for the whole block |
| [`StabilisingDose`](@ref)`(; m_max = 1.0)` | after every event, the level at which the population's expected net growth is zero, capped at `m_max` |
| [`SizeDose`](@ref)`(f)` | `m = f(N)`, re-evaluated after every event |
| [`TimeDose`](@ref)`(f; m_max = 1.0)` | `m = f(t_rel)`, where `t_rel` is the time since the block start, which may change *between* events |

The **stabilising level** solves

```
D(m) = [ B(N)·Σₖ Nₖβₖ − Σₖ Nₖδₖ − C(N)·N ] / Σₖ Nₖαₖ
```

with `m = 0` when the numerator is not positive and `m = m_max` when the target exceeds
`D(m_max)` or no drug-sensitive cells are left (`Σₖ Nₖαₖ = 0`: never `NaN`). The three sums
are running sums updated in O(1) at every birth and death. Because mutation never changes
`N`, this is exact for the whole population.

A [`TimeDose`](@ref) is simulated exactly by thinning against the bound `(max α)·D(m_max)`
and evaluating `m(t)` at each tick time; the tighter the declared `m_max`, the fewer
rejected ticks. The bound's `D(m_max)` is raised by a few ulps, so that floating-point
rounding of a monotone response just below `m_max` (a Hill response can round a level
`m < m_max` one ulp above `D(m_max)`) cannot exceed it; a response that is not monotone
still exceeds it and is an error. A level outside `[0, m_max]` is an error. A block's log
entry records the level at its end time. Piecewise-constant doses are better
written as a sequence of blocks (exact stops at the change points, no wasted ticks), which is
what [`Metronomic`](@ref) does.

```@example treatment
base = BirthDeathBlock(
    birth = s -> 1.0, death = s -> 0.5,
    response = LinearResponse(4.0),
    mutations = [Resistance(ν = 0.005, α = 0.1)],
    label = "block",
)
# 700 drug-sensitive cells and 300 cells of a strongly resistant clone
grown = initialize_population([(n = 700, s = 0.0, α = 1.0), (n = 300, s = 0.0, α = 0.1)];
                              mode = :cells)

# a fixed dose for a duration of 3 time units (tmax is absolute)
treated = copy(grown)
simulate!(treated, remake(base; dose = FixedDose(1.0), tmax = grown.t + 3.0, label = "drug"),
          Xoshiro(2))
(popsize(grown), popsize(treated), treated.m)
```

## Schedules

A **schedule** produces the next block from the current population and its history until it
is done. It is built from a **base block** (rates, mutations, `μ`) and varies only the dose
and the stop conditions. Run it with `simulate!(pop, schedule, rng; accumulator = ...)`; one
[`BlockSummary`](@ref) per block lands in `pop.log`. The built-in schedules end on
progression (`N ≥ N_prog`), extinction, a `stop` of the base block, or when the clock reaches
`t0 + t_max`, where `t0` is the time at which the schedule started.

| schedule | behaviour |
|:--|:--|
| [`MTD`](@ref)`(base; N_prog, t_max, m)` | one block with `FixedDose(m)` (default 1) until progression, extinction or `t_max` |
| [`Metronomic`](@ref)`(base; intervals, N_prog, t_max, m)` | `m` on each `(t_on, t_off)` of `intervals` (relative to the schedule start), none otherwise |
| [`AdaptiveOnOff`](@ref)`(base; N_low, N_up, N_prog, t_max, m)` | `FixedDose(m)` (default 1) until `N ≤ N_low`, no drug until `N ≥ N_up`, repeat; needs `1 ≤ N_low < N_up < N_prog` and a finite `t_max` |
| [`Stabilising`](@ref)`(base; N_stable, N_prog, t_max, m_max)` | `FixedDose(m_max)` until `N ≤ N_stable`, then a `StabilisingDose(m_max = m_max)` until progression or `t_max`; finite `t_max` |
| [`CustomDose`](@ref)`(base, dose; N_prog, t_max)` | one block with a [`SizeDose`](@ref) `m(N)` or [`TimeDose`](@ref) `m(t)` |

Note that `t_max` of a schedule is measured from the schedule start (the blocks get
`tmax = t0 + t_max`), for every schedule including `Stabilising`; the `tmax` of a block is
absolute. The `Stabilising` convention differs from BirthDeathMutation_AT, which measured the
maintenance `t_max` from the start of the maintenance phase. `MTD`, `Metronomic` and
`CustomDose` default to `t_max = Inf`: a run that can neither reach `N_prog` nor die out (for
example one held below `N_prog` by density dependence) then never ends, so give a finite
`t_max` unless progression or extinction is certain.

The pre-treatment growth to a detection size is an ordinary block (here we start directly from
a mixed population), so one grown population serves several schedules (copy it first):

```@example treatment
mtd = simulate!(copy(grown), MTD(base; N_prog = 2_000, t_max = 40.0), Xoshiro(3))
(schedule_outcome(mtd), popsize(mtd), round(mtd.t - grown.t, digits = 2))
```

`schedule_outcome(pop)` reports how the last block ended: `:progressed`, `:extinct`, `:tmax`,
`:below` or `:stopped`. A population that dies out in the middle of a schedule (for example
under MTD with no resistance) ends the schedule cleanly with `:extinct`.

```@example treatment
at = AdaptiveOnOff(base; N_low = 500, N_up = 1_000, N_prog = 2_000, t_max = 12.0)
cyc = simulate!(copy(grown), at, Xoshiro(4))
(schedule_outcome(cyc), [b.label for b in cyc.log])
```

```@example treatment
met = Metronomic(base; intervals = [(0.0, 1.0), (3.0, 4.0)], N_prog = 5_000, t_max = 6.0)
mpop = simulate!(copy(grown), met, Xoshiro(5))
[(b.label, round(b.t_start - grown.t, digits = 2), round(b.t_end - grown.t, digits = 2)) for b in mpop.log]
```

Treatment switches are block ends, so an [`AtBlockEnd`](@ref) trigger (see
[Measurements and output](@ref)) records every toggle.

### Stabilising therapy

```@example treatment
stab = Stabilising(base; N_stable = 700, N_prog = 2_000, t_max = 10.0)
spop = simulate!(copy(grown), stab, Xoshiro(6))
(schedule_outcome(spop), [b.label for b in spop.log], round(spop.m, digits = 3))
```

### Custom doses

A dose that depends on the population size or on time uses [`CustomDose`](@ref):

```@example treatment
dose = SizeDose(N -> clamp((N - 500) / 500, 0.0, 1.0))      # more drug at larger N
cpop = simulate!(copy(grown), CustomDose(base, dose; N_prog = 2_000, t_max = 8.0), Xoshiro(7))
(schedule_outcome(cpop), popsize(cpop))
```

```@example treatment
ramp = TimeDose(t -> min(t / 4, 1.0))                       # ramp up over 4 time units
tpop = simulate!(copy(grown), CustomDose(base, ramp; N_prog = 2_000, t_max = 5.0), Xoshiro(8))
(schedule_outcome(tpop), round(tpop.m, digits = 3))
```

### Writing your own schedule

Implement [`next_block`](@ref)`(s, pop, i, t0)` for a subtype of [`AbstractSchedule`](@ref), or
wrap a function in [`FunctionSchedule`](@ref): `f(pop, i, t0)` returns the next block, or
`nothing` when done. `i` is the number of blocks the schedule has already run. A function
schedule's blocks are used as returned (no built-in end rules), and a block's `tmax` is
absolute: write `remake(base; tmax = t0 + 5.0)`.

```@example treatment
# three alternating blocks: 2 time units on, 1 off, 2 on
function plan(pop, i, t0)
    i == 3 && return nothing
    on = iseven(i)
    return remake(base; dose = FixedDose(on ? 1.0 : 0.0), tmax = pop.t + (on ? 2.0 : 1.0),
                  label = on ? "on" : "off")
end
fpop = simulate!(copy(grown), FunctionSchedule(plan), Xoshiro(9))
[b.label for b in fpop.log]
```
