# The simulation block

A [`BirthDeathBlock`](@ref) is one stretch of birth–death dynamics with fixed rules. It is
immutable; [`remake`](@ref) returns a copy with some fields changed, which is how one base
block is reused for growth, treatment and everything else. [`simulate!`](@ref) runs a block
on a population in place and returns the same population.

```@example blocks
using ClonalEvolutionTx, Random

base = BirthDeathBlock(
    birth = s -> 1.0 + s,        # β(s)
    death = s -> 0.5,            # δ(s)
    μ     = 1.0,                 # neutral mutations per daughter
    label = "growth",
)
pop = initialize_population(1; mode = :cells)
simulate!(pop, remake(base; stop_above = 500, restart_on_extinction = true), Xoshiro(1))
(popsize(pop), pop.log[end].end_reason)
```

## Rates and density dependence

The rates are those of [Concepts](@ref): `β(s)·B(N)` for birth and `δ(s) + C(N) + α·D(m)`
for death. The keywords are

| keyword | meaning | default |
|:--|:--|:--|
| `birth`, `death` | functions `s -> β(s)`, `s -> δ(s)`, evaluated once per clone | required |
| `birth_density` | `N -> B(N)`, multiplies every birth rate | `Returns(1.0)` |
| `death_density` | `N -> C(N)`, added to every death rate | `Returns(0.0)` |
| `response`, `dose` | drug effect, see [Treatment and schedules](@ref) | no drug |

Logistic growth is `death_density = N -> N / (2K)` with `death = s -> 0.5`, settling at `N ≈ K`;
a birth-limited version is `birth_density = N -> max(0, 1 - N/K)`. Cost of resistance is
`birth = s -> b0 * (1 + s)` together with a negative `Δs` in [`Resistance`](@ref).

A negative or non-finite rate is an `ArgumentError`, never silently clamped: at clone
creation for `β` and `δ`, at evaluation for `B` and `C`.

```@example blocks
K = 400
logistic = BirthDeathBlock(birth = s -> 1.0, death = s -> 0.5,
                           death_density = N -> 0.5 * N / K, tmax = 25.0)
pop = initialize_population(50; mode = :cells)
simulate!(pop, logistic, Xoshiro(2))
(popsize(pop), pop.t)       # fluctuates around K = 400 (b - d = 0.5 N/K gives N = K)
```

## Stop conditions

At the start of a block and after every event, the block ends when one of these holds, checked
in this order:

1. the population is extinct (reason `:extinct`);
2. `N ≥ stop_above` (`:above`);
3. `N ≤ stop_below` (`:below`);
4. the user function `stop(pop)` returns `true` (`:stopped`).

`tmax` is the **absolute** simulation time at which the block ends (`:tmax`), not a
duration; `Inf` means none. A block whose `tmax` lies before the population's current time
throws. To run for a duration `T` from the current time, pass `tmax = pop.t + T`. (The
schedules of [Treatment and schedules](@ref) have a `t_max` that *is* relative to the
schedule start.)

Stops are exact. The stop condition is checked before every event, and if the next event
would fall after `tmax` the clock is set to `tmax` exactly and the event is discarded,
which is exact because waiting times are memoryless. Nothing overshoots by one event.

```@example blocks
pop = initialize_population(200; mode = :cells)
simulate!(pop, remake(base; tmax = pop.t + 2.0), Xoshiro(3))
(pop.t, pop.log[end].end_reason)
```

If no event can ever happen (all rates zero), the clock jumps to `tmax`; with `tmax = Inf`
and no stop that can still fire, `simulate!` throws an `ErrorException` instead of hanging.

## Restarts on extinction

With `restart_on_extinction = true`, a population that dies out is restored to its state at
the start of the block and simulated again with the continuing random stream, until a stop
other than extinction ends the block. This is how growth from a single founder to a
detection size is usually run. The restarts are counted in the log. A population that is
certain to die out restarts forever, so use it only with a viable growth rate. A restart
discards what an attached [`MeasurementAccumulator`](@ref) recorded since the block start.

## The log

Every `simulate!` call appends a [`BlockSummary`](@ref) to `pop.log`: label, start and end
time, end reason, number of ticks and accepted events (their ratio is the acceptance
fraction of the thinning), restarts, a description of the dose and the final drug level.

```@example blocks
b = pop.log[end]
(b.label, b.end_reason, b.ticks >= b.accepted)
```

## Chaining and copying

Because `simulate!` continues the population's clock, cells and clone table, blocks can be
chained: grow to a detection size, then continue with a different block. A grown population
can be reused for several scenarios by copying it first. `copy(pop)`
is deep, so simulating the copy leaves the original untouched.

```@example blocks
grown = initialize_population(1; mode = :cells)
simulate!(grown, remake(base; stop_above = 300, restart_on_extinction = true), Xoshiro(4))

a = copy(grown); b = copy(grown)
simulate!(a, remake(base; stop_above = 600), Xoshiro(5))
simulate!(b, remake(base; death = s -> 2.0, stop_below = 150), Xoshiro(6))   # declining
(popsize(grown), popsize(a), popsize(b))
```

## Limiting the clone count

`max_clones` stops the creation of new clones once that many exist (a daughter that would have
mutated simply stays in its clone). The default is unlimited.

## Errors

Invalid inputs throw `ArgumentError` (a negative `μ`, a `NaN` `tmax`, a block whose `tmax`
lies in the past, a bad rate); an impossible state throws `ErrorException`. The package
never replaces a bad value by a clamped one.
