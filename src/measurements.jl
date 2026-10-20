# Measurements (SPEC §7): snapshot triggers, statistics, the measurement spec, the
# recorded output and the accumulator that `simulate!` writes into. The accumulator's
# hook methods live with the loop in simulate.jl; the machinery they use is here.

# ── Triggers ──────────────────────────────────────────────────────────────────

"""
    AbstractTrigger

Supertype of the snapshot triggers [`AtTime`](@ref), [`AtSize`](@ref), [`Every`](@ref),
[`When`](@ref) and [`AtBlockEnd`](@ref). Triggers compare with `==` by value.
"""
abstract type AbstractTrigger end

"""
    AtTime(t)

Snapshot of the state exactly at the absolute time `t`: after every event at or before `t`
and before any later one. Fires once per accumulator; `t` equal to a block's `tmax` fires
at the end of that block. An `AtTime` that lies before the population's time when the
accumulator starts a block throws `ArgumentError` there (that state is gone). `t` must be
finite.
"""
struct AtTime <: AbstractTrigger
    t::Float64
    function AtTime(t::Real)
        isfinite(t) || throw(ArgumentError("AtTime: t must be finite, got $t"))
        return new(t)
    end
end

"""
    AtSize(N; direction = :up, every = false)

Snapshot when an accepted event brings the population size to `N` from `direction`: `:up`
(arrival from below), `:down` (arrival from above) or `:either`. A population that already
sits at `N` when the accumulator first sees it has not arrived there. Fires once per
accumulator, or at every arrival with `every = true`. If the population was simulated
without the accumulator between two of its blocks, the change from the size it last saw
to the size at the next block start counts as one step. Requires `N ≥ 0` (`N ≥ 1` for
`:up`).
"""
struct AtSize <: AbstractTrigger
    N::Int
    direction::Symbol
    every::Bool
    function AtSize(N::Integer; direction::Symbol = :up, every::Bool = false)
        direction in (:up, :down, :either) || throw(ArgumentError(
            "AtSize: direction must be :up, :down or :either, got :$direction"))
        N >= (direction === :up ? 1 : 0) || throw(ArgumentError(
            "AtSize: N = $N can never be reached from direction :$direction"))
        return new(N, direction, every)
    end
end

"""
    Every(Δt)

Snapshots on the time grid `t₀ + k·Δt`, `k = 0, 1, …`, where `t₀` is the population's time
when the accumulator first starts a block. Grid times that pass while the population is
simulated without the accumulator are skipped. `Δt` must be finite and positive.
"""
struct Every <: AbstractTrigger
    Δt::Float64
    function Every(Δt::Real)
        (isfinite(Δt) && Δt > 0) || throw(ArgumentError("Every: Δt must be finite and > 0, got $Δt"))
        return new(Δt)
    end
end

"""
    When(f; every = false)

Snapshot when `f(pop)` returns `true`; `f` is checked at block start and after each
accepted event, so keep it cheap. Fires once per accumulator, or at every check that
returns `true` with `every = true`. The start of a block chained onto one the accumulator
saw end is not checked (that state was already checked).
"""
struct When{F} <: AbstractTrigger
    f::F
    every::Bool
end
When(f; every::Bool = false) = When(f, every)

"""
    AtBlockEnd()

Snapshot at the end of every block, whatever ended it (a stop, `tmax` or extinction).
Treatment switches are block ends, so in a schedule this marks every toggle.
"""
struct AtBlockEnd <: AbstractTrigger end

const _KnownTrigger = Union{AtTime,AtSize,Every,When,AtBlockEnd}

# ── Statistics ────────────────────────────────────────────────────────────────

"""
    AbstractStatistic

Supertype of statistic types that snapshots can evaluate. A new statistic defines
[`measure`](@ref)`(stat, pop)` and, optionally, [`statistic_name`](@ref)`(stat)`:

```julia
struct MeanBurden <: AbstractStatistic end
ClonalEvolutionTx.measure(::MeanBurden, pop) = mean(mutations_per_cell(pop))
ClonalEvolutionTx.statistic_name(::MeanBurden) = :mean_burden
```
"""
abstract type AbstractStatistic end

"""
    measure(stat, pop)
    measure(pop::AbstractPopulation, stats::AbstractVector) -> NamedTuple

The value of the statistic `stat` on `pop`. A statistic is an [`AbstractStatistic`](@ref)
(each type defines its `measure` method), a plain function `pop -> value`, or a pair
`name => function`.

The second form evaluates every statistic in `stats` and returns a `NamedTuple` keyed by
[`statistic_name`](@ref), in the given order — between `simulate!` calls, or as snapshots
do. It never changes `pop` (as long as the statistics themselves do not). Throws
`ArgumentError` for an entry that is not a statistic or for two statistics with the same
name.
"""
function measure end

measure(f::Function, pop) = f(pop)
measure(p::Pair{Symbol,<:Function}, pop) = last(p)(pop)

function measure(pop::AbstractPopulation, stats::AbstractVector)
    names = _statistic_names(stats)
    return NamedTuple{names}(Tuple(measure(s, pop) for s in stats))
end

"""
    statistic_name(stat) -> Symbol

The key of `stat`'s value in a snapshot's `stats`: for an [`AbstractStatistic`](@ref) its
type name unless the type defines its own method, for a named function its name
(`popsize` → `:popsize`), for a pair `name => f` the `name`. An anonymous function has no
usable name and throws `ArgumentError`: pass it as `:name => f`.
"""
statistic_name(stat::AbstractStatistic) = Symbol(nameof(typeof(stat)))
function statistic_name(f::Function)
    name = nameof(f)
    startswith(String(name), '#') && throw(ArgumentError(
        "an anonymous function has no usable statistic name; pass it as `:name => f`"))
    return name
end
statistic_name(p::Pair{Symbol,<:Function}) = first(p)
statistic_name(x) = throw(ArgumentError(
    "a statistic must be an AbstractStatistic, a function or a `name => function` pair, got $(repr(x))"))

function _statistic_names(stats)
    names = Tuple(statistic_name(s) for s in stats)
    allunique(names) || throw(ArgumentError("two statistics share a name in $names"))
    return names
end

# ── Specification ─────────────────────────────────────────────────────────────

"""
    MeasurementSpec(; trajectory_dt = nothing, trajectory_clone_sizes = false,
                    triggers = AbstractTrigger[], stats = Any[], samples = Int[],
                    sample_seed = 0)

What a [`MeasurementAccumulator`](@ref) records during `simulate!`:

- `trajectory_dt`: a [`TrajectoryPoint`](@ref) on the time grid `t₀ + k·trajectory_dt`,
  where `t₀` is the time at which the accumulator first starts a block (`nothing`: no
  trajectory). `trajectory_clone_sizes = true` adds the living clones' sizes to each point.
- `triggers`: when to take a [`Snapshot`](@ref): [`AtTime`](@ref), [`AtSize`](@ref),
  [`Every`](@ref), [`When`](@ref), [`AtBlockEnd`](@ref).
- `stats`: the statistics each snapshot evaluates, with distinct names (see
  [`measure`](@ref)).
- `samples`, `sample_seed`: sizes of random samples on which every snapshot also evaluates
  its statistics, and the base seed of their draws. For the `n` in `samples` and a snapshot
  that will be the `i`-th of the accumulator, the sample is
  [`sample_population`](@ref)`(pop, n; seed = derive_seed(sample_seed, i, n))`. A snapshot of a
  population with fewer than `n` cells (for example at extinction) records
  `(n = n, stats = nothing)` for that size: no sample is drawn, and the other sizes keep
  their seeds.

Throws `ArgumentError` for a `trajectory_dt` that is not finite and positive,
`trajectory_clone_sizes` without `trajectory_dt`, an entry of `triggers` or `stats` of the
wrong kind, two statistics with the same name, or a sample size below 1.
"""
struct MeasurementSpec
    trajectory_dt::Union{Nothing,Float64}
    trajectory_clone_sizes::Bool
    triggers::Vector{AbstractTrigger}
    stats::Vector{Any}
    samples::Vector{Int}
    sample_seed::UInt64
    function MeasurementSpec(; trajectory_dt = nothing, trajectory_clone_sizes::Bool = false,
                             triggers = AbstractTrigger[], stats = Any[], samples = Int[],
                             sample_seed::Integer = 0)
        if trajectory_dt !== nothing
            (trajectory_dt isa Real && isfinite(trajectory_dt) && trajectory_dt > 0) ||
                throw(ArgumentError("trajectory_dt must be finite and > 0 (or nothing), got $trajectory_dt"))
        elseif trajectory_clone_sizes
            throw(ArgumentError("trajectory_clone_sizes = true needs a trajectory_dt"))
        end
        for tr in triggers
            tr isa _KnownTrigger || throw(ArgumentError(
                "triggers must be AtTime, AtSize, Every, When or AtBlockEnd, got $(repr(tr))"))
        end
        _statistic_names(stats)
        for n in samples
            n >= 1 || throw(ArgumentError("sample sizes must be ≥ 1, got $n"))
        end
        dt = trajectory_dt === nothing ? nothing : Float64(trajectory_dt)
        return new(dt, trajectory_clone_sizes, AbstractTrigger[tr for tr in triggers],
                   Any[s for s in stats], Int[n for n in samples], sample_seed % UInt64)
    end
end

# ── Output ────────────────────────────────────────────────────────────────────

"""
    TrajectoryPoint

The state at grid time `t`: population size `N`, drug level `m`, the numbers of cells in
fully drug-sensitive clones (`n_sensitive`, `α == 1`) and in all other clones
(`n_resistant`), and `clone_sizes`, the living clones' sizes in
[`living_clones`](@ref) order (`nothing` unless `trajectory_clone_sizes = true`).
"""
struct TrajectoryPoint
    t::Float64
    N::Int
    m::Float64
    n_sensitive::Int
    n_resistant::Int
    clone_sizes::Union{Nothing,Vector{Int}}
end

"""
    Snapshot

The state when a trigger fired: time `t`, population size `N`, drug level `m`, the
`trigger`, `block` (the index the block has in the population's log) and its `label`,
the `stats` as a `NamedTuple` (see [`measure`](@ref)), and `samples`, one
`(n = …, stats = NamedTuple)` per sample size of the spec, in its order (see
[`MeasurementSpec`](@ref)), with `stats = nothing` for a size larger than the population;
empty when the spec has no `samples`.
"""
struct Snapshot
    t::Float64
    N::Int
    m::Float64
    trigger::AbstractTrigger
    block::Int
    label::String
    stats::NamedTuple
    samples::Vector{NamedTuple}
end

"""
    Measurements

What an accumulator recorded: `trajectory::Vector{TrajectoryPoint}` and
`snapshots::Vector{Snapshot}`, both in time order. Returned by
[`finalize_measurements`](@ref).
"""
struct Measurements
    trajectory::Vector{TrajectoryPoint}
    snapshots::Vector{Snapshot}
end

# ── Accumulator ───────────────────────────────────────────────────────────────

"""
    MeasurementAccumulator(spec::MeasurementSpec)

Records what `spec` asks for while `simulate!` runs: pass it as `accumulator` to
`simulate!` (of a block or a schedule) and read the result with
[`finalize_measurements`](@ref).

- Time items (trajectory points, [`AtTime`](@ref), [`Every`](@ref)) record the exact state
  at their time, the clock included (statistics see `pop.t` equal to the item's time):
  before each event, every item due at or before the event's time is recorded with the
  state before the event, and at a `tmax` end every item up to `tmax`.
- [`AtSize`](@ref) and [`When`](@ref) are checked at block start and after each accepted
  event; [`AtBlockEnd`](@ref) fires at every block end.
- At one instant, time items come first (in time order; at equal times the trajectory
  point, then the triggers in list order), then `AtSize`/`When` in list order, then
  `AtBlockEnd`.
- The recorded drug level `m` is the dose's level at the measurement time.
- One accumulator can be carried across blocks and schedules: its time grids continue from
  its first start, each one-shot trigger fires at most once, and the start of a block
  chained onto one it saw end does not fire again for that state. Grid times that pass
  while the population is simulated without it are skipped; an `AtTime` that passes so
  throws `ArgumentError` at the next block start.
- An extinction restart discards what was recorded since the block start.
"""
mutable struct MeasurementAccumulator
    spec::MeasurementSpec
    trajectory::Vector{TrajectoryPoint}
    snapshots::Vector{Snapshot}
    # Indices into spec.triggers by kind, in list order (AtTime sorted by time, then index).
    every::Vector{Int}
    attime::Vector{Int}
    state::Vector{Int}              # AtSize and When
    blockend::Vector{Int}
    # Time items: grid origin t_acc, next grid indices, next AtTime, earliest pending time.
    started::Bool
    t_acc::Float64
    dt::Float64                     # trajectory step; Inf for no trajectory
    traj_k::Int
    every_k::Vector{Int}
    attime_next::Int
    next_due::Float64
    # State triggers: fired one-shot triggers (by trigger index), the last size seen.
    fired::Vector{Bool}
    last_N::Int
    # The current block, and the state seen at the end of the previous one.
    dose::AbstractDose
    t_block::Float64
    label::String
    end_t::Float64
    end_N::Int
    end_log::Int
end

function MeasurementAccumulator(spec::MeasurementSpec)
    trig = spec.triggers
    kind(T) = Int[i for i in eachindex(trig) if trig[i] isa T]
    every = kind(Every)
    attime = sort(kind(AtTime); by = i -> (trig[i].t, i))
    dt = spec.trajectory_dt === nothing ? Inf : spec.trajectory_dt
    return MeasurementAccumulator(spec, TrajectoryPoint[], Snapshot[], every, attime,
                                  kind(Union{AtSize,When}), kind(AtBlockEnd), false, 0.0, dt,
                                  0, zeros(Int, length(every)), 1, Inf,
                                  fill(false, length(trig)), 0, FixedDose(0.0), 0.0, "",
                                  NaN, -1, -1)
end

"""
    finalize_measurements(acc::MeasurementAccumulator) -> Measurements

Copy what `acc` has recorded into a [`Measurements`](@ref). The accumulator stays usable.
"""
finalize_measurements(acc::MeasurementAccumulator) =
    Measurements(copy(acc.trajectory), copy(acc.snapshots))

# ── Recording machinery (used by the hooks in simulate.jl) ────────────────────

@inline _grid_time(t0::Float64, k::Int, Δt::Float64) = t0 + k * Δt   # by multiplication: no drift

# The first grid index ≥ k whose time is at or after t.
function _first_index_from(t0::Float64, k::Int, Δt::Float64, t::Float64)
    j = max(k, floor(Int, (t - t0) / Δt) - 1)
    while _grid_time(t0, j, Δt) < t
        j += 1
    end
    return j
end

# Throws if the next AtTime lies before t (AtTime is sorted, so the next one is the earliest).
function _check_pending_times(acc::MeasurementAccumulator, t::Float64)
    a = acc.attime_next
    a <= length(acc.attime) || return nothing
    tr = acc.spec.triggers[acc.attime[a]]::AtTime
    tr.t < t && throw(ArgumentError(
        "$tr lies before the population's time t = $t: the state at that time is gone"))
    return nothing
end

# First start: the grids start at t. Later starts: skip grid times before t.
function _align_grids!(acc::MeasurementAccumulator, t::Float64)
    if !acc.started
        acc.started = true
        acc.t_acc = t
        return nothing
    end
    isfinite(acc.dt) && (acc.traj_k = _first_index_from(acc.t_acc, acc.traj_k, acc.dt, t))
    for q in eachindex(acc.every)
        Δt = (acc.spec.triggers[acc.every[q]]::Every).Δt
        acc.every_k[q] = _first_index_from(acc.t_acc, acc.every_k[q], Δt, t)
    end
    return nothing
end

# Whether this block start continues directly from the block end the accumulator saw last.
_chained(acc::MeasurementAccumulator, pop) =
    pop.t == acc.end_t && popsize(pop) == acc.end_N && length(pop.log) == acc.end_log

# The earliest pending time item as (time, kind, j): kind 1 is the next trajectory point,
# 2 the next grid time of the j-th Every trigger, 3 the j-th AtTime of the sorted list;
# (Inf, 0, 0) if none. Ties: the trajectory point first, then the trigger listed first.
function _earliest(acc::MeasurementAccumulator)
    τ, kind, j, rank = Inf, 0, 0, typemax(Int)
    if isfinite(acc.dt)
        τ, kind, rank = _grid_time(acc.t_acc, acc.traj_k, acc.dt), 1, 0
    end
    trig = acc.spec.triggers
    for q in eachindex(acc.every)
        i = acc.every[q]
        τq = _grid_time(acc.t_acc, acc.every_k[q], (trig[i]::Every).Δt)
        if τq < τ || (τq == τ && i < rank)
            τ, kind, j, rank = τq, 2, q, i
        end
    end
    a = acc.attime_next
    if a <= length(acc.attime)
        i = acc.attime[a]
        τq = (trig[i]::AtTime).t
        if τq < τ || (τq == τ && i < rank)
            τ, kind, j, rank = τq, 3, a, i
        end
    end
    return τ, kind, j
end

# Record every time item due at or before t with the current state, in order, and cache
# the earliest time still pending. The clock is part of the state, so it is moved to each
# item's time τ first: no event lies between `pop.t` and τ (pending items are never
# earlier than `pop.t`, and τ ≤ t, the time of the next tick, which the loop sets next).
@noinline function _record_due!(acc::MeasurementAccumulator, pop, t::Float64)
    while true
        τ, kind, j = _earliest(acc)
        if kind == 0 || τ > t
            acc.next_due = τ
            return nothing
        end
        pop.t = τ
        if kind == 1
            push!(acc.trajectory, _trajectory_point(acc, pop, τ))
            acc.traj_k += 1
        elseif kind == 2
            _snapshot!(acc, pop, τ, acc.spec.triggers[acc.every[j]])
            acc.every_k[j] += 1
        else
            _snapshot!(acc, pop, τ, acc.spec.triggers[acc.attime[j]])
            acc.attime_next += 1
        end
    end
end

# AtSize and When triggers, in list order, for the size change prev → N into the current
# state; When is skipped if `check_when` is false.
function _check_state!(acc::MeasurementAccumulator, pop, prev::Int, N::Int, check_when::Bool)
    trig = acc.spec.triggers
    for i in acc.state
        tr = trig[i]
        if tr isa AtSize
            (tr.every || !acc.fired[i]) && _arrives(tr, prev, N) && _fire!(acc, pop, i, tr)
        elseif check_when
            w = tr::When
            (w.every || !acc.fired[i]) && _holds(w, pop) && _fire!(acc, pop, i, w)
        end
    end
    return nothing
end

function _arrives(tr::AtSize, prev::Int, N::Int)
    up = prev < tr.N <= N
    down = prev > tr.N >= N
    return tr.direction === :up ? up : tr.direction === :down ? down : (up || down)
end

function _holds(w::When, pop)
    r = w.f(pop)
    r isa Bool || throw(ArgumentError("When: f(pop) must return a Bool, got $(repr(r))"))
    return r
end

function _fire!(acc::MeasurementAccumulator, pop, i::Int, tr::Union{AtSize,When})
    tr.every || (acc.fired[i] = true)
    _snapshot!(acc, pop, pop.t, tr)
    return true
end

# The drug level at time τ of the current block: a time dose's own level there, otherwise
# the level of the current state (constant between events).
_level(acc::MeasurementAccumulator, pop, τ::Float64) = _level_at(acc.dose, acc.t_block, pop, τ)
_level_at(d::TimeDose, t0::Float64, pop, τ::Float64) = _time_level(d, τ - t0)
_level_at(::AbstractDose, t0::Float64, pop, τ::Float64) = pop.m

function _trajectory_point(acc::MeasurementAccumulator, pop, τ::Float64)
    ct = pop.clones
    ns = nr = 0
    for k in ct.alive
        ct.α[k] == 1.0 ? (ns += ct.size[k]) : (nr += ct.size[k])
    end
    sizes = acc.spec.trajectory_clone_sizes ? Int[ct.size[k] for k in ct.alive] : nothing
    return TrajectoryPoint(τ, popsize(pop), _level(acc, pop, τ), ns, nr, sizes)
end

# A sample size larger than the population (e.g. at extinction) is recorded as
# `stats = nothing`: no sample is drawn, and the seeds of the other sizes do not change.
function _snapshot!(acc::MeasurementAccumulator, pop, τ::Float64, trigger::AbstractTrigger)
    index = length(acc.snapshots) + 1
    samples = NamedTuple[]
    for n in acc.spec.samples
        if n > popsize(pop)
            push!(samples, (n = n, stats = nothing))
            continue
        end
        sample = sample_population(pop, n; seed = derive_seed(acc.spec.sample_seed, index, n))
        push!(samples, (n = n, stats = measure(sample, acc.spec.stats)))
    end
    push!(acc.snapshots, Snapshot(τ, popsize(pop), _level(acc, pop, τ), trigger,
                                  length(pop.log) + 1, acc.label,
                                  measure(pop, acc.spec.stats), samples))
    return nothing
end
