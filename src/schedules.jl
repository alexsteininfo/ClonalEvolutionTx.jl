# Treatment schedules (SPEC §4.6): sequences of blocks built from one base block.

"""
    AbstractSchedule

Supertype of treatment schedules. A schedule produces its next block from the population
and its history; implement `next_block(s, pop, i, t0)` for a new design (or wrap a function
in [`FunctionSchedule`](@ref)). Run a schedule with `simulate!(pop, schedule, rng)`.
"""
abstract type AbstractSchedule end

"""
    next_block(s::AbstractSchedule, pop, i::Int, t0::Float64) -> Union{BirthDeathBlock, Nothing}

The next block of schedule `s`, or `nothing` when the schedule is done. `i` is the number of
blocks this schedule has already run (the last one is `pop.log[end]` when `i ≥ 1`) and `t0`
the time at which the schedule started.
"""
function next_block end

"""
    simulate!(pop, schedule::AbstractSchedule, rng; accumulator = nothing) -> pop

Run a treatment schedule: with `t0 = pop.t`, repeatedly ask the schedule for its next block
and simulate it, until the schedule returns `nothing`. One [`BlockSummary`](@ref) per block is
appended to `pop.log`. The built-in schedules end on progression (`N ≥ N_prog`), extinction,
a user `stop` in the base block, or when the clock reaches `t0 + t_max`.
"""
function simulate!(pop::AbstractPopulation, s::AbstractSchedule, rng::AbstractRNG = Xoshiro();
                   accumulator = nothing)
    t0 = pop.t
    i = 0
    while true
        block = next_block(s, pop, i, t0)
        block === nothing && break
        simulate!(pop, block, rng; accumulator)
        i += 1
    end
    return pop
end

"""
    schedule_outcome(pop) -> Symbol

How the last simulated block of `pop` ended: `:progressed` if it ended at or above its upper size
threshold (`N ≥` the threshold), otherwise its end reason (`:extinct`, `:tmax`, `:below`, `:stopped`).
"""
function schedule_outcome(pop::AbstractPopulation)
    isempty(pop.log) && throw(ArgumentError("the population's log is empty: nothing has been simulated"))
    r = pop.log[end].end_reason
    return r === :above ? :progressed : r
end

# Common end-of-schedule rule of the built-ins (Ruling 5).
function _schedule_done(pop, i::Int, t0::Float64, N_prog::Int, t_max::Float64)
    i >= 1 || return false
    r = pop.log[end].end_reason
    (r === :extinct || r === :stopped) && return true
    r === :above && popsize(pop) >= N_prog && return true
    return pop.t >= t0 + t_max
end

function _check_common(N_prog, t_max)
    N_prog >= 1 || throw(ArgumentError("N_prog must be ≥ 1, got $N_prog"))
    (isnan(t_max) || t_max < 0) && throw(ArgumentError("t_max must be ≥ 0 and not NaN, got $t_max"))
    return Int(N_prog), Float64(t_max)
end

# Block of a built-in schedule: base rules, new dose/label/stops, never restarting.
_sched_block(base::BirthDeathBlock, label, dose, tmax; stop_above, stop_below = -1) =
    remake(base; dose, label, tmax, stop_above, stop_below, restart_on_extinction = false)

"""
    MTD(base; N_prog, t_max = Inf, m = 1.0)

Maximum tolerated dose: one block labelled `"MTD"` with `FixedDose(m)` until progression
(`N ≥ N_prog`), extinction or `t_max`. Counterpart of `manuscript_v1`'s MTD.
`base` is a [`BirthDeathBlock`](@ref) whose rates, mutations and `stop` are kept.
`t_max` is measured from the schedule start (the block gets `tmax = t0 + t_max`).
With the default `t_max = Inf` a run that can neither reach `N_prog` nor die out (for
example one held below `N_prog` by its density dependence) never ends: give a finite
`t_max` unless progression or extinction is certain.
"""
struct MTD{B<:BirthDeathBlock} <: AbstractSchedule
    base::B
    N_prog::Int
    t_max::Float64
    m::Float64
    function MTD(base::B; N_prog, t_max = Inf, m = 1.0) where {B<:BirthDeathBlock}
        Np, tm = _check_common(N_prog, t_max)
        return new{B}(base, Np, tm, FixedDose(m).m)
    end
end

function next_block(s::MTD, pop, i::Int, t0::Float64)
    _schedule_done(pop, i, t0, s.N_prog, s.t_max) && return nothing
    i == 0 || return nothing
    return _sched_block(s.base, "MTD", FixedDose(s.m), t0 + s.t_max; stop_above = s.N_prog)
end

"""
    Metronomic(base; intervals, N_prog, t_max = Inf, m = 1.0)

Metronomic therapy: `FixedDose(m)` during each `(t_on, t_off)` of `intervals` and no drug
otherwise, until progression, extinction or `t_max`. `intervals` is a sorted, non-overlapping
vector of `Tuple{Float64,Float64}` with `0 ≤ t_on < t_off`, **relative to the schedule start**.
Blocks alternate `"metronomic off"` and `"metronomic on"` (zero-length gaps are skipped),
followed by a last off block. Counterpart of `manuscript_v1`'s metronomic schedule.
`t_max` is measured from the schedule start. Throws `ArgumentError` for invalid intervals.
With the default `t_max = Inf` a run that can neither reach `N_prog` nor die out (for
example one held below `N_prog` by its density dependence) never ends: give a finite
`t_max` unless progression or extinction is certain.
"""
struct Metronomic{B<:BirthDeathBlock} <: AbstractSchedule
    base::B
    intervals::Vector{Tuple{Float64,Float64}}
    N_prog::Int
    t_max::Float64
    m::Float64
    function Metronomic(base::B; intervals, N_prog, t_max = Inf, m = 1.0) where {B<:BirthDeathBlock}
        Np, tm = _check_common(N_prog, t_max)
        iv = Tuple{Float64,Float64}[(Float64(a), Float64(b)) for (a, b) in intervals]
        prev = 0.0
        for (a, b) in iv
            (0 <= a < b) || throw(ArgumentError("each interval needs 0 ≤ t_on < t_off, got ($a, $b)"))
            a >= prev || throw(ArgumentError("intervals must be sorted and non-overlapping, got ($a, $b) after an interval ending at $prev"))
            prev = b
        end
        return new{B}(base, iv, Np, tm, FixedDose(m).m)
    end
end

# The (on?, end offset) segments of the schedule, zero-length gaps skipped.
function _segments(s::Metronomic)
    segs = Tuple{Bool,Float64}[]
    cursor = 0.0
    for (a, b) in s.intervals
        a > cursor && push!(segs, (false, a))
        push!(segs, (true, b))
        cursor = b
    end
    push!(segs, (false, Inf))
    return segs
end

function next_block(s::Metronomic, pop, i::Int, t0::Float64)
    _schedule_done(pop, i, t0, s.N_prog, s.t_max) && return nothing
    segs = _segments(s)
    i < length(segs) || return nothing
    on, off_at = segs[i + 1]
    tmax = t0 + min(off_at, s.t_max)
    return _sched_block(s.base, on ? "metronomic on" : "metronomic off",
                        FixedDose(on ? s.m : 0.0), tmax; stop_above = s.N_prog)
end

"""
    AdaptiveOnOff(base; N_low, N_up, N_prog, t_max, m = 1.0)

Adaptive therapy: `"AT on"` (`FixedDose(m)`) until `N` falls to `N_low` (`N ≤ N_low`), then
`"AT off"` (no drug) until `N` reaches `N_up` (`N ≥ N_up`), repeated until progression (`N_prog`), extinction or
`t_max`. Requires `1 ≤ N_low < N_up < N_prog` and a finite `t_max` (`ArgumentError`).
Counterpart of `manuscript_v1`'s adaptive therapy (`BirthDeathMutation_AT`'s
`adaptiveprocess!`). `t_max` is measured from the schedule start.
"""
struct AdaptiveOnOff{B<:BirthDeathBlock} <: AbstractSchedule
    base::B
    N_low::Int
    N_up::Int
    N_prog::Int
    t_max::Float64
    m::Float64
    function AdaptiveOnOff(base::B; N_low, N_up, N_prog, t_max, m = 1.0) where {B<:BirthDeathBlock}
        Np, tm = _check_common(N_prog, t_max)
        isfinite(tm) || throw(ArgumentError("AdaptiveOnOff needs a finite t_max, got $t_max"))
        (1 <= N_low < N_up < Np) || throw(ArgumentError(
            "AdaptiveOnOff needs 1 ≤ N_low < N_up < N_prog, got $N_low, $N_up, $Np"))
        return new{B}(base, Int(N_low), Int(N_up), Np, tm, FixedDose(m).m)
    end
end

function next_block(s::AdaptiveOnOff, pop, i::Int, t0::Float64)
    _schedule_done(pop, i, t0, s.N_prog, s.t_max) && return nothing
    tmax = t0 + s.t_max
    if iseven(i)
        return _sched_block(s.base, "AT on", FixedDose(s.m), tmax;
                            stop_above = s.N_prog, stop_below = s.N_low)
    else
        return _sched_block(s.base, "AT off", FixedDose(0.0), tmax; stop_above = s.N_up)
    end
end

"""
    Stabilising(base; N_stable, N_prog, t_max, m_max = 1.0)

Stabilising therapy: `"stabilising: reduce"` (`FixedDose(m_max)`) until `N` falls to
`N_stable` (`N ≤ N_stable`), then `"stabilising: hold"` with [`StabilisingDose`](@ref)`(m_max = m_max)` until
progression, extinction or `t_max`. Requires `1 ≤ N_stable < N_prog` and a finite `t_max`
(`ArgumentError`). Counterpart of `manuscript_v1`'s stabilising therapy
(`BirthDeathMutation_AT`'s `maintenanceprocess!`). Unlike `BirthDeathMutation_AT`, which
measured the maintenance `t_max` from the start of the maintenance phase, `t_max` here is
measured from the schedule start.
"""
struct Stabilising{B<:BirthDeathBlock} <: AbstractSchedule
    base::B
    N_stable::Int
    N_prog::Int
    t_max::Float64
    m_max::Float64
    function Stabilising(base::B; N_stable, N_prog, t_max, m_max = 1.0) where {B<:BirthDeathBlock}
        Np, tm = _check_common(N_prog, t_max)
        isfinite(tm) || throw(ArgumentError("Stabilising needs a finite t_max, got $t_max"))
        (1 <= N_stable < Np) || throw(ArgumentError(
            "Stabilising needs 1 ≤ N_stable < N_prog, got $N_stable, $Np"))
        return new{B}(base, Int(N_stable), Np, tm, StabilisingDose(; m_max).m_max)
    end
end

function next_block(s::Stabilising, pop, i::Int, t0::Float64)
    _schedule_done(pop, i, t0, s.N_prog, s.t_max) && return nothing
    tmax = t0 + s.t_max
    if i == 0
        return _sched_block(s.base, "stabilising: reduce", FixedDose(s.m_max), tmax;
                            stop_above = s.N_prog, stop_below = s.N_stable)
    elseif i == 1
        return _sched_block(s.base, "stabilising: hold", StabilisingDose(; m_max = s.m_max), tmax;
                            stop_above = s.N_prog)
    end
    return nothing
end

"""
    CustomDose(base, dose::AbstractDose; N_prog, t_max = Inf)

One block labelled `"custom dose"` run with any dose (typically a [`SizeDose`](@ref) `m(N)` or
[`TimeDose`](@ref) `m(t)`) until progression, extinction or `t_max`. The block starts at the
schedule start, so a `TimeDose` measures time from the schedule start. No `manuscript_v1`
counterpart (an extension). `t_max` is measured from the schedule start.
With the default `t_max = Inf` a run that can neither reach `N_prog` nor die out (for
example one held below `N_prog` by its density dependence) never ends: give a finite
`t_max` unless progression or extinction is certain.
"""
struct CustomDose{B<:BirthDeathBlock,Z<:AbstractDose} <: AbstractSchedule
    base::B
    dose::Z
    N_prog::Int
    t_max::Float64
    function CustomDose(base::B, dose::Z; N_prog, t_max = Inf) where {B<:BirthDeathBlock,Z<:AbstractDose}
        Np, tm = _check_common(N_prog, t_max)
        return new{B,Z}(base, dose, Np, tm)
    end
end

function next_block(s::CustomDose, pop, i::Int, t0::Float64)
    _schedule_done(pop, i, t0, s.N_prog, s.t_max) && return nothing
    i == 0 || return nothing
    return _sched_block(s.base, "custom dose", s.dose, t0 + s.t_max; stop_above = s.N_prog)
end

"""
    FunctionSchedule(f)

A schedule defined by a function `f(pop, i, t0) -> Union{BirthDeathBlock, Nothing}`, called
with the population, the number of blocks already run and the schedule start time. The blocks
are used as returned (no validation, no built-in end rules): return `nothing` to finish.
Remember that a block's `tmax` is absolute, e.g. `remake(base; tmax = t0 + 5.0)`.
"""
struct FunctionSchedule{F} <: AbstractSchedule
    f::F
end

next_block(s::FunctionSchedule, pop, i::Int, t0::Float64) = s.f(pop, i, t0)
