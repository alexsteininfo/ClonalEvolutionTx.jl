# The simulation engine: exact thinning (uniformisation) of the clone-structured
# birth–death process (SPEC §5). One exponential clock at rate R̄·N drives the whole
# population, where R̄ bounds the total rate b + d of every living clone. At each tick one
# living cell is picked uniformly; it divides with probability b/R̄, dies with probability
# d/R̄, and otherwise the tick is rejected. R̄ comes in O(1) per event from the clone
# table's cached maxima, so a tick costs O(1).

# ── Accumulator hooks ─────────────────────────────────────────────────────────
# Measurements hook into the loop here; `accumulator = nothing` records nothing.
# `_acc_start!`: once at block start, before the first stop check.
# `_acc_until!`: the state is about to advance to time `t` (every tick, and `t = tmax`
#     at a `:tmax` end); the population still holds the state before `t`. A measurement
#     may move the clock `pop.t` forward to any time ≤ `t` (the state there is the same),
#     so the loop must compute `t` before the call and set `pop.t = t` right after it.
# `_acc_after_event!`: after each accepted event (and the dose update that follows it).
# `_acc_block_end!`: once at block end, before the log entry is written.
# `_acc_mark` / `_acc_rollback!`: a restart rolls the accumulator back to the mark taken
#     right after `_acc_start!`.
_acc_start!(::Nothing, pop, block) = nothing
_acc_until!(::Nothing, pop, t) = nothing
_acc_after_event!(::Nothing, pop) = nothing
_acc_block_end!(::Nothing, pop) = nothing
_acc_mark(::Nothing) = nothing
_acc_rollback!(::Nothing, mark) = nothing

# A `MeasurementAccumulator` (measurements.jl; its docstring states the semantics). It
# never keeps a reference to a field of `pop`: a restart replaces them.
function _acc_start!(acc::MeasurementAccumulator, pop, block)
    t = pop.t
    _check_pending_times(acc, t)                 # throws before the accumulator changes
    first_start = !acc.started
    _align_grids!(acc, t)
    acc.dose, acc.t_block, acc.label = block.dose, t, block.label
    _record_due!(acc, pop, t)                    # time items at exactly the start time
    N = popsize(pop)
    prev = first_start ? N : acc.last_N          # a first sight is no arrival
    _check_state!(acc, pop, prev, N, !_chained(acc, pop))
    acc.last_N = N
    return nothing
end

@inline function _acc_until!(acc::MeasurementAccumulator, pop, t)
    t >= acc.next_due && _record_due!(acc, pop, t)   # one comparison when nothing is due
    return nothing
end

@inline function _acc_after_event!(acc::MeasurementAccumulator, pop)
    N = popsize(pop)
    isempty(acc.state) || _check_state!(acc, pop, acc.last_N, N, true)
    acc.last_N = N
    return nothing
end

function _acc_block_end!(acc::MeasurementAccumulator, pop)
    for i in acc.blockend
        _snapshot!(acc, pop, pop.t, acc.spec.triggers[i])
    end
    acc.end_t, acc.end_N, acc.end_log = pop.t, popsize(pop), length(pop.log) + 1   # after the push
    return nothing
end

_acc_mark(acc::MeasurementAccumulator) =
    (n_trajectory = length(acc.trajectory), n_snapshots = length(acc.snapshots),
     traj_k = acc.traj_k, every_k = copy(acc.every_k), attime_next = acc.attime_next,
     next_due = acc.next_due, fired = copy(acc.fired), last_N = acc.last_N)

function _acc_rollback!(acc::MeasurementAccumulator, mark)
    resize!(acc.trajectory, mark.n_trajectory)
    resize!(acc.snapshots, mark.n_snapshots)
    acc.traj_k = mark.traj_k
    copyto!(acc.every_k, mark.every_k)
    acc.attime_next = mark.attime_next
    acc.next_due = mark.next_due
    copyto!(acc.fired, mark.fired)
    acc.last_N = mark.last_N
    return nothing
end

# ── Dose hooks ────────────────────────────────────────────────────────────────
# `_drug_bound`: the drug death D the bound must cover (its term in R̄ is maxα·D), at
#     block start and after a restart.
# `_drug_at_tick`: the drug death used in the acceptance step of a tick at time `t`
#     (`t0` is the block start); it must not exceed the current bound term `D`.
# `_dose_after_event!`: update `pop.m` after an accepted event; returns the new bound term.
# `_dose_at_end!`: set `pop.m` to the level at the block's end time (before the log entry).
_drug_bound(dose::FixedDose, block::BirthDeathBlock, pop) = _checked_drug_death(block.response, dose.m)
_drug_at_tick(::FixedDose, block, pop, t, t0, D) = D
_dose_after_event!(::FixedDose, block, pop, D) = D

# Stabilising and size doses: `pop.m` holds the level of the current state (set at block
# start, restored by a restart), and is refreshed after every accepted event. D is computed
# once per level and used both in the bound and in the per-cell death, so b + d ≤ R̄ holds
# exactly.
const _StateDose = Union{StabilisingDose,SizeDose}
_drug_bound(::_StateDose, block::BirthDeathBlock, pop) = _checked_drug_death(block.response, pop.m)
_drug_at_tick(::_StateDose, block, pop, t, t0, D) = D
function _dose_after_event!(dose::_StateDose, block, pop, D)
    popsize(pop) == 0 && return D            # the block ends at the next stop check
    pop.m = _dose_level(dose, block, pop)
    return _checked_drug_death(block.response, pop.m)
end

# Time dose: the bound covers the largest admissible level, m_max; each tick evaluates the
# level at its own time. In floating point a monotone response can still give D(m) a few
# ulps above D(m_max) for m just below m_max (HillResponse does, through mⁿ), so the bound's
# drug term is D(m_max) raised by 4 ulps: a larger bound keeps the thinning exact and only
# wastes a negligible share of ticks. A D(m(t)) above even that raised bound means the
# response is not monotone, which would make the thinning inexact, so it is an error.
# D(m_max) = 0 stays 0, so that a block in which no event can happen is still detected.
function _drug_bound(d::TimeDose, block::BirthDeathBlock, pop)
    D = _checked_drug_death(block.response, d.m_max)
    return D == 0.0 ? D : nextfloat(D, 4)
end
function _drug_at_tick(d::TimeDose, block, pop, t, t0, D)
    pop.m = _time_level(d, t - t0)
    Dt = _checked_drug_death(block.response, pop.m)
    Dt <= D || throw(ArgumentError(
        "the drug response $(block.response) is not monotone: D(m = $(pop.m)) = $Dt exceeds " *
        "the bound $D (D(m_max = $(d.m_max)) plus 4 ulps), so the thinning would be inexact"))
    return Dt
end
_dose_after_event!(::TimeDose, block, pop, D) = D
# A time dose's level changes between events: at block end it is set to the level at the
# end time (the other doses already hold the level of the final state).
_dose_at_end!(d::TimeDose, pop, t0) = (pop.m = _time_level(d, pop.t - t0); nothing)
_dose_at_end!(::AbstractDose, pop, t0) = nothing

function _checked_drug_death(r::AbstractResponse, m::Float64)
    D = Float64(drug_death(r, m))
    (D >= 0 && isfinite(D)) ||
        throw(ArgumentError("the drug response must give a finite D(m) ≥ 0, got $D at m = $m"))
    return D
end

# ── Loop helpers ──────────────────────────────────────────────────────────────

# B(N) or C(N), checked: a negative or non-finite value is an error, never clamped.
@inline function _density(f, N::Int, name::String)
    x = Float64(f(N))
    (x >= 0 && isfinite(x)) || _bad_density(name, x, N)
    return x
end
@noinline _bad_density(name, x, N) =
    throw(ArgumentError("$name must be finite and ≥ 0, got $x at N = $N"))

# Clone of the k-th living cell.
@inline _clone_of(pop::TreePopulation, k::Int) = @inbounds pop.tree.clone[pop.alive[k]]
@inline _clone_of(pop::CellPopulation, k::Int) = @inbounds pop.cells[k].clone

# The stop checks, in order; `:none` when the block goes on.
@inline function _stop_reason(pop::AbstractPopulation, block::BirthDeathBlock, N::Int)
    N == 0 && return :extinct
    N >= block.stop_above && return :above
    N <= block.stop_below && return :below
    block.stop(pop) && return :stopped
    return :none
end

# Restore every field of `pop` from a fresh copy of `snapshot` (kept for later restarts).
function _restore!(pop::P, snapshot::P) where {P<:AbstractPopulation}
    fresh = copy(snapshot)
    for name in fieldnames(P)
        setfield!(pop, name, getfield(fresh, name))
    end
    return pop
end

"""
    simulate!(pop, block::BirthDeathBlock, rng = Xoshiro(); accumulator = nothing) -> pop

Simulate `block` on the population `pop` in place, continuing its clock, cells and clone
table, and append one [`BlockSummary`](@ref) (label, start and end time, end reason, ticks,
accepted events, restarts, dose, final drug level) to `pop.log`.

The process is simulated exactly, by thinning against a bound on the total rate of the
living clones. At the start and after every event the block ends, in this order of checks,
when the population is extinct (`:extinct`), `N ≥ stop_above` (`:above`),
`N ≤ stop_below` (`:below`) or `stop(pop)` returns `true` (`:stopped`). `tmax` is an
**absolute** time: if the next event would fall after it, the clock is set to exactly
`tmax` and the block ends (`:tmax`); discarding that event is exact because waiting times
are memoryless. If no event can occur (all rates zero), the clock jumps to `tmax`, or an
`ErrorException` is thrown when `tmax = Inf`.

With `restart_on_extinction = true`, a population that dies out is restored to its state
at the start of the block and simulated again with the continuing random stream, until a
stop other than extinction ends the block (a population certain to die out restarts
forever). The restarts are counted, and the ticks and accepted events of all attempts
are added up.

All randomness of the dynamics comes from `rng`; the neutral-mutation stream is seeded
from it at the block start. `accumulator` (a [`MeasurementAccumulator`](@ref)) collects
measurements during the run; `nothing` records nothing.

Throws `ArgumentError` if `tmax` lies before `pop.t`, if a rate is negative or not finite
(`β(s)`, `δ(s)`, `B(N)`, `C(N)` or `D(m)`), or if `restart_on_extinction` is set for an
empty population.
"""
function simulate!(pop::AbstractPopulation, block::BirthDeathBlock, rng::AbstractRNG = Xoshiro();
                   accumulator = nothing)
    t_start = pop.t
    block.tmax < t_start && throw(ArgumentError(
        "tmax = $(block.tmax) lies before the population's time t = $t_start (tmax is absolute)"))
    block.restart_on_extinction && popsize(pop) == 0 && throw(ArgumentError(
        "restart_on_extinction is set but the population is already empty, so no restart can succeed"))
    _start_block!(pop, block, rng)
    _acc_start!(accumulator, pop, block)
    reason, ticks, accepted, restarts = _run_block!(pop, block, rng, accumulator, t_start)
    _dose_at_end!(block.dose, pop, t_start)
    _acc_block_end!(accumulator, pop)
    push!(pop.log, BlockSummary(block.label, t_start, pop.t, reason, ticks, accepted,
                                restarts, _describe(block.dose), pop.m))
    return pop
end

# The loop of one block, specialised on the population, block (hence dose) and
# accumulator types. Returns (end reason, ticks, accepted events, restarts).
function _run_block!(pop::AbstractPopulation, block::BirthDeathBlock, rng::AbstractRNG,
                     acc, t0::Float64)
    snapshot = block.restart_on_extinction ? copy(pop) : nothing
    mark = _acc_mark(acc)
    dose = block.dose
    tmax = block.tmax
    ticks = accepted = restarts = 0
    D = _drug_bound(dose, block, pop)
    while true
        N = popsize(pop)
        reason = _stop_reason(pop, block, N)
        if reason === :extinct && snapshot !== nothing
            _restore!(pop, snapshot)
            _acc_rollback!(acc, mark)
            restarts += 1
            D = _drug_bound(dose, block, pop)
            continue
        end
        reason === :none || return reason, ticks, accepted, restarts

        # The bound, with the same expression tree as b + d below, so that b + d ≤ R̄ holds
        # exactly in floating point too.
        ct = pop.clones
        B = _density(block.birth_density, N, "birth_density B(N)")
        C = _density(block.death_density, N, "death_density C(N)")
        R̄ = ct.maxβ * B + ((ct.maxδ + C) + ct.maxα * D)
        if R̄ == 0.0                           # nothing can ever happen again
            isfinite(tmax) || error("no event can occur and the block has no tmax")
            _acc_until!(acc, pop, tmax)
            pop.t = tmax
            return :tmax, ticks, accepted, restarts
        end

        # Ticks until one is accepted; a rejected tick only advances the clock.
        while true
            t = pop.t + randexp(rng) / (R̄ * N)
            if t > tmax
                _acc_until!(acc, pop, tmax)
                pop.t = tmax
                return :tmax, ticks, accepted, restarts
            end
            ticks += 1
            _acc_until!(acc, pop, t)
            pop.t = t
            k = rand(rng, 1:N)
            u = rand(rng) * R̄
            c = _clone_of(pop, k)
            Dt = _drug_at_tick(dose, block, pop, t, t0, D)
            b = @inbounds ct.β[c] * B
            d = @inbounds (ct.δ[c] + C) + ct.α[c] * Dt
            if u < b
                _divide!(pop, k, t, block, rng)
                break
            elseif u < b + d
                _die!(pop, k)
                break
            end
        end
        accepted += 1
        D = _dose_after_event!(dose, block, pop, D)
        _acc_after_event!(acc, pop)
    end
end
