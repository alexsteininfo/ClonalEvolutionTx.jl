"""
    AbstractDose

How the treatment level `m ∈ [0, 1]` of a block is set: [`FixedDose`](@ref),
[`StabilisingDose`](@ref), [`SizeDose`](@ref) or [`TimeDose`](@ref) (SPEC §4.6).
"""
abstract type AbstractDose end

function _check_level(name, m)
    (isfinite(m) && 0.0 <= m <= 1.0) || throw(ArgumentError("$name: level must lie in [0, 1], got $m"))
    return Float64(m)
end

function _check_m_max(name, m_max)
    (isfinite(m_max) && 0.0 < m_max <= 1.0) || throw(ArgumentError("$name: m_max must lie in (0, 1], got $m_max"))
    return Float64(m_max)
end

"""
    FixedDose(m)

Constant treatment level `m ∈ [0, 1]` for the whole block.
"""
struct FixedDose <: AbstractDose
    m::Float64
    FixedDose(m::Real) = new(_check_level("FixedDose", m))
end

"""
    StabilisingDose(; m_max = 1.0)

Treatment level set after every event so that the expected net growth of the whole
population is zero, capped at `m_max ∈ (0, 1]` (SPEC §4.6).
"""
struct StabilisingDose <: AbstractDose
    m_max::Float64
    StabilisingDose(; m_max::Real = 1.0) = new(_check_m_max("StabilisingDose", m_max))
end

"""
    SizeDose(f)

Treatment level `m = f(N)`, a user function of the population size `N::Int`, re-evaluated
after every event. A level outside `[0, 1]` is an error.
"""
struct SizeDose{F} <: AbstractDose
    f::F
end

"""
    TimeDose(f; m_max = 1.0)

Treatment level `m = f(t_rel)`, a user function of the time since the block (or schedule)
start, which may vary between events. `m_max ∈ (0, 1]` bounds the level (a tighter value
only raises the acceptance rate of the exact thinning). An evaluated level outside
`[0, m_max]` is an error.
"""
struct TimeDose{F} <: AbstractDose
    f::F
    m_max::Float64
    TimeDose(f::F; m_max::Real = 1.0) where {F} = new{F}(f, _check_m_max("TimeDose", m_max))
end

# Stabilising level (SPEC §4.6): D(m) = (B·ΣNβ − ΣNδ − C·N) / ΣNα, clamped to [0, m_max].
# Order of the checks avoids any division by zero and any NaN.
function _stabilising_level(r::AbstractResponse, m_max::Real; B::Real, C::Real, N::Real,
                            sumNβ::Real, sumNδ::Real, sumNα::Real)
    num = B * sumNβ - sumNδ - C * N
    num <= 0 && return 0.0
    sumNα <= 0 && return Float64(m_max)   # ≤: a running sum may hold tiny negative round-off
    x = num / sumNα
    x >= drug_death(r, m_max) && return Float64(m_max)
    return inverse_drug_death(r, x)
end

function _size_level(d::SizeDose, N::Integer)
    m = d.f(N)
    (m isa Real && isfinite(m) && 0.0 <= m <= 1.0) ||
        throw(ArgumentError("SizeDose: m(N) must lie in [0, 1], got $m at N = $N"))
    return Float64(m)
end

function _time_level(d::TimeDose, t_rel::Real)
    m = d.f(Float64(t_rel))
    (m isa Real && isfinite(m) && 0.0 <= m <= d.m_max) ||
        throw(ArgumentError("TimeDose: m(t) must lie in [0, m_max = $(d.m_max)], got $m at t = $t_rel"))
    return Float64(m)
end

_describe(d::FixedDose) = "fixed m = $(d.m)"
_describe(::StabilisingDose) = "stabilising"
_describe(::SizeDose) = "m(N)"
_describe(::TimeDose) = "m(t)"
