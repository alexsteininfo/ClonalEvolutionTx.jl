"""
    AbstractMutation

Supertype of mutation kinds (`Drivers`, `Resistance`, `CustomMutation`). Each has a
probability `ν` per daughter per division and a method `mutate(m, parent, rng)`.
"""
abstract type AbstractMutation end

function _check_ν(ν)
    0 <= ν <= 1 || throw(ArgumentError("ν must lie in [0, 1], got $ν"))
    return Float64(ν)
end

_check_Δs(Δs::Real) = Float64(Δs)
_check_Δs(Δs::Sampleable{Univariate}) = Δs
_check_Δs(Δs) = throw(ArgumentError("Δs must be a number or a univariate distribution, got $(typeof(Δs))"))

function _check_α(α::Real)
    0 <= α <= 1 || throw(ArgumentError("α must lie in [0, 1], got $α"))
    return Float64(α)
end
_check_α(α::Sampleable{Univariate}) = α
_check_α(α) = throw(ArgumentError("α must be a number in [0, 1] or a univariate distribution, got $(typeof(α))"))

"""
    Drivers(; ν, Δs, update = :additive)

Driver mutations: with probability `ν` per daughter per division the child's `s` changes by
`Δs` (a number or univariate distribution, drawn per mutation); `α` is untouched.
`update` is `:additive` (`s′ = s + Δs`), `:multiplicative` (`1 + s′ = (1 + s)(1 + Δs)`) or a
function `(s, Δs) -> s′`.
"""
struct Drivers{D,U} <: AbstractMutation
    ν::Float64
    Δs::D
    update::U
end

function Drivers(; ν, Δs, update = :additive)
    ν = _check_ν(ν)
    Δs = _check_Δs(Δs)
    rule = _update_rule(update)
    return Drivers{typeof(Δs),typeof(rule)}(ν, Δs, rule)
end

"""
    Resistance(; ν, α, Δs = 0.0, update = :additive)

Resistance mutations: with probability `ν` per daughter per division a new `α` is drawn
(`α` a number in [0, 1] or a univariate distribution; a draw outside [0, 1] throws
`ArgumentError`) and the child gets `min(α_drawn, α_parent)`, so cells only become more
resistant. `Δs` (cost or benefit, either sign) is a number, a distribution, or a function of
the drawn `α`; it is applied to `s` with the `update` rule of [`Drivers`](@ref), even if the
min rule leaves `α` unchanged.
"""
struct Resistance{A,D,U} <: AbstractMutation
    ν::Float64
    α::A
    Δs::D
    update::U
end

_check_resistance_Δs(Δs::Function) = Δs
_check_resistance_Δs(Δs) = _check_Δs(Δs)

function Resistance(; ν, α, Δs = 0.0, update = :additive)
    ν = _check_ν(ν)
    α = _check_α(α)
    Δs = _check_resistance_Δs(Δs)
    rule = _update_rule(update)
    return Resistance{typeof(α),typeof(Δs),typeof(rule)}(ν, α, Δs, rule)
end

"""
    CustomMutation(; ν, f)

Mutation with probability `ν` per daughter per division whose effect is
`f(rng, traits) -> Traits`, for anything not covered by [`Drivers`](@ref) or [`Resistance`](@ref).
"""
struct CustomMutation{F} <: AbstractMutation
    ν::Float64
    f::F
end

CustomMutation(; ν, f) = CustomMutation(_check_ν(ν), f)

"""
    mutate(m::AbstractMutation, parent::Traits, rng::AbstractRNG) -> Traits

Traits of a daughter in which mutation `m` occurred. All randomness comes from `rng`.
"""
function mutate end

function mutate(m::Drivers, parent::Traits, rng::AbstractRNG)
    Δs = _draw(m.Δs, rng)
    return Traits(_update_s(m.update, parent.s, Δs), parent.α)
end

function mutate(m::Resistance, parent::Traits, rng::AbstractRNG)
    α = Float64(_draw(m.α, rng))
    0.0 <= α <= 1.0 || throw(ArgumentError("drawn resistance α must lie in [0, 1], got $α"))
    Δs = _resistance_Δs(m.Δs, α, rng)
    return Traits(_update_s(m.update, parent.s, Δs), min(α, parent.α))
end

_resistance_Δs(Δs::Function, α, rng) = Δs(α)
_resistance_Δs(Δs, α, rng) = _draw(Δs, rng)

mutate(m::CustomMutation, parent::Traits, rng::AbstractRNG) = m.f(rng, parent)::Traits

"""
    mutation_kind(m::AbstractMutation) -> Symbol

`:driver`, `:resistance` or `:custom`.
"""
function mutation_kind end
mutation_kind(::Drivers) = :driver
mutation_kind(::Resistance) = :resistance
mutation_kind(::CustomMutation) = :custom
