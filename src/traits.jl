"""
    Traits(s, α)

Heritable traits of a clone: net fitness effect `s` (finite, either sign) and drug
resistance `α ∈ [0, 1]` (`α = 1`: fully drug-sensitive start, smaller `α` is more resistant).
Throws `ArgumentError` for non-finite `s` or `α` outside `[0, 1]`.
"""
struct Traits
    s::Float64
    α::Float64
    function Traits(s::Float64, α::Float64)
        isfinite(s) || throw(ArgumentError("Traits: s must be finite, got $s"))
        0.0 <= α <= 1.0 || throw(ArgumentError("Traits: α must lie in [0, 1], got $α"))
        new(s, α)
    end
end

Traits(s::Real, α::Real) = Traits(Float64(s), Float64(α))

_additive(s, Δs) = s + Δs
_multiplicative(s, Δs) = (1 + s) * (1 + Δs) - 1

# Resolve an `update` keyword (:additive, :multiplicative or a function) to a function (s, Δs) -> s′.
_update_rule(f::Function) = f
function _update_rule(u::Symbol)
    u === :additive && return _additive
    u === :multiplicative && return _multiplicative
    throw(ArgumentError("update must be :additive, :multiplicative or a function (s, Δs) -> s′, got :$u"))
end
_update_rule(u) = throw(ArgumentError("update must be :additive, :multiplicative or a function (s, Δs) -> s′, got $u"))

_update_s(update::Function, s::Float64, Δs::Real) = Float64(update(s, Δs))
_update_s(update::Symbol, s::Float64, Δs::Real) = _update_s(_update_rule(update), s, Δs)

_draw(x::Real, rng::AbstractRNG) = x
_draw(d::Sampleable{Univariate}, rng::AbstractRNG) = rand(rng, d)
