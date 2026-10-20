"""
    AbstractResponse

Drug response `D(m)`: the extra death rate per unit of resistance trait `α` at treatment
level `m ∈ [0, 1]` (SPEC §4.2). A subtype must define [`drug_death`](@ref) and may define
[`inverse_drug_death`](@ref); `D` must be monotone non-decreasing on `[0, 1]` with
`D(0) = 0`. Without an inverse, a bisection fallback (which relies on exactly that
assumption) is used.
"""
abstract type AbstractResponse end

"""
    LinearResponse(d_S)

Linear drug response `D(m) = d_S·m`. `d_S ≥ 0` (`d_S = 0` means no drug effect); throws
`ArgumentError` otherwise.
"""
struct LinearResponse <: AbstractResponse
    d_S::Float64
    function LinearResponse(d_S::Real)
        (isfinite(d_S) && d_S >= 0) || throw(ArgumentError("LinearResponse: d_S must be finite and ≥ 0, got $d_S"))
        new(Float64(d_S))
    end
end

"""
    HillResponse(d_max, n, EC50)

Hill drug response `D(m) = d_max·mⁿ / (EC50ⁿ + mⁿ)`. All parameters must be finite and
positive; throws `ArgumentError` otherwise.
"""
struct HillResponse <: AbstractResponse
    d_max::Float64
    n::Float64
    EC50::Float64
    function HillResponse(d_max::Real, n::Real, EC50::Real)
        for (name, v) in (("d_max", d_max), ("n", n), ("EC50", EC50))
            (isfinite(v) && v > 0) || throw(ArgumentError("HillResponse: $name must be finite and > 0, got $v"))
        end
        new(Float64(d_max), Float64(n), Float64(EC50))
    end
end

"""
    drug_death(r::AbstractResponse, m::Real) -> Float64

The drug-induced death rate `D(m)` of response `r` at treatment level `m`.
"""
function drug_death end

drug_death(r::LinearResponse, m::Real) = r.d_S * Float64(m)
function drug_death(r::HillResponse, m::Real)
    mn = Float64(m)^r.n
    return r.d_max * mn / (r.EC50^r.n + mn)
end

"""
    inverse_drug_death(r::AbstractResponse, x::Real) -> Float64

The treatment level `m ∈ [0, 1]` with `D(m) = x`. Closed form for [`LinearResponse`](@ref)
and [`HillResponse`](@ref); for any other subtype that only defines [`drug_death`](@ref),
bisection on `[0, 1]` to `1e-12`, assuming `D` is monotone non-decreasing with `D(0) = 0`.
Throws `ArgumentError` if `x` lies outside `[0, D(1)]` (nothing is clamped).
"""
function inverse_drug_death(r::AbstractResponse, x::Real)
    xf = Float64(x)
    _check_inverse_range(r, xf)
    lo, hi = 0.0, 1.0
    while hi - lo > 1e-12
        mid = (lo + hi) / 2
        drug_death(r, mid) < xf ? (lo = mid) : (hi = mid)
    end
    return (lo + hi) / 2
end

function _check_inverse_range(r::AbstractResponse, x::Float64)
    top = drug_death(r, 1.0)
    (0.0 <= x <= top) ||
        throw(ArgumentError("inverse_drug_death: x = $x is outside the range [0, D(1)] = [0, $top]"))
    return nothing
end

function inverse_drug_death(r::LinearResponse, x::Real)
    xf = Float64(x)
    _check_inverse_range(r, xf)
    return xf == 0.0 ? 0.0 : xf / r.d_S      # xf > 0 implies d_S > 0 after the range check
end

function inverse_drug_death(r::HillResponse, x::Real)
    xf = Float64(x)
    _check_inverse_range(r, xf)
    return r.EC50 * (xf / (r.d_max - xf))^(1 / r.n)
end
