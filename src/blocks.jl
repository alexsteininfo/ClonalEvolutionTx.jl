"""
    BirthDeathBlock(; birth, death, birth_density = Returns(1.0), death_density = Returns(0.0),
                    response = LinearResponse(0.0), dose = FixedDose(0.0), μ = 0.0,
                    mutations = (), stop = Returns(false), stop_above = typemax(Int),
                    stop_below = -1, tmax = Inf, restart_on_extinction = false,
                    max_clones = typemax(Int), label = "block")

One stretch of birth–death dynamics with fixed rules, simulated by `simulate!`.

- `birth`, `death`: functions `s -> β(s)` and `s -> δ(s)` of a clone's selection trait.
- `birth_density`, `death_density`: functions `N -> B(N)` and `N -> C(N)` of the population
  size, common to all clones: `B(N)` multiplies every birth rate and `C(N)` is added to
  every death rate, so a cell of a clone with traits `(s, α)` has birth rate
  `b = β(s)·B(N)` and death rate `d = δ(s) + C(N) + α·D(m)`.
- `response`: how the drug level raises death (see [`LinearResponse`](@ref), [`HillResponse`](@ref)).
- `dose`: how the drug level `m` is set (see [`FixedDose`](@ref) and the other doses).
- `μ`: mean number of neutral mutations per daughter per division (Poisson), `μ ≥ 0`.
- `mutations`: the mutation kinds ([`Drivers`](@ref), [`Resistance`](@ref),
  [`CustomMutation`](@ref)) as a tuple or vector, applied to each daughter in list order.
- `stop`: a function of the population; the block ends as soon as it returns `true`.
- `stop_above`, `stop_below`: the block ends when the population size reaches
  `stop_above` (`N ≥ stop_above`) or falls to `stop_below` (`N ≤ stop_below`).
- `tmax`: the **absolute** simulation time at which the block ends (not a duration);
  `Inf` for none, never `NaN`.
- `restart_on_extinction`: restart from the state at the start of the block if the
  population dies out.
- `max_clones`: no new clone is created once this many clones exist (`≥ 1`).
- `label`: name recorded in the population's log.

Throws `ArgumentError` for `μ < 0`, `isnan(tmax)`, an entry of `mutations` that is not an
[`AbstractMutation`](@ref), or `max_clones < 1`. Use [`remake`](@ref) to change fields.
"""
struct BirthDeathBlock{B,D,BD,CD,R<:AbstractResponse,Z<:AbstractDose,M<:Tuple,S}
    birth::B
    death::D
    birth_density::BD
    death_density::CD
    response::R
    dose::Z
    μ::Float64
    mutations::M
    stop::S
    stop_above::Int
    stop_below::Int
    tmax::Float64
    restart_on_extinction::Bool
    max_clones::Int
    label::String
end

_as_mutation_tuple(m::Tuple) = _check_mutations(m)
_as_mutation_tuple(m::AbstractVector) = _check_mutations(Tuple(m))
_as_mutation_tuple(m) = throw(ArgumentError(
    "mutations must be a tuple or vector of AbstractMutation, got $(typeof(m))"))

function _check_mutations(m::Tuple)
    for x in m
        x isa AbstractMutation || throw(ArgumentError(
            "mutations must all be AbstractMutation (Drivers, Resistance, CustomMutation), got $(typeof(x))"))
    end
    return m
end

function BirthDeathBlock(; birth, death, birth_density = Returns(1.0),
                         death_density = Returns(0.0), response = LinearResponse(0.0),
                         dose = FixedDose(0.0), μ = 0.0, mutations = (),
                         stop = Returns(false), stop_above = typemax(Int), stop_below = -1,
                         tmax = Inf, restart_on_extinction = false,
                         max_clones = typemax(Int), label = "block")
    μ >= 0 || throw(ArgumentError("μ must be ≥ 0, got $μ"))
    isnan(tmax) && throw(ArgumentError("tmax must not be NaN"))
    max_clones >= 1 || throw(ArgumentError("max_clones must be ≥ 1, got $max_clones"))
    muts = _as_mutation_tuple(mutations)
    return BirthDeathBlock(birth, death, birth_density, death_density, response, dose,
                           Float64(μ), muts, stop, Int(stop_above), Int(stop_below),
                           Float64(tmax), Bool(restart_on_extinction), Int(max_clones),
                           String(label))
end

"""
    remake(block::BirthDeathBlock; kwargs...) -> BirthDeathBlock

A copy of `block` with the given fields replaced, e.g. `remake(block; μ = 2.0, label = "growth")`.
The result is validated like a freshly constructed block.
"""
function remake(block::BirthDeathBlock; kwargs...)
    names = fieldnames(BirthDeathBlock)
    for k in keys(kwargs)
        k in names || throw(ArgumentError("BirthDeathBlock has no field $k"))
    end
    return BirthDeathBlock(; (n => getfield(block, n) for n in names)..., kwargs...)
end
