"""
    CloneTable()

Clone-level state of a population, stored as struct-of-arrays. Clones are numbered in
order of creation (`parent == 0` marks a founder clone) and are never removed, so the
history of extinct clones is kept. Per clone: `s`, `α`, `kind`
(`:founder | :driver | :resistance | :custom`), `origin_time`, `origin_cell`, `size`, and
the cached rates `β`, `δ`. The living list `alive` (with `alive_pos`, 0 when not living),
the maxima `maxβ`, `maxδ`, `maxα` over living clones (0.0 if none) and the running sums
`sumNβ = Σ Nₖβₖ`, `sumNδ`, `sumNα` give the simulation loop its exact O(1) rate bound and
dose level.
"""
mutable struct CloneTable
    parent::Vector{Int32}
    s::Vector{Float64}
    α::Vector{Float64}
    kind::Vector{Symbol}
    origin_time::Vector{Float64}
    origin_cell::Vector{Int64}
    size::Vector{Int64}
    β::Vector{Float64}
    δ::Vector{Float64}
    alive::Vector{Int32}
    alive_pos::Vector{Int32}
    maxβ::Float64
    maxδ::Float64
    maxα::Float64
    sumNβ::Float64
    sumNδ::Float64
    sumNα::Float64
    ever_alive::Vector{Bool}        # internal: enforces "no revival"
    argmaxβ::Int32                  # internal: clone holding each maximum (0 if none)
    argmaxδ::Int32
    argmaxα::Int32
    birth::Function                 # bound rate functions, called only at creation/binding
    death::Function
end

CloneTable() = CloneTable(Int32[], Float64[], Float64[], Symbol[], Float64[], Int64[], Int64[],
                          Float64[], Float64[], Int32[], Int32[], 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
                          Bool[], Int32(0), Int32(0), Int32(0), s -> 1.0, s -> 0.0)

"""
    n_clones(ct::CloneTable) -> Int

Number of clones ever created (living or extinct).
"""
n_clones(ct::CloneTable) = length(ct.s)

"""
    living_clones(ct::CloneTable) -> Vector{Int32}

Copy of the ids of the clones with `size > 0`, in living-list order.
"""
living_clones(ct::CloneTable) = copy(ct.alive)

"""
    traits(ct::CloneTable, k) -> Traits

Traits `(s, α)` of clone `k`.
"""
traits(ct::CloneTable, k::Integer) = Traits(ct.s[k], ct.α[k])

function _checked_rates(birth, death, s::Float64, k)
    β = Float64(birth(s))
    δ = Float64(death(s))
    β >= 0 || throw(ArgumentError("negative or non-finite birth rate β = $β for clone $k (s = $s)"))
    δ >= 0 || throw(ArgumentError("negative or non-finite death rate δ = $δ for clone $k (s = $s)"))
    isfinite(β) && isfinite(δ) ||
        throw(ArgumentError("non-finite rate (β = $β, δ = $δ) for clone $k (s = $s)"))
    return β, δ
end

function _add_clone!(ct::CloneTable, parent::Integer, t::Traits, kind::Symbol, time::Float64, cell::Int64)
    k = n_clones(ct) + 1
    β, δ = _checked_rates(ct.birth, ct.death, t.s, k)
    push!(ct.parent, parent)
    push!(ct.s, t.s)
    push!(ct.α, t.α)
    push!(ct.kind, kind)
    push!(ct.origin_time, time)
    push!(ct.origin_cell, cell)
    push!(ct.size, 0)
    push!(ct.β, β)
    push!(ct.δ, δ)
    push!(ct.alive_pos, 0)
    push!(ct.ever_alive, false)
    return Int32(k)
end

# Recompute the three maxima (and their holders) over the living clones.
function _recompute_maxima!(ct::CloneTable)
    ct.maxβ = ct.maxδ = ct.maxα = 0.0
    ct.argmaxβ = ct.argmaxδ = ct.argmaxα = 0
    for k in ct.alive
        _raise_maxima!(ct, k)
    end
    return ct
end

# O(1): fold clone k into the maxima.
function _raise_maxima!(ct::CloneTable, k::Integer)
    if ct.argmaxβ == 0 || ct.β[k] > ct.maxβ
        ct.maxβ = ct.β[k]; ct.argmaxβ = k
    end
    if ct.argmaxδ == 0 || ct.δ[k] > ct.maxδ
        ct.maxδ = ct.δ[k]; ct.argmaxδ = k
    end
    if ct.argmaxα == 0 || ct.α[k] > ct.maxα
        ct.maxα = ct.α[k]; ct.argmaxα = k
    end
    return nothing
end

function _change_size!(ct::CloneTable, k::Integer, Δ::Integer)
    old = ct.size[k]
    new = old + Δ
    new >= 0 || error("clone $k would get negative size $new (size $old, change $Δ)")
    if old == 0 && new > 0 && ct.ever_alive[k]
        error("clone $k went extinct and cannot be revived")
    end
    ct.size[k] = new
    ct.sumNβ += Δ * ct.β[k]
    ct.sumNδ += Δ * ct.δ[k]
    ct.sumNα += Δ * ct.α[k]
    if old == 0 && new > 0
        push!(ct.alive, k)
        ct.alive_pos[k] = length(ct.alive)
        ct.ever_alive[k] = true
        _raise_maxima!(ct, k)
    elseif old > 0 && new == 0
        pos = ct.alive_pos[k]
        last = ct.alive[end]
        ct.alive[pos] = last
        ct.alive_pos[last] = pos
        pop!(ct.alive)
        ct.alive_pos[k] = 0
        if k == ct.argmaxβ || k == ct.argmaxδ || k == ct.argmaxα
            _recompute_maxima!(ct)
        end
        isempty(ct.alive) && (ct.sumNβ = ct.sumNδ = ct.sumNα = 0.0)   # no drift once empty
    end
    return ct
end

function _bind_rates!(ct::CloneTable, birth, death)
    n = n_clones(ct)
    β = Vector{Float64}(undef, n)
    δ = Vector{Float64}(undef, n)
    for k in 1:n
        β[k], δ[k] = _checked_rates(birth, death, ct.s[k], k)
    end
    copyto!(ct.β, β)
    copyto!(ct.δ, δ)
    ct.birth = birth
    ct.death = death
    ct.sumNβ = ct.sumNδ = ct.sumNα = 0.0
    for k in ct.alive
        ct.sumNβ += ct.size[k] * ct.β[k]
        ct.sumNδ += ct.size[k] * ct.δ[k]
        ct.sumNα += ct.size[k] * ct.α[k]
    end
    _recompute_maxima!(ct)
    return ct
end

# Deep copy: every column is copied; the bound rate functions are shared by reference.
function Base.copy(ct::CloneTable)
    return CloneTable(copy(ct.parent), copy(ct.s), copy(ct.α), copy(ct.kind),
                      copy(ct.origin_time), copy(ct.origin_cell), copy(ct.size),
                      copy(ct.β), copy(ct.δ), copy(ct.alive), copy(ct.alive_pos),
                      ct.maxβ, ct.maxδ, ct.maxα, ct.sumNβ, ct.sumNδ, ct.sumNα,
                      copy(ct.ever_alive), ct.argmaxβ, ct.argmaxδ, ct.argmaxα,
                      ct.birth, ct.death)
end

# Rebuild everything derived from `size` (and the rate columns): the living list and its
# positions, the running sums, the maxima with their holders, and `ever_alive` (a clone that
# is alive now has been alive). Used where sizes are set wholesale (sampling, loading), which
# the one-step `_change_size!` cannot do for a clone that was alive before.
function _rebuild_derived!(ct::CloneTable)
    empty!(ct.alive)
    fill!(ct.alive_pos, 0)
    ct.sumNβ = ct.sumNδ = ct.sumNα = 0.0
    for k in 1:n_clones(ct)
        ct.size[k] > 0 || continue
        push!(ct.alive, k)
        ct.alive_pos[k] = length(ct.alive)
        ct.ever_alive[k] = true
        ct.sumNβ += ct.size[k] * ct.β[k]
        ct.sumNδ += ct.size[k] * ct.δ[k]
        ct.sumNα += ct.size[k] * ct.α[k]
    end
    _recompute_maxima!(ct)
    return ct
end
