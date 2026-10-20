# Statistics (SPEC §8): tree statistics from neutral mutations, clone statistics for both
# storage modes, Hill numbers and the statistic types that snapshots evaluate.

function _need_tree(pop::AbstractPopulation, what::String)
    pop isa TreePopulation || throw(ArgumentError(
        "$what needs a population in tree mode (initialize_population(...; mode = :tree)), got a cell-mode population"))
    return pop
end

# ── Tree statistics ───────────────────────────────────────────────────────────

"""
    site_frequency_spectrum(pop) -> Vector{Int}

Neutral-mutation site-frequency spectrum of a tree-mode population: `sfs[k]` is the number
of mutations carried by exactly `k` living cells; the vector has length `popsize(pop)`
(empty for an empty population). `pop.clonal_mutations` are counted at `k = N`. One
reverse sweep over the tree rows; a forest contributes every tree. Throws `ArgumentError`
in cell mode.
"""
function site_frequency_spectrum(pop::TreePopulation)
    tree = pop.tree
    N = popsize(pop)
    sfs = zeros(Int, N)
    N == 0 && return sfs
    counts = zeros(Int32, length(tree))
    @inbounds for a in pop.alive
        counts[a] = 1
    end
    parent, neutral = tree.parent, tree.neutral
    @inbounds for i in length(tree):-1:1
        _isdead(tree, i) && continue
        c = counts[i]
        c == 0 && continue
        sfs[c] += neutral[i]
        p = parent[i]
        p > NOPARENT && (counts[p] += c)
    end
    sfs[N] += pop.clonal_mutations
    return sfs
end
site_frequency_spectrum(pop::AbstractPopulation) = _need_tree(pop, "site_frequency_spectrum")

"""
    total_mutational_burden(pop) -> Int

Number of distinct neutral mutations in the living cells, `sum(site_frequency_spectrum(pop))`.
Tree mode only.
"""
total_mutational_burden(pop::AbstractPopulation) = sum(site_frequency_spectrum(pop))

# Index of the most recent common ancestor of nodes `i` and `j` (0 if in different trees).
# Children have larger indices than parents, so the larger index is never the ancestor.
function _mrca_idx(tree::LineageTree, i::Integer, j::Integer)
    i, j = Int32(i), Int32(j)
    @inbounds while i != j
        if i > j
            i = tree.parent[i]
        else
            j = tree.parent[j]
        end
        (i == NOPARENT || j == NOPARENT) && return Int32(0)
    end
    return i
end

"""
    clonal_mutations(pop) -> Int

Neutral mutations carried by every living cell: `pop.clonal_mutations` plus the path total
of the most recent common ancestor of all living cells. In a forest whose living cells sit
under more than one root there is no common ancestor and the result is
`pop.clonal_mutations`; the same holds for an empty population. Tree mode only.
"""
function clonal_mutations(pop::TreePopulation)
    isempty(pop.alive) && return pop.clonal_mutations
    tree = pop.tree
    m = pop.alive[1]
    for a in pop.alive
        m = _mrca_idx(tree, m, a)
        m == NOPARENT && return pop.clonal_mutations
    end
    return pop.clonal_mutations + Int(tree.total[m])
end
clonal_mutations(pop::AbstractPopulation) = _need_tree(pop, "clonal_mutations")

"""
    pairwise_distances(pop) -> Vector{Int}

Number of differing neutral mutations for every pair of living cells `i < j` (row-major
order, in living-cell order): the mutations on the paths below their most recent common
ancestor; cells in different trees differ by both full burdens (without the shared
`clonal_mutations`). Quadratic in the population size. Tree mode only.
"""
function pairwise_distances(pop::TreePopulation)
    tree, alive = pop.tree, pop.alive
    n = length(alive)
    out = Vector{Int}(undef, n * (n - 1) ÷ 2)
    k = 0
    for x in 1:n, y in (x + 1):n
        a, b = alive[x], alive[y]
        m = _mrca_idx(tree, a, b)
        shared = m == NOPARENT ? 0 : Int(tree.total[m])
        out[k += 1] = Int(tree.total[a]) + Int(tree.total[b]) - 2shared
    end
    return out
end
pairwise_distances(pop::AbstractPopulation) = _need_tree(pop, "pairwise_distances")

"""
    coalescence_times(pop) -> Vector{Float64}

For every pair `i < j` of living cells (as in [`pairwise_distances`](@ref)) the time from
`pop.t` back to the division of their most recent common ancestor. Cells in different
trees coalesce, by convention, at the earlier of their roots' birth times. Tree mode only.
"""
function coalescence_times(pop::TreePopulation)
    tree, alive, t = pop.tree, pop.alive, pop.t
    n = length(alive)
    out = Vector{Float64}(undef, n * (n - 1) ÷ 2)
    k = 0
    for x in 1:n, y in (x + 1):n
        a, b = alive[x], alive[y]
        m = _mrca_idx(tree, a, b)
        out[k += 1] = if m == NOPARENT
            t - min(tree.birthtime[_rootidx(tree, a)], tree.birthtime[_rootidx(tree, b)])
        else
            c = tree.left[m] != NOPARENT ? tree.left[m] : tree.right[m]   # a daughter's birth is the division
            c == NOPARENT ? 0.0 : t - tree.birthtime[c]
        end
    end
    return out
end
coalescence_times(pop::AbstractPopulation) = _need_tree(pop, "coalescence_times")

"""
    repaint_neutral!(pop::TreePopulation, μ, rng) -> pop

Redraw the neutral mutations of a tree-mode population (SPEC §4.5): every non-root node
gets a fresh `Poisson(μ)` count from `rng`, roots keep theirs, and the path totals are
recomputed. Dynamics, clones and `clonal_mutations` are untouched. Throws `ArgumentError`
for `μ < 0` and in cell mode.
"""
function repaint_neutral!(pop::TreePopulation, μ::Real, rng::AbstractRNG)
    μ >= 0 || throw(ArgumentError("repaint_neutral!: μ must be ≥ 0, got $μ"))
    tree = pop.tree
    d = Poisson(Float64(μ))
    @inbounds for i in 1:length(tree)
        _isdead(tree, i) && continue
        p = tree.parent[i]
        if p == NOPARENT
            tree.total[i] = tree.neutral[i]
        else
            tree.neutral[i] = rand(rng, d)
            tree.total[i] = tree.total[p] + tree.neutral[i]
        end
    end
    return pop
end
repaint_neutral!(pop::AbstractPopulation, μ::Real, rng::AbstractRNG) =
    _need_tree(pop, "repaint_neutral!")

# ── Statistics of both modes ──────────────────────────────────────────────────

"""
    mean_burden(pop) -> Float64

Mean number of neutral mutations per living cell (`NaN` for an empty population).
"""
mean_burden(pop::AbstractPopulation) =
    popsize(pop) == 0 ? NaN : mean(mutations_per_cell(pop))

"""
    clone_sizes(pop) -> Vector{Int}

Sizes of the living clones, in [`living_clones`](@ref) order.
"""
clone_sizes(pop::AbstractPopulation) = Int[pop.clones.size[k] for k in pop.clones.alive]

"""
    mutation_frequencies(pop) -> Vector{Int}

For every non-founder clone, in id order, its size plus the sizes of all its descendant
clones: the number of cells carrying that driver or resistance mutation.
"""
function mutation_frequencies(pop::AbstractPopulation)
    ct = pop.clones
    acc = Int[s for s in ct.size]
    for k in length(acc):-1:1            # parents have smaller ids than their children
        p = ct.parent[k]
        p > 0 && (acc[p] += acc[k])
    end
    return Int[acc[k] for k in eachindex(acc) if ct.kind[k] !== :founder]
end

"""
    n_sensitive(pop) -> Int

Cells in living clones with `α == 1` (fully drug-sensitive).
"""
function n_sensitive(pop::AbstractPopulation)
    ct = pop.clones
    return sum((ct.size[k] for k in ct.alive if ct.α[k] == 1.0); init = 0)
end

"""
    n_resistant(pop) -> Int

Cells in living clones with `α < 1`.
"""
function n_resistant(pop::AbstractPopulation)
    ct = pop.clones
    return sum((ct.size[k] for k in ct.alive if ct.α[k] != 1.0); init = 0)
end

# ── Hill numbers ──────────────────────────────────────────────────────────────

function _check_q(q::Real)
    isnan(q) && throw(ArgumentError("the order q must not be NaN"))
    return nothing
end

# Hill number of order q of the distribution with weights `w > 0` over probabilities `p`:
# `(Σ wᵢ pᵢ^q)^(1/(1−q))`, with the limits `1/max p` at q = Inf and `1/min p` at q = -Inf
# (the weights drop out of both). The power sum is taken in log space (log-sum-exp), so
# that `pᵢ^q` cannot overflow or underflow at large |q| and the result tends to those limits.
function _hill(p::AbstractVector{Float64}, w::AbstractVector{Float64}, q::Real)
    if q == 1
        return exp(-sum(w[i] * p[i] * log(p[i]) for i in eachindex(p)))
    elseif q == Inf
        return 1 / maximum(p)
    elseif q == -Inf
        return 1 / minimum(p)
    else
        x(i) = log(w[i]) + q * log(p[i])
        L = maximum(x(i) for i in eachindex(p))
        return exp((L + log(sum(exp(x(i) - L) for i in eachindex(p)))) / (1 - q))
    end
end

"""
    hill_number(counts, q) -> Float64

Hill number (effective number of types) of order `q` of the distribution `p = counts/sum`:
`q = 0` the number of nonzero entries, `q = 1` `exp(Shannon)`, `q = Inf` `1/max p`,
`q = -Inf` `1/min p` over the nonzero entries (the limits of `q → ±∞`), otherwise
`(Σ pᵢ^q)^(1/(1−q))` over the nonzero entries, evaluated in log space so that a large `|q|`
tends to those limits instead of overflowing. All counts zero (an empty distribution)
gives `0.0` for every `q`. Throws `ArgumentError` for `q = NaN` or negative counts.
"""
function hill_number(counts::AbstractVector{<:Real}, q::Real)
    _check_q(q)
    any(<(0), counts) && throw(ArgumentError("counts must be non-negative"))
    nz = Float64[c for c in counts if c > 0]
    isempty(nz) && return 0.0
    q == 0 && return Float64(length(nz))
    p = nz ./ sum(nz)
    return _hill(p, ones(length(p)), q)
end

"""
    renyi_entropy(counts, q) -> Float64

Rényi entropy `log(hill_number(counts, q))` (`-Inf` for an empty distribution).
"""
renyi_entropy(counts::AbstractVector{<:Real}, q::Real) = log(hill_number(counts, q))

function _check_frequencies(frequencies::Symbol)
    frequencies in (:normalised, :cell_fraction) || throw(ArgumentError(
        "frequencies must be :normalised or :cell_fraction, got :$frequencies"))
    return nothing
end

"""
    sfs_hill_number(sfs, q; frequencies = :normalised) -> Float64

Hill number of order `q` of the mutations in a site-frequency spectrum (`sfs[k]` mutations
carried by `k` cells, `N = length(sfs)`), under one of two definitions of the frequencies:

- `:normalised` (default): `p = k/A` with `A = Σ k·Sₖ`, a proper distribution; `∞D = A/k_max`
  and `₋∞D = A/k_min`.
- `:cell_fraction`: `p = k/N`, as in the v1 papers; the `p` sum to the mean burden, so this
  is not a proper Hill family. `q = 1` is `exp(−Σ Sₖ (k/N) ln(k/N))`, `∞D = N/k_max` and
  `₋∞D = N/k_min`.

`k_max` and `k_min` are the largest and smallest `k` with `Sₖ > 0`; `q = ±Inf` are the limits
of `q → ±∞`, which large finite `|q|` approach (the power sum is evaluated in log space).

`q = 0` is the number of mutations `ΣSₖ` in both. With `m̄ = A/N`: `²D_norm = m̄²·²D_cell`
and `H₁,norm = H₁,cell/m̄ + ln m̄`. An SFS without mutations gives `0.0` for every `q`.
Throws `ArgumentError` for another `frequencies` or `q = NaN`.
"""
function sfs_hill_number(sfs::AbstractVector{<:Integer}, q::Real; frequencies::Symbol = :normalised)
    _check_frequencies(frequencies)
    _check_q(q)
    any(<(0), sfs) && throw(ArgumentError("an SFS must be non-negative"))
    ks = [k for k in eachindex(sfs) if sfs[k] > 0]
    isempty(ks) && return 0.0
    q == 0 && return Float64(sum(sfs))
    w = Float64[sfs[k] for k in ks]
    kmax = maximum(ks)
    denom = frequencies === :normalised ? Float64(sum(k * sfs[k] for k in ks)) : Float64(length(sfs))
    p = Float64[k / denom for k in ks]
    q == Inf && return denom / kmax
    q == -Inf && return denom / minimum(ks)
    return _hill(p, w, q)
end

"""
    sfs_renyi_entropy(sfs, q; frequencies = :normalised) -> Float64

`log(sfs_hill_number(sfs, q; frequencies))`; see [`sfs_hill_number`](@ref).
"""
sfs_renyi_entropy(sfs::AbstractVector{<:Integer}, q::Real; frequencies::Symbol = :normalised) =
    log(sfs_hill_number(sfs, q; frequencies))

# ── Statistic types ───────────────────────────────────────────────────────────

"""
    SFS()

Statistic `:sfs`: [`site_frequency_spectrum`](@ref) (tree mode).
"""
struct SFS <: AbstractStatistic end
measure(::SFS, pop) = site_frequency_spectrum(pop)
statistic_name(::SFS) = :sfs

"""
    MutationsPerCell()

Statistic `:mutations_per_cell`: [`mutations_per_cell`](@ref).
"""
struct MutationsPerCell <: AbstractStatistic end
measure(::MutationsPerCell, pop) = mutations_per_cell(pop)
statistic_name(::MutationsPerCell) = :mutations_per_cell

"""
    MeanBurden()

Statistic `:mean_burden`: [`mean_burden`](@ref).
"""
struct MeanBurden <: AbstractStatistic end
measure(::MeanBurden, pop) = mean_burden(pop)
statistic_name(::MeanBurden) = :mean_burden

"""
    CloneSizes()

Statistic `:clone_sizes`: [`clone_sizes`](@ref).
"""
struct CloneSizes <: AbstractStatistic end
measure(::CloneSizes, pop) = clone_sizes(pop)
statistic_name(::CloneSizes) = :clone_sizes

"""
    MutationFrequencies()

Statistic `:mutation_frequencies`: [`mutation_frequencies`](@ref).
"""
struct MutationFrequencies <: AbstractStatistic end
measure(::MutationFrequencies, pop) = mutation_frequencies(pop)
statistic_name(::MutationFrequencies) = :mutation_frequencies

"""
    HillNumbers(; q = [0, 1, 2, Inf])

Statistic `:hill_numbers`: the Hill numbers of the living clone sizes ([`hill_number`](@ref)),
a `Vector{Float64}` in the order of `q`.
"""
struct HillNumbers <: AbstractStatistic
    q::Vector{Float64}
    function HillNumbers(; q = [0, 1, 2, Inf])
        foreach(_check_q, q)
        return new(Float64[q...])
    end
end
measure(h::HillNumbers, pop) = (c = clone_sizes(pop); Float64[hill_number(c, q) for q in h.q])
statistic_name(::HillNumbers) = :hill_numbers

"""
    SFSHillNumbers(; q = [0, 1, 2, Inf], frequencies = :normalised)

Statistic `:sfs_hill_numbers`: the Hill numbers of the site-frequency spectrum
([`sfs_hill_number`](@ref)), a `Vector{Float64}` in the order of `q` (tree mode).
"""
struct SFSHillNumbers <: AbstractStatistic
    q::Vector{Float64}
    frequencies::Symbol
    function SFSHillNumbers(; q = [0, 1, 2, Inf], frequencies::Symbol = :normalised)
        _check_frequencies(frequencies)
        foreach(_check_q, q)
        return new(Float64[q...], frequencies)
    end
end
measure(h::SFSHillNumbers, pop) =
    (sfs = site_frequency_spectrum(pop);
     Float64[sfs_hill_number(sfs, q; frequencies = h.frequencies) for q in h.q])
statistic_name(::SFSHillNumbers) = :sfs_hill_numbers

"""
    Divisions()

Statistic `:divisions`: [`divisions_per_cell`](@ref).
"""
struct Divisions <: AbstractStatistic end
measure(::Divisions, pop) = divisions_per_cell(pop)
statistic_name(::Divisions) = :divisions

"""
    SensitiveResistant()

Statistic `:sensitive_resistant`: the tuple `(n_sensitive(pop), n_resistant(pop))`.
"""
struct SensitiveResistant <: AbstractStatistic end
measure(::SensitiveResistant, pop) = (n_sensitive(pop), n_resistant(pop))
statistic_name(::SensitiveResistant) = :sensitive_resistant
