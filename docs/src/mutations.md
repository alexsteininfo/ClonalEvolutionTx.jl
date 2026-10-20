# Mutations

When a cell divides, each daughter independently

1. acquires `Poisson(μ)` **neutral** mutations (the block's `μ`; see below), and
2. for each mutation kind in the block's `mutations` list, **in list order**, with that
   kind's probability `ν` creates a new clone whose traits are the traits of its current
   clone changed by the mutation. If several kinds fire in one daughter, the clones nest:
   each firing is its own clone, so clone statistics count every mutation.

Mutation kinds are passed to the block as a tuple or vector. `ν ∈ [0, 1]` is the probability
per daughter per division. All randomness comes from the dynamics stream, so runs are
reproducible from one seed.

## Drivers

[`Drivers`](@ref)`(; ν, Δs, update = :additive)` changes `s` only. `Δs` is a number or a
univariate distribution (drawn per mutation). The update rule is

- `:additive` (default): `s′ = s + Δs`;
- `:multiplicative`: `1 + s′ = (1 + s)(1 + Δs)`;
- or a function `(s, Δs) -> s′`.

```@example mutations
using ClonalEvolutionTx, Distributions, Random

block = BirthDeathBlock(
    birth = s -> 1.0 + s,
    death = s -> 0.5,
    mutations = [Drivers(ν = 0.02, Δs = Exponential(0.1))],
    stop_above = 2_000,
    restart_on_extinction = true,
)
pop = initialize_population(1; mode = :cells)
simulate!(pop, block, Xoshiro(1))

ct = pop.clones
(n_clones(ct), length(living_clones(ct)), round(maximum(ct.s), digits = 3))
```

## Resistance

[`Resistance`](@ref)`(; ν, α, Δs = 0.0, update = :additive)` changes the drug resistance and,
optionally, the fitness.

- `α` is a number in `[0, 1]` or a distribution on `[0, 1]`; a draw outside the interval is
  an `ArgumentError`.
- The child's resistance is `min(α_drawn, α_parent)`: cells only become more resistant.
- `Δs` is the **cost or benefit** of resistance, of either sign. It is a number, a
  distribution or a function of the drawn `α` (for example `α -> -0.1 * (1 - α)`), applied
  with the same `update` rule as drivers, and it is applied even if the min rule leaves `α`
  unchanged.

A cost of resistance is a negative `Δs` together with a fitness-dependent birth rate:

```@example mutations
costly = BirthDeathBlock(
    birth = s -> 1.0 + s,
    death = s -> 0.5,
    response = LinearResponse(4.0),
    mutations = [Resistance(ν = 0.01, α = Uniform(0.0, 0.9), Δs = α -> -0.1 * (1 - α))],
    stop_above = 2_000,
    restart_on_extinction = true,
)
pop = initialize_population(1; mode = :cells)
simulate!(pop, costly, Xoshiro(2))
resistant = [k for k in living_clones(pop.clones) if pop.clones.α[k] < 1]
(length(resistant), n_resistant(pop) > 0)
```

## Custom mutations

[`CustomMutation`](@ref)`(; ν, f)` with `f(rng, traits) -> Traits` covers anything else. It
must draw only from the `rng` it is given.

```@example mutations
flip = CustomMutation(ν = 0.01, f = (rng, t) -> Traits(t.s + 0.2 * rand(rng), t.α))
custom = BirthDeathBlock(birth = s -> 1.0 + s, death = s -> 0.5, mutations = [flip],
                         stop_above = 500, restart_on_extinction = true)
pop = initialize_population(1; mode = :cells)
simulate!(pop, custom, Xoshiro(3))
n_clones(pop.clones)
```

[`mutate`](@ref)`(m, parent_traits, rng)` applies one mutation kind to a [`Traits`](@ref) value
directly, which is handy for testing a custom rule, and [`mutation_kind`](@ref) names the kind
of a mutation (`:driver`, `:resistance`, `:custom`) as it is stored in the clone table.

## Neutral mutations

Neutral mutations follow an infinite-sites model on the **division clock**: each daughter
acquires `Poisson(μ)` new mutations at each division. They do not affect the dynamics. They
are drawn from the separate neutral stream (see [Concepts](@ref)), so the same seed gives the
same living cells for any `μ`.

- In **tree mode** each node stores the mutations acquired at that cell's birth and the running
  total along the path from its root.
- In **cell mode** each cell stores its running total.

After a tree-mode run, [`repaint_neutral!`](@ref)`(pop, μ′, rng)` redraws the neutral
mutations with a new rate (new `Poisson(μ′)` per node), so a sweep over `μ` can reuse one tree.
`initialize_population(...; clonal_mutations = c)` adds `c` mutations carried by every cell,
counted at `k = N` in the site frequency spectrum, identically in both modes. There are no
per-founder mutations: a clone specification sets sizes and traits only (see
[Limitations and open questions](@ref)).

```@example mutations
tree = initialize_population(1; mode = :tree)
simulate!(tree, BirthDeathBlock(birth = s -> 1.0, death = s -> 0.5, μ = 1.0,
                                stop_above = 500, restart_on_extinction = true), Xoshiro(4))
before = mean_burden(tree)
repaint_neutral!(tree, 3.0, Xoshiro(5))
(round(before, digits = 1), round(mean_burden(tree), digits = 1))   # larger after μ′ = 3
```
