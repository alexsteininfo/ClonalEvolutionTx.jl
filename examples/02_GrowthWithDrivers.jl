##
## Example 2: Growth with driver mutations
##
## Exponential growth (b = 1, d = 0.5) from one cell to 10⁴ cells. Each daughter acquires a
## driver with probability ν = 0.01; a driver adds Δs ~ Exponential(0.1) to the clone's
## selection trait s, and the birth rate is β(s) = 1 + s, so drivers make clones fitter.
## Neutral mutations (μ = 1) accumulate in parallel. We report the clone structure and the
## fittest clones.
##
## Inputs : none.  Outputs: summary statistics printed to stdout.
## Run    : julia --project=. examples/02_GrowthWithDrivers.jl
##

using ClonalEvolutionTx
using Random
using Distributions
using Statistics: mean

rng = Xoshiro(2)

block = BirthDeathBlock(
    birth = s -> 1.0 + s,
    death = s -> 0.5,
    mutations = [Drivers(ν = 0.01, Δs = Exponential(0.1))],
    μ = 1.0,
    stop_above = 10_000,
    restart_on_extinction = true,
    label = "growth with drivers",
)

pop = initialize_population(1)
t_wall = @elapsed simulate!(pop, block, rng)

ct = pop.clones
alive = living_clones(ct)
sizes = clone_sizes(pop)

println("=== Growth with drivers ===")
println("Population size   : ", popsize(pop), "   (t = ", round(pop.t, digits = 2), ")")
println("Clones created    : ", n_clones(ct))
println("Clones alive      : ", length(alive))
println("Largest clone     : ", maximum(sizes), " cells")
println("Wall time         : ", round(t_wall, digits = 2), " s")

println("\nTen largest clones (id, cells, s):")
for j in first(sortperm(sizes, rev = true), 10)
    k = alive[j]
    println("  clone ", lpad(k, 5), "  ", lpad(ct.size[k], 6), " cells   s = ", round(ct.s[k], digits = 3))
end

println("\nFittest living clones (id, s, cells):")
for j in first(sortperm(ct.s[alive], rev = true), 5)
    k = alive[j]
    println("  clone ", lpad(k, 5), "  s = ", round(ct.s[k], digits = 3), "   ", ct.size[k], " cells")
end

freqs = mutation_frequencies(pop)
println("\nDriver carriers (cells carrying each driver): max ", maximum(freqs),
        ", median ", round(Int, median(freqs)))
println("Mean neutral burden: ", round(mean_burden(pop), digits = 2))
