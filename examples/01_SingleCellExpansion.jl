##
## Example 1: Single-cell expansion with neutral mutation accumulation
##
## One founder cell grows logistically (birth rate 1, death rate 0.5 + l·N, carrying capacity
## K = 20 000) until the population reaches 10⁴ cells. `restart_on_extinction` retries the
## block if the founder lineage dies out early. Each daughter acquires Poisson(μ = 1)
## neutral mutations; the lineage tree gives the burden and the site-frequency spectrum.
##
## Inputs : none.  Outputs: summary statistics printed to stdout.
## Run    : julia --project=. examples/01_SingleCellExpansion.jl
##

using ClonalEvolutionTx
using Random
using Statistics: mean

rng = Xoshiro(1)

N_target = 10_000
K = 20_000
l = 0.5 / K                       # death rate d₀ + l·N has b = d at N = K

block = BirthDeathBlock(
    birth = s -> 1.0,
    death = s -> 0.5,
    death_density = N -> l * N,   # C(N) = l·N
    μ = 1.0,
    stop_above = N_target,        # the block ends as soon as N ≥ N_target
    restart_on_extinction = true,
    label = "growth",
)

pop = initialize_population(1)    # tree mode, one founder
t_wall = @elapsed simulate!(pop, block, rng)

println("=== Single-cell expansion ===")
println("Population size : ", popsize(pop))
println("Simulation time : ", round(pop.t, digits = 2))
println("Restarts        : ", pop.log[end].restarts)
println("Wall time       : ", round(t_wall, digits = 2), " s")

burden = mutations_per_cell(pop)
println("\nNeutral mutations per cell:")
println("  Mean   : ", round(mean(burden), digits = 2))
println("  Min    : ", minimum(burden))
println("  Max    : ", maximum(burden))
println("  Clonal : ", clonal_mutations(pop), " (shared by all cells)")
println("  Total distinct mutations : ", total_mutational_burden(pop))

sfs = site_frequency_spectrum(pop)
println("\nSFS head (mutations carried by exactly k cells):")
for k in 1:8
    println("  k = $k : ", sfs[k])
end
