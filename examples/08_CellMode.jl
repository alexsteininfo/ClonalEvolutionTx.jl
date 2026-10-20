##
## Example 8: Cell mode versus tree mode
##
## Tree mode keeps the whole lineage tree (needed for the SFS, genealogy and Newick
## export); cell mode keeps only the living cells and is much lighter. Neutral mutations
## are drawn from a stream separate from the dynamics, so one seed gives identical living
## cells in both modes. First this is checked at 10⁴ cells, then cell mode is timed on a
## 10⁶-cell expansion.
##
## Inputs : none.  Outputs: comparison and timings printed to stdout.
## Run    : julia --project=. examples/08_CellMode.jl
##

using ClonalEvolutionTx
using Random
using Distributions
using Statistics: mean

block(N; kw...) = BirthDeathBlock(
    birth = s -> 1.0 + s,
    death = s -> 0.5,
    mutations = [Drivers(ν = 0.001, Δs = Exponential(0.05))],
    μ = 1.0,
    stop_above = N,
    restart_on_extinction = true,
    kw...,
)

println("=== Cell mode ===")
seed = 8
tree = initialize_population(1; mode = :tree)
cells = initialize_population(1; mode = :cells)
t_tree = @elapsed simulate!(tree, block(10_000), Xoshiro(seed))
t_cells = @elapsed simulate!(cells, block(10_000), Xoshiro(seed))
println("10⁴ cells, seed ", seed, ":")
println("  tree mode : ", round(t_tree, digits = 2), " s,  cell mode : ", round(t_cells, digits = 2), " s")
println("  same cell ids    : ", cell_ids(tree) == cell_ids(cells))
println("  same clones      : ", cell_clones(tree) == cell_clones(cells))
println("  same burdens     : ", mutations_per_cell(tree) == mutations_per_cell(cells))
println("  same final time  : ", tree.t == cells.t)

big = initialize_population(1; mode = :cells)
t_big = @elapsed simulate!(big, block(1_000_000), Xoshiro(seed))
println("\nCell mode at 10⁶ cells:")
println("  N = ", popsize(big), ", t = ", round(big.t, digits = 2), ", wall ", round(t_big, digits = 2), " s")
println("  mean neutral burden = ", round(mean_burden(big), digits = 1),
        ", living clones = ", length(living_clones(big.clones)))
println("  events: ", big.log[end].accepted, " accepted of ", big.log[end].ticks, " ticks")
println("  (the SFS and the genealogy need tree mode: ",
        try site_frequency_spectrum(big); "available" catch e; e isa ArgumentError ? "throws ArgumentError" : "error" end, ")")
