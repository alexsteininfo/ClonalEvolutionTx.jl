##
## Example 3: Maximum tolerated dose against a growing tumour (manuscript_v1 parameters)
##
## Birth b = 1, death d₀ = 0.5, linear drug response d_S = 4, so a fully sensitive cell has
## net growth 1 - 0.5 - 4 < 0 under full dose. Resistance arises with probability
## ν = 0.002 per daughter and sets α = 0.8 (a resistant cell still feels 80 % of the drug,
## so MTD is hard on it too). The tumour grows to 10⁴ cells, then `MTD` is applied until
## progression (N ≥ 2·10⁴), extinction, or t_max = 40. The grown tumour is copied, so
## one growth phase serves several treatment replicates.
##
## Inputs : none.  Outputs: outcome and treatment log printed to stdout.
## Run    : julia --project=. examples/03_MaximumToleratedDose.jl
##

using ClonalEvolutionTx
using Random

rng = Xoshiro(3)

base = BirthDeathBlock(
    birth = s -> 1.0,
    death = s -> 0.5,
    response = LinearResponse(4.0),            # D(m) = d_S·m
    mutations = [Resistance(ν = 0.002, α = 0.8)],
    label = "growth",
)

growth = remake(base; stop_above = 10_000, restart_on_extinction = true)
tumour = initialize_population(1; mode = :cells)
simulate!(tumour, growth, rng)
println("=== Maximum tolerated dose ===")
println("Tumour before treatment: N = ", popsize(tumour), ", t = ", round(tumour.t, digits = 2),
        ", resistant cells = ", n_resistant(tumour))

mtd = MTD(base; N_prog = 20_000, t_max = 40.0)

pop = copy(tumour)
t_wall = @elapsed simulate!(pop, mtd, rng)
println("\nOutcome of one MTD run: ", schedule_outcome(pop), "   (N = ", popsize(pop),
        ", t = ", round(pop.t, digits = 2), ", wall ", round(t_wall, digits = 2), " s)")
println("Log:")
for b in pop.log
    println("  ", rpad(b.label, 8), "t ", round(b.t_start, digits = 2), " -> ", round(b.t_end, digits = 2),
            "   ended: ", b.end_reason, "   dose: ", b.dose, "   m_end = ", b.m_end)
end
println("Resistant cells at the end: ", n_resistant(pop), " of ", popsize(pop))

println("\nTwenty replicates from the same grown tumour:")
outcomes = Symbol[]
for r in 1:20
    p = copy(tumour)
    simulate!(p, mtd, Xoshiro(100 + r))
    push!(outcomes, schedule_outcome(p))
end
for o in unique(outcomes)
    println("  ", rpad(o, 12), count(==(o), outcomes))
end
println("Resistant clones here have α = 0.8, so full dose still kills them (net growth 1 - 0.5 - 3.2 < 0):\nMTD cures this tumour. Lower d_S or raise α for a regime where MTD selects resistance.")
println("(the original tumour is untouched: N = ", popsize(tumour), ")")
