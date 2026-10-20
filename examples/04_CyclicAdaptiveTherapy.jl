##
## Example 4: Cyclic adaptive therapy
##
## The same model and growth phase as example 3 (b = 1, d₀ = 0.5, d_S = 4, resistance
## ν = 0.002 with α = 0.8), then `AdaptiveOnOff`: full dose until N < N_low = 3000, no drug
## until N ≥ N_up = 10 000, repeated until progression (N ≥ 20 000) or t_max = 14.5 (about six
## cycles). An `AtBlockEnd` trigger records, at every toggle, the population size and
## the number of resistant clones (living clones with α < 1).
##
## Inputs : none.  Outputs: per-cycle table printed to stdout.
## Run    : julia --project=. examples/04_CyclicAdaptiveTherapy.jl
##

using ClonalEvolutionTx
using Random

# number of living clones with α < 1 (a named function, so it can serve as a statistic)
resistant_clones(pop) = count(k -> pop.clones.α[k] < 1.0, living_clones(pop.clones))

rng = Xoshiro(4)

base = BirthDeathBlock(
    birth = s -> 1.0,
    death = s -> 0.5,
    response = LinearResponse(4.0),
    mutations = [Resistance(ν = 0.002, α = 0.8)],
    label = "growth",
)

pop = initialize_population(1; mode = :cells)
simulate!(pop, remake(base; stop_above = 10_000, restart_on_extinction = true), rng)
println("=== Cyclic adaptive therapy ===")
println("Tumour before treatment: N = ", popsize(pop), ", t = ", round(pop.t, digits = 2))

at = AdaptiveOnOff(base; N_low = 3_000, N_up = 10_000, N_prog = 20_000, t_max = 14.5)
spec = MeasurementSpec(triggers = [AtBlockEnd()],
                       stats = [popsize, resistant_clones, n_resistant])
acc = MeasurementAccumulator(spec)
t_wall = @elapsed simulate!(pop, at, rng; accumulator = acc)
meas = finalize_measurements(acc)

println("\nOutcome: ", schedule_outcome(pop), "  (t = ", round(pop.t, digits = 2),
        ", wall ", round(t_wall, digits = 2), " s)")
println("\nAt each toggle (end of a block):")
println("  ", rpad("block", 8), rpad("time", 9), rpad("N", 8), rpad("resistant clones", 18), "resistant cells")
for sn in meas.snapshots
    println("  ", rpad(sn.label, 8), rpad(round(sn.t, digits = 2), 9), rpad(sn.stats.popsize, 8),
            rpad(sn.stats.resistant_clones, 18), sn.stats.n_resistant)
end
n_cycles = count(b -> b.label == "AT on", pop.log)
println("\nCycles (AT on blocks): ", n_cycles)
