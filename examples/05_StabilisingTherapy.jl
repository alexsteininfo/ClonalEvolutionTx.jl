##
## Example 5: Stabilising therapy
##
## Same model and growth as example 3. `Stabilising` first applies the full dose until
## N falls below N_stable = 5000, then holds the tumour with the dose m(t) that makes
## the expected net growth of the drug-sensitive part cancel the growth of the resistant
## part (dose rises as resistance spreads, capped at m_max = 1). A trajectory recorded
## every 1.0 time units shows N and m(t).
##
## Inputs : none.  Outputs: trajectory table printed to stdout.
## Run    : julia --project=. examples/05_StabilisingTherapy.jl
##

using ClonalEvolutionTx
using Random

rng = Xoshiro(5)

base = BirthDeathBlock(
    birth = s -> 1.0,
    death = s -> 0.5,
    response = LinearResponse(4.0),
    mutations = [Resistance(ν = 0.002, α = 0.8)],
    label = "growth",
)

pop = initialize_population(1; mode = :cells)
simulate!(pop, remake(base; stop_above = 10_000, restart_on_extinction = true), rng)
println("=== Stabilising therapy ===")
println("Tumour before treatment: N = ", popsize(pop), ", t = ", round(pop.t, digits = 2))

stab = Stabilising(base; N_stable = 5_000, N_prog = 20_000, t_max = 40.0)
acc = MeasurementAccumulator(MeasurementSpec(trajectory_dt = 1.0))
t_wall = @elapsed simulate!(pop, stab, rng; accumulator = acc)
traj = finalize_measurements(acc).trajectory

println("\nOutcome: ", schedule_outcome(pop), "  (t = ", round(pop.t, digits = 2),
        ", wall ", round(t_wall, digits = 2), " s)")
println("Phases: ", join(("$(b.label) [$(round(b.t_start, digits = 1)), $(round(b.t_end, digits = 1))]" for b in pop.log[2:end]), "; "))

println("\n  ", rpad("t", 8), rpad("N", 8), rpad("m(t)", 9), "resistant cells")
for p in traj[1:max(1, length(traj) ÷ 12):end]       # about a dozen rows
    println("  ", rpad(round(p.t, digits = 1), 8), rpad(p.N, 8), rpad(round(p.m, digits = 3), 9), p.n_resistant)
end
