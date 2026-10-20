##
## Example 6: Cost of resistance under adaptive therapy
##
## Resistance arises with probability ν = 0.002 per daughter and draws α ~ Uniform(0, 1)
## (α is the fraction of the drug a cell feels; smaller means more resistant). Two
## scenarios share the same model otherwise (b = β(s) = 1 + s, d₀ = 0.5, d_S = 4):
##   - no cost : resistance leaves s unchanged;
##   - cost    : resistance changes s by Δs = -0.1 (1 - α), so more resistant cells divide slower.
## Each grown tumour (10⁴ cells) is treated with `AdaptiveOnOff` (N_low = 3000,
## N_up = 10 000, N_prog = 20 000, t_max = 60). We compare the time to progression over
## 30 replicates each (the replicate's growth phase is part of the replicate).
##
## Inputs : none.  Outputs: summary table printed to stdout.
## Run    : julia --project=. examples/06_CostOfResistance.jl
##

using ClonalEvolutionTx
using Random
using Distributions
using Statistics: mean, median

function replicate(seed; cost)
    rng = Xoshiro(seed)
    Δs = cost ? (α -> -0.1 * (1 - α)) : 0.0
    base = BirthDeathBlock(
        birth = s -> 1.0 + s,
        death = s -> 0.5,
        response = LinearResponse(4.0),
        mutations = [Resistance(ν = 0.002, α = Uniform(0, 1), Δs = Δs)],
    )
    pop = initialize_population(1; mode = :cells)
    simulate!(pop, remake(base; stop_above = 10_000, restart_on_extinction = true, label = "growth"), rng)
    t0 = pop.t
    simulate!(pop, AdaptiveOnOff(base; N_low = 3_000, N_up = 10_000, N_prog = 20_000, t_max = 60.0), rng)
    return (outcome = schedule_outcome(pop), duration = pop.t - t0)
end

println("=== Cost of resistance under adaptive therapy ===")
nrep = 30
t_wall = @elapsed begin
    results = Dict(c => [replicate(1000 * i + r; cost = c) for r in 1:nrep] for (i, c) in enumerate((false, true)))
end

println("Wall time: ", round(t_wall, digits = 2), " s\n")
println("  ", rpad("scenario", 10), rpad("progressed", 12), rpad("extinct", 9), rpad("time-out", 10), "time to progression (mean / median)")
for c in (false, true)
    rs = results[c]
    prog = [r.duration for r in rs if r.outcome == :progressed]
    ttp = isempty(prog) ? "n/a" : string(round(mean(prog), digits = 1), " / ", round(median(prog), digits = 1))
    println("  ", rpad(c ? "cost" : "no cost", 10),
            rpad(length(prog), 12),
            rpad(count(r -> r.outcome == :extinct, rs), 9),
            rpad(count(r -> r.outcome == :tmax, rs), 10), ttp)
end
println("\nReplicates that did not progress ran until extinction or t_max = 60.")
