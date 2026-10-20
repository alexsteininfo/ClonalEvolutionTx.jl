##
## Example 7: Measurements, sampling and storage
##
## A `MeasurementSpec` is attached to a simulation through a `MeasurementAccumulator`. It
## records a trajectory on a time grid and snapshots at `AtTime`, `AtSize` and
## `AtBlockEnd` triggers. Each snapshot evaluates statistics (the SFS, mean burden, clone
## sizes) on the whole population and on random samples of 100 and 1000 cells. At the
## end, the population is written to a temporary file with `save_population` and read back
## with `load_population`.
##
## Inputs : none.  Outputs: summary printed to stdout; one file in `tempdir()` (removed).
## Run    : julia --project=. examples/07_MeasurementsAndSampling.jl
##

using ClonalEvolutionTx
using Random
using Distributions
using Statistics: mean

rng = Xoshiro(7)

block = BirthDeathBlock(
    birth = s -> 1.0 + s,
    death = s -> 0.5,
    mutations = [Drivers(ν = 0.01, Δs = Exponential(0.05))],
    μ = 1.0,
    stop_above = 5_000,
    restart_on_extinction = true,
    label = "growth",
)

spec = MeasurementSpec(
    trajectory_dt = 3.0,
    triggers = [AtTime(15.0), AtSize(2_000), AtBlockEnd()],
    stats = [SFS(), MeanBurden(), CloneSizes(), popsize],
    samples = [100, 1000],
    sample_seed = 7,
)
acc = MeasurementAccumulator(spec)
pop = initialize_population(1)
simulate!(pop, block, rng; accumulator = acc)
meas = finalize_measurements(acc)

println("=== Measurements and sampling ===")
println("Final N = ", popsize(pop), " at t = ", round(pop.t, digits = 2))

println("\nTrajectory (every 3 time units):")
for p in meas.trajectory
    println("  t = ", rpad(p.t, 6), "N = ", p.N)
end

println("\nSnapshots:")
for sn in meas.snapshots
    println("  ", rpad(string(sn.trigger), 26), "t = ", rpad(round(sn.t, digits = 2), 7), "N = ", rpad(sn.N, 6),
            "mean burden = ", rpad(round(sn.stats.mean_burden, digits = 1), 7),
            "singletons = ", sn.stats.sfs[1], "   clones = ", length(sn.stats.clone_sizes))
    for s in sn.samples
        println("      sample of ", rpad(s.n, 5), "mean burden = ", rpad(round(s.stats.mean_burden, digits = 1), 7),
                "singletons = ", s.stats.sfs[1], "   clones = ", length(s.stats.clone_sizes))
    end
end

println("\nDirect sampling of the final population:")
sub = sample_population(pop, 200; seed = 1)
println("  sample of ", popsize(sub), " cells, mean burden ", round(mean_burden(sub), digits = 1),
        " (population: ", round(mean_burden(pop), digits = 1), ")")

path = joinpath(tempdir(), "clonalevolutiontx_example07.cepop")
try
    save_population(path, pop)
    back = load_population(path)
    println("\nSave/load round trip (", filesize(path), " bytes):")
    println("  same cells      : ", cell_ids(back) == cell_ids(pop))
    println("  same SFS        : ", site_frequency_spectrum(back) == site_frequency_spectrum(pop))
    println("  same clone sizes: ", sort(clone_sizes(back)) == sort(clone_sizes(pop)), "  (living-clone order may differ)")
    println("  same time       : ", back.t == pop.t)
finally
    rm(path; force = true)
end
