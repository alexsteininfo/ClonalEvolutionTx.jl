# Benchmark of ClonalEvolutionTx: growth to N cells in tree and cell mode (neutral, and with
# resistance and neutral mutations), treatment schedules (MTD, cyclic adaptive therapy,
# stabilising therapy), statistics on one 10^6 tree population, and a parallel throughput mode.
#
# Inputs : population sizes on the command line (optional), or --parallel N.
# Outputs: timing and memory tables printed to stdout; nothing is written to disk.
# Run:
#   julia --project=. benchmark/benchmark.jl                     # N = 10^3, 10^4, 10^5, 10^6, 10^7
#   julia --project=. benchmark/benchmark.jl 1000 10000          # chosen sizes
#   julia -t 8 --project=. benchmark/benchmark.jl --parallel 100000   # throughput at N = 10^5
#
# Growth modes (birth b = 1, death d = 0.5, one seeded RNG per run from derive_seed):
#   neutral/tree, neutral/cells   μ = 0, no mutations
#   resist/tree,  resist/cells    Resistance(ν = 0.002, α = 0.8), μ = 1 neutral mutations
# Treatment modes (linear response d_S = 4, resistance as above; the tumour is grown first
# in cell mode, untimed):
#   mtd α=0, mtd α=0.8   MTD from NTREAT cells to progression (2·NTREAT), extinction or t = 40
#   cyclic               exactly 6 on/off cycles (a FunctionSchedule of the two blocks of
#                        AdaptiveOnOff, N_low = N/2, N_up = N; ends early on extinction/progression)
#   stabilising          Stabilising, N_stable = N/2, reports ns per event (tick)
# Statistics are timed on one tree population of min(10^6, largest size) cells.
#
# Sizes run in increasing order. Each run measures the live memory per cell (the GC's
# live-bytes counter, after GC.gc()); a larger run whose estimated peak memory
# (1.5 × bytes per cell × N) exceeds 60% of the free RAM, read again before every size, is
# skipped. The first size of a mode uses a stated default estimate (DEFAULT_BYTES), later
# sizes the measured value. The 10^6 treatment and statistics setups and --parallel
# (nthreads concurrent runs) are guarded in the same way.
# --parallel N: 4 × nthreads runs to N cells, first serially, then with Threads.@threads.
# ns/cell is wall time divided by the final N (more cells are created than survive, so the
# cost per created cell is lower); ticks/N is thinning ticks (accepted or not) per final cell.
using ClonalEvolutionTx, Random, Printf
using Distributions: Exponential

const PARALLEL = "--parallel" in ARGS
const PARALLEL_N = let i = findfirst(==("--parallel"), ARGS)
    isnothing(i) ? 0 :
    i < length(ARGS) ? parse(Int, ARGS[i + 1]) :
    error("--parallel needs a population size, e.g. --parallel 100000")
end
const SIZES = sort(PARALLEL ? [PARALLEL_N] :
                   length(ARGS) >= 1 ? [parse(Int, a) for a in ARGS] :
                   [10^3, 10^4, 10^5, 10^6, 10^7])
const PEAK_FACTOR = 1.5
# conservative bytes per living cell by storage mode (measured here: tree ≈ 175, cells ≈ 47)
const DEFAULT_BYTES = Dict(:tree => 200.0, :cells => 60.0)
memory_limit() = 0.6 * Sys.free_memory()          # re-read before every guarded run

# true if `copies` runs of `N` cells at `bytes` per cell fit; otherwise prints a skip row
function fits(name, N, bytes; copies = 1)
    est = PEAK_FACTOR * bytes * N * copies
    lim = memory_limit()
    est <= lim && return true
    @printf("%-14s %10d  skipped (estimated %.2f GB > 60%% of free %.2f GB)\n", name, N,
            gb(est), gb(lim / 0.6))
    return false
end
const BASE_SEED = 20260101

gb(x) = x / 2^30

neutral_block(N) = BirthDeathBlock(birth = s -> 1.0, death = s -> 0.5, stop_above = N,
                                   restart_on_extinction = true)
resist_block(N) = BirthDeathBlock(birth = s -> 1.0, death = s -> 0.5,
                                  mutations = [Resistance(ν = 0.002, α = 0.8)], μ = 1.0,
                                  stop_above = N, restart_on_extinction = true)

# mode name => (block constructor, storage mode)
const MODES = [
    "neutral/tree"  => (neutral_block, :tree),
    "neutral/cells" => (neutral_block, :cells),
    "resist/tree"   => (resist_block, :tree),
    "resist/cells"  => (resist_block, :cells),
]

treat_base(α) = BirthDeathBlock(birth = s -> 1.0, death = s -> 0.5,
                                response = LinearResponse(4.0),
                                mutations = [Resistance(ν = 0.002, α = α)],
                                label = "growth")
function grown(N, α, seed)
    pop = initialize_population(1; mode = :cells)
    simulate!(pop, remake(treat_base(α); stop_above = N, restart_on_extinction = true),
              Xoshiro(seed))
    return pop
end

# warm-up with the same block types as the measured runs (compilation is not timed)
for (_, (mk, m)) in MODES
    simulate!(initialize_population(1; mode = m), mk(200), Xoshiro(0))
end
let w = grown(500, 0.8, 0)
    simulate!(copy(w), MTD(treat_base(0.8); N_prog = 1000, t_max = 40.0), Xoshiro(0))
    simulate!(copy(w), AdaptiveOnOff(treat_base(0.8); N_low = 250, N_up = 500, N_prog = 5000,
                                     t_max = 10.0), Xoshiro(0))
    simulate!(copy(w), Stabilising(treat_base(0.8); N_stable = 250, N_prog = 1000, t_max = 5.0),
              Xoshiro(0))
end

# One measured growth run, inside a function so that nothing survives it.
function measure(mk, mode, N, seed)
    GC.gc()
    baseline = Base.gc_live_bytes()
    pop = initialize_population(1; mode = mode)
    stats = @timed simulate!(pop, mk(N), Xoshiro(seed))
    GC.gc()
    live = Base.gc_live_bytes() - baseline
    return (time = stats.time, gctime = stats.gctime, alloc = stats.bytes, live = live,
            n = popsize(pop), ticks = sum(b.ticks for b in pop.log),
            accepted = sum(b.accepted for b in pop.log))
end

if PARALLEL
    run_one(N, seed) = popsize(simulate!(initialize_population(1; mode = :cells),
                                         resist_block(N), Xoshiro(seed)))
    nt = Threads.nthreads()
    nruns = 4 * nt
    seeds = [derive_seed(BASE_SEED, r) for r in 1:nruns]
    fits("parallel", PARALLEL_N, DEFAULT_BYTES[:cells]; copies = nt) || exit(0)
    run_one(PARALLEL_N, 0)
    GC.gc()
    tserial = @elapsed for r in 1:nruns
        run_one(PARALLEL_N, seeds[r])
    end
    GC.gc()
    tpar = @elapsed Threads.@threads for r in 1:nruns
        run_one(PARALLEL_N, seeds[r])
    end
    @printf("parallel throughput (resist/cells): N = %d, %d runs, %d threads (CPU threads: %d)\n",
            PARALLEL_N, nruns, nt, Sys.CPU_THREADS)
    @printf("  serial:   %8.2f s  -> %8.2f runs/s\n", tserial, nruns / tserial)
    @printf("  parallel: %8.2f s  -> %8.2f runs/s  (speedup %.1fx on %d threads, efficiency %.0f%%)\n",
            tpar, nruns / tpar, tserial / tpar, nt, 100 * tserial / tpar / nt)
    nt == 1 && println("  (start Julia with -t N to use N threads)")
else

@printf("Julia %s, %s, %d CPU threads\n", VERSION, Sys.cpu_info()[1].model, Sys.CPU_THREADS)
@printf("RAM: %.1f GB total, %.1f GB free, limit per run: %.1f GB\n\n",
        gb(Sys.total_memory()), gb(Sys.free_memory()), gb(memory_limit()))
@printf("%-14s %10s %10s %7s %11s %11s %9s %9s\n", "mode", "N", "time (s)", "GC (%)",
        "alloc (GB)", "live (MB)", "ns/cell", "ticks/N")

for (mode, (mk, storage)) in MODES
    bytes_per_cell = DEFAULT_BYTES[storage]
    for N in SIZES
        fits(mode, N, bytes_per_cell) || continue
        r = measure(mk, storage, N, derive_seed(BASE_SEED, N))
        bytes_per_cell = r.live / r.n
        @printf("%-14s %10d %10.3f %7.1f %11.3f %11.3f %9.0f %9.2f\n", mode, N, r.time,
                100 * r.gctime / r.time, gb(r.alloc), r.live / 2^20, 1e9 * r.time / r.n,
                r.ticks / r.n)
    end
end

# treatment ----------------------------------------------------------------------------
const NTREAT = min(10^6, maximum(SIZES))
const NCYC = min(10^5, maximum(SIZES))
println("\ntreatment (tumour grown in cell mode, untimed; linear response d_S = 4):")
@printf("%-26s %10s %10s %12s %12s %10s\n", "scenario", "N0", "time (s)", "ticks",
        "accepted", "accepted/ticks")
function report(name, N0, pop, t)
    ticks = sum(b.ticks for b in pop.log)
    acc = sum(b.accepted for b in pop.log)
    @printf("%-26s %10d %10.3f %12d %12d %10.3f   outcome %s\n", name, N0, t, ticks, acc,
            acc / max(ticks, 1), schedule_outcome(pop))
    return ticks, acc
end
for α in (0.0, 0.8)
    # the grown tumour and its copy are alive together
    fits("mtd α=$α", NTREAT, DEFAULT_BYTES[:cells]; copies = 2) || continue
    tumour = grown(NTREAT, α, derive_seed(BASE_SEED, 1, round(Int, 10α)))
    pop = copy(tumour)
    n0 = length(pop.log)
    t = @elapsed simulate!(pop, MTD(treat_base(α); N_prog = 2NTREAT, t_max = 40.0),
                           Xoshiro(derive_seed(BASE_SEED, 2, round(Int, 10α))))
    # count only the treatment blocks (the log also holds the growth phase)
    treat_log = pop.log[(n0 + 1):end]
    ticks = sum(b.ticks for b in treat_log); acc = sum(b.accepted for b in treat_log)
    @printf("%-26s %10d %10.3f %12d %12d %10.3f   outcome %s\n", "mtd α=$α", NTREAT, t, ticks,
            acc, acc / max(ticks, 1), schedule_outcome(pop))
end
# Exactly 6 complete on/off cycles: the two blocks of AdaptiveOnOff (on until N ≤ N_low, off
# until N ≥ N_up), 12 blocks, ending early only on extinction or progression.
function six_cycles(base, N_low, N_up, N_prog)
    return FunctionSchedule() do pop, i, t0
        i >= 12 && return nothing
        i >= 1 && pop.log[end].end_reason in (:extinct, :stopped) && return nothing
        i >= 1 && pop.log[end].end_reason === :above && popsize(pop) >= N_prog && return nothing
        if iseven(i)
            return remake(base; dose = FixedDose(1.0), label = "AT on", tmax = Inf,
                          stop_above = N_prog, stop_below = N_low, restart_on_extinction = false)
        else
            return remake(base; dose = FixedDose(0.0), label = "AT off", tmax = Inf,
                          stop_above = N_up, stop_below = -1, restart_on_extinction = false)
        end
    end
end
if fits("cyclic", NCYC, DEFAULT_BYTES[:cells]; copies = 2)
  let tumour = grown(NCYC, 0.8, derive_seed(BASE_SEED, 3)), pop = copy(tumour), n0 = length(tumour.log)
    t = @elapsed simulate!(pop, six_cycles(treat_base(0.8), NCYC ÷ 2, NCYC, 10NCYC),
                           Xoshiro(derive_seed(BASE_SEED, 4)))
    tl = pop.log[(n0 + 1):end]
    ticks = sum(b.ticks for b in tl); acc = sum(b.accepted for b in tl)
    @printf("%-26s %10d %10.3f %12d %12d %10.3f   ended after %d blocks, %d complete cycles\n", "cyclic", NCYC, t,
            ticks, acc, acc / max(ticks, 1), length(tl),
            count(b -> b.label == "AT off", tl))
  end
end
if fits("stabilising", NCYC, DEFAULT_BYTES[:cells]; copies = 2)
  let tumour = grown(NCYC, 0.8, derive_seed(BASE_SEED, 5)), pop = copy(tumour), n0 = length(tumour.log)
    t = @elapsed simulate!(pop, Stabilising(treat_base(0.8); N_stable = NCYC ÷ 2,
                                            N_prog = 2NCYC, t_max = 40.0),
                           Xoshiro(derive_seed(BASE_SEED, 6)))
    tl = pop.log[(n0 + 1):end]
    ticks = sum(b.ticks for b in tl); acc = sum(b.accepted for b in tl)
    @printf("%-26s %10d %10.3f %12d %12d %10.3f   outcome %s, %.0f ns per tick, %.0f ns per accepted event\n",
            "stabilising", NCYC, t, ticks, acc, acc / max(ticks, 1), schedule_outcome(pop),
            1e9 * t / max(ticks, 1), 1e9 * t / max(acc, 1))
  end
end

# statistics ---------------------------------------------------------------------------
function time_statistics(pop)
    sfs = site_frequency_spectrum(pop)
    tasks = Pair{String, Function}[
        "SFS"                    => () -> site_frequency_spectrum(pop),
        "mutations_per_cell"     => () -> mutations_per_cell(pop),
        "mutation_frequencies"   => () -> mutation_frequencies(pop),
        "sfs_hill_number q=0"    => () -> sfs_hill_number(sfs, 0),
        "sfs_hill_number q=1"    => () -> sfs_hill_number(sfs, 1),
        "sfs_hill_number q=2"    => () -> sfs_hill_number(sfs, 2),
        "sfs_hill_number q=Inf"  => () -> sfs_hill_number(sfs, Inf),
        "sample_population 1000" => () -> sample_population(pop, 1000; seed = 1),
    ]
    results = Pair{String, Float64}[]
    for (name, f) in tasks
        f()
        GC.gc()
        push!(results, name => @elapsed f())
    end
    return results
end
GC.gc()
const NSTAT = min(10^6, maximum(SIZES))
if fits("statistics", NSTAT, DEFAULT_BYTES[:tree])
  let pop = initialize_population(1; mode = :tree)
    simulate!(pop, resist_block(NSTAT), Xoshiro(derive_seed(BASE_SEED, 7)))
    println("\nstatistics on one tree population of N = $NSTAT (resist, μ = 1):")
    for (name, t) in time_statistics(pop)
        @printf("  %-24s %9.4f s  (%6.0f ns per cell)\n", name, t, 1e9 * t / NSTAT)
    end
end
end
end
