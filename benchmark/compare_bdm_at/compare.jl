# Statistical and speed comparison of ClonalEvolutionTx with its predecessor
# BirthDeathMutation_AT.jl (the reference implementation of the manuscript_v1 model): the same
# v1 model and treatment schedules are run in both packages, and the distributions of a fixed
# set of statistics are compared by their means and two-sample Kolmogorov–Smirnov tests.
#
# Inputs : optional number of runs per package and configuration (default 200).
# Outputs: one table per configuration and a timing table, printed to stdout; nothing is
#          written to disk.
# Run    : julia --project=benchmark/compare_bdm_at benchmark/compare_bdm_at/compare.jl [runs]
#          (the first instantiation of this environment downloads and precompiles the
#          predecessor's heavy dependencies, e.g. CairoMakie; that takes several minutes)
#
# ── The manuscript_v1 model in both packages ────────────────────────────────────────────────
# Each mapping was checked against the predecessor's code (paths under
# BirthDeatMutation_AT.jl/src; the predecessor is read-only and used as it is).
#
#   v1 quantity               ClonalEvolutionTx                      BirthDeathMutation_AT
#   birth b = 1               birth = s -> 1.0                     birthrate = 1.0
#   death d₀ = 0.5            death = s -> 0.5                     deathrate = 0.5
#   competition l·N,          death_density = N -> l*N             competition = l (death gains
#     l = 0.5/10⁵ (K = 10⁵)                                          l·N, simulations.jl:155-159)
#   drug response d_S·m,      response = LinearResponse(4.0)       LinearSensitivity(4.0)
#     d_S = 4                                                        (drug_sensitivity.jl:48)
#   resistance factor α       Resistance(ν = 0.002, α): the clone  SelectionFixedAlpha(α, 0.002,
#     (death α·d_S·m),          kills at α·D(m)                      typemax(Int)): the subclone gets
#     ν = 0.002 per daughter                                         LinearSensitivity(α·d_S)
#                                                                    (drug_sensitivity.jl:64; one
#                                                                    draw per daughter,
#                                                                    simulations.jl:166-179)
#   min rule min(α, α_parent) applied by Resistance                SelectionFixedAlpha uses the
#                                                                    drawn α as it is
#                                                                    (selection.jl:136). Identical
#                                                                    here: α is fixed and parents
#                                                                    have α_parent ∈ {1, α}, so
#                                                                    min(α, α_parent) = α.
#   neutral μ = 1             μ = 1.0: Poisson(μ) per daughter     μ = 1.0, :poisson: Poisson(μ)
#                                                                    per daughter
#                                                                    (cellupdates.jl:57-60). Same
#                                                                    convention, no conversion.
#   new clone per mutation    every resistance mutation is a       every mutation is a new
#                               new clone                            subclone (cellmutation!)
#   max. number of clones     unlimited (max_clones default)       max_subclones = typemax(Int)
#   detection N = 10⁴         growth block, stop_above = 10⁴,      phase 1 of the protocol
#                               restart_on_extinction = true         (until N ≥ N_detect); a run
#                                                                    that dies out before detection
#                                                                    is discarded and rerun on the
#                                                                    continuing RNG (see below)
#   MTD                       MTD(base; N_prog, t_max = Inf)       SimpleTreatment(N_prog,
#                                                                    t_max = Inf, N_min = 0)
#   cyclic adaptive therapy   AdaptiveOnOff(base; N_low, N_up,     AdaptiveTreatment(N_low, N_up,
#                               N_prog, t_max)                       N_prog, t_max)
#   on → off switch           N ≤ N_low                            N ≤ N_low (loop N_low < N)
#   off → on switch           N ≥ N_up                             N ≥ N_up (loop N < N_up)
#   t_max origin              the schedule start = detection       detection (simulations.jl:223,
#                                                                    358)
#   thinning bound            max over living clones               max over all subclones ever
#                                                                    created (only wastes ticks)
#
# Thresholds (the manuscript leaves them open; these are the SPEC §9 illustration
# N_low = N_detect/2, N_up = N_detect, N_prog = 2·N_detect, t_max = 100, at N_detect = 10⁴):
# N_prog = 2·10⁴ for every configuration; adaptive N_low = 5000, N_up = 10⁴, t_max = 100.
#
# Configurations: (A) MTD with α = 0.8 and (B) AdaptiveOnOff with α = 0.8 are the plan's.
# With α = 0.8 and d_S = 4 a resistant cell still has net growth 1 - 0.5 - 3.2 < 0 under full
# dose, so (A) cures every tumour and (B) can never progress (the on phase always shrinks the
# tumour): time to progression is undefined in both, and (A) is compared on time to extinction
# and on the state at fixed times after detection, (B) on the number of cycles and the state at
# fixed times and at t_max. (C) MTD with α = 0 is an extra configuration (not in the plan) in
# which resistant cells are immune, so MTD selects them and the tumour progresses: it exercises
# time to progression and the state at progression.
#
# Predecessor protocol use: `runsimulation(input, selection, rng; returnextinct = true,
# measurement_spec)`. Without `returnextinct` the predecessor reruns any run whose *final*
# population is empty (run.jl:36), which would discard every MTD cure; with it, a run that died
# out before detection (no AtDetection snapshot) is discarded here and rerun with the same,
# continuing RNG, which matches ClonalEvolutionTx's restart_on_extinction and the manuscript
# ("discard simulations that lead to extinction during the initial pre-detection growth
# phase"). Fixed-time states come from the predecessor's AtTime(T) snapshots (T relative to
# detection; taken before the first event after T, i.e. the state at T) and from
# ClonalEvolutionTx's AtTime(t_detect + T). A run that died out before T counts as N(T) = 0 and
# 0 resistant clones (its burden is undefined); a run that ended otherwise before T is left out.
#
# Seeds: run r of configuration c uses Xoshiro(derive_seed(SEED_CE, c, r)) in ClonalEvolutionTx
# and Xoshiro(derive_seed(SEED_BDM, c, r)) passed as the predecessor's `rng` argument (its own
# seeding mechanism; its examples seed the default Xoshiro generator); the two bases differ, so
# the samples are independent.
#
# Statistics and pass criterion (from the plan): for every statistic, the means differ by at most
# 5% (relative to the predecessor's mean) and the two-sample KS test (HypothesisTests'
# ApproximateTwoSampleKSTest) gives p > 0.01. The table also shows the standard error of the
# difference of the means relative to the predecessor's mean ("SE") and z = difference / SE, to
# tell a systematic difference from sampling noise; KS on integer-valued statistics (ties) is
# conservative. Outcome frequencies are compared with Fisher's exact test (most common outcome
# versus the rest).
#
# Timing: each run is simulated twice from the same seed: once with the measurements (for the
# statistics) and once without them (the timing pass: ClonalEvolutionTx without an accumulator,
# the predecessor with an AtDetection trigger and no statistics, the minimum needed to recognise
# a run that died before detection). Measurements draw no random numbers in either package, so
# both passes simulate the same runs. Wall time per run is the timing pass's total over all runs
# (growth to detection, restarts and treatment) divided by the number of runs, after one warm-up
# run per package and configuration. ClonalEvolutionTx runs in tree mode, which keeps the
# genealogy as the predecessor does.

using ClonalEvolutionTx, HypothesisTests, Statistics, Random, Printf
import BirthDeathMutation_AT as BDM

const RUNS = length(ARGS) >= 1 ? parse(Int, ARGS[1]) : 200
const SEED_CE = 20261019
const SEED_BDM = 19102026

# v1 parameters
const B0 = 1.0
const D0 = 0.5
const L = 0.5 / 1e5
const DS = 4.0
const NU = 0.002
const MU = 1.0
const N_DETECT = 10_000
const N_PROG = 20_000
const N_LOW = 5_000
const N_UP = 10_000

const CONFIGS = [
    (id = 1, name = "A: MTD, α = 0.8 (plan)", kind = :mtd, α = 0.8, t_max = Inf,
     times = [0.5, 1.0, 1.5],
     stats = [:t_detect, :det_rclones, :det_rcells, :det_burden, :det_seg, :t_extinct,
              Symbol("N@0.5"), Symbol("N@1.0"), Symbol("N@1.5"), Symbol("rclones@0.5"),
              Symbol("rclones@1.0"), Symbol("burden@1.0")]),
    (id = 2, name = "B: AdaptiveOnOff, α = 0.8 (plan)", kind = :adaptive, α = 0.8, t_max = 100.0,
     times = [25.0, 50.0, 75.0],
     stats = [:cycles, Symbol("N@25.0"), Symbol("N@50.0"), Symbol("N@75.0"),
              Symbol("rcells@50.0"), :end_N, :end_rclones, :end_maxrclone, :end_burden,
              :end_seg]),
    (id = 3, name = "C: MTD, α = 0 (extra: exercises progression)", kind = :mtd, α = 0.0,
     t_max = Inf, times = [1.0, 3.0],
     stats = [:t_progressed, Symbol("N@1.0"), Symbol("N@3.0"), :end_rclones, :end_maxrclone,
              :end_burden, :end_seg]),
]

const LABELS = Dict(
    :t_detect => "time to detection", :det_rclones => "resistant clones at detection",
    :det_rcells => "resistant cells at detection", :det_burden => "mean burden at detection",
    :det_seg => "segregating mutations at det.", :t_extinct => "time to extinction",
    :t_progressed => "time to progression", :cycles => "AT on phases by t_max",
    :end_N => "N at the end", :end_rclones => "resistant clones at the end",
    :end_maxrclone => "largest resistant clone, end", :end_burden => "mean burden at the end",
    :end_seg => "segregating mutations, end")
function label(k::Symbol)
    haskey(LABELS, k) && return LABELS[k]
    name, T = split(String(k), '@')
    what = Dict("N" => "N", "rclones" => "resistant clones", "rcells" => "resistant cells",
                "burden" => "mean burden")[name]
    return "$what at t = $T"
end

# ── State statistics (same definitions in both packages) ─────────────────────────────────────
# resistant: clones with α < 1 (with a fixed α every non-founder clone is resistant);
# burden: mean neutral mutations per living cell; segregating: distinct neutral mutations in
# the living cells (the sum of the site-frequency spectrum).

resistant_clones(pop) = count(k -> pop.clones.α[k] < 1.0, living_clones(pop.clones))
largest_resistant_clone(pop) =
    maximum((pop.clones.size[k] for k in living_clones(pop.clones) if pop.clones.α[k] < 1.0);
            init = 0)
const CE_STATS = [popsize, resistant_clones, n_resistant, largest_resistant_clone, mean_burden,
                  total_mutational_burden]

ce_state(st) = (N = st.popsize, rclones = st.resistant_clones, rcells = st.n_resistant,
                maxrclone = st.largest_resistant_clone, burden = st.mean_burden,
                seg = st.total_mutational_burden)

# The predecessor's snapshot: clone_size_dist is [wildtype cells, size of every mutant subclone].
function bdm_state(s)
    csd = s.clone_size_dist
    res = csd[2:end]
    m = s.neutral_muts_per_cell
    return (N = sum(csd), rclones = count(>(0), res), rcells = sum(res; init = 0),
            maxrclone = maximum(res; init = 0), burden = isempty(m) ? NaN : mean(m),
            seg = sum(s.sfs; init = 0))
end

function put_state!(rec, prefix, st)
    for (k, v) in pairs(st)
        rec[Symbol(prefix, k)] = Float64(v)
    end
    return rec
end

# Fixed-time entries: the state at T, or N = 0 and no resistant clones if the run died out before.
function put_time!(rec, T, st, outcome)
    if st !== nothing
        for (k, v) in pairs(st)
            rec[Symbol(k, "@", T)] = Float64(v)
        end
    elseif outcome === :extinct
        rec[Symbol("N@", T)] = 0.0
        rec[Symbol("rclones@", T)] = 0.0
        rec[Symbol("rcells@", T)] = 0.0
    end
    return rec
end

function put_end!(rec, outcome, t_end, st)
    rec[:t_extinct] = outcome === :extinct ? t_end : NaN
    rec[:t_progressed] = outcome === :progressed ? t_end : NaN
    st === nothing || put_state!(rec, "end_", st)
    return rec
end

# ── ClonalEvolutionTx ──────────────────────────────────────────────────────────────────────────

ce_base(α) = BirthDeathBlock(birth = s -> B0, death = s -> D0, death_density = N -> L * N,
                             response = LinearResponse(DS), μ = MU,
                             mutations = [Resistance(ν = NU, α = α)], label = "growth")

ce_schedule(cfg, base) = cfg.kind === :mtd ? MTD(base; N_prog = N_PROG, t_max = cfg.t_max) :
    AdaptiveOnOff(base; N_low = N_LOW, N_up = N_UP, N_prog = N_PROG, t_max = cfg.t_max)

# One run; with `stats = false` nothing is measured (timing pass) and `nothing` is returned.
function ce_run(cfg, seed; stats::Bool = true)
    rng = Xoshiro(seed)
    base = ce_base(cfg.α)
    pop = initialize_population(1; mode = :tree)
    simulate!(pop, remake(base; stop_above = N_DETECT, restart_on_extinction = true), rng)
    sched = ce_schedule(cfg, base)
    stats || (simulate!(pop, sched, rng); return nothing)
    t_det = pop.t
    rec = Dict{Symbol,Float64}(:t_detect => t_det, :restarts => pop.log[1].restarts)
    put_state!(rec, "det_", ce_state(measure(pop, CE_STATS)))
    acc = MeasurementAccumulator(MeasurementSpec(triggers = [AtTime(t_det + T) for T in cfg.times],
                                                 stats = CE_STATS))
    simulate!(pop, sched, rng; accumulator = acc)
    snaps = finalize_measurements(acc).snapshots
    outcome = schedule_outcome(pop)
    for T in cfg.times
        i = findfirst(s -> s.trigger == AtTime(t_det + T), snaps)
        put_time!(rec, T, i === nothing ? nothing : ce_state(snaps[i].stats), outcome)
    end
    rec[:cycles] = count(b -> b.label == "AT on", pop.log)
    put_end!(rec, outcome, pop.t - t_det, popsize(pop) > 0 ? ce_state(measure(pop, CE_STATS)) : nothing)
    return rec, outcome
end

# ── BirthDeathMutation_AT ────────────────────────────────────────────────────────────────────

function bdm_input(cfg)
    common = (birthrate = B0, deathrate = D0, competition = L,
              wildtype_sensitivity = BDM.LinearSensitivity(DS), N_detect = N_DETECT,
              N_prog = N_PROG, t_max = cfg.t_max, μ = MU)
    return cfg.kind === :mtd ? BDM.SimpleTreatment(; common...) :
           BDM.AdaptiveTreatment(; common..., N_low = N_LOW, N_up = N_UP)
end

function bdm_spec(cfg, stats::Bool)
    stats || return BDM.MeasurementSpec(snapshot_triggers = [BDM.AtDetection()])
    triggers = [BDM.AtDetection(), BDM.AtTreatmentToggle(), [BDM.AtTime(T) for T in cfg.times]...,
                BDM.AtEnd()]
    return BDM.MeasurementSpec(snapshot_triggers = triggers,
                               snapshot_stats = [BDM.SFS(), BDM.CloneSizeDist(), BDM.NeutralMutPerCell()])
end

# One run, rerun on the continuing RNG while the population dies out before detection.
function bdm_run(cfg, seed; stats::Bool = true)
    rng = Xoshiro(seed)
    input = bdm_input(cfg)
    selection = BDM.SelectionFixedAlpha(cfg.α, NU, typemax(Int))
    spec = bdm_spec(cfg, stats)
    restarts = 0
    while true
        sim = BDM.runsimulation(input, selection, rng; returnextinct = true, measurement_spec = spec)
        snaps = sim.measurements.snapshots
        if !any(s -> s.trigger isa BDM.AtDetection, snaps)
            restarts += 1
            continue
        end
        stats || return nothing
        return bdm_record(cfg, sim, restarts)
    end
end

function bdm_record(cfg, sim, restarts)
    snaps = sim.measurements.snapshots
    pop = sim.output
    N_end = length(BDM.allcells(pop))
    outcome = N_end == 0 ? :extinct : N_end >= N_PROG ? :progressed : :tmax
    det = snaps[findfirst(s -> s.trigger isa BDM.AtDetection, snaps)]
    fin = snaps[findfirst(s -> s.trigger isa BDM.AtEnd, snaps)]
    t_det = BDM.age(pop) - fin.t                    # snapshot times are relative to detection
    rec = Dict{Symbol,Float64}(:t_detect => t_det, :restarts => restarts)
    put_state!(rec, "det_", bdm_state(det))
    for T in cfg.times
        i = findfirst(s -> s.trigger isa BDM.AtTime && s.trigger.t == T, snaps)
        put_time!(rec, T, i === nothing ? nothing : bdm_state(snaps[i]), outcome)
    end
    toggles = count(s -> s.trigger isa BDM.AtTreatmentToggle, snaps)   # before every on and off phase
    rec[:cycles] = cld(toggles, 2)
    put_end!(rec, outcome, fin.t, N_end > 0 ? bdm_state(fin) : nothing)
    return rec, outcome
end

# ── Comparison ───────────────────────────────────────────────────────────────────────────────

function compare_stat(k, recs_ce, recs_bdm)
    x = filter(!isnan, [get(r, k, NaN) for r in recs_ce])
    y = filter(!isnan, [get(r, k, NaN) for r in recs_bdm])
    if length(x) < 2 || length(y) < 2
        @printf("  %-32s %5d %5d   too few values to compare\n", label(k), length(x), length(y))
        return (stat = k, pass = false)
    end
    mx, my = mean(x), mean(y)
    rel = my == 0 ? (mx == 0 ? 0.0 : Inf) : (mx - my) / abs(my)
    se = sqrt(var(x) / length(x) + var(y) / length(y))
    z = se == 0 ? 0.0 : (mx - my) / se
    serel = my == 0 ? NaN : se / abs(my)
    p = allequal([x; y]) ? 1.0 : pvalue(ApproximateTwoSampleKSTest(x, y))
    pass = abs(rel) <= 0.05 && p > 0.01
    @printf("  %-32s %5d %5d %11.4g %11.4g %+8.2f%% %6.2f%% %6.2f %8.3f   %s\n", label(k),
            length(x), length(y), mx, my, 100rel, 100serel, z, p, pass ? "pass" : "FAIL")
    return (stat = k, pass = pass, rel = rel, p = p, z = z)
end

function compare_outcomes(out_ce, out_bdm)
    kinds = sort(unique([out_ce; out_bdm]))
    @printf("  outcomes   ClonalEvolutionTx: %s\n", join(["$o $(count(==(o), out_ce))" for o in kinds], ", "))
    @printf("             BirthDeathMutation_AT: %s\n", join(["$o $(count(==(o), out_bdm))" for o in kinds], ", "))
    length(kinds) == 1 && (println("             identical (one outcome only)"); return 1.0)
    top = argmax(o -> count(==(o), [out_ce; out_bdm]), kinds)
    a, c = count(==(top), out_ce), count(==(top), out_bdm)
    p = pvalue(FisherExactTest(a, length(out_ce) - a, c, length(out_bdm) - c))
    @printf("             Fisher exact test (%s versus the rest): p = %.3f\n", top, p)
    return p
end

function timed_runs(f, cfg, base, n)
    f(cfg, derive_seed(base, 0, 0); stats = false)       # warm-up (compilation), not timed
    GC.gc()
    return @elapsed for r in 1:n
        f(cfg, derive_seed(base, cfg.id, r); stats = false)
    end
end

function main()
    t_script = time()
    println("Comparison of ClonalEvolutionTx with BirthDeathMutation_AT (manuscript_v1 model)")
    @printf("Julia %s, %s, %d threads; %s\n", VERSION, Sys.cpu_info()[1].model, Threads.nthreads(),
            Libc.strftime("%Y-%m-%d %H:%M", time()))
    @printf("ClonalEvolutionTx %s, BirthDeathMutation_AT %s; %d runs per package and configuration\n",
            pkgversion(ClonalEvolutionTx), pkgversion(BDM), RUNS)
    @printf("b = %g, d0 = %g, l = %g, d_S = %g, nu = %g, mu = %g, N_detect = %d, N_prog = %d, ",
            B0, D0, L, DS, NU, MU, N_DETECT, N_PROG)
    @printf("adaptive N_low = %d, N_up = %d\n", N_LOW, N_UP)
    println("Pass: |rel. diff of means| <= 5% and KS p > 0.01. SE: standard error of the difference ",
            "of the means, relative to the BDM_AT mean; z = difference / SE.")
    results, timings = [], []
    for cfg in CONFIGS
        @printf("\n=== %s, t_max = %g, fixed times after detection %s ===\n", cfg.name, cfg.t_max,
                join(cfg.times, ", "))
        ce_run(cfg, derive_seed(SEED_CE, 0, 0)); bdm_run(cfg, derive_seed(SEED_BDM, 0, 0))  # warm-up
        runs_ce = [ce_run(cfg, derive_seed(SEED_CE, cfg.id, r)) for r in 1:RUNS]
        runs_bdm = [bdm_run(cfg, derive_seed(SEED_BDM, cfg.id, r)) for r in 1:RUNS]
        recs_ce, recs_bdm = first.(runs_ce), first.(runs_bdm)
        compare_outcomes(last.(runs_ce), last.(runs_bdm))
        @printf("  restarts before detection per run (mean): CE %.3f, BDM_AT %.3f\n",
                mean(r[:restarts] for r in recs_ce), mean(r[:restarts] for r in recs_bdm))
        @printf("  %-32s %5s %5s %11s %11s %9s %7s %6s %8s\n", "statistic", "n CE", "n BDM",
                "mean CE", "mean BDM", "rel.diff", "SE", "z", "KS p")
        for k in cfg.stats
            push!(results, merge((cfg = cfg.name,), compare_stat(k, recs_ce, recs_bdm)))
        end
        t_ce = timed_runs(ce_run, cfg, SEED_CE, RUNS) / RUNS
        t_bdm = timed_runs(bdm_run, cfg, SEED_BDM, RUNS) / RUNS
        @printf("  wall time per run (timing pass): ClonalEvolutionTx %.4f s, BirthDeathMutation_AT %.4f s, speed-up %.1f×\n",
                t_ce, t_bdm, t_bdm / t_ce)
        push!(timings, (cfg = cfg.name, ce = t_ce, bdm = t_bdm))
    end
    npass = count(r -> r.pass, results)
    @printf("\nSummary: %d of %d statistics pass.\n", npass, length(results))
    for r in results
        r.pass || @printf("  FAIL  %s: %s\n", r.cfg, label(r.stat))
    end
    println("\nWall time per run (s):")
    for t in timings
        @printf("  %-48s CE %.4f   BDM_AT %.4f   speed-up %.1f×\n", t.cfg, t.ce, t.bdm, t.bdm / t.ce)
    end
    @printf("\nTotal script time %.0f s\n", time() - t_script)
end

main()
