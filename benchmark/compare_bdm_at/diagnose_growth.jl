# Diagnostic of the growth phase of the comparison with BirthDeathMutation_AT.jl (compare.jl;
# see "Comparison with BirthDeathMutation_AT.jl" in docs/src/performance.md). It is not part
# of compare.jl, whose run count stays at 200.
#
# What it shows: compare.jl's one failing statistic, the mean number of resistant cells at
# detection, has a heavy-tailed Luria–Delbrück distribution, so 200 runs cannot resolve a 5%
# difference of the means. This script runs only the growth phase (one cell to N = 10⁴, the
# growth block of compare.jl, no treatment) many times in both packages and compares per run:
#   R         resistant cells at detection;
#   R_rb      Σᵢ (1 - (1-ν)^Gᵢ), the expectation of R given the genealogy (Gᵢ = divisions on
#             cell i's lineage; resistance is neutral at m = 0, so the tree is independent of
#             the resistance marks and the stopping time depends on N only), which is free of
#             the Luria–Delbrück noise;
#   burden    mean neutral burden (E[burden | tree] = μ·mean(G));
#   restarts  extinctions before detection, against the analytic mean of the birth–death
#             chain with b = 1, d(N) = 0.5 + l·N (1.00004).
# Output: a table of means, relative difference, SE, z, KS p-value and medians; the raw and
# Rao–Blackwellised means of R; the tails of R. Printed to stdout, nothing written to disk.
#
# Inputs : optional runs per package (default 20 000, the run on the Performance page) and
#          the two base seeds (default 7001 7002).
# Run    : julia --project=benchmark/compare_bdm_at benchmark/compare_bdm_at/diagnose_growth.jl [runs [seed_ce seed_bdm]]
# Runtime: 20 000 runs per package took 62 s (ClonalEvolutionTx) and 190 s (BirthDeathMutation_AT)
#          on an Apple M4 Pro with Julia 1.12.1, plus start-up and compilation (about 5 minutes
#          in all); pass a smaller run count for a quick check.
using ClonalEvolutionTx, HypothesisTests, Statistics, Random, Printf
import BirthDeathMutation_AT as BDM

const RUNS = length(ARGS) >= 1 ? parse(Int, ARGS[1]) : 20_000
const L = 0.5 / 1e5
const NU = 0.002
const ND = 10_000

function ce_growth(seed)
    rng = Xoshiro(seed)
    base = BirthDeathBlock(birth = s -> 1.0, death = s -> 0.5, death_density = N -> L * N,
                           response = LinearResponse(4.0), μ = 1.0,
                           mutations = [Resistance(ν = NU, α = 0.8)],
                           stop_above = ND, restart_on_extinction = true)
    pop = initialize_population(1; mode = :tree)
    simulate!(pop, base, rng)
    G = divisions_per_cell(pop)
    return (R = n_resistant(pop), R_rb = sum(g -> 1 - (1 - NU)^g, G), burden = mean_burden(pop),
            G = mean(G), restarts = pop.log[1].restarts, t = pop.t)
end

# Phase 1 of the predecessor's protocols is branchingprocess! with m = 0 (simulations.jl:208-211);
# SingleCellExpansion runs the same function with the same arguments (simulations.jl:24-38).
function bdm_growth(seed)
    rng = Xoshiro(seed)
    input = BDM.SingleCellExpansion(Nmax = ND, birthrate = 1.0, deathrate = 0.5, competition = L,
                                    μ = 1.0)
    sel = BDM.SelectionFixedAlpha(0.8, NU, typemax(Int))
    restarts = 0
    while true
        sim = BDM.runsimulation(input, sel, rng; returnextinct = true)
        pop = sim.output
        cells = BDM.allcells(pop)
        if isempty(cells)
            restarts += 1
            continue
        end
        G = Int[]
        for c in cells
            g, node = 0, c
            while !isnothing(node.parent)
                node = node.parent
                g += 1
            end
            push!(G, g)
        end
        return (R = sum(s.size for s in pop.subclones if s.parentid > 0; init = 0),
                R_rb = sum(g -> 1 - (1 - NU)^g, G), burden = mean(BDM.mutations_per_cell(pop)),
                G = mean(G), restarts = restarts, t = BDM.age(pop))
    end
end

function row(name, x, y; ref = nothing)
    mx, my = mean(x), mean(y)
    se = sqrt(var(x) / length(x) + var(y) / length(y))
    p = pvalue(ApproximateTwoSampleKSTest(Float64.(x), Float64.(y)))
    @printf("  %-34s %11.5g %11.5g %+8.2f%% %7.2f%% %6.2f %8.3f  medians %g / %g", name, mx, my,
            100(mx - my) / my, 100se / my, (mx - my) / se, p, median(x), median(y))
    ref === nothing || @printf("   analytic %.5f", ref)
    println()
end

ce_growth(1); bdm_growth(1)
const SEED_CE_D, SEED_BDM_D = (length(ARGS) >= 3 ? parse.(Int, ARGS[2:3]) : (7001, 7002))
t_ce = @elapsed ce = [ce_growth(derive_seed(SEED_CE_D, r)) for r in 1:RUNS]
t_bdm = @elapsed bd = [bdm_growth(derive_seed(SEED_BDM_D, r)) for r in 1:RUNS]
@printf("Growth to N = %d, %d runs per package (CE %.1f s, BDM_AT %.1f s)\n", ND, RUNS, t_ce, t_bdm)
@printf("  %-34s %11s %11s %9s %8s %6s %8s\n", "statistic", "mean CE", "mean BDM", "rel.diff", "SE", "z", "KS p")
# analytic restart mean: P(reach N_D before 0 | 1 cell) for b = 1, d(N) = 0.5 + l·N
let s = 0.0, p = 1.0
    for k in 1:ND-1
        p *= 0.5 + L * k
        s += p
    end
    global Q = 1 - 1 / (1 + s)
end
row("restarts before detection", getfield.(ce, :restarts), getfield.(bd, :restarts); ref = Q / (1 - Q))
row("time to detection", getfield.(ce, :t), getfield.(bd, :t))
row("resistant cells R", getfield.(ce, :R), getfield.(bd, :R))
row("E[R | genealogy] (Rao-Blackwell)", getfield.(ce, :R_rb), getfield.(bd, :R_rb))
row("mean lineage divisions G", getfield.(ce, :G), getfield.(bd, :G))
row("mean neutral burden", getfield.(ce, :burden), getfield.(bd, :burden))
for (name, v) in (("CE", ce), ("BDM_AT", bd))
    R, Rrb = getfield.(v, :R), getfield.(v, :R_rb)
    @printf("  %-6s mean R %.1f ± %.1f (SE) versus mean E[R | genealogy] %.1f ± %.1f; first 200 runs: mean R %.1f\n",
            name, mean(R), std(R) / sqrt(length(R)), mean(Rrb), std(Rrb) / sqrt(length(Rrb)), mean(R[1:200]))
end
# Tail: where does the difference of the means come from?
for (name, v) in (("CE", ce), ("BDM_AT", bd))
    R = sort(Float64.(getfield.(v, :R)))
    n = length(R)
    tm(f) = mean(R[1:floor(Int, (1 - f) * n)])
    @printf("  %-6s trimmed means (drop top 1%% / 0.1%%): %.1f / %.1f; runs with R ≥ 2000: %d, ≥ 5000: %d, ≥ 9000: %d; contribution of runs with R ≥ 2000 to the mean: %.1f\n",
            name, tm(0.01), tm(0.001), count(>=(2000), R), count(>=(5000), R), count(>=(9000), R),
            sum(r for r in R if r >= 2000) / n)
end
