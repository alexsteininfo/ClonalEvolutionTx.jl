# Performance

The benchmark is a local script, `benchmark/benchmark.jl`; nothing is run on the docs build and
nothing is written to disk. Sizes are given on the command line in increasing order, warm-up runs
are excluded from the timings, live memory comes from the GC counter, and a size whose
estimated peak memory exceeds 60% of the free RAM is skipped.

## Running the benchmark

```bash
# one or more population sizes; --parallel N for the threaded mode
julia --project=. benchmark/benchmark.jl 1000 10000 100000 [--parallel N]

# the full local set (Julia from PATH; THREADS env var for the parallel mode, default half the
# CPU threads): serial sizes 10³ … 10⁷, then --parallel 100000
bash benchmark/run_benchmarks.sh
```

The log of a full run is kept in `.development/benchmark/` (not tracked).

## Modes

1. Neutral growth to `N` (10³ … 10⁷), tree mode versus cell mode.
2. Growth with resistance (`ν = 0.002`) and `μ = 1`, tree versus cell mode.
3. MTD from 10⁶ cells with `α ∈ {0, 0.8}`; reports the thinning acceptance fraction.
4. Cyclic adaptive therapy (the two blocks of `AdaptiveOnOff`, `N_low = N/2`, `N_up = N`, in a
   `FunctionSchedule` that stops after exactly six on/off cycles), on a tumour of 10⁵ cells.
5. Stabilising therapy (the dose is recomputed after every event), 10⁵ cells; reported per event.
6. Statistics on one 10⁶ population: SFS, burdens, mutation frequencies, Hill numbers and a
   sample of 1000.
7. `--parallel N`: `4 × nthreads` runs serially, then on threads; runs per second.

## Results

**What was measured here, and what was not.** Machine: Apple M4 Pro (14 cores, 10 CPU threads
reported by Julia), 48 GB RAM, macOS; Julia 1.12.1; ClonalEvolutionTx 0.1.0-DEV; 2026-10-19; Julia
started with `-t 4`. Measured for this page, one invocation each:
`julia --project=. -t 4 benchmark/benchmark.jl 1000 10000 100000 1000000` (serial modes up to 10⁶,
treatment and statistics) and `julia --project=. -t 4 benchmark/benchmark.jl --parallel 10000`.
**Not measured here:** the 10⁷ rows and the full `--parallel 100000` run of
`bash benchmark/run_benchmarks.sh`; they are left to the author and marked "not run" below.
Every value is a single seeded run (seeds from `derive_seed`), warm-up excluded, with no repeats;
timings at 10⁶ vary by roughly 10-40% between invocations (the tree-mode 10⁶ rows of three
invocations ranged 0.36-0.50 s), so read the figures as orders of magnitude. Time per cell is the
wall time divided by the final `N`; more cells are created than survive (with `b = 1`, `d = 0.5`,
about 3 thinning ticks are spent per final cell), so the cost per *created* cell is lower than the
figure shown. The script skips any run whose estimated peak memory (1.5 × bytes per cell × N)
exceeds 60% of the free RAM; nothing was skipped here.

### Neutral growth, tree versus cell mode

Birth 1, death 0.5, `μ = 0`, no mutations.

| N | tree: time (s) | tree: µs per cell | tree: live memory | cells: time (s) | cells: µs per cell | cells: live memory |
|---:|---:|---:|---:|---:|---:|---:|
| 10³ | 0.000 | 0.16 | 0.3 MB | 0.000 | 0.08 | 0.05 MB |
| 10⁴ | 0.001 | 0.12 | 2.3 MB | 0.001 | 0.06 | 0.4 MB |
| 10⁵ | 0.015 | 0.15 | 20.7 MB | 0.006 | 0.06 | 3.8 MB |
| 10⁶ | 0.362 | 0.36 | 172 MB | 0.116 | 0.12 | 45.9 MB |
| 10⁷ | not run | | | not run | | |

### Growth with resistance and neutral mutations

`Resistance(ν = 0.002, α = 0.8)` and `μ = 1` neutral mutations per daughter.

| N | tree: time (s) | tree: µs per cell | tree: live memory | cells: time (s) | cells: µs per cell | cells: live memory |
|---:|---:|---:|---:|---:|---:|---:|
| 10³ | 0.000 | 0.26 | 0.3 MB | 0.000 | 0.19 | 0.05 MB |
| 10⁴ | 0.002 | 0.22 | 2.3 MB | 0.001 | 0.15 | 0.4 MB |
| 10⁵ | 0.028 | 0.28 | 20.8 MB | 0.016 | 0.16 | 3.9 MB |
| 10⁶ | 0.504 | 0.50 | 173 MB | 0.234 | 0.23 | 47.0 MB |
| 10⁷ | not run | | | not run | | |

In tree mode the GC share is 6-12% of the time at 10⁶; in cell mode it is 0.3% or less.

### Treatment

Linear response `d_S = 4`, `Resistance(ν = 0.002, α)`. The tumour is grown first in cell mode
(untimed); the time is the treatment phase only. Ticks and accepted events are summed over the
treatment blocks. The cyclic scenario is a `FunctionSchedule` built from the two blocks of
`AdaptiveOnOff` (`N_low = N/2`, `N_up = N`) that stops after exactly 6 complete on/off cycles
(12 blocks).

| scenario | N at start | time (s) | accepted / ticks | outcome |
|:--|---:|---:|---:|:--|
| MTD, `α = 0` | 10⁶ | 0.498 | 0.560 (7.27·10⁶ of 12.98·10⁶) | progressed |
| MTD, `α = 0.8` | 10⁶ | 0.065 | 0.988 (1.58·10⁶ of 1.60·10⁶) | extinct |
| cyclic adaptive therapy, 6 cycles | 10⁵ | 0.041 | 0.994 | 12 blocks, 6 complete cycles |
| stabilising therapy | 10⁵ | 0.119 | 0.971 | `t_max` reached |

Stabilising therapy costs about 29 ns per tick (30 ns per accepted event), so recomputing the dose
after every event does not slow the loop noticeably. With `α = 0` the resistant cells are immune
and the bound `(max αₖ)·D(m)` stays at the full drug term while sensitive cells die out, which
is why only 56% of the ticks are accepted there; with `α = 0.8` the bound is tight.

### Statistics on a 10⁶ population

One tree population grown to 10⁶ cells (resistance, `μ = 1`).

| statistic | time (s) |
|:--|---:|
| SFS | 0.0122 |
| `mutations_per_cell` | 0.0014 |
| `mutation_frequencies` | < 0.0001 |
| `sfs_hill_number`, q = 0 / 1 / 2 / ∞ | 0.0007 / 0.0005 / 0.0006 / 0.0006 |
| `sample_population(pop, 1000)` | 0.0031 |

### Parallel runs

Resistance model in cell mode, `4 × nthreads` = 16 runs, serial then `Threads.@threads`.

| N | threads | runs per second (serial) | runs per second (threaded) | speedup |
|---:|---:|---:|---:|---:|
| 10⁴ | 4 | 547 | 246 | 0.4 |
| 10⁵ (`run_benchmarks.sh`, more threads) | | not run | not run | |

The 10⁴ run is far too short to measure scaling (16 runs take 0.03 s serially; thread start and
allocation dominate, so threading is slower), and it is not representative. The scaling of the
full `--parallel 100000` run is left to the author; a one-off look at 10⁵ on 4 threads during
development suggested a modest speedup (below the thread count), which the author's run should
confirm or refute.

## Comparison with the targets of the specification

| target (SPEC §5) | measured (10⁶ cells, single runs) | verdict |
|:--|:--|:--|
| cell mode ≤ ~0.5 µs per created cell | 0.12 µs (neutral), 0.23 µs (resistance, `μ = 1`) per final cell | met, by a factor of 2 or more (cost per created cell is lower still) |
| tree mode within ~1.5× of NonMarkovEvolution's thinning path (1.2 µs per cell at 10⁶, from its performance page) | 0.36 µs (neutral), 0.50 µs (resistance, `μ = 1`) per final cell | met: more than 2× faster, not slower |

The comparison with NonMarkovEvolution is indicative only: the models differ (that benchmark
uses Exponential division times with fitness-dependent rates and `ν = 0.5` driver mutations per
daughter, this one a constant-rate birth-death process), and the machines differ (its figures
come from a cluster node). The 10⁷ rows would show whether the cost per cell stays flat at the
largest size; they are not part of this page yet.

## Comparison with BirthDeathMutation_AT.jl

`benchmark/compare_bdm_at/` (its own environment, with BirthDeathMutation_AT.jl added by path)
runs the `manuscript_v1` model in both packages and compares the speed and the distributions of a
fixed set of statistics. BirthDeathMutation_AT is the predecessor and the reference
implementation of v1. The comparison is run once per release and is not part of the test suite.

```bash
julia --project=benchmark/compare_bdm_at benchmark/compare_bdm_at/compare.jl        # 200 runs
# growth-phase diagnostic (see the end of the results below): 20 000 runs, about 5 minutes
julia --project=benchmark/compare_bdm_at benchmark/compare_bdm_at/diagnose_growth.jl
```

The first instantiation of this environment downloads and precompiles the predecessor's
dependencies (CairoMakie, DataFrames, GLM, …), which takes a few minutes. The full run took 626 s.

**Model and mapping.** `b = 1`, `d₀ = 0.5`, competition `l = 0.5/10⁵` (`death_density = N -> l*N`
↔ `competition = l`), `LinearResponse(4.0)` ↔ `LinearSensitivity(4.0)`, `Resistance(ν = 0.002, α)`
↔ `SelectionFixedAlpha(α, 0.002, typemax(Int))`, `μ = 1` neutral mutations per daughter in both
(the predecessor also draws `Poisson(μ)` per daughter, so no conversion is needed), detection at
`N = 10⁴`. Thresholds: `N_prog = 2·10⁴`; adaptive `N_low = 5000`, `N_up = 10⁴`, `t_max = 100`,
measured from detection in both packages (the SPEC §9 illustration scaled to `N_detect = 10⁴`;
the manuscript does not fix them). The header of `compare.jl` lists every mapping, with the
lines of the predecessor's code it was checked against. ClonalEvolutionTx runs in tree mode, which
keeps the genealogy as the predecessor does.

**Configurations.** (A) MTD and (B) cyclic adaptive therapy (`MTD` ↔ `SimpleTreatment`,
`AdaptiveOnOff` ↔ `AdaptiveTreatment`) with `α = 0.8` are the planned comparisons. With `α = 0.8` a
resistant cell still has net growth `1 − 0.5 − 3.2 < 0` under full dose. So (A) cures every
tumour (200 of 200 extinct in both packages), and (B) cannot progress (200 of 200 reach `t_max`
in both). Time to progression is therefore undefined in both. (A) is compared on time to
extinction and on the state at fixed times after detection; (B) on the number of on phases and
on the state at fixed times and at `t_max`. (C) MTD with `α = 0` is an **extra configuration, not
part of the plan**: resistant cells are immune, so the tumour progresses (200 of 200 in both).
It exercises time to progression and the state at progression.

**Method.** There are 200 runs per package and configuration, each with its own seed:
`Xoshiro(derive_seed(base, configuration, run))`, with different bases for the two packages, so
the samples are independent. The predecessor is called through
`runsimulation(…; returnextinct = true)`; a run that dies out before detection is discarded and
rerun on the continuing random stream, as `restart_on_extinction` does. Fixed-time states come
from `AtTime` snapshots in both packages; a run that died out before that time counts as `N = 0`
with no resistant clones. "Resistant" means `α < 1`, "burden" is the mean number of neutral
mutations per living cell, and "segregating" is the number of distinct neutral mutations in the
living cells (the sum of the SFS).

The pass criterion is from the plan: the means differ by at most 5% (relative to the
predecessor's mean) and the two-sample KS test (`ApproximateTwoSampleKSTest`) gives
`p > 0.01`. SE is the standard error of the difference of the means, relative to the
predecessor's mean. KS is conservative for integer-valued statistics (ties). Machine:
Apple M4 Pro, macOS, Julia 1.12.1, one thread; ClonalEvolutionTx 0.1.0-DEV, BirthDeathMutation_AT
0.2.1; 2026-10-19.

### (A) MTD, `α = 0.8`

| statistic | mean, ClonalEvolutionTx | mean, BDM_AT | difference | SE | KS p | verdict |
|:--|---:|---:|---:|---:|---:|:--|
| time to detection | 18.25 | 18.42 | −0.96% | 1.52% | 0.711 | pass |
| resistant clones at detection | 56.02 | 56.37 | −0.62% | 1.29% | 0.923 | pass |
| resistant cells at detection | 766.7 | 710.8 | +7.87% | 17.03% | 0.393 | **fail** (mean; see below) |
| mean burden at detection | 34.54 | 34.82 | −0.82% | 1.25% | 0.923 | pass |
| segregating mutations at detection | 28 400 | 28 390 | +0.01% | 0.08% | 0.997 | pass |
| time to extinction | 2.870 | 2.822 | +1.73% | 1.44% | 0.465 | pass |
| `N` at `t = 0.5` | 1783 | 1778 | +0.26% | 0.62% | 0.711 | pass |
| `N` at `t = 1.0` | 327.6 | 322.7 | +1.54% | 1.46% | 0.923 | pass |
| `N` at `t = 1.5` | 60.41 | 60.02 | +0.66% | 2.82% | 0.270 | pass |
| resistant clones at `t = 0.5` | 28.70 | 28.83 | −0.45% | 1.81% | 0.964 | pass |
| resistant clones at `t = 1.0` | 11.88 | 11.88 | 0.00% | 2.98% | 1.000 | pass |
| mean burden at `t = 1.0` | 36.60 | 36.91 | −0.85% | 1.21% | 0.627 | pass |

### (B) `AdaptiveOnOff`, `α = 0.8`, `t_max = 100`

| statistic | mean, ClonalEvolutionTx | mean, BDM_AT | difference | SE | KS p | verdict |
|:--|---:|---:|---:|---:|---:|:--|
| on phases by `t_max` | 58.05 | 58.06 | −0.03% | 0.04% | 1.000 | pass |
| `N` at `t = 25` | 7778 | 7782 | −0.06% | 0.83% | 0.792 | pass |
| `N` at `t = 50` | 6714 | 7013 | −4.27% | 2.52% | 0.270 | pass |
| `N` at `t = 75` | 6468 | 6434 | +0.53% | 1.45% | 0.220 | pass |
| resistant cells at `t = 50` | 6391 | 6677 | −4.29% | 2.60% | 0.220 | pass |
| `N` at the end | 7840 | 7875 | −0.45% | 1.40% | 0.964 | pass |
| resistant clones at the end | 152.4 | 151.7 | +0.44% | 0.95% | 0.544 | pass |
| largest resistant clone at the end | 1389 | 1413 | −1.73% | 6.17% | 0.220 | pass |
| mean burden at the end | 236.0 | 235.5 | +0.18% | 0.23% | 0.711 | pass |
| segregating mutations at the end | 65 430 | 65 540 | −0.16% | 0.40% | 0.792 | pass |

### (C) MTD, `α = 0` (extra configuration)

| statistic | mean, ClonalEvolutionTx | mean, BDM_AT | difference | SE | KS p | verdict |
|:--|---:|---:|---:|---:|---:|:--|
| time to progression | 8.026 | 7.913 | +1.42% | 1.92% | 0.393 | pass |
| `N` at `t = 1.0` | 1407 | 1429 | −1.53% | 12.39% | 0.178 | pass |
| `N` at `t = 3.0` (196 runs per package not yet progressed) | 2377 | 2492 | −4.59% | 9.01% | 0.135 | pass |
| resistant clones at progression | 160.7 | 159.2 | +0.93% | 0.86% | 0.792 | pass |
| largest resistant clone at progression | 6409 | 6618 | −3.16% | 6.12% | 0.864 | pass |
| mean burden at progression | 51.41 | 50.98 | +0.84% | 1.02% | 0.465 | pass |
| segregating mutations at progression | 59 040 | 58 990 | +0.08% | 0.06% | 0.068 | pass |

**Result: 28 of 29 statistics pass**, and the outcome frequencies are identical. The failure is
the mean number of resistant cells at detection: +7.87%, while KS gives `p = 0.39`. This is
sampling noise, not a difference between the models. The statistic has a Luria–Delbrück
distribution with a heavy right tail: rare runs with an early resistance mutation reach thousands
of resistant cells. Its standard error at 200 runs is 17% of the mean, three times the
tolerance, so the 5% criterion cannot be resolved for it at this run count. A separate
diagnostic of the growth phase, `benchmark/compare_bdm_at/diagnose_growth.jl` (not part of
`compare.jl`, whose run count is unchanged; its header lists what it computes), checked this
with more runs:

- **20 000 runs per package, new seeds:** the means are 667.1 ± 7.2 (ClonalEvolutionTx) and
  677.0 ± 7.3 (BDM_AT), a difference of −1.46% (z = −0.97) with KS `p = 0.46`. The direction of
  the difference reversed. The tails agree: runs with `R ≥ 2000`, `≥ 5000` and `≥ 9000` number
  900 / 267 / 93 and 907 / 266 / 105.
- **Expectation given the genealogy:** at `m = 0` resistance is neutral, so the genealogy does
  not depend on the resistance marks. A cell whose lineage has `G` divisions is then resistant
  with probability `1 − (1 − ν)^G`, which makes `E[R | genealogy] = Σᵢ (1 − (1 − ν)^{Gᵢ})`. This
  quantity is free of the Luria–Delbrück noise. It gives 672.5 ± 0.4 and 673.0 ± 0.4, which agree
  to 0.07%, and both raw means match it within one standard error.
- **Neutral burden:** the mean neutral burden at detection (34.79 and 34.82) equals `μ` times the
  mean number of lineage divisions (34.82 and 34.84), as v1's `Poisson(μ)` per daughter requires.
- **Restarts:** the mean number of restarts before detection (0.989 and 0.994) matches the
  birth–death value of 1.00004 (the script prints it with five decimals).

**Speed.** Wall time per run is the growth from one cell (restarts included) plus treatment, from
a second pass over the same seeds without measurements. After one warm-up run, it is the total
divided by 200.

| configuration | ClonalEvolutionTx (tree mode) | BirthDeathMutation_AT | speed-up |
|:--|---:|---:|---:|
| (A) MTD, `α = 0.8` | 4.1 ms | 13.0 ms | 3.2× |
| (B) `AdaptiveOnOff`, `α = 0.8`, `t_max = 100` (≈ 58 cycles) | 148 ms | 624 ms | 4.2× |
| (C) MTD, `α = 0` | 9.7 ms | 45.7 ms | 4.7× |

At these sizes (`N ≤ 2·10⁴`) the predecessor's O(N) `deleteat!` per death is still cheap, so the
speed-up is modest. It is largest in (C). There, after the sensitive cells die out, the
predecessor's bound still includes the extinct wild type's drug death
(`1 + 0.5 + 4 + lN` against the true `1 + 0.5 + lN`), so, from the rates, at most about 27% of
its ticks can be accepted. ClonalEvolutionTx's bound over the living clones is exact there and
accepts every tick (these acceptance rates were not measured).

### Known semantic differences

None of these differences changes the distribution of the statistics above. They matter when
porting analyses.

- **Stabilising `t_max` origin:** BirthDeathMutation_AT's `MaintenanceTreatment` measures `t_max`
  from the start of the maintenance phase; `Stabilising` measures it from the schedule start, like
  every other schedule. This is not exercised here.
- **Extinct clones in the bound:** the predecessor's thinning bound is a maximum over every
  subclone ever created, including extinct ones. This wastes ticks (see Speed) but does not bias
  the results.
- **Clock at a `t_max` end:** both packages discard the event that would cross `t_max`, so neither
  overshoots. The predecessor leaves the clock at the last tick before `t_max`; ClonalEvolutionTx
  sets it to exactly `t_max`. The state is the same; the reported end time differs by less than one
  inter-tick interval (of order 10⁻⁵ in (B)).
- **Restarts:** by default, the predecessor's `runsimulation` reruns every run whose *final*
  population is empty, including a cure under MTD. Unless `returnextinct = true`, its treatment
  outcomes are therefore conditioned on non-extinction. `restart_on_extinction` restarts only the
  block it is set on (here the growth block).
- **Min rule:** `SelectionFixedAlpha` and `SelectionPredefined` use the drawn `α` as it is; only
  `SelectionDistribution` applies `min(α, α_parent)`. `Resistance` always applies the min rule.
  The two are identical for a fixed `α`, as here.
- **Measurement timing:** the predecessor's `AtTime(T)` holds the state at `T`, but it is stamped
  with the time of the first tick after `T`; ClonalEvolutionTx's is stamped exactly `T`. The
  predecessor's trajectory grid starts at absolute time 0, not at detection, and its
  `AtTreatmentToggle` fires at the start of every on and off phase. ClonalEvolutionTx's
  `AtBlockEnd` fires at the end of every block.
