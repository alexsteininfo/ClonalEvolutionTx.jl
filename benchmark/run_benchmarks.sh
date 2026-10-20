#!/usr/bin/env bash
# Run the full benchmark (serial sizes 10^3-10^7, treatment, statistics, and parallel
# throughput at N = 10^5) and keep the output for docs/src/performance.md.
#
# Inputs : THREADS (optional env var, default half the CPU threads) - threads for the
#          parallel mode. Julia is taken from PATH.
# Outputs: .development/benchmark/benchmark-<host>-<date>.txt (gitignored), also on stdout.
# Run    : bash benchmark/run_benchmarks.sh      (from anywhere; may take tens of minutes)
set -euo pipefail

REPO=$(cd "$(dirname "$0")/.." && pwd)
THREADS=${THREADS:-$(( $(sysctl -n hw.ncpu 2>/dev/null || nproc) / 2 ))}
OUTDIR="$REPO/.development/benchmark"
LOG="$OUTDIR/benchmark-$(hostname -s)-$(date +%Y-%m-%d).txt"
mkdir -p "$OUTDIR"

cd "$REPO"
{
    echo "# $(date) on $(hostname), $THREADS threads for the parallel mode"
    echo "## serial growth, treatment and statistics"
    julia --project=. benchmark/benchmark.jl
    echo
    echo "## parallel throughput"
    julia -t "$THREADS" --project=. benchmark/benchmark.jl --parallel 100000
} 2>&1 | tee "$LOG"
echo "saved to $LOG"
