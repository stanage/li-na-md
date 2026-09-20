#!/bin/bash
# Run a list of prepared run dirs on the local machine's GPU, N at a time.
#
#   run_local_batch.sh <keyfile> [concurrency] [omp_threads]
#
# Used when we already hold an interactive GPU allocation and do not want to
# wait behind the gpu partition queue. run_md.sh is stage-resumable, so a run
# interrupted here can be finished later either locally or via sbatch.
set -uo pipefail
KEYFILE="${1:?usage: run_local_batch.sh <keyfile> [concurrency] [omp]}"
CONC="${2:-3}"
OMP="${3:-8}"
HT="$(cd "$(dirname "${BASH_SOURCE[0]}")/.." && pwd)"

export SLURM_CPUS_PER_TASK="$OMP"
mkdir -p "$HT/logs"

run_one () {
    local key="$1"
    if [[ -f "$HT/runs/$key/solvation.json" ]]; then
        echo "[skip] $key already complete"; return 0
    fi
    echo "[start] $key $(date +%H:%M:%S)"
    bash "$HT/scripts/run_md.sh" "$HT/runs/$key" > "$HT/logs/local_$key.log" 2>&1
    local rc=$?
    if [[ -f "$HT/runs/$key/solvation.json" ]]; then
        echo "[done]  $key $(date +%H:%M:%S)"
    else
        echo "[FAIL]  $key rc=$rc $(date +%H:%M:%S) -- see logs/local_$key.log"
    fi
    return 0
}
export -f run_one
export HT OMP

grep -v '^\s*$' "$KEYFILE" | xargs -P "$CONC" -I{} bash -c 'run_one "$@"' _ {}
echo "LOCAL_BATCH_COMPLETE $(date +%H:%M:%S)"
