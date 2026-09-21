#!/bin/bash
# One lane of the campaign: build, run and analyse a strided slice of the runs.
#
#   campaign_worker.sh <lane> <n_lanes> [omp_threads]
#
# Lanes are the unit of parallelism. Each owns one GPU and walks the rows of
# index/run_manifest.csv where (row_index mod n_lanes) == lane, so the 1000
# runs spread evenly over every lane with no coordination and no shared work
# queue to corrupt.
#
# ht1000 runs both of its lanes inside a single 2-GPU job rather than as
# separate array tasks, so the caller (submit_campaign.sbatch) has already
# narrowed CUDA_VISIBLE_DEVICES to this lane's one device before invoking us.
# Nothing here touches that variable.
#
# Everything here is restartable. A lane skips any run that already has
# clusters.json, and run_md.sh itself skips finished stages, so a job killed
# at the wall clock resumes mid-run on resubmission rather than starting over.
set -uo pipefail

LANE="${1:?usage: campaign_worker.sh <lane> <n_lanes> [omp]}"
NLANES="${2:?usage: campaign_worker.sh <lane> <n_lanes> [omp]}"
OMP="${3:-12}"

HT="$(cd "$(dirname "${BASH_SOURCE[0]}")/.." && pwd)"
PY=/scratch/midway3/eshiemogie/moleng/bin/python
MANIFEST="$HT/index/run_manifest.csv"
PROGRESS="$HT/logs/campaign/lane_$(printf '%02d' "$LANE").csv"

# run_md.sh reads its OpenMP width from this; the lanes share a node, so each
# must take a slice of the cores rather than the whole allocation.
export SLURM_CPUS_PER_TASK="$OMP"

mkdir -p "$HT/logs/campaign"
[[ -f "$PROGRESS" ]] || echo "utc,row_index,key,status,seconds,stage" > "$PROGRESS"

# Stop before SLURM kills us mid-stage. Graceful exit keeps the logs readable
# and leaves the run dir in a state the next submission picks straight up.
RESERVE_S="${HT_RESERVE_S:-2700}"          # 45 min
seconds_left () {
    [[ -n "${SLURM_JOB_ID:-}" ]] || { echo 999999; return; }
    local end
    end="$(scontrol show job "$SLURM_JOB_ID" -o 2>/dev/null \
           | grep -o 'EndTime=[^ ]*' | head -1 | cut -d= -f2)"
    [[ -n "$end" ]] || { echo 999999; return; }
    echo $(( $(date -d "$end" +%s) - $(date +%s) ))
}

log () { echo "[lane $LANE $(date -u +%H:%M:%S)] $*"; }

mapfile -t ROWS < <($PY - "$MANIFEST" "$LANE" "$NLANES" <<'PYEOF'
import csv, sys
path, lane, n = sys.argv[1], int(sys.argv[2]), int(sys.argv[3])
for r in csv.DictReader(open(path)):
    if int(r["row_index"]) % n == lane:
        print("\t".join((r["row_index"], r["key"], r["salt"],
                         r["solvent"], r["concentration"], r["solvent_ff"])))
PYEOF
)

log "assigned ${#ROWS[@]} runs (lane $LANE of $NLANES), omp=$OMP, GPU=${CUDA_VISIBLE_DEVICES:-all}"

done_n=0; fail_n=0; skip_n=0
for row in "${ROWS[@]}"; do
    IFS=$'\t' read -r IDX KEY SALT SOLV CONC HAVE_FF <<< "$row"
    RUNDIR="$HT/runs/$KEY"
    LOG="$RUNDIR/pipeline.log"

    if [[ -f "$RUNDIR/clusters.json" ]]; then
        skip_n=$((skip_n+1)); continue
    fi
    if [[ "$HAVE_FF" != "True" ]]; then
        echo "$(date -u +%FT%TZ),$IDX,$KEY,no_force_field,0," >> "$PROGRESS"
        skip_n=$((skip_n+1)); continue
    fi

    left="$(seconds_left)"
    if (( left < RESERVE_S )); then
        log "only ${left}s of wall left -- stopping cleanly before row $IDX"
        break
    fi

    t0=$SECONDS
    mkdir -p "$RUNDIR"
    {
        echo "========================================================="
        echo "row $IDX   key $KEY"
        echo "salt $SALT"
        echo "solvent $SOLV"
        echo "concentration $CONC M"
        echo "node $(hostname)   lane $LANE   GPU ${CUDA_VISIBLE_DEVICES:-?}   omp $OMP"
        echo "start $(date)"
        echo "========================================================="
    } >> "$LOG"

    if [[ ! -f "$RUNDIR/mixture.pdb" ]]; then
        log "setup $KEY (row $IDX)"
        echo "---> [SETUP] packmol + topology" >> "$LOG"
        $PY "$HT/scripts/setup_run.py" --salt "$SALT" --solvent "$SOLV" \
            --molarity "$CONC" >> "$LOG" 2>&1
    fi

    if [[ -f "$RUNDIR/mixture.pdb" ]]; then
        log "run $KEY (row $IDX)"
        bash "$HT/scripts/run_md.sh" "$RUNDIR" >> "$LOG" 2>&1
    fi

    dt=$((SECONDS-t0))
    if [[ -f "$RUNDIR/clusters.json" ]]; then
        status=ok; done_n=$((done_n+1))
    else
        status=FAILED; fail_n=$((fail_n+1))
    fi
    # last stage whose .gro exists, so a failure says how far it got
    stage=none
    for s in em nvt_heat npt prod; do [[ -f "$RUNDIR/$s.gro" ]] && stage=$s; done

    echo "PIPELINE_$status  stage=$stage  ${dt}s  end $(date)" >> "$LOG"
    echo "$(date -u +%FT%TZ),$IDX,$KEY,$status,$dt,$stage" >> "$PROGRESS"
    log "$status $KEY stage=$stage ${dt}s"
done

log "lane finished: ok=$done_n failed=$fail_n skipped=$skip_n"
