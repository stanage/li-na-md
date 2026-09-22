#!/bin/bash
# The whole 2000-run campaign as ONE submission.
#
#   sbatch submit_campaign.sbatch
#
# A 3-task job array, one GPU each -- 3 runs in flight. One sbatch instead of
# 2000, using 3 of the 4 permitted job slots.
#
# Why one GPU per task rather than a whole node: every GPU in the partition is
# held by single-GPU jobs from other users, so a --gres=gpu:4 request has to
# wait for some node's last job to end -- measured at 25 h minimum, 84 h on one
# node. A gpu:1 request slots into the first GPU that frees anywhere. It also
# leaves the rest of each node available to everyone else.
#
# The cost is throughput: 3 concurrent runs instead of 12, so roughly 330 h of
# wall rather than 85. Raise LANES and --gres together if the partition frees
# up (the QOS ceiling is 4 nodes / 16 GPUs / 12 jobs).
#
# Wall clock is capped at 36 h by the QOS, so the campaign runs as a chain of
# segments. Each lane RESUBMITS ITSELF when it finishes -- immediately, not
# after the other two -- as a one-task job for its own slice only:
#
#   sbatch --array=<lane> --dependency=afterany:<this job> submit_campaign.sbatch
#
# That is safe without coordination because lane k only ever owns rows with
# row_index % 3 == k, and only ever submits lane k. Lanes skip finished runs
# and run_md.sh skips finished stages / resumes from checkpoints, so every
# segment picks up exactly where the last one stopped.
#
# The chain stops by itself when any of these holds (reason is logged):
#   - the lane has no runs left
#   - the segment completed nothing (crash-loop guard: stops a lane that is
#     failing every run from burning the allocation 36 h at a time)
#   - CHAIN_MAX segments have run (default 20; about 2x what 2000 runs need)
#   - the file  MD/ht2000/STOP  exists  (touch it to stop every lane cleanly
#     after its current segment; delete it and resubmit to carry on)
#   - another job for the same lane is already queued (e.g. a manual sbatch)
#
#   sbatch submit_campaign.sbatch          # start all three lanes once
#   touch  MD/ht2000/STOP                  # stop the chain after this segment
#   python campaign_status.py              # progress, failures, ETA
#
# NB: sbatch copies THIS file at submit time, so edits here reach the next
# segment only. campaign_worker.sh and run_md.sh are read from disk while
# they run -- never edit those while any lane is running.
#
#SBATCH --job-name=ht2000
#SBATCH --account=pi-chibueze
#SBATCH --partition=gpu
#SBATCH --qos=gpu
#SBATCH --array=0-2
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=12
#SBATCH --gres=gpu:1
# Host memory only; GPU memory is not separately requestable on Midway3, where
# --gres=gpu:1 hands over one whole 24 GB device. A single GROMACS run of these
# box sizes needs ~2 GB, so 32 G is generous -- and asking for a share rather
# than --mem=0 leaves the rest of the node usable by other jobs.
#SBATCH --mem=32G
#SBATCH --time=36:00:00
#SBATCH --output=/project/chibueze/stanley/projects/li-na-md/MD/ht2000/logs/slurm/ht2000_%A_%a.out

set -uo pipefail

# Absolute, because sbatch copies the script and $BASH_SOURCE is then a
# path in SLURM's spool directory rather than the repo.
HT=/project/chibueze/stanley/projects/li-na-md/MD/ht2000
SCRIPTS="$HT/scripts"
mkdir -p "$HT/logs/slurm" "$HT/logs/campaign"

# One lane per array task, one GPU per lane.
NLANES=3
OMP=${SLURM_CPUS_PER_TASK:-12}

echo "========================================================="
echo "ht2000 campaign -- lane ${SLURM_ARRAY_TASK_ID} of 0-2, segment ${CHAIN_N:-1}"
echo "node    : $(hostname)"
echo "lane    : ${SLURM_ARRAY_TASK_ID} of 0-$((NLANES-1)), omp=$OMP"
# Nodes in this partition do not isolate GPUs per job, so the devices visible
# here are the node's, not ours. Trust SLURM's CUDA_VISIBLE_DEVICES and say so
# loudly if it is missing -- picking a device ourselves would risk landing on
# one allocated to another user.
echo "visible : ${CUDA_VISIBLE_DEVICES:-UNSET}  (node has $(nvidia-smi -L | wc -l) devices)"
if [[ -z "${CUDA_VISIBLE_DEVICES:-}" ]]; then
    echo "WARNING: SLURM did not set CUDA_VISIBLE_DEVICES; gromacs may pick a"
    echo "         GPU belonging to another job on this node."
fi
echo "start   : $(date)"
echo "========================================================="

LANE="$SLURM_ARRAY_TASK_ID"
PROGRESS="$HT/logs/campaign/lane_$(printf '%02d' "$LANE").csv"
# grep -c prints 0 AND exits 1 on no match, so "|| echo 0" would print it twice
ok_count () { [[ -f "$PROGRESS" ]] || { echo 0; return; }; grep -c ',ok,' "$PROGRESS" || true; }
OK_BEFORE="$(ok_count)"

# CUDA_VISIBLE_DEVICES is left exactly as SLURM set it.
bash "$SCRIPTS/campaign_worker.sh" "$LANE" "$NLANES" "$OMP"

echo "========================================================="
echo "array task ${LANE} done: $(date)"

# ------------------------------------------------------------ self-chaining
PY=/scratch/midway3/eshiemogie/moleng/bin/python
CHAIN_N="${CHAIN_N:-1}"
CHAIN_MAX="${CHAIN_MAX:-20}"
DONE_HERE=$(( $(ok_count) - OK_BEFORE ))
LEFT="$($PY "$SCRIPTS/campaign_status.py" --remaining --lane "$LANE" --nlanes "$NLANES" 2>/dev/null)"
SELF="${SLURM_ARRAY_JOB_ID}_${SLURM_ARRAY_TASK_ID}"
# any OTHER queued/running job already covering this lane?
OTHER="$(squeue -u "$USER" -h -r -n ht2000 -o '%i %K' 2>/dev/null \
         | awk -v l="$LANE" -v me="$SELF" '$2==l && $1!=me' | head -1)"

echo "segment $CHAIN_N/$CHAIN_MAX: completed $DONE_HERE run(s); lane $LANE has ${LEFT:-?} left"
CHAINLOG="$HT/logs/campaign/chain.log"
reason=""
if   [[ -z "$LEFT" ]];                 then reason="could not count remaining runs"
elif (( LEFT == 0 ));                  then reason="lane finished -- nothing left"
elif [[ -f "$HT/STOP" ]];              then reason="STOP file present"
elif (( DONE_HERE == 0 ));             then reason="segment completed nothing (crash-loop guard)"
elif (( CHAIN_N >= CHAIN_MAX ));       then reason="reached CHAIN_MAX=$CHAIN_MAX"
elif [[ -n "$OTHER" ]];                then reason="lane already queued as ${OTHER%% *}"
fi

if [[ -n "$reason" ]]; then
    echo "chain stopped: $reason"
    echo "$(date -u +%FT%TZ) lane=$LANE job=$SELF seg=$CHAIN_N done=$DONE_HERE left=${LEFT:-?} STOP: $reason" >> "$CHAINLOG"
else
    # Absolute path: inside a running job $0 is SLURM's spool copy.
    # afterany: the successor cannot start until this job has fully exited.
    NEXT="$(sbatch --parsable --array="$LANE" \
                  --dependency=afterany:"$SLURM_JOB_ID" \
                  --export=ALL,CHAIN_N=$((CHAIN_N+1)),CHAIN_MAX="$CHAIN_MAX" \
                  "$SCRIPTS/submit_campaign.sbatch" 2>&1)"
    echo "chain: submitted lane $LANE segment $((CHAIN_N+1)) as job $NEXT"
    echo "$(date -u +%FT%TZ) lane=$LANE job=$SELF seg=$CHAIN_N done=$DONE_HERE left=$LEFT NEXT: $NEXT" >> "$CHAINLOG"
fi
echo "========================================================="
