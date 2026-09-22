#!/bin/bash
# The whole 2000-run campaign as ONE submission.
#
#   sbatch submit_campaign.sbatch
#
# A 3-task job array. Each task is one NODE holding 3 GPUs and 36 cores, and
# runs 3 lanes side by side -- one lane per GPU, 12 cores each -- so 9 runs are
# in flight at once.
#
# Why 3 GPUs per node rather than 9 single-GPU jobs: the gpu QOS caps a user at
# 4 nodes (and 16 GPUs, 192 cores, 12 jobs). Single-GPU jobs land wherever a
# GPU frees up and spread across nodes, so they hit the node cap long before
# the GPU cap. Packing 3 per node is the only way to get 9 GPUs inside it.
# The cost is scheduling: a task has to wait for a node with 3 free GPUs AND 36
# free cores, rather than any one free GPU. A 4-GPU request once waited 25-84 h.
#
# Lanes. There are NLANES = 3 nodes x 3 GPUs = 9. The lane on GPU g of node n
# is  lane = 3n + g  and owns exactly the rows with  row_index % 9 == lane.
# Which runs are finished is decided on disk (clusters.json), so re-striding
# from the earlier 3-lane layout loses nothing: any lane can resume any row.
#
# Wall clock is capped at 36 h by the QOS, so the campaign runs as a chain of
# segments. Each NODE resubmits itself when its three lanes are done --
# immediately, not after the other nodes -- as a one-task job for itself only:
#
#   sbatch --array=<node> --dependency=afterany:<this job> submit_campaign.sbatch
#
# Safe without coordination because node n only ever owns lanes 3n..3n+2 and
# only ever submits node n. afterany means the successor cannot start until
# this job has fully exited, so a node's lanes never run twice at once.
#
# The chain stops by itself when any of these holds (reason is logged):
#   - the node's three lanes have no runs left
#   - the segment completed nothing on any of its lanes (crash-loop guard:
#     stops a node failing every run from burning 36 h of allocation at a time)
#   - CHAIN_MAX segments have run (default 20)
#   - the file  MD/ht2000/STOP  exists  (touch it to stop every node cleanly
#     after its current segment; delete it and resubmit to carry on)
#   - another job for the same node is already queued (e.g. a manual sbatch)
#   - SLURM did not hand over exactly 3 GPUs (refuses to guess devices)
#
#   sbatch submit_campaign.sbatch          # start all three nodes once
#   sbatch --array=1 submit_campaign.sbatch  # restart just node 1
#   touch  MD/ht2000/STOP                  # stop the chain after this segment
#   cat    MD/ht2000/logs/campaign/chain.log
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
#SBATCH --cpus-per-task=36
#SBATCH --gres=gpu:3
# Host memory only; GPU memory is not separately requestable on Midway3. One
# GROMACS run of these box sizes needs ~2 GB, so 64G for three is generous
# while leaving most of a node's 180 GB to other jobs.
#SBATCH --mem=64G
#SBATCH --time=36:00:00
#SBATCH --output=/project/chibueze/stanley/projects/li-na-md/MD/ht2000/logs/slurm/ht2000_%A_%a.out

set -uo pipefail

# Absolute, because sbatch copies the script and $BASH_SOURCE is then a
# path in SLURM's spool directory rather than the repo.
HT=/project/chibueze/stanley/projects/li-na-md/MD/ht2000
SCRIPTS="$HT/scripts"
PY=/scratch/midway3/eshiemogie/moleng/bin/python
mkdir -p "$HT/logs/slurm" "$HT/logs/campaign"
CHAINLOG="$HT/logs/campaign/chain.log"

GPN=3                          # GPUs (= lanes) per node
NNODES=3                       # array size
NLANES=$(( GPN * NNODES ))     # 9
NODE="$SLURM_ARRAY_TASK_ID"
CHAIN_N="${CHAIN_N:-1}"
CHAIN_MAX="${CHAIN_MAX:-20}"
SELF="${SLURM_ARRAY_JOB_ID}_${SLURM_ARRAY_TASK_ID}"
LANES=()
for (( g = 0; g < GPN; g++ )); do LANES+=( $(( GPN * NODE + g )) ); done

echo "========================================================="
echo "ht2000 campaign -- node ${NODE} of 0-$((NNODES-1)), segment ${CHAIN_N}"
echo "host    : $(hostname)"
echo "lanes   : ${LANES[*]} of 0-$((NLANES-1))"
echo "visible : ${CUDA_VISIBLE_DEVICES:-UNSET}  (node has $(nvidia-smi -L | wc -l) devices)"
echo "cpus    : $(grep Cpus_allowed_list /proc/self/status | awk '{print $2}')"
echo "start   : $(date)"
echo "========================================================="

# ------------------------------------------------------------ GPUs
# Nodes in this partition do not isolate GPUs per job: every device on the
# node is visible, including other users'. So the ONLY safe device ids are the
# ones SLURM put in CUDA_VISIBLE_DEVICES. Never pick ids ourselves.
IFS=, read -ra GPUS <<< "${CUDA_VISIBLE_DEVICES:-}"
if (( ${#GPUS[@]} != GPN )); then
    msg="expected $GPN GPUs from SLURM, got '${CUDA_VISIBLE_DEVICES:-UNSET}'"
    echo "FATAL: $msg -- not running, not chaining"
    echo "$(date -u +%FT%TZ) node=$NODE job=$SELF seg=$CHAIN_N STOP: $msg" >> "$CHAINLOG"
    exit 1
fi

# ------------------------------------------------------------ cores
# Split the job's cores into GPN equal contiguous blocks, one per lane, and
# confine each lane to its block with taskset. Without this the three mdruns
# would each start 12 threads over all 36 cores and contend for them.
expand_cpus () {        # "0-3,8,10-11" -> "0 1 2 3 8 10 11"
    local part out=()
    IFS=, read -ra parts <<< "$1"
    for part in "${parts[@]}"; do
        if [[ $part == *-* ]]; then
            for (( c = ${part%-*}; c <= ${part#*-}; c++ )); do out+=( "$c" ); done
        else
            out+=( "$part" )
        fi
    done
    echo "${out[@]}"
}
read -ra CPUS <<< "$(expand_cpus "$(grep Cpus_allowed_list /proc/self/status | awk '{print $2}')")"
PER=$(( ${#CPUS[@]} / GPN ))
if (( PER < 1 )); then
    echo "FATAL: only ${#CPUS[@]} cores for $GPN lanes"; exit 1
fi

# ------------------------------------------------------------ run the lanes
progress () { echo "$HT/logs/campaign/lane_$(printf '%02d' "$1").csv"; }
# grep -c prints 0 AND exits 1 on no match, so "|| echo 0" would print it twice
ok_count () { local f; f="$(progress "$1")"
              [[ -f "$f" ]] || { echo 0; return; }; grep -c ',ok,' "$f" || true; }

declare -A OK_BEFORE
for (( g = 0; g < GPN; g++ )); do
    lane=${LANES[$g]}
    OK_BEFORE[$lane]="$(ok_count "$lane")"
    block="$(IFS=,; echo "${CPUS[*]:$(( g * PER )):$PER}")"
    echo "lane $lane -> GPU ${GPUS[$g]}, cores $block"
    CUDA_VISIBLE_DEVICES="${GPUS[$g]}" \
        taskset -c "$block" \
        bash "$SCRIPTS/campaign_worker.sh" "$lane" "$NLANES" "$PER" &
done
wait

echo "========================================================="
echo "node ${NODE} done: $(date)"

# ------------------------------------------------------------ self-chaining
DONE_HERE=0; LEFT=0; summary=""
for lane in "${LANES[@]}"; do
    d=$(( $(ok_count "$lane") - ${OK_BEFORE[$lane]} ))
    l="$($PY "$SCRIPTS/campaign_status.py" --remaining --lane "$lane" --nlanes "$NLANES" 2>/dev/null)"
    [[ "$l" =~ ^[0-9]+$ ]] || { LEFT=""; break; }
    DONE_HERE=$(( DONE_HERE + d )); LEFT=$(( LEFT + l ))
    summary+=" lane$lane:+$d/${l}left"
done
# any OTHER queued/running job already covering this node?
OTHER="$(squeue -u "$USER" -h -r -n ht2000 -o '%i %K' 2>/dev/null \
         | awk -v n="$NODE" -v me="$SELF" '$2==n && $1!=me' | head -1)"

echo "segment $CHAIN_N/$CHAIN_MAX: completed $DONE_HERE run(s); node $NODE has ${LEFT:-?} left ($summary )"
reason=""
if   [[ -z "$LEFT" ]];                 then reason="could not count remaining runs"
elif (( LEFT == 0 ));                  then reason="node finished -- nothing left"
elif [[ -f "$HT/STOP" ]];              then reason="STOP file present"
elif (( DONE_HERE == 0 ));             then reason="segment completed nothing (crash-loop guard)"
elif (( CHAIN_N >= CHAIN_MAX ));       then reason="reached CHAIN_MAX=$CHAIN_MAX"
elif [[ -n "$OTHER" ]];                then reason="node already queued as ${OTHER%% *}"
fi

if [[ -n "$reason" ]]; then
    echo "chain stopped: $reason"
    echo "$(date -u +%FT%TZ) node=$NODE job=$SELF seg=$CHAIN_N done=$DONE_HERE left=${LEFT:-?}$summary STOP: $reason" >> "$CHAINLOG"
else
    # Absolute path: inside a running job $0 is SLURM's spool copy.
    # afterany: the successor cannot start until this job has fully exited.
    NEXT="$(sbatch --parsable --array="$NODE" \
                  --dependency=afterany:"$SLURM_JOB_ID" \
                  --export=ALL,CHAIN_N=$((CHAIN_N+1)),CHAIN_MAX="$CHAIN_MAX" \
                  "$SCRIPTS/submit_campaign.sbatch" 2>&1)"
    echo "chain: submitted node $NODE segment $((CHAIN_N+1)) as job $NEXT"
    echo "$(date -u +%FT%TZ) node=$NODE job=$SELF seg=$CHAIN_N done=$DONE_HERE left=$LEFT$summary NEXT: $NEXT" >> "$CHAINLOG"
fi
echo "========================================================="
