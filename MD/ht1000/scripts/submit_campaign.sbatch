#!/bin/bash
# The whole 1000-run ht1000 campaign as ONE submission.
#
#   sbatch submit_campaign.sbatch
#
# A single job on ONE node holding TWO GPUs, running two lanes -- 2 runs in
# flight. This is the fourth node slot: the `gpu` QOS caps this user at
# node=4, and ht2000 is already occupying three nodes with its 3-task array.
# Taking one node with gpu:2 here uses the last permitted node without
# displacing that campaign.
#
# Why one node x 2 GPUs rather than a 2-task array of gpu:1 (which is what
# ht2000 does): an array of two gpu:1 tasks can land on two *different* nodes
# and would then consume two of the four node slots, leaving none spare while
# ht2000 holds three. One node with two GPUs is the same throughput for half
# the node budget.
#
# Wall clock is capped at 36 h by the QOS, so expect to resubmit several
# times. That is by design: lanes skip finished runs and run_md.sh skips
# finished stages, so the exact same command starts the campaign and
# continues it.
#
#   sbatch submit_campaign.sbatch          # start, or continue
#   python campaign_status.py --lanes      # progress, failures, ETA
#
#SBATCH --job-name=ht1000
#SBATCH --account=pi-chibueze
#SBATCH --partition=gpu
#SBATCH --qos=gpu
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=24
#SBATCH --gres=gpu:2
# Host memory only; GPU memory is not separately requestable on Midway3, where
# each --gres=gpu unit hands over one whole 24 GB device. Two GROMACS runs of
# these box sizes need ~2 GB of device memory each, so the 64 G here is host
# RAM for two lanes' packmol and analysis -- generous, and asking for a share
# rather than --mem=0 leaves the rest of the 48-core node usable by others.
#SBATCH --mem=64G
#SBATCH --time=36:00:00
#SBATCH --output=/project/chibueze/stanley/projects/li-na-md/MD/ht1000/logs/slurm/ht1000_%j.out

set -uo pipefail

# Absolute, because sbatch copies the script and $BASH_SOURCE is then a path
# in SLURM's spool directory rather than the repo.
HT=/project/chibueze/stanley/projects/li-na-md/MD/ht1000
SCRIPTS="$HT/scripts"
mkdir -p "$HT/logs/slurm" "$HT/logs/campaign"

NLANES=2
# 24 cores over 2 lanes. run_md.sh reads SLURM_CPUS_PER_TASK for -ntomp, and
# campaign_worker.sh re-exports it per lane, so each lane must be told its
# own slice rather than the whole allocation.
OMP=$(( ${SLURM_CPUS_PER_TASK:-24} / NLANES ))

echo "========================================================="
echo "ht1000 campaign -- $NLANES lanes on one node, one GPU each"
echo "job     : ${SLURM_JOB_ID:-none}"
echo "node    : $(hostname)"
echo "visible : ${CUDA_VISIBLE_DEVICES:-UNSET}  (node has $(nvidia-smi -L | wc -l) devices)"
echo "omp     : $OMP per lane"
echo "start   : $(date)"
echo "========================================================="

# The gpu nodes do NOT isolate GPUs per job -- nvidia-smi lists every device on
# the node, including other people's, and CUDA_VISIBLE_DEVICES is the only
# thing telling us which two are ours. ht2000 learned this by trespassing on
# three GPUs belonging to another user (PLAN.md 6, gotcha 4).
#
# Two lanes in one job do have to be split across our two devices, or both
# would land on the first one and halve each other's throughput. The rule that
# keeps that safe: only ever hand a lane a device SLURM already named. We
# subset the allocated list, never extend it, and refuse to run at all if the
# list is missing or shorter than the number of lanes.
if [[ -z "${CUDA_VISIBLE_DEVICES:-}" ]]; then
    echo "FATAL: SLURM did not set CUDA_VISIBLE_DEVICES. Refusing to guess a"
    echo "       device -- on this partition that risks running on a GPU"
    echo "       allocated to another user."
    exit 1
fi
IFS=',' read -r -a DEVS <<< "$CUDA_VISIBLE_DEVICES"
if (( ${#DEVS[@]} < NLANES )); then
    echo "FATAL: asked for $NLANES lanes but SLURM allocated ${#DEVS[@]} GPU(s)"
    echo "       (CUDA_VISIBLE_DEVICES=$CUDA_VISIBLE_DEVICES)."
    exit 1
fi

pids=()
for (( L=0; L<NLANES; L++ )); do
    # Each lane sees exactly one device, taken from our own allocation.
    echo "lane $L -> GPU ${DEVS[$L]}"
    CUDA_VISIBLE_DEVICES="${DEVS[$L]}" \
        bash "$SCRIPTS/campaign_worker.sh" "$L" "$NLANES" "$OMP" \
        > "$HT/logs/campaign/lane_$(printf '%02d' "$L").log" 2>&1 &
    pids+=($!)
done

echo "lanes started: ${pids[*]}"
echo "  tail -f $HT/logs/campaign/lane_00.log"

rc=0
for p in "${pids[@]}"; do
    wait "$p" || rc=1
done

echo "========================================================="
echo "job ${SLURM_JOB_ID:-none} done: $(date)   lane exit status: $rc"
echo "resubmit the same script to continue; finished runs are skipped"
echo "========================================================="
exit $rc
