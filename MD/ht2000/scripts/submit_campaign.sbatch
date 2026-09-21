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
# Wall clock is capped at 36 h by the QOS, so expect to resubmit several times.
# That is fine and by design:
# lanes skip finished runs and run_md.sh skips finished stages, so the exact
# same command picks up where the last one stopped.
#
#   sbatch submit_campaign.sbatch          # start, or continue
#   python campaign_status.py              # progress, failures, ETA
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
echo "ht2000 campaign -- array task ${SLURM_ARRAY_TASK_ID} of 0-2"
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

# CUDA_VISIBLE_DEVICES is left exactly as SLURM set it.
bash "$SCRIPTS/campaign_worker.sh" "$SLURM_ARRAY_TASK_ID" "$NLANES" "$OMP"

echo "========================================================="
echo "array task ${SLURM_ARRAY_TASK_ID} done: $(date)"
echo "resubmit the same script to continue; finished runs are skipped"
echo "========================================================="
