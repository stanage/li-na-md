#!/bin/bash
# The whole 2000-run campaign as ONE submission.
#
#   sbatch submit_campaign.sbatch
#
# A 3-task job array, one task per node, four lanes per task pinned to the
# node's four GPUs -- 12 runs in flight. That is one sbatch instead of 2000,
# and it sits inside the gpu QOS ceiling (4 nodes / 16 GPUs / 12 jobs) while
# leaving a node free for interactive work.
#
# Wall clock is capped at 36 h by the QOS and the campaign needs roughly 80
# GPU-hours per lane, so expect to resubmit. That is fine and by design:
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
#SBATCH --cpus-per-task=48
#SBATCH --gres=gpu:4
#SBATCH --time=36:00:00
#SBATCH --output=/project/chibueze/stanley/projects/li-na-md/MD/ht2000/logs/slurm/ht2000_%A_%a.out

set -uo pipefail

# Absolute, because sbatch copies the script and $BASH_SOURCE is then a
# path in SLURM's spool directory rather than the repo.
HT=/project/chibueze/stanley/projects/li-na-md/MD/ht2000
SCRIPTS="$HT/scripts"
mkdir -p "$HT/logs/slurm" "$HT/logs/campaign"

LANES_PER_NODE=4
NODES=3
NLANES=$(( LANES_PER_NODE * NODES ))
OMP=$(( ${SLURM_CPUS_PER_TASK:-48} / LANES_PER_NODE ))

echo "========================================================="
echo "ht2000 campaign -- array task ${SLURM_ARRAY_TASK_ID} of 0-2"
echo "node    : $(hostname)"
echo "gpus    : $(nvidia-smi -L | wc -l)"
echo "lanes   : $LANES_PER_NODE here, $NLANES total, omp=$OMP each"
echo "start   : $(date)"
echo "========================================================="

for L in $(seq 0 $((LANES_PER_NODE - 1))); do
    GLOBAL_LANE=$(( SLURM_ARRAY_TASK_ID * LANES_PER_NODE + L ))
    # One GPU per lane. gmx then sees a single device and needs no -gpu_id.
    CUDA_VISIBLE_DEVICES="$L" \
        bash "$SCRIPTS/campaign_worker.sh" "$GLOBAL_LANE" "$NLANES" "$OMP" &
done
wait

echo "========================================================="
echo "array task ${SLURM_ARRAY_TASK_ID} done: $(date)"
echo "resubmit the same script to continue; finished runs are skipped"
echo "========================================================="
