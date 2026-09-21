#!/bin/bash
# Run the four-stage GROMACS protocol for one prepared run directory.
#
#   run_md.sh <run_dir>
#
# Shared by the per-run submit.sh scripts and by the job-array driver, so the
# protocol lives in exactly one place. Resumable: a stage whose .gro already
# exists is skipped, so a job killed by a wall-clock limit picks up where it
# stopped when resubmitted.
set -euo pipefail

# Resolve our own location to an ABSOLUTE path before doing anything else.
# We cd into the run directory below, so a relative $BASH_SOURCE (which is what
# you get from `bash scripts/run_md.sh ...`) would stop resolving afterwards --
# that silently broke the post-run solvation analysis on the first pilot run.
SELF="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)/$(basename "${BASH_SOURCE[0]}")"
SCRIPTS="$(dirname "$SELF")"

RUNDIR="$(cd "${1:?usage: run_md.sh <run_dir>}" && pwd)"
cd "$RUNDIR"

# `module` is a shell function set up by the login profile, and Midway3's
# gromacs modulefile additionally needs site variables (SOFTPATH, CPUTYPE) that
# only the profile defines -- sourcing modules/init/bash by hand is not enough.
# So if we were not started from a login shell, re-exec as one.
if ! type module >/dev/null 2>&1; then
    exec bash -l "$SELF" "$RUNDIR"
fi

# NB: do NOT `module purge` first. Midway3's site module defines SOFTPATH, and
# the gromacs modulefile reads it -- purging makes the load fail with
# "eval set [array get env SOFTPATH] / wrong # args". Verified in a real batch job.
module load gromacs/2025.3
command -v gmx_mpi >/dev/null || { echo "FATAL: gmx_mpi not on PATH"; exit 1; }
export OMP_NUM_THREADS="${SLURM_CPUS_PER_TASK:-8}"
GMX=gmx_mpi
NT="${SLURM_CPUS_PER_TASK:-8}"

echo "run:  $RUNDIR"
echo "node: $(hostname)  GPU: $(nvidia-smi --query-gpu=name --format=csv,noheader 2>/dev/null | head -1)"
echo "start: $(date)"

# Hours of wall clock left, for -maxh. mdrun then checkpoints and exits
# cleanly at 99% of it instead of being killed mid-write. Outside SLURM there
# is no deadline, so pass a value large enough to never bind.
maxh () {
    local end left
    [[ -n "${SLURM_JOB_ID:-}" ]] || { echo 999; return; }
    end="$(scontrol show job "$SLURM_JOB_ID" -o 2>/dev/null \
           | grep -o 'EndTime=[^ ]*' | head -1 | cut -d= -f2)"
    [[ -n "$end" ]] || { echo 999; return; }
    left=$(( $(date -d "$end" +%s) - $(date +%s) ))
    # 120 s of margin for the checkpoint write and the analysis that follows
    awk -v s="$((left - 120))" 'BEGIN{ printf "%.4f", (s>0 ? s : 0)/3600 }'
}

# NB: no -ntmpi. Midway3's gromacs is built against real MPI, not thread-MPI,
# and `-ntmpi` is a hard error there ("Setting the number of thread-MPI ranks is
# only supported with thread-MPI"). The old campaign's scripts used -ntmpi 1
# because that build was thread-MPI. Running gmx_mpi without mpirun already
# gives the single rank we want; -ntomp sets the OpenMP width.
#
# -maxwarn 0 is deliberate. Measured across 7 runs x 4 stages spanning Beti,
# TFSI, OTs, BNZ, PF6, SCN and NO3, grompp emits zero warnings, so any warning
# here is a real problem -- a non-neutral system under Ewald, a timestep too
# long for a bond, an atom-type override. The pipeline is stage-resumable, so
# stopping costs one rebuild and nothing else. Raise this only with evidence.
#
# -cpi/-maxh make a wall-clock kill recoverable. Without them mdrun is killed
# mid-production and the stage restarts from zero on the next submission,
# which for a 47k-atom run is hours of GPU time lost every time.
stage () {           # stage <name> <input.gro> <extra mdrun flags...>
    local name="$1"; shift
    local input="$1"; shift
    if [[ -f "${name}.gro" ]]; then
        echo "== ${name}: already done, skipping"
        return 0
    fi

    local -a resume=()
    if [[ -f "${name}.cpt" && -f "${name}.tpr" ]]; then
        # Interrupted part-way: continue from the checkpoint and keep the
        # existing .tpr, so the trajectory stays one continuous run.
        echo "== ${name}: resuming from checkpoint"
        resume=(-cpi "${name}.cpt" -append)
    else
        echo "== ${name}"
        $GMX grompp -f "${name}.mdp" -c "${input}" -p system.top \
             -o "${name}.tpr" -maxwarn 0
    fi
    $GMX mdrun -deffnm "${name}" -ntomp "${NT}" -maxh "$(maxh)" "${resume[@]}" "$@"

    if [[ ! -f "${name}.gro" ]]; then
        echo "== ${name}: stopped before finishing (wall clock); checkpoint kept"
        exit 0        # not a failure -- resubmitting resumes from ${name}.cpt
    fi
}

stage em       mixture.gro  -nb gpu -bonded cpu -pme cpu -v
stage nvt_heat em.gro       -nb gpu -bonded gpu -pme gpu -update gpu
stage npt      nvt_heat.gro -nb gpu -bonded gpu -pme gpu -update gpu
stage prod     npt.gro      -nb gpu -bonded gpu -pme gpu -update gpu

set +e
echo "=== NPT density ==="
echo Density | $GMX energy -f npt.edr -o density.xvg 2>&1 | tail -3
echo "=== Performance ==="
grep Performance prod.log

PYBIN=/scratch/midway3/eshiemogie/moleng/bin/python

echo "=== solvation analysis ==="
"$PYBIN" "$SCRIPTS/analyze_solvation.py" --run-dir "$RUNDIR" \
    || echo "solvation analysis FAILED"

# Ion clustering is a separate pass because it uses its own cutoff: the first
# minimum of the cation-ANION rdf, not the cation-solvent-O one above.
echo "=== ion cluster analysis ==="
"$PYBIN" "$SCRIPTS/cluster_analysis.py" --run-dir "$RUNDIR" \
    || echo "cluster analysis FAILED"

echo "end: $(date)"
exit 0
