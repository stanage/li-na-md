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

# NB: no -ntmpi. Midway3's gromacs is built against real MPI, not thread-MPI,
# and `-ntmpi` is a hard error there ("Setting the number of thread-MPI ranks is
# only supported with thread-MPI"). The old campaign's scripts used -ntmpi 1
# because that build was thread-MPI. Running gmx_mpi without mpirun already
# gives the single rank we want; -ntomp sets the OpenMP width.
stage () {           # stage <name> <input.gro> <extra mdrun flags...>
    local name="$1"; shift
    local input="$1"; shift
    if [[ -f "${name}.gro" ]]; then
        echo "== ${name}: already done, skipping"
        return 0
    fi
    echo "== ${name}"
    $GMX grompp -f "${name}.mdp" -c "${input}" -p system.top -o "${name}.tpr" -maxwarn 5
    $GMX mdrun -deffnm "${name}" -ntomp "${NT}" "$@"
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

echo "=== solvation analysis ==="
PYBIN=/scratch/midway3/eshiemogie/moleng/bin/python
"$PYBIN" "$SCRIPTS/analyze_solvation.py" --run-dir "$RUNDIR" \
    || echo "solvation analysis FAILED"

echo "end: $(date)"
exit 0
