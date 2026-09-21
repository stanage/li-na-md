#!/bin/bash
# Everything the campaign does NOT do for itself, in one command.
#
#   ./maintain.sh                 refresh indexes, show status, upload to Box
#   ./maintain.sh --prune         ... and free the uploaded trajectories
#   ./maintain.sh --no-box        local bookkeeping only (safe on a compute node)
#   ./maintain.sh --quiet         machine-readable-ish, skip the status report
#
# Run it from a LOGIN NODE. The Box step needs outbound internet, which compute
# nodes do not have -- the same constraint that keeps LigParGen fetching off
# them. Without --no-box on a compute node the archive step reports that and
# the rest still runs.
#
# Safe to run at any time, including while the campaign is going: every step
# only reads run directories or writes files the lanes never touch.
set -uo pipefail

HERE="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
HT="$(dirname "$HERE")"
PY=/scratch/midway3/eshiemogie/moleng/bin/python

PRUNE=0; BOX=1; QUIET=0
for a in "$@"; do
    case "$a" in
        --prune)   PRUNE=1 ;;
        --no-box)  BOX=0 ;;
        --quiet)   QUIET=1 ;;
        -h|--help) sed -n '2,16p' "${BASH_SOURCE[0]}"; exit 0 ;;
        *) echo "unknown option: $a"; exit 2 ;;
    esac
done

# One at a time. Two concurrent runs would have two collect_results.py
# processes writing index/results.csv, which is how you get a torn CSV.
LOCK="$HT/logs/.maintain.lock"
mkdir -p "$HT/logs"
exec 9>"$LOCK"
if ! flock -n 9; then
    echo "another maintain.sh is running (lock: $LOCK) -- nothing to do"
    exit 0
fi

step () { echo; echo "=== $* ==="; }
rc=0

step "1/4  refresh run manifest (stage column)"
$PY "$HERE/build_run_manifest.py" | tail -6 || rc=1

step "2/4  refresh results.csv"
$PY "$HERE/collect_results.py" | tail -20 || rc=1

if (( QUIET == 0 )); then
    step "3/4  campaign status"
    $PY "$HERE/campaign_status.py" || rc=1
else
    step "3/4  campaign status (skipped, --quiet)"
fi

if (( BOX == 1 )); then
    if (( PRUNE == 1 )); then
        step "4/4  archive to Box and free local trajectories"
        $PY "$HERE/archive_to_box.py" --prune || rc=1
    else
        step "4/4  archive to Box (upload only; add --prune to free space)"
        $PY "$HERE/archive_to_box.py" || rc=1
    fi
else
    step "4/4  Box archive skipped (--no-box)"
fi

step "done"
echo "scratch: $(du -sh "$(readlink -f "$HT/runs")" 2>/dev/null | cut -f1) in $(ls -d "$HT"/runs/el* 2>/dev/null | wc -l) run dirs"
echo "project: $(quota -s 2>/dev/null | awk '/^chibueze/ && /blocks/ {print $4" of "$5; exit}')"
exit $rc
