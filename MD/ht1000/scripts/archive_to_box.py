#!/usr/bin/env python
"""Archive finished trajectories to UChicago Box, then free the local copy.

Why this exists. prod.xtc is ~95% of the campaign's footprint -- for ht1000's
1000 runs at a median 11.6k atoms that is roughly 50 GB, against maybe 40 GB
of headroom in the pi-chibueze /project quota, and ht2000 is competing for the
same space. Runs therefore live on scratch (2 TB limit) and their trajectories
are pushed to Box, which is effectively unlimited.

The remote below is ht1000's own folder. ht2000 writes to a sibling path and
the two never share a destination directory, so neither campaign can overwrite
or prune the other's trajectories.

Why it cannot be part of run_md.sh. Compute nodes have no outbound internet --
the same constraint that forces the LigParGen fetches onto a login node. So
this is a separate pass, run from a LOGIN NODE while the campaign continues on
the GPU nodes. It only touches runs that are already fully analysed, so it is
safe to run repeatedly at any point.

What gets archived: prod.xtc and prod.tpr. That pair is what re-running
analyse/cluster needs; everything else a run produces is small enough to keep.

Deleting is opt-in and verified. --prune removes the local file only after
`rclone check` confirms the remote copy matches by hash. Without --prune this
uploads and leaves everything in place.

Usage (login node):
  python archive_to_box.py                 # upload only
  python archive_to_box.py --prune         # upload, verify, then free space
  python archive_to_box.py --dry-run
"""
from __future__ import annotations

import argparse
import subprocess
import sys
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parent))

from common import DONE_MARKER, RUNS  # noqa: E402

RCLONE = Path("/scratch/midway3/eshiemogie/myrcc/bin/rclone")
#: Shared-lab destination. Spaces are fine -- rclone is invoked through
#: subprocess with a list, so nothing is shell-split.
REMOTE = ("box:Amanchukwu Lab Shared Box Folder/1. Lab Members/"
          "1. Grad Students/17. Stanley/Projects/li-na-md/ht1000/runs")
#: the pair needed to redo any analysis; everything else is KB-scale
ARCHIVE = ("prod.xtc", "prod.tpr")


def rclone(*args, check: bool = True) -> subprocess.CompletedProcess:
    return subprocess.run([str(RCLONE), *args], capture_output=True,
                          text=True, check=check)


def human(n: int) -> str:
    for unit in ("B", "KB", "MB", "GB", "TB"):
        if n < 1024 or unit == "TB":
            return f"{n:.1f} {unit}"
        n /= 1024.0
    return str(n)


def main() -> int:
    ap = argparse.ArgumentParser()
    ap.add_argument("--remote", default=REMOTE)
    ap.add_argument("--prune", action="store_true",
                    help="delete the local file after the remote copy verifies")
    ap.add_argument("--dry-run", action="store_true")
    ap.add_argument("--limit", type=int, help="stop after N runs")
    ap.add_argument("--transfers", type=int, default=4)
    args = ap.parse_args()

    if not RCLONE.exists():
        print(f"FATAL: rclone not found at {RCLONE}")
        return 1
    # A compute node will fail here rather than silently doing nothing.
    probe = rclone("lsd", args.remote, check=False)
    if probe.returncode != 0:
        print("FATAL: cannot reach Box. Run this from a LOGIN node -- compute "
              "nodes have no outbound internet.")
        print(probe.stderr.strip()[:400])
        return 1

    todo, skipped, total_bytes = [], 0, 0
    for d in sorted(RUNS.glob("el*")):
        if not (d / DONE_MARKER).exists():
            skipped += 1
            continue
        files = [f for f in ARCHIVE if (d / f).exists()]
        if not files:
            continue
        todo.append((d, files))
        total_bytes += sum((d / f).stat().st_size for f in files)
        if args.limit and len(todo) >= args.limit:
            break

    print(f"{len(todo)} run(s) ready to archive ({human(total_bytes)}); "
          f"{skipped} not yet analysed")
    if not todo:
        return 0
    if args.dry_run:
        for d, files in todo[:10]:
            print(f"  would upload {d.name}: {', '.join(files)}")
        if len(todo) > 10:
            print(f"  ... and {len(todo) - 10} more")
        return 0

    freed = failed = 0
    for i, (d, files) in enumerate(todo, 1):
        dest = f"{args.remote}/{d.name}"
        ok = True
        for f in files:
            up = rclone("copy", str(d / f), dest,
                        "--transfers", str(args.transfers), check=False)
            if up.returncode != 0:
                print(f"[{i}/{len(todo)}] {d.name}: upload FAILED for {f}")
                print("   ", up.stderr.strip()[:300])
                ok = False
                break
        if not ok:
            failed += 1
            continue

        # verify before anything is deleted; --size-only would not catch
        # a truncated transfer, so let rclone hash both sides
        chk = rclone("check", str(d), dest, "--include",
                     "{" + ",".join(files) + "}", check=False)
        if chk.returncode != 0:
            print(f"[{i}/{len(todo)}] {d.name}: VERIFY FAILED, keeping local copy")
            failed += 1
            continue

        msg = "uploaded+verified"
        if args.prune:
            for f in files:
                if f == "prod.tpr":
                    continue          # 280 KB, and grompp inputs are gone
                p = d / f
                freed += p.stat().st_size
                p.unlink()
            msg += ", local prod.xtc removed"
        print(f"[{i}/{len(todo)}] {d.name}: {msg}")

    print(f"\ndone. {len(todo) - failed} archived, {failed} failed"
          + (f", {human(freed)} freed" if args.prune else ""))
    if not args.prune and todo:
        print("run again with --prune to free the local trajectories")
    # NOT shutil.disk_usage: that reports the whole shared scratch filesystem
    # (hundreds of TB), which says nothing about this campaign or the 2 TB
    # per-user quota it actually has to fit inside.
    local = sum(f.stat().st_size for d in RUNS.glob("el*")
                for f in d.iterdir() if f.is_file())
    print(f"campaign on scratch: {human(local)} in "
          f"{len(list(RUNS.glob('el*')))} run dirs")
    return 1 if failed else 0


if __name__ == "__main__":
    raise SystemExit(main())
