#!/usr/bin/env python
"""Give ht1000 its own copy of the solvent force fields it shares with ht2000.

ht1000's solvent pool is the oxygen-bearing subset of ht2000's library, so 272
of its 278 solvents are molecules ht2000 has already fetched from LigParGen and
composition-verified. Asking the Yale server for them a second time would take
~45 minutes and would re-run the one failure mode that scares us most: under
load LigParGen has been observed returning *another molecule's* parameters
(PLAN.md 3.1). A verified local copy has neither cost nor risk.

What this does NOT do is let the two campaigns share files. Every byte is
copied into ht1000/ff_solvents and re-checked against its SMILES here; after
this runs, nothing in ht1000 reads anything under ht2000 at any later point,
and deleting ht2000 entirely would not affect a single run. The source is
opened read-only and never written to.

Anything the source does not have is left for fetch_solvent_ff.py, which pulls
it fresh from LigParGen on a login node.

Usage (safe to re-run; already-present files are left alone):
  python import_solvent_ff.py
  python import_solvent_ff.py --dry-run
  python import_solvent_ff.py --source /path/to/other/ff_solvents
"""
from __future__ import annotations

import argparse
import shutil
import sys
from pathlib import Path

import pandas as pd

sys.path.insert(0, str(Path(__file__).resolve().parent))

from common import (FF_SOLV, MD, SOLVENT_POOL, canon,  # noqa: E402
                    check_solvent_ff, solvent_key)

#: Read-only source. Deliberately not importable as a writable path anywhere
#: else in ht1000.
DEFAULT_SOURCE = MD / "ht2000" / "ff_solvents"

#: A solvent is only usable if all four are present. .lmp is LAMMPS output we
#: do not use, so it is copied when available but never required.
REQUIRED = ("itp", "gro", "pdb", "smi")
OPTIONAL = ("lmp",)


def main() -> int:
    ap = argparse.ArgumentParser()
    ap.add_argument("--source", type=Path, default=DEFAULT_SOURCE)
    ap.add_argument("--pool", type=Path, default=SOLVENT_POOL)
    ap.add_argument("--dry-run", action="store_true")
    args = ap.parse_args()

    if not args.source.is_dir():
        print(f"FATAL: source not found: {args.source}")
        return 1
    if args.source.resolve() == FF_SOLV.resolve():
        print("FATAL: source and destination are the same directory")
        return 1

    pool = pd.read_csv(args.pool)
    FF_SOLV.mkdir(parents=True, exist_ok=True)

    copied = present = absent = bad = 0
    missing: list[str] = []

    for smi in pool["smiles"]:
        c = canon(smi)
        key = solvent_key(c)
        dst_itp = FF_SOLV / f"{key}.itp"

        if all((FF_SOLV / f"{key}.{e}").exists() for e in REQUIRED):
            present += 1
            continue

        src_files = {e: args.source / f"{key}.{e}" for e in REQUIRED + OPTIONAL}
        if not all(src_files[e].exists() for e in REQUIRED):
            absent += 1
            missing.append(c)
            continue

        # Verify at the source, before copying: an .itp that does not describe
        # this SMILES is exactly the LigParGen mix-up we are avoiding, and
        # copying it would launder the error into ht1000.
        ok, why = check_solvent_ff(src_files["itp"].read_text(), c)
        if not ok:
            print(f"  REJECT {key}  {c}\n         {why}")
            bad += 1
            missing.append(c)
            continue

        if not args.dry_run:
            for e, src in src_files.items():
                if src.exists():
                    shutil.copy2(src, FF_SOLV / f"{key}.{e}")
            # copy2 preserves mtime, which would make these look older than the
            # campaign that owns them; stamp them as written now instead.
            for e in REQUIRED + OPTIONAL:
                p = FF_SOLV / f"{key}.{e}"
                if p.exists():
                    p.touch()
        copied += 1

    verb = "would copy" if args.dry_run else "copied"
    print(f"pool           : {len(pool)} solvents")
    print(f"already in ht1000: {present}")
    print(f"{verb:<16} : {copied}   (from {args.source})")
    print(f"not in source  : {absent}")
    if bad:
        print(f"rejected at source (composition mismatch): {bad}")
    if missing:
        out = args.pool.parent / "solvents_to_fetch.txt"
        if not args.dry_run:
            out.write_text("\n".join(missing) + "\n")
            print(f"\n{len(missing)} solvent(s) still need LigParGen -> {out}")
            print(f"  python fetch_solvent_ff.py --smiles-file {out}")
        for m in missing:
            print(f"    {m}")
    # A dry run copies nothing, so the destination is legitimately incomplete.
    if not missing and not args.dry_run:
        print(f"\nht1000/ff_solvents is complete for the pool "
              f"({present + copied}/{len(pool)})")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
