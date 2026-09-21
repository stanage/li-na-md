#!/usr/bin/env python
"""Audit ff_solvents/ -- does every cached .itp actually match its SMILES?

The LigParGen server can return another molecule's files when requests overlap
(observed: a 49-atom solvent came back with a 58-atom topology). Such a file is
perfectly valid GROMACS input, so nothing downstream would complain -- it would
just silently simulate the wrong chemistry. This script re-checks the whole
library by element composition.

    python verify_solvent_ff.py            # report only
    python verify_solvent_ff.py --fix      # delete bad entries and re-fetch

Run --fix on a login node (needs internet).
"""
from __future__ import annotations

import argparse
import sys
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parent))

from common import FF_SOLV, INDEX, check_solvent_ff, solvent_key  # noqa: E402


def main() -> int:
    ap = argparse.ArgumentParser()
    ap.add_argument("--fix", action="store_true",
                    help="delete mismatched entries and re-fetch them serially")
    args = ap.parse_args()

    entries = sorted(FF_SOLV.glob("*.smi"))
    if not entries:
        print("ff_solvents/ is empty")
        return 0

    bad: list[tuple[str, str, str]] = []
    for p in entries:
        smi = p.read_text().strip()
        itp = p.with_suffix(".itp")
        if not itp.exists():
            bad.append((p.stem, smi, "missing .itp"))
            continue
        ok, why = check_solvent_ff(itp.read_text(), smi)
        if not ok:
            bad.append((p.stem, smi, why))

    print(f"checked {len(entries)} solvent force fields: "
          f"{len(entries) - len(bad)} ok, {len(bad)} bad")
    for key, smi, why in bad:
        print(f"  BAD {key}  {smi}\n      {why}")

    if bad and args.fix:
        from ligpargen import fetch, LigParGenError
        print(f"\nre-fetching {len(bad)} entries serially...")
        still_bad = []
        for key, smi, _ in bad:
            prefix = FF_SOLV / key
            for e in ("lmp", "pdb", "itp", "gro", "smi"):
                prefix.with_suffix(f".{e}").unlink(missing_ok=True)
            try:
                fetch(smi, prefix, net_charge=0, verbose=False)
                ok, why = check_solvent_ff(prefix.with_suffix(".itp").read_text(), smi)
            except LigParGenError as e:
                ok, why = False, str(e)[:120]
            if ok:
                prefix.with_suffix(".smi").write_text(smi + "\n")
                print(f"  fixed {key}  {smi}")
            else:
                still_bad.append(key)
                print(f"  STILL BAD {key}: {why}")
        if still_bad:
            print(f"\nunresolved: {still_bad}")
            return 1
        print("\nall entries now verified")
        return 0

    if bad:
        print("\nre-run with --fix (on a login node) to repair these")
        return 1
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
