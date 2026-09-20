#!/usr/bin/env python
"""Smoke-test every ion in ff_ions/ by running it through gmx grompp.

For each anion, packs one cation + one anion into a small box and asks GROMACS
to build a .tpr. This catches malformed [atomtypes] rows, missing parameters,
bad section ordering and atom-count mismatches before any of it reaches the
2000-run campaign.

Usage:  python validate_ions.py [--cation Li|Na] [--keep]
"""
from __future__ import annotations

import argparse
import shutil
import subprocess
import sys
import tempfile
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parent))

import common as C  # noqa: E402
from common import FF_IONS, extract_atomtypes  # noqa: E402

EM_MDP = """integrator = steep
nsteps     = 0
nstlist    = 10
cutoff-scheme = Verlet
coulombtype   = PME
rcoulomb   = 0.9
rvdw       = 0.9
pbc        = xyz
"""

PACKMOL_INP = """tolerance 2.5
filetype pdb
output mix.pdb
seed 12345
structure cat.pdb
  number 1
  inside box 1.0 1.0 1.0 29.0 29.0 29.0
end structure
structure an.pdb
  number 1
  inside box 1.0 1.0 1.0 29.0 29.0 29.0
end structure
"""


def run(cmd, cwd, **kw):
    return subprocess.run(cmd, cwd=cwd, capture_output=True, text=True, **kw)


def check(anion: str, cation: str, workroot: Path) -> tuple[bool, str]:
    w = workroot / anion
    w.mkdir(parents=True, exist_ok=True)

    an_itp = (FF_IONS / f"{anion}.itp").read_text()
    cat_itp = (FF_IONS / f"{cation}.itp").read_text()
    an_types, an_body = extract_atomtypes(an_itp)
    cat_types, cat_body = extract_atomtypes(cat_itp)

    (w / "an.itp").write_text(an_body)
    (w / "cat.itp").write_text(cat_body)
    shutil.copy(FF_IONS / f"{anion}.pdb", w / "an.pdb")
    shutil.copy(FF_IONS / f"{cation}.pdb", w / "cat.pdb")

    (w / "system.top").write_text(
        "[ defaults ]\n1  3  yes  0.5  0.5\n\n"
        "[ atomtypes ]\n" + "\n".join(cat_types + an_types) + "\n\n"
        '#include "cat.itp"\n#include "an.itp"\n\n'
        f"[ system ]\n{cation}{anion} test\n\n"
        f"[ molecules ]\n{cation}+   1\nan    1\n")
    (w / "em.mdp").write_text(EM_MDP)
    (w / "pack.inp").write_text(PACKMOL_INP)

    p = run([str(C.PACKMOL)], w, stdin=(w / "pack.inp").open())
    if not (w / "mix.pdb").exists():
        return False, "packmol produced no output"

    # packmol writes a PDB with no box; grompp takes the box from -c, so set one
    g = run(["bash", "-lc",
             f"module load {C.GMX_MODULE} >/dev/null 2>&1; "
             f"{C.GMX} editconf -f mix.pdb -o mix.gro -box 3 3 3 >/dev/null 2>&1 && "
             f"{C.GMX} grompp -f em.mdp -c mix.gro -p system.top -o em.tpr -maxwarn 2"],
            w)
    if (w / "em.tpr").exists():
        return True, ""
    err = [l for l in (g.stdout + g.stderr).splitlines()
           if "ERROR" in l or "Fatal" in l or "not found" in l.lower()]
    return False, " | ".join(err[:3]) or (g.stderr.strip().splitlines() or ["?"])[-1]


def main() -> int:
    ap = argparse.ArgumentParser()
    ap.add_argument("--cation", default="Li", choices=["Li", "Na"])
    ap.add_argument("--keep", action="store_true")
    args = ap.parse_args()

    anions = sorted(p.stem for p in FF_IONS.glob("*.itp")
                    if p.stem not in ("Li", "Na"))
    tmp = Path(tempfile.mkdtemp(prefix="ionval_"))
    bad = []
    for a in anions:
        ok, msg = check(a, args.cation, tmp)
        print(f"  {'PASS' if ok else 'FAIL'}  {args.cation}{a:<6s} {msg}")
        if not ok:
            bad.append(a)
    print(f"\n{len(anions) - len(bad)}/{len(anions)} passed grompp with {args.cation}+")
    if args.keep:
        print(f"work dir kept: {tmp}")
    else:
        shutil.rmtree(tmp, ignore_errors=True)
    if bad:
        print(f"FAILED: {bad}")
        return 1
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
