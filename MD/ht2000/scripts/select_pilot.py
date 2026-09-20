#!/usr/bin/env python
"""Choose the pilot formulations, deterministically.

The pilot is picked to *stress the pipeline*, not at random: it spans both
cations, all four force-field provenances in ff_ions/, monatomic through
21-atom anions, and both ends of the concentration range -- so that anything
structurally broken shows up now rather than at run 1500.

Writes index/pilot_rows.csv and index/pilot_solvents.txt.

Usage:  python select_pilot.py [-n 10]
"""
from __future__ import annotations

import argparse
import sys
from pathlib import Path

import pandas as pd

sys.path.insert(0, str(Path(__file__).resolve().parent))

from common import (DATASET, INDEX, anion_name, canon, cation_name,  # noqa: E402
                    fix_salt, run_key, split_salt)

#: (cation, anion) pairs to cover, in priority order, with why each is here.
PILOT_TARGETS = [
    ("Li", "FSI",  "baseline -- the salt the old 800-run campaign used"),
    ("Na", "FSI",  "same anion, Na+ -- isolates the new cation"),
    ("Li", "TFSI", "large fluorinated imide, 2009IL source"),
    ("Na", "PF6",  "hypervalent inorganic, 2009IL source"),
    ("Li", "Cl",   "monatomic anion -- degenerate case, no bonded terms"),
    ("Na", "OAc",  "carboxylate -- exercises the generated [pairs]"),
    ("Li", "TCM",  "nitrile-rich planar anion"),
    ("Na", "BF4",  "small inorganic, very common salt"),
    ("Li", "Beti", "largest anion (21 atoms), fftool/CL&P source"),
    ("Na", "OTs",  "aromatic sulfonate, LigParGen source"),
]


def main() -> int:
    ap = argparse.ArgumentParser()
    ap.add_argument("-n", type=int, default=10)
    args = ap.parse_args()

    df = pd.read_csv(DATASET)
    df["salt_fixed"] = df["salt"].map(fix_salt)
    cat_an = df["salt_fixed"].map(lambda s: split_salt(s))
    df["cation"] = [cation_name(c) for c, _ in cat_an]
    df["anion"] = [anion_name(a) for _, a in cat_an]

    picked = []
    for i, (cat, an, why) in enumerate(PILOT_TARGETS[:args.n]):
        sub = df[(df.cation == cat) & (df.anion == an)]
        if sub.empty:
            print(f"  !! no row for {cat}{an}")
            continue
        # alternate between the lowest and highest concentration available,
        # so the pilot covers both the big-box and small-box extremes
        sub = sub.sort_values("concentration")
        row = (sub.iloc[0] if i % 2 == 0 else sub.iloc[-1])
        picked.append({
            "key": run_key(row.salt, row.solvent, row.concentration),
            "salt": fix_salt(row.salt),
            "salt_orig": row.salt,
            "solvent": canon(row.solvent),
            "concentration": row.concentration,
            "cation": cat,
            "anion": an,
            "reason": why,
        })

    out = pd.DataFrame(picked)
    INDEX.mkdir(parents=True, exist_ok=True)
    out.to_csv(INDEX / "pilot_rows.csv", index=False)
    (INDEX / "pilot_solvents.txt").write_text(
        "\n".join(sorted(set(out.solvent))) + "\n")

    print(out[["key", "cation", "anion", "concentration", "solvent"]].to_string(index=False))
    print(f"\n{len(out)} formulations, {out.solvent.nunique()} unique solvents")
    print(f"-> {INDEX/'pilot_rows.csv'}")
    print(f"-> {INDEX/'pilot_solvents.txt'}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
