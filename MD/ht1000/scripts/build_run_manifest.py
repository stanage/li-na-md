#!/usr/bin/env python
"""Index every row of the dataset up front, before anything is built.

results.csv only knows about runs that exist on disk, so it answers
"which dataset row is el9f2c...?" but not "was row 742 ever built?". This
writes the complete row -> key mapping for all 1000 rows, together with the
composition each row will produce and whether its force fields are in hand,
so campaign coverage and cost can be checked without launching anything.

Usage:  python build_run_manifest.py [--out ../index/run_manifest.csv]
"""
from __future__ import annotations

import argparse
import sys
from pathlib import Path

import pandas as pd

sys.path.insert(0, str(Path(__file__).resolve().parent))

from common import (DATASET, FF_SOLV, INDEX, N_AVOGADRO, RUNS,  # noqa: E402
                    anion_name, cation_name, compute_counts, fix_salt,
                    mol_weight, n_atoms_with_h, run_key, solvent_key,
                    split_salt)

STAGES = ["em", "nvt_heat", "npt", "prod"]


def stage_reached(d: Path) -> str:
    last = "not_built"
    if d.is_dir():
        last = "built"
        for s in STAGES:
            if (d / f"{s}.gro").exists():
                last = s
    return last


def main() -> int:
    ap = argparse.ArgumentParser()
    ap.add_argument("--out", type=Path, default=INDEX / "run_manifest.csv")
    args = ap.parse_args()

    df = pd.read_csv(DATASET)
    mw_cache: dict[str, float] = {}
    nat_cache: dict[str, int] = {}

    def mw(smi: str) -> float:
        if smi not in mw_cache:
            mw_cache[smi] = mol_weight(smi)
        return mw_cache[smi]

    def nat(smi: str) -> int:
        if smi not in nat_cache:
            nat_cache[smi] = n_atoms_with_h(smi)
        return nat_cache[smi]

    rows = []
    for i, r in df.iterrows():
        salt, solv, conc = r["salt"], r["solvent"], float(r["concentration"])
        cat_smi, an_smi = split_salt(fix_salt(salt))
        skey = solvent_key(solv)
        c = compute_counts(mw(solv), mw(fix_salt(salt)), conc)

        # atoms the packed box will hold, so the campaign's cost is knowable now
        n_atoms = (c["n_solvent"] * nat(solv)
                   + c["n_salt"] * (nat(cat_smi) + nat(an_smi)))

        key = run_key(salt, solv, conc)
        rows.append({
            "row_index": int(i),
            "csv_line": int(i) + 2,          # 1-based, header included
            "key": key,
            "salt": salt,
            "solvent": solv,
            "concentration": conc,
            "cation": cation_name(cat_smi),
            "anion": anion_name(an_smi),
            "solvent_key": skey,
            "mw_solvent": round(mw(solv), 3),
            "mw_salt": round(mw(fix_salt(salt)), 3),
            "n_solvent": c["n_solvent"],
            "n_salt": c["n_salt"],
            "solv_per_ion_pair": round(c["n_solvent"] / c["n_salt"], 3),
            "box_build_nm": round(c["box_build_nm"], 4),
            "n_atoms_est": n_atoms,
            "solvent_ff": (FF_SOLV / f"{skey}.itp").exists(),
            "clamped": bool(c["clamped"]),
            "stage": stage_reached(RUNS / key),
        })

    out = pd.DataFrame(rows)
    if out.key.duplicated().any():
        dup = out[out.key.duplicated(keep=False)]
        print(f"WARNING: {len(dup)} rows share a run key -- they would collide "
              f"in the same directory:\n{dup[['row_index', 'key']]}")

    args.out.parent.mkdir(parents=True, exist_ok=True)
    out.to_csv(args.out, index=False)

    runnable = out[out.solvent_ff]
    print(f"{len(out)} rows, {out.key.nunique()} unique keys")
    print(f"solvent force field present: {len(runnable)}  "
          f"missing: {len(out) - len(runnable)}")
    print(f"\nstage:\n{out.stage.value_counts().to_string()}")
    print(f"\nsolvent per ion pair: median {out.solv_per_ion_pair.median():.1f}, "
          f"<4 in {(out.solv_per_ion_pair < 4).sum()} rows")
    print(f"atoms per box: median {int(out.n_atoms_est.median())}, "
          f"max {int(out.n_atoms_est.max())}, total {int(out.n_atoms_est.sum()):,}")
    print(f"\n-> {args.out}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
