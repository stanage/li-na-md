#!/usr/bin/env python
"""Build ht1000's dataset: 1000 Li/Na electrolytes over the oxygen-only pool.

Same construction as ht2000's dataset builder, with two differences that are
the whole point of this campaign:

  * solvents come from `datasets/solvents_oxygen_pool.csv` -- 278 molecules
    that all contain at least one oxygen, so the cation-solvent-O RDF that
    `analyze_solvation.py` measures is defined for every single row. In
    ht2000, 222 of 500 solvents had no oxygen at all and 106 were pure
    hydrocarbons, which made those rows uninformative about solvation
    structure (see select_solvent_pool.py).
  * 1000 rows rather than 2000, matching the two GPUs this campaign runs on.

Everything else is deliberately identical to ht2000 so the two datasets are
comparable and can be pooled: 40 salts (20 anions x Li/Na), the same anion
SMILES strings character for character (including the malformed nitrate that
`fix_salt` repairs, so both campaigns key the same ion force fields), the
0.3-1.5 M / 28-level concentration grid, and one solvent per row.

1000 / 40 = 25 rows per salt, drawn without replacement, so no (salt, solvent)
pair repeats. With 278 solvents in the pool each solvent appears about 3.6
times across *different* salts -- close to ht2000's 4.1, so neither campaign
leans harder on any one molecule.

The draw uses its own seed (1000, not ht2000's 42). Sharing the seed would be
harmless but misleading: these are different pools, so identical seeds would
not produce corresponding rows anyway.

Usage:  python build_dataset.py
"""
from __future__ import annotations

import argparse
import sys
from pathlib import Path

import numpy as np
import pandas as pd

sys.path.insert(0, str(Path(__file__).resolve().parent))

from common import (ANION_NAMES, DATASET, N_AVOGADRO, RHO_SALT,  # noqa: E402
                    RHO_SOLVENT, SOLVENT_POOL, fix_salt, mol_weight)

LI, NA = "[Li+]", "[Na+]"

#: Anion SMILES exactly as ht2000's builder writes them, so the two datasets
#: share ion force fields and the salt column is directly comparable.
ANIONS = {v: k for k, v in ANION_NAMES.items()}
ANIONS["NO3"] = "[O-][N+](=O)=O"          # malformed on purpose; fix_salt() repairs it


def solv_per_pair(molarity: float, mw_solvent: float, mw_salt: float,
                  n_salt: int = 64) -> float:
    """Solvent molecules per ion pair the box builder will produce.

    Same arithmetic as common.compute_counts, without the integer rounding or
    the MIN_SOLVENT floor -- this is a dataset-design diagnostic, not the
    composition that gets packed.
    """
    v_target = n_salt / (molarity * N_AVOGADRO / 1e24)
    v_salt = n_salt * mw_salt / N_AVOGADRO * 1e21 / RHO_SALT
    return (v_target - v_salt) * RHO_SOLVENT / 1e21 * N_AVOGADRO / mw_solvent / n_salt


def main() -> int:
    ap = argparse.ArgumentParser()
    ap.add_argument("--solvents", type=Path, default=SOLVENT_POOL)
    ap.add_argument("--out", type=Path, default=DATASET)
    ap.add_argument("--n-rows", type=int, default=1000)
    ap.add_argument("--conc-min", type=float, default=0.3)
    ap.add_argument("--conc-max", type=float, default=1.5)
    ap.add_argument("--conc-steps", type=int, default=28)
    ap.add_argument("--seed", type=int, default=1000)
    args = ap.parse_args()

    salts = ([f"{LI}.{smi}" for smi in ANIONS.values()]
             + [f"{NA}.{smi}" for smi in ANIONS.values()])
    per_salt, rem = divmod(args.n_rows, len(salts))
    if rem:
        raise SystemExit(f"--n-rows must divide by {len(salts)} salts")

    solvents = pd.read_csv(args.solvents)["smiles"].tolist()
    if len(solvents) < per_salt:
        raise SystemExit(f"need >= {per_salt} solvents, have {len(solvents)}")

    rng = np.random.default_rng(args.seed)
    conc_options = np.round(np.linspace(args.conc_min, args.conc_max,
                                        args.conc_steps), 2)

    mwc: dict[str, float] = {}

    def mw(smi: str) -> float:
        if smi not in mwc:
            mwc[smi] = mol_weight(smi)
        return mwc[smi]

    rows = []
    for salt in salts:
        chosen = rng.choice(solvents, size=per_salt, replace=False)
        concs = rng.choice(conc_options, size=per_salt, replace=True)
        mw_salt = mw(fix_salt(salt))
        for solv, conc in zip(chosen, concs):
            rows.append({
                "salt": salt,
                "solvent": solv,
                "concentration": float(conc),
                "solv_per_ion_pair": round(
                    solv_per_pair(float(conc), mw(solv), mw_salt), 2),
            })

    df = pd.DataFrame(rows)
    n_uniq = df[["salt", "solvent"]].drop_duplicates().shape[0]
    if n_uniq != len(df):
        raise SystemExit(f"{len(df) - n_uniq} duplicate (salt, solvent) pairs")

    args.out.parent.mkdir(parents=True, exist_ok=True)
    df.to_csv(args.out, index=False)

    smw = np.array([mw(s) for s in df.solvent])
    r = df.solv_per_ion_pair
    reuse = df.solvent.value_counts()
    print(f"{len(df)} rows, {df.solvent.nunique()} unique solvents "
          f"(of {len(solvents)} in the pool), {df.salt.nunique()} salts")
    print(f"concentration : {df.concentration.min()}-{df.concentration.max()} M "
          f"over {len(conc_options)} levels")
    print(f"solvent MW    : {smw.min():.0f}-{smw.max():.0f} "
          f"(median {np.median(smw):.0f})")
    print(f"solvent reuse : {reuse.min()}-{reuse.max()} rows each "
          f"(mean {reuse.mean():.1f})")
    print(f"solv per pair : {r.min():.1f}-{r.max():.1f} (median {r.median():.1f})")
    print(f"  below 4     : {(r < 4).sum()} / {len(df)}")
    print(f"  below 2     : {(r < 2).sum()} / {len(df)}")
    print(f"\n-> {args.out}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
