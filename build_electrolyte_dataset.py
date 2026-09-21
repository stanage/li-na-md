"""Build the electrolyte dataset: 2000 rows over the generated solvent library.

Two choices differ from the first version of this dataset:

  * solvents come from the MW <= 200, aldehyde-free, aprotic library
  * the concentration grid tops out at 1.5 M instead of 3.0 M

Together those keep solvent molecules per ion pair between roughly 16 and 3
across the whole set, which is the range EC spans over 0.3-3.0 M. The v1
pairing put 1463 of 2000 rows below 4 -- solvate melts whose salt solubility
was never established.

The concentration is still drawn independently of the solvent, as in v1. That
is safe here only because the MW window is narrow enough that every row lands
in a sane regime; it is the MW cap, not the sampling, that does the work. A
`solv_per_ion_pair` column is written so the assumption stays checkable.

Usage:  python build_electrolyte_dataset.py
"""
from __future__ import annotations

import argparse
import sys
from pathlib import Path

import numpy as np
import pandas as pd

sys.path.insert(0, str(Path(__file__).resolve().parent / "MD" / "ht2000" / "scripts"))

from common import (N_AVOGADRO, RHO_SALT, RHO_SOLVENT,  # noqa: E402
                    ANION_NAMES, fix_salt, mol_weight)

HERE = Path(__file__).resolve().parent
LI, NA = "[Li+]", "[Na+]"

#: SMILES as build_electrolyte_dataset.py writes them -- kept verbatim so the
#: two datasets use identical salt strings and share the ff_ions cache.
ANIONS = {v: k for k, v in ANION_NAMES.items()}
ANIONS["NO3"] = "[O-][N+](=O)=O"          # v1 spelling; fix_salt() repairs it


def solv_per_pair(molarity: float, mw_solvent: float, mw_salt: float,
                  n_salt: int = 64) -> float:
    v_target = n_salt / (molarity * N_AVOGADRO / 1e24)
    v_salt = n_salt * mw_salt / N_AVOGADRO * 1e21 / RHO_SALT
    return (v_target - v_salt) * RHO_SOLVENT / 1e21 * N_AVOGADRO / mw_solvent / n_salt


def main() -> int:
    ap = argparse.ArgumentParser()
    ap.add_argument("--solvents", type=Path,
                    default=HERE / "datasets" / "generated_500_organic_solvents.csv")
    ap.add_argument("--out", type=Path,
                    default=HERE / "datasets" / "electrolyte_dataset_2000.csv")
    ap.add_argument("--n-rows", type=int, default=2000)
    ap.add_argument("--conc-min", type=float, default=0.3)
    ap.add_argument("--conc-max", type=float, default=1.5)
    ap.add_argument("--conc-steps", type=int, default=28)
    ap.add_argument("--seed", type=int, default=42)
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
    print(f"{len(df)} rows, {df.solvent.nunique()} unique solvents, "
          f"{df.salt.nunique()} salts")
    print(f"concentration : {df.concentration.min()}-{df.concentration.max()} M "
          f"over {len(conc_options)} levels")
    print(f"solvent MW    : {smw.min():.0f}-{smw.max():.0f} "
          f"(median {np.median(smw):.0f})")
    print(f"solv per pair : {r.min():.1f}-{r.max():.1f} (median {r.median():.1f})")
    print(f"  below 4     : {(r < 4).sum()} / {len(df)}   "
          f"[the 0.3-3.0 M / MW-242 version had 1463/2000]")
    print(f"  below 2     : {(r < 2).sum()} / {len(df)}")
    print(f"\n-> {args.out}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
