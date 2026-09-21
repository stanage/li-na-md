#!/usr/bin/env python
"""Select ht1000's solvent pool: every oxygen-bearing solvent in the library.

Why this exists rather than a fresh generation run. ht2000's 500-solvent
library has a hole in it that only shows up in the analysis: **222 of the 500
contain no oxygen at all, and 106 are pure hydrocarbons.** The whole point of
`analyze_solvation.py` is the cation-solvent-oxygen RDF -- first peak position,
first minimum, coordination number. For an O-free solvent that quantity is not
small, it is *undefined*: the run produces a solvation.json whose `cn_solv_O`
is 0 and whose `rdf_first_peak_A` is noise, and the row contributes nothing to
the structure-property relationship the campaign exists to measure.

So ht1000 keeps the chemistry filters ht2000 already applies (no boron, no
charged or multi-fragment species, no metals, no aldehydes, no protic species,
no radicals, no phosphate/phosphite esters, MW <= 200) and adds one more:

    every solvent must contain at least one oxygen atom.

That single rule removes all 106 pure hydrocarbons and all 222 O-free
molecules by construction, because a pure hydrocarbon has no oxygen.

The pool is a *subset of ht2000's library*, taken rather than regenerated at
the user's direction (generating needs a GPU node, and the existing library
already contains 278 qualifying molecules with force fields almost all in
hand). Sharing the *molecules* is not sharing *files*: ht1000 writes its own
copies of every .itp/.gro/.pdb under ht1000/ff_solvents and never reads
ht2000's at run time -- see import_solvent_ff.py.

Usage:
  python select_solvent_pool.py                 # -> datasets/solvents_oxygen_pool.csv
  python select_solvent_pool.py --report        # composition census, no write
"""
from __future__ import annotations

import argparse
import sys
from collections import Counter
from pathlib import Path

import pandas as pd
from rdkit import Chem, RDLogger
from rdkit.Chem import Descriptors

sys.path.insert(0, str(Path(__file__).resolve().parent))

from common import REPO, SOLVENT_POOL, solvent_key  # noqa: E402

RDLogger.DisableLog("rdApp.*")

#: ht2000's generated library. Read-only: ht1000 never writes here.
SOURCE_LIBRARY = REPO / "datasets" / "generated_500_organic_solvents.csv"

#: Same element set generate_500_solvents.py allows, restated here so this
#: script is a complete statement of what the pool may contain rather than a
#: reference to a file in another campaign's history.
ALLOWED_Z = {1, 6, 7, 8, 9, 15, 16, 17, 35, 53}

ALDEHYDE = Chem.MolFromSmarts("[CX3H1](=O)[#6]")
PROTIC = Chem.MolFromSmarts("[OX2H,NX3;H1,H2]")
BAD_P = Chem.MolFromSmarts("[#15;!$([#15]~[#6])]")

#: Functional families, for the census. A molecule can match several; these
#: are reported to document what the pool actually spans, not to gate it.
FAMILIES = {
    "carbonate": "[OX2][CX3](=O)[OX2]",
    "ester": "[CX3](=O)[OX2H0][#6]",
    "ether": "[OD2]([#6])[#6]",
    "ketone": "[#6][CX3](=O)[#6]",
    "sulfone": "[SX4](=O)(=O)",
    "sulfoxide": "[SX3](=O)[#6]",
    "amide": "[NX3][CX3](=O)",
    "nitrile": "C#N",
    "fluorinated": "[F]",
    "chlorinated": "[Cl]",
    "phosphoryl": "[PX4]=O",
}


def rejects(mol, mw_max: float) -> str:
    """'' if the molecule belongs in the pool, else why it does not."""
    if mol is None:
        return "unparseable"
    has_c = has_o = False
    for atom in mol.GetAtoms():
        an = atom.GetAtomicNum()
        if atom.GetFormalCharge() != 0:
            return "charged"
        if an not in ALLOWED_Z:
            return "element_boss_cannot_type"
        if an == 6:
            has_c = True
        elif an == 8:
            has_o = True
        if atom.GetNumRadicalElectrons():
            return "radical"
    if not has_c:
        return "no_carbon"
    # The ht1000 rule. Everything above this line is inherited from ht2000.
    if not has_o:
        return "no_oxygen"
    if "." in Chem.MolToSmiles(mol):
        return "multi_fragment"
    if Descriptors.MolWt(mol) > mw_max:
        return "mw_high"
    if mol.HasSubstructMatch(ALDEHYDE):
        return "aldehyde"
    if mol.HasSubstructMatch(PROTIC):
        return "protic"
    if mol.HasSubstructMatch(BAD_P):
        return "phosphate_ester"
    return ""


def main() -> int:
    ap = argparse.ArgumentParser()
    ap.add_argument("--source", type=Path, default=SOURCE_LIBRARY)
    ap.add_argument("--out", type=Path, default=SOLVENT_POOL)
    ap.add_argument("--mw-max", type=float, default=200.0)
    ap.add_argument("--report", action="store_true",
                    help="print the census and exit without writing")
    args = ap.parse_args()

    lib = pd.read_csv(args.source)
    kept, why = [], Counter()
    for smi in lib["smiles"]:
        mol = Chem.MolFromSmiles(smi)
        reason = rejects(mol, args.mw_max)
        if reason:
            why[reason] += 1
            continue
        canon = Chem.MolToSmiles(mol)
        kept.append({
            "smiles": canon,
            "mol_weight": round(Descriptors.MolWt(mol), 3),
            "n_oxygen": sum(1 for a in mol.GetAtoms() if a.GetSymbol() == "O"),
            "n_heavy": mol.GetNumHeavyAtoms(),
            "solvent_key": solvent_key(canon),
        })

    pool = pd.DataFrame(kept).drop_duplicates(subset="smiles")
    print(f"source library : {len(lib)} solvents ({args.source})")
    print(f"ht1000 pool    : {len(pool)} solvents (all oxygen-bearing)")
    print("rejected:")
    for k, v in why.most_common():
        print(f"   {k:26s} {v:4d}")

    mols = [Chem.MolFromSmiles(s) for s in pool.smiles]
    print(f"\nMW {pool.mol_weight.min():.0f}-{pool.mol_weight.max():.0f} "
          f"(mean {pool.mol_weight.mean():.0f})")
    print(f"oxygens per molecule: "
          f"{dict(sorted(Counter(pool.n_oxygen).items()))}")
    print("functional families (molecules may match more than one):")
    for name, sma in FAMILIES.items():
        q = Chem.MolFromSmarts(sma)
        print(f"   {name:13s} {sum(1 for m in mols if m.HasSubstructMatch(q)):4d}")
    els = Counter()
    for m in mols:
        for e in {a.GetSymbol() for a in m.GetAtoms()}:
            els[e] += 1
    print(f"element presence: {dict(els.most_common())}")

    if args.report:
        return 0
    args.out.parent.mkdir(parents=True, exist_ok=True)
    pool.to_csv(args.out, index=False)
    print(f"\n-> {args.out}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
