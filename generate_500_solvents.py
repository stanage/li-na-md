"""Generate a light, aldehyde-free solvent library with Electrolyte-GPT.

Why the filters are what they are. The first version of this library had
median MW 242 and a 113-359 range, which is 2.7x EC/DMC. Paired with a
0.3-3.0 M grid that gave under 4 solvent molecules per ion pair in 1463 of
2000 rows -- solvate melts rather than solutions, for solvents whose salt
solubility nobody has checked. It also contains 108 aldehydes, which reduce
on alkali metal and are not used as electrolyte solvents.

So this run tightens two filters on top of v1's (valid / neutral / single
fragment / has carbon / no metals):

  * MW <= 200  -- with a 0.3-1.5 M grid this keeps solvent-per-ion-pair
    between 16.1 and 2.7, matching what EC spans over 0.3-3.0 M.
  * no aldehyde -- SMARTS [CX3H1](=O)[#6]

Roughly 8% of v1's accepted molecules would pass both, so expect ~12x the
generation of v1. Resumable: point --out at a partial CSV and it tops it up.

Usage (GPU node):
  python generate_500_solvents.py --target 500
"""

import argparse
import json
import os
import re
import sys

import pandas as pd
import torch
from rdkit import Chem, RDLogger
from rdkit.Chem.Descriptors import ExactMolWt, MolWt
from tqdm import tqdm

RDLogger.DisableLog("rdApp.*")

ORIG = "/project2/chibueze/jaemink/genMolGPT/Electrolyte-GPT"
MODEL_DIR = f"{ORIG}/pretrained_models"
JSON_DIR = f"{ORIG}/json"

sys.path.insert(0, ORIG)
from model import GPT, GPTConfig          # noqa: E402
from utils import sample                  # noqa: E402
from get_mol import get_mol               # noqa: E402

METAL_ATOMIC_NUMS = (set(range(3, 5)) | set(range(11, 15)) | set(range(19, 35))
                     | set(range(37, 53)) | set(range(55, 85))
                     | set(range(87, 119)))

#: R-CHO. Formate esters (HC(=O)O-) are deliberately not matched: the carbonyl
#: there is an ester, which is ordinary electrolyte chemistry.
ALDEHYDE = Chem.MolFromSmarts("[CX3H1](=O)[#6]")

#: Alcohols/acids and primary/secondary amines. Protic solvents react with
#: alkali metal, so they are not usable here whatever else they look like.
PROTIC = Chem.MolFromSmarts("[OX2H,NX3;H1,H2]")

#: Elements BOSS/OPLS can type. Boron is absent deliberately: both boron
#: solvents in the v1 campaign failed at LigParGen and none succeeded.
#: Silicon is absent too, but for a different reason -- METAL_ATOMIC_NUMS
#: (inherited from v1) covers Z=11..14, so Si is already rejected as a
#: "metal" and never reaches this check. Listing it here would be dead code.
#: LigParGen does in fact handle Si; add 14 here AND drop it from
#: METAL_ATOMIC_NUMS if silanes/siloxanes are ever wanted.
ALLOWED_Z = {1, 6, 7, 8, 9, 15, 16, 17, 35, 53}

#: Phosphorus with no P-C bond, i.e. phosphate and phosphite esters. Measured
#: from the v1 campaign: all 6 P-containing solvents LigParGen handled have at
#: least one P-C bond, and the one that failed had P bonded to four oxygens.
#: So this rejects the failing motif without excluding phosphorus wholesale.
BAD_P = Chem.MolFromSmarts("[#15;!$([#15]~[#6])]")


def classify(mol, mw_min: float, mw_max: float) -> str:
    """'' if the molecule is kept, else the name of the filter that rejected it.

    Returning the reason rather than a bool is what lets the run report where
    the yield is going, so the oversampling factor for a future run is known
    rather than guessed.
    """
    if mol is None:
        return "unparseable"
    has_c = False
    for atom in mol.GetAtoms():
        an = atom.GetAtomicNum()
        if an in METAL_ATOMIC_NUMS:
            return "metal"
        if atom.GetFormalCharge() != 0:
            return "charged"
        if an not in ALLOWED_Z:
            return "element_boss_cannot_type"
        if an == 6:
            has_c = True
    if not has_c:
        return "no_carbon"
    if "." in Chem.MolToSmiles(mol):
        return "multi_fragment"
    mw = MolWt(mol)
    if mw < mw_min:
        return "mw_low"
    if mw > mw_max:
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
    ap.add_argument("--model_weight", default="Jaemin_unconditioned_20epoch.pt")
    ap.add_argument("--target", type=int, default=500)
    ap.add_argument("--mw-min", type=float, default=30.0)
    ap.add_argument("--mw-max", type=float, default=200.0)
    ap.add_argument("--batch_size", type=int, default=192)
    ap.add_argument("--max_iters", type=int, default=4000)
    ap.add_argument("--vocab_size", type=int, default=82)
    ap.add_argument("--block_size", type=int, default=190)
    ap.add_argument("--n_layer", type=int, default=8)
    ap.add_argument("--n_head", type=int, default=8)
    ap.add_argument("--n_embd", type=int, default=256)
    ap.add_argument("--temperature", type=float, default=1.0)
    ap.add_argument("--save-every", type=int, default=25,
                    help="checkpoint the CSV every N batches")
    ap.add_argument("--out", default=os.path.join(
        os.path.dirname(os.path.abspath(__file__)),
        "datasets", "generated_500_organic_solvents.csv"))
    args = ap.parse_args()

    device = "cuda" if torch.cuda.is_available() else "cpu"
    print(f"device: {device}   MW window: {args.mw_min}-{args.mw_max}   "
          f"aldehydes: rejected")
    if device == "cpu":
        print("WARNING: no GPU visible; this will be very slow")

    stoi = json.load(open(os.path.join(
        JSON_DIR, args.model_weight.replace(".pt", "_stoi.json"))))
    itos = {i: ch for ch, i in stoi.items()}

    mconf = GPTConfig(args.vocab_size, args.block_size, num_props=0,
                      n_layer=args.n_layer, n_head=args.n_head,
                      n_embd=args.n_embd)
    model = GPT(mconf)
    model.load_state_dict(torch.load(
        os.path.join(MODEL_DIR, args.model_weight), map_location=device))
    model.to(device).eval()
    print("model loaded")

    regex = re.compile(r"(\[[^\]]+]|<|Br?|Cl?|N|O|S|P|F|I|b|c|n|o|s|p|\(|\)|\.|"
                       r"=|#|-|\+|\\\\|\/|:|~|@|\?|>|\*|\$|\%[0-9]{2}|[0-9])")
    ctx = torch.tensor([stoi[s] for s in regex.findall("C")],
                       dtype=torch.long)[None, ...].repeat(args.batch_size, 1)

    # resume from a partial run rather than starting over
    keep: dict[str, object] = {}
    out_path = os.path.abspath(args.out)
    n_start = 0
    if os.path.exists(out_path):
        prev = pd.read_csv(out_path)
        dropped: dict[str, int] = {}
        for smi in prev["smiles"]:
            m = Chem.MolFromSmiles(smi)
            why = classify(m, args.mw_min, args.mw_max)
            if why:
                dropped[why] = dropped.get(why, 0) + 1
                continue
            keep[Chem.MolToSmiles(m)] = m
        print(f"resuming from {len(keep)} molecules already in {out_path}")
        if dropped:
            print(f"  dropped {sum(dropped.values())} that fail the current "
                  f"filters: {dropped}")
    n_start = len(keep)

    def flush():
        os.makedirs(os.path.dirname(out_path), exist_ok=True)
        pd.DataFrame([{"smiles": s, "mol_weight": MolWt(m),
                       "exact_mol_weight": ExactMolWt(m)}
                      for s, m in keep.items()]).to_csv(out_path, index=False)

    reasons: dict[str, int] = {}
    total = dupes = 0
    pbar = tqdm(total=args.target, initial=len(keep), desc="solvents")

    for it in range(args.max_iters):
        if len(keep) >= args.target:
            break
        with torch.no_grad():
            y = sample(model, ctx.clone().to(device), args.block_size,
                       temperature=args.temperature, sample=True,
                       top_k=None, prop=None)
        for gen in y:
            total += 1
            smi = "".join(itos[int(i)] for i in gen).replace("<", "")
            mol = get_mol(smi)
            why = classify(mol, args.mw_min, args.mw_max)
            if why:
                reasons[why] = reasons.get(why, 0) + 1
                continue
            canon = Chem.MolToSmiles(mol)
            if canon in keep:
                dupes += 1
                continue
            keep[canon] = mol
            pbar.update(1)
            if len(keep) >= args.target:
                break
        if (it + 1) % args.save_every == 0:
            flush()
            pbar.set_postfix(kept=len(keep), tried=total)

    pbar.close()
    flush()

    new = len(keep) - n_start
    print(f"\n{len(keep)} unique solvents total "
          f"({new} new this session from {total} samples, "
          f"{100*new/max(total,1):.2f}% yield)")
    print(f"duplicates discarded: {dupes}")
    print("rejected by filter:")
    for k, v in sorted(reasons.items(), key=lambda kv: -kv[1]):
        print(f"   {k:16s} {v:7d}  ({100*v/max(total,1):5.2f}%)")
    print(f"\n-> {out_path}")
    if len(keep) < args.target:
        print(f"NOTE: short of --target {args.target}; rerun the same command "
              f"to top up (it resumes from the CSV)")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
