"""
Generate 500 unique organic solvents using the ElectrolyteGPT unconditioned model.

Organic solvent filters applied:
  - Valid SMILES (parseable by RDKit)
  - Contains at least one carbon
  - No metals
  - Neutral (no formal charge on any atom)
  - Single connected fragment (no salts/mixtures)
  - Molecular weight between 30 and 600 g/mol
"""

import sys
import os
import math
import re
import json
import argparse

import pandas as pd
import numpy as np
import torch
from tqdm import tqdm
from rdkit import Chem
from rdkit.Chem.Descriptors import ExactMolWt

ORIG = '/project2/chibueze/jaemink/genMolGPT/MolGPT_DrugDesign'
MODEL_DIR = f'{ORIG}/pretrained_models'
JSON_DIR = f'{ORIG}/json'

# Use the original model/utils that match the checkpoint format
sys.path.insert(0, ORIG)
from model import GPT, GPTConfig
from utils import sample, canonic_smiles
from get_mol import get_mol

METAL_ATOMIC_NUMS = set(range(3, 5)) | set(range(11, 15)) | set(range(19, 35)) | \
                    set(range(37, 53)) | set(range(55, 85)) | set(range(87, 119))

def is_organic_solvent(mol):
    if mol is None:
        return False
    atoms = mol.GetAtoms()
    has_carbon = False
    for atom in atoms:
        an = atom.GetAtomicNum()
        if an in METAL_ATOMIC_NUMS:
            return False
        if atom.GetFormalCharge() != 0:
            return False
        if an == 6:
            has_carbon = True
    if not has_carbon:
        return False
    smi = Chem.MolToSmiles(mol)
    if '.' in smi:
        return False
    mw = ExactMolWt(mol)
    if not (30 <= mw <= 600):
        return False
    return True


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument('--model_weight', default='Jaemin_unconditioned_range_20epoch.pt')
    parser.add_argument('--data_name', default='Jaemin_trainingdataset_7properties2')
    parser.add_argument('--target', type=int, default=500)
    parser.add_argument('--batch_size', type=int, default=192)
    parser.add_argument('--max_iters', type=int, default=200,
                        help='max generation batches (safety cap)')
    parser.add_argument('--vocab_size', type=int, default=82)
    parser.add_argument('--block_size', type=int, default=190)
    parser.add_argument('--n_layer', type=int, default=8)
    parser.add_argument('--n_head', type=int, default=8)
    parser.add_argument('--n_embd', type=int, default=256)
    parser.add_argument('--out', default='/project/chibueze/stanley/projects/li-na-md/datasets/generated_500_organic_solvents.csv')
    args = parser.parse_args()

    device = 'cuda' if torch.cuda.is_available() else 'cpu'
    print(f'Using device: {device}')

    # Load vocabulary
    stoi_path = os.path.join(JSON_DIR, args.model_weight.replace('.pt', '_stoi.json'))
    stoi = json.load(open(stoi_path))
    itos = {i: ch for ch, i in stoi.items()}

    # Load model — scaffold_maxlen=1 matches how the checkpoint was saved
    mconf = GPTConfig(args.vocab_size, args.block_size, num_props=0,
                      n_layer=args.n_layer, n_head=args.n_head, n_embd=args.n_embd,
                      scaffold=False, scaffold_maxlen=1, lstm=False, lstm_layers=0)
    model = GPT(mconf)
    model.load_state_dict(torch.load(os.path.join(MODEL_DIR, args.model_weight),
                                     map_location=device))
    model.to(device)
    model.eval()
    print('Model loaded.')

    pattern = r"(\[[^\]]+]|<|Br?|Cl?|N|O|S|P|F|I|b|c|n|o|s|p|\(|\)|\.|=|#|-|\+|\\\\|\/|:|~|@|\?|>|\*|\$|\%[0-9]{2}|[0-9])"
    regex = re.compile(pattern)
    context = "C"
    context_tensor = torch.tensor(
        [stoi[s] for s in regex.findall(context)],
        dtype=torch.long
    )[None, ...].repeat(args.batch_size, 1)

    unique_solvents = {}  # canonical SMILES -> RDKit mol
    total_generated = 0

    pbar = tqdm(total=args.target, desc='Unique organic solvents')

    for iteration in range(args.max_iters):
        if len(unique_solvents) >= args.target:
            break

        x = context_tensor.clone().to(device)
        with torch.no_grad():
            y = sample(model, x, args.block_size, temperature=1.0,
                       sample=True, top_k=None, prop=None)

        for gen_mol in y:
            total_generated += 1
            completion = ''.join([itos[int(i)] for i in gen_mol]).replace('<', '')
            mol = get_mol(completion)
            if not is_organic_solvent(mol):
                continue
            canon = Chem.MolToSmiles(mol)
            if canon not in unique_solvents:
                unique_solvents[canon] = mol
                pbar.update(1)
                if len(unique_solvents) >= args.target:
                    break

    pbar.close()

    records = [{'smiles': smi, 'mol_weight': ExactMolWt(mol)}
               for smi, mol in unique_solvents.items()]
    df = pd.DataFrame(records)

    out_path = args.out
    os.makedirs(os.path.dirname(out_path), exist_ok=True)
    df.to_csv(out_path, index=False)

    print(f'\nGenerated {len(df)} unique organic solvents from {total_generated} total attempts.')
    print(f'Results saved to: {out_path}')


if __name__ == '__main__':
    main()
