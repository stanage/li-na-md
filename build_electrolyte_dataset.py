"""
Build a 2000-entry electrolyte dataset with columns: salt, solvent, concentration.
- Li and Na salts (SMILES), split evenly (1000 each)
- Solvents drawn from datasets/generated_500_organic_solvents.csv
- Concentrations 0.3–3.0 M
- All (salt, solvent) pairs are unique
"""

import pandas as pd
import numpy as np
import os

ANIONS = {
    "FSI":   "[N-](S(=O)(=O)F)S(=O)(=O)F",
    "TFSI":  "[N-](S(=O)(=O)C(F)(F)F)S(=O)(=O)C(F)(F)F",
    "Beti":  "[N-](S(=O)(=O)C(F)(F)C(F)(F)F)S(=O)(=O)C(F)(F)C(F)(F)F",
    "BF4":   "[B-](F)(F)(F)F",
    "PF6":   "F[P-](F)(F)(F)(F)F",
    "ClO4":  "[O-]Cl(=O)(=O)=O",
    "OTf":   "[O-]S(=O)(=O)C(F)(F)F",
    "DCA":   "[N-](C#N)C#N",
    "TCM":   "[C-](C#N)(C#N)C#N",
    "NO3":   "[O-][N+](=O)=O",
    "OAc":   "CC(=O)[O-]",
    "TFA":   "[O-]C(=O)C(F)(F)F",
    "PROP":  "CCC(=O)[O-]",
    "BNZ":   "[O-]C(=O)c1ccccc1",
    "HCOO":  "[O-]C=O",
    "MS":    "COS(=O)(=O)[O-]",
    "SCN":   "[S-]C#N",
    "Br":    "[Br-]",
    "Cl":    "[Cl-]",
    "OTs":   "[O-]S(=O)(=O)c1ccc(C)cc1",
}

def make_salt_smiles(cation_smiles, anion_smiles):
    return f"{cation_smiles}.{anion_smiles}"

LI = "[Li+]"
NA = "[Na+]"

li_salts = {f"Li{name}": make_salt_smiles(LI, smi) for name, smi in ANIONS.items()}
na_salts = {f"Na{name}": make_salt_smiles(NA, smi) for name, smi in ANIONS.items()}

# 20 Li salts, 20 Na salts
all_salts = list(li_salts.values()) + list(na_salts.values())  # 40 total
n_salts = len(all_salts)  # 40
entries_per_salt = 2000 // n_salts  # 50 each

solvents_df = pd.read_csv(
    '/project/chibueze/stanley/projects/li-na-md/datasets/generated_500_organic_solvents.csv'
)
solvents = solvents_df['smiles'].tolist()  # 500 solvents

rng = np.random.default_rng(42)

# Concentration grid: 28 evenly spaced values from 0.3 to 3.0
conc_options = np.round(np.linspace(0.3, 3.0, 28), 2)

rows = []
for salt_smiles in all_salts:
    # Pick 50 unique solvents for this salt (no repeats within the salt)
    chosen_solvents = rng.choice(solvents, size=entries_per_salt, replace=False)
    chosen_concs = rng.choice(conc_options, size=entries_per_salt, replace=True)
    for solv, conc in zip(chosen_solvents, chosen_concs):
        rows.append({
            'salt': salt_smiles,
            'solvent': solv,
            'concentration': conc,
        })

df = pd.DataFrame(rows)

# Verify uniqueness of (salt, solvent) pairs
n_unique = df[['salt', 'solvent']].drop_duplicates().shape[0]
print(f"Total entries     : {len(df)}")
print(f"Unique Li salts   : {df['salt'].str.startswith('[Li+]').sum()}")
print(f"Unique Na salts   : {df['salt'].str.startswith('[Na+]').sum()}")
print(f"Unique (salt, solvent) pairs: {n_unique}")
print(f"Concentration range: {df['concentration'].min()} – {df['concentration'].max()} M")

out_path = '/project/chibueze/stanley/projects/li-na-md/datasets/electrolyte_dataset_2000.csv'
df.to_csv(out_path, index=False)
print(f"Saved to: {out_path}")
