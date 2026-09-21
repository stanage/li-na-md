"""Shared helpers for the 1000-electrolyte high-throughput MD campaign.

Canonical SMILES handling, deterministic run/solvent keys, the box-and-count
math reverse-engineered from the ce_solvation_md archive, and paths.
"""
from __future__ import annotations

import hashlib
import os
import re
from pathlib import Path

from rdkit import Chem, RDLogger
from rdkit.Chem import Descriptors

RDLogger.DisableLog("rdApp.*")

# ---------------------------------------------------------------- paths
HT = Path(__file__).resolve().parent.parent          # MD/ht1000
MD = HT.parent                                        # MD
REPO = MD.parent                                      # project root

FF_SOLV = HT / "ff_solvents"
FF_IONS = HT / "ff_ions"
RUNS = HT / "runs"
INDEX = HT / "index"
LOGS = HT / "logs"

#: ht1000 keeps its dataset inside the campaign directory rather than at the
#: repo root (where ht2000's lives). That is deliberate: the two campaigns must
#: share no writable file anywhere, and a dataset path that resolves through
#: HT cannot accidentally point at the other one.
DATASET = HT / "datasets" / "electrolyte_dataset_1000.csv"
#: The oxygen-filtered solvent pool this dataset is drawn from.
SOLVENT_POOL = HT / "datasets" / "solvents_oxygen_pool.csv"
OLD_RUNS = MD / "ce_solvation_md" / "solvation" / "runs"
OLD_SCRIPTS = (MD / "ce_solvation_md" / "simulation_workflow-main" /
               "simulation_workflow-main" / "scripts")
ANION_ITP_DIR = MD / "anions"
CLANDP = MD / "Force Fields" / "clandp-master"

# Midway3 environment
CONDA_ENV = Path("/scratch/midway3/eshiemogie/moleng")
PYTHON = CONDA_ENV / "bin" / "python"
PACKMOL = CONDA_ENV / "bin" / "packmol"
GMX = "gmx_mpi"                       # Midway3 ships gmx_mpi only, no plain `gmx`
GMX_MODULE = "gromacs/2025.3"
SLURM_ACCOUNT = "pi-chibueze"
SLURM_PARTITION = "gpu"

# ---------------------------------------------------------------- constants
N_AVOGADRO = 6.02214076e23
RHO_SOLVENT = 1.00     # g/cm^3, assumed for every solvent (as the old campaign did)
RHO_SALT = 1.70        # g/cm^3, recovered by fitting the archived summary.json files
PACK_FRACTION = 0.65   # packmol builds the box at 65% of the target density
N_SALT_DEFAULT = 64

#: The single definition of "this run is finished". run_md.sh writes
#: solvation.json before clusters.json, so keying on the earlier one would
#: call a run done whose cluster analysis had failed -- and did, because the
#: drivers disagreed: campaign_worker.sh used clusters.json while
#: launch_batch.py and run_local_batch.sh used solvation.json. Anything that
#: needs the test imports this; the two shell drivers hardcode the same name
#: with a pointer back here.
DONE_MARKER = "clusters.json"
MIN_SOLVENT = 10
Q_SCALE_DEFAULT = 0.8  # ECC charge scaling, ions only

T_K = 298.15
P_BAR = 1.0
PROD_NS = 4.0

#: The dataset's nitrate SMILES is invalid (N valence 5). Corrected here rather
#: than editing the source CSV, which lives outside MD/.
SALT_FIXUPS = {
    "[Li+].[O-][N+](=O)=O": "[Li+].[O-][N+](=O)[O-]",
    "[Na+].[O-][N+](=O)=O": "[Na+].[O-][N+](=O)[O-]",
}

#: Anion SMILES -> short name, mirroring build_electrolyte_dataset.py.
ANION_NAMES = {
    "[N-](S(=O)(=O)F)S(=O)(=O)F": "FSI",
    "[N-](S(=O)(=O)C(F)(F)F)S(=O)(=O)C(F)(F)F": "TFSI",
    "[N-](S(=O)(=O)C(F)(F)C(F)(F)F)S(=O)(=O)C(F)(F)C(F)(F)F": "Beti",
    "[B-](F)(F)(F)F": "BF4",
    "F[P-](F)(F)(F)(F)F": "PF6",
    "[O-]Cl(=O)(=O)=O": "ClO4",
    "[O-]S(=O)(=O)C(F)(F)F": "OTf",
    "[N-](C#N)C#N": "DCA",
    "[C-](C#N)(C#N)C#N": "TCM",
    "[O-][N+](=O)[O-]": "NO3",
    "CC(=O)[O-]": "OAc",
    "[O-]C(=O)C(F)(F)F": "TFA",
    "CCC(=O)[O-]": "PROP",
    "[O-]C(=O)c1ccccc1": "BNZ",
    "[O-]C=O": "HCOO",
    "COS(=O)(=O)[O-]": "MS",
    "[S-]C#N": "SCN",
    "[Br-]": "Br",
    "[Cl-]": "Cl",
    "[O-]S(=O)(=O)c1ccc(C)cc1": "OTs",
}


# ---------------------------------------------------------------- SMILES
def canon(smi: str) -> str:
    """Canonical SMILES, or '' if RDKit cannot parse it."""
    m = Chem.MolFromSmiles(smi)
    return Chem.MolToSmiles(m) if m is not None else ""


def fix_salt(smi: str) -> str:
    """Apply known dataset SMILES corrections."""
    return SALT_FIXUPS.get(smi, smi)


def split_salt(salt_smi: str) -> tuple[str, str]:
    """Split a dotted salt SMILES into (cation, anion) canonical SMILES."""
    parts = [p for p in fix_salt(salt_smi).split(".") if p]
    cation = anion = ""
    for p in parts:
        m = Chem.MolFromSmiles(p)
        if m is None:
            continue
        q = Chem.GetFormalCharge(m)
        if q > 0:
            cation = Chem.MolToSmiles(m)
        elif q < 0:
            anion = Chem.MolToSmiles(m)
    return cation, anion


#: canonical anion SMILES -> name (built once, keyed on canonical form)
ANION_BY_CANON = {canon(k): v for k, v in ANION_NAMES.items()}


def anion_name(anion_smi: str) -> str:
    """Short name for an anion, falling back to a hash-based label."""
    c = canon(anion_smi)
    return ANION_BY_CANON.get(c, "AN" + _h(c, 6).upper())


def cation_name(cation_smi: str) -> str:
    c = canon(cation_smi)
    return {"[Li+]": "Li", "[Na+]": "Na"}.get(c, "CAT")


def mol_weight(smi: str) -> float:
    m = Chem.MolFromSmiles(smi)
    if m is None:
        raise ValueError(f"unparseable SMILES: {smi!r}")
    return Descriptors.MolWt(m)


def n_atoms_with_h(smi: str) -> int:
    m = Chem.MolFromSmiles(smi)
    if m is None:
        raise ValueError(f"unparseable SMILES: {smi!r}")
    return Chem.AddHs(m).GetNumAtoms()


# ---------------------------------------------------------------- keys
def _h(s: str, n: int = 12) -> str:
    return hashlib.sha1(s.encode()).hexdigest()[:n]


def solvent_key(smi: str) -> str:
    """Stable force-field cache key for a solvent."""
    return "S" + _h(canon(smi))


def run_key(salt_smi: str, solvent_smi: str, molarity: float) -> str:
    """Stable run-directory name for one dataset row."""
    tag = f"{canon(fix_salt(salt_smi))}|{canon(solvent_smi)}|{float(molarity):.4f}"
    return "el" + _h(tag)


# ---------------------------------------------------------------- box math
def compute_counts(mw_solvent: float, mw_salt: float, molarity: float,
                   n_salt: int = N_SALT_DEFAULT) -> dict:
    """Solvent count and packmol box edge for one formulation.

    Reproduces the ce_solvation_md sizing rule (see PLAN.md 5.1):
      - molarity + a fixed ion count set the *target* volume
      - packmol builds at 65% of that density; NPT compresses to the true one
      - the salt's own volume is excluded at an assumed 1.70 g/cm^3

    Verified against the archive: box edge matches exactly on all single-solvent
    runs, solvent count on 48/53 (misses are 5-10 M cases outside our range).
    """
    if molarity <= 0:
        raise ValueError("molarity must be > 0")
    # nm^3 implied by the fixed ion count at the requested molarity
    v_target = n_salt / (molarity * N_AVOGADRO / 1e24)
    v_salt = n_salt * mw_salt / N_AVOGADRO * 1e21 / RHO_SALT
    n_solvent = round((v_target - v_salt) * RHO_SOLVENT / 1e21 * N_AVOGADRO / mw_solvent)
    clamped = n_solvent < MIN_SOLVENT
    n_solvent = max(MIN_SOLVENT, n_solvent)
    return {
        "n_solvent": int(n_solvent),
        "n_salt": int(n_salt),
        "box_build_nm": (v_target / PACK_FRACTION) ** (1.0 / 3.0),
        "v_target_nm3": v_target,
        "clamped": clamped,
    }


# ---------------------------------------------------------------- itp utils
_SECTION = re.compile(r"^\s*\[\s*([a-zA-Z_]+)\s*\]")


def split_sections(itp_text: str) -> list[tuple[str, list[str]]]:
    """Split an .itp into [(section_name, lines), ...]; preamble is ('', lines)."""
    out: list[tuple[str, list[str]]] = []
    cur_name, cur_lines = "", []
    for ln in itp_text.splitlines():
        m = _SECTION.match(ln)
        if m:
            out.append((cur_name, cur_lines))
            cur_name, cur_lines = m.group(1).lower(), []
        else:
            cur_lines.append(ln)
    out.append((cur_name, cur_lines))
    return out


def strip_comment(ln: str) -> str:
    return ln.split(";", 1)[0].rstrip()


def extract_atomtypes(itp_text: str) -> tuple[list[str], str]:
    """Pull [ atomtypes ] out of an .itp.

    Returns (atomtype_lines, itp_without_atomtypes). LigParGen emits a local
    [ atomtypes ] block inside each molecule .itp, but GROMACS only allows it
    once, at the top of the master topology.
    """
    types: list[str] = []
    kept: list[str] = []
    for name, lines in split_sections(itp_text):
        if name == "atomtypes":
            types.extend(l for l in lines if strip_comment(l).strip())
            continue
        if name:
            kept.append(f"[ {name} ]")
        kept.extend(lines)
    return types, "\n".join(kept)


def total_charge(itp_text: str) -> float:
    """Sum the charge column of an .itp's [ atoms ] section."""
    q = 0.0
    for name, lines in split_sections(itp_text):
        if name != "atoms":
            continue
        for ln in lines:
            f = strip_comment(ln).split()
            if len(f) >= 7:
                try:
                    q += float(f[6])
                except ValueError:
                    pass
    return q


# ---------------------------------------------------------------- FF integrity
#: atomic mass -> element symbol, for reading composition back out of an .itp
_MASS_ELEMENT = [
    (1.008, "H"), (6.941, "Li"), (10.811, "B"), (12.011, "C"), (14.007, "N"),
    (15.999, "O"), (18.998, "F"), (22.990, "Na"), (28.086, "Si"),
    (30.974, "P"), (32.066, "S"), (35.453, "Cl"), (79.904, "Br"),
]


def element_from_mass(mass: float) -> str:
    return min(_MASS_ELEMENT, key=lambda me: abs(me[0] - mass))[1]


def itp_formula(itp_text: str) -> dict[str, int]:
    """Element counts of the molecule described by an .itp, read from masses."""
    counts: dict[str, int] = {}
    for name, lines in split_sections(itp_text):
        if name != "atoms":
            continue
        for ln in lines:
            f = strip_comment(ln).split()
            if len(f) < 8:
                continue
            try:
                el = element_from_mass(float(f[7]))
            except ValueError:
                continue
            counts[el] = counts.get(el, 0) + 1
    return counts


def smiles_formula(smi: str) -> dict[str, int]:
    """Element counts implied by a SMILES, hydrogens included."""
    m = Chem.MolFromSmiles(smi)
    if m is None:
        raise ValueError(f"unparseable SMILES: {smi!r}")
    counts: dict[str, int] = {}
    for a in Chem.AddHs(m).GetAtoms():
        s = a.GetSymbol()
        counts[s] = counts.get(s, 0) + 1
    return counts


def check_solvent_ff(itp_text: str, smi: str) -> tuple[bool, str]:
    """Does this .itp actually describe `smi`?

    Necessary because the LigParGen server can return *another molecule's*
    files when requests overlap -- observed in practice, not hypothetical: a
    49-atom solvent came back with a 58-atom topology. An atom-count check
    alone would miss an isomer swap, so compare the full element composition.
    """
    got, want = itp_formula(itp_text), smiles_formula(smi)
    if got == want:
        return True, ""
    fmt = lambda d: "".join(f"{k}{v}" for k, v in sorted(d.items()))  # noqa: E731
    return False, f"composition mismatch: itp={fmt(got)} expected={fmt(want)}"

