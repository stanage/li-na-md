#!/usr/bin/env python
"""Assemble the ion force-field library for the 2000-electrolyte campaign.

Produces, once, a uniform self-contained .itp + .pdb per ion in ff_ions/:

    ff_ions/<NAME>.itp   [atomtypes] (namespaced) + [moleculetype] + bonded terms
    ff_ions/<NAME>.pdb   single-molecule geometry for packmol
    ff_ions/ions_manifest.json

Charges are stored UNSCALED (anion sums to -1, cation to +1). ECC scaling
(q = 0.8, as in the old campaign) is applied later by setup_run.py, so q_scale
stays a run-time knob.

Why four different sources
--------------------------
LigParGen cannot parameterise hypervalent inorganic anions (BF4, PF6, FSI, ...)
-- the Yale server rejects them outright -- so no single source covers all 20
anions in this dataset. Each anion is taken from the best available source:

  2009IL     Acevedo OPLS-2009IL (16 anions). Validated IL parameters, but ships
             atomtypes in a separate file and with NO [pairs] section -- both
             fixed here (see below).
  oldrun     FSI, lifted from the ce_solvation_md runs exactly as simulated.
  fftool     TFA and Beti, generated from the CL&P database via Padua's fftool.
  ligpargen  OTs (tosylate), which LigParGen handles fine as an organic anion.
  builtin    Li+ and Na+, Aqvist parameters as listed in CL&P il.ff. Li+ matches
             the old campaign's Li.itp exactly.

Two corrections applied to the 2009IL files
-------------------------------------------
1. [pairs] are generated from the bond graph (all atom pairs exactly 3 bonds
   apart). The Acevedo distribution omits them, which would silently drop the
   OPLS 1-4 scaled LJ/Coulomb terms for acetate, propionate, benzoate,
   mesylate, TFSI and triflate.
2. [atomtypes] rows are normalised to the 7-column GROMACS form
   (name at.nr mass charge ptype sigma epsilon) so that ion and LigParGen
   solvent atomtypes can share one block without the parser guessing formats.

Usage:  python build_ion_ff.py [--only NAME[,NAME...]] [--force]
"""
from __future__ import annotations

import argparse
import collections
import json
import shutil
import subprocess
import sys
import tempfile
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parent))

import common as C  # noqa: E402
from common import (ANION_NAMES, FF_IONS, MD, canon, split_sections,  # noqa: E402
                    strip_comment)

IL2009 = MD / "Force Fields" / "2009IL FF Orlando" / "2009IL"
IL_ITP = IL2009 / "ITP"
IL_PDB = MD / "Force Fields" / "2009IL FF Orlando" / "PDB"
CLANDP = MD / "Force Fields" / "clandp-master"
FFTOOL = (MD / "ce_solvation_md" / "simulation_workflow-main" /
          "simulation_workflow-main" / "tools" / "fftool" / "fftool")

ANION_PREFIX = "an_"       # namespace for anion atom types, as in the old runs
ANION_MOLNAME = "an"       # moleculetype name, as in the old runs

# mass -> atomic number, for normalising atomtype rows
_MASS_Z = [(1.008, 1), (6.941, 3), (10.811, 5), (12.011, 6), (14.007, 7),
           (15.999, 8), (18.998, 9), (22.990, 11), (30.974, 15), (32.066, 16),
           (35.453, 17), (79.904, 35)]


def atomic_number(mass: float) -> int:
    return min(_MASS_Z, key=lambda mz: abs(mz[0] - mass))[1]


def resname_of(name: str) -> str:
    """Residue name used in BOTH the .itp and the .pdb for an ion.

    Three characters, because that is all the PDB resName field holds; keeping
    the two in sync avoids grompp name-mismatch notes later. All 20 anion names
    stay distinct under this truncation.
    """
    return name.lower()[:3]


# --------------------------------------------------------------- sourcing map
# name -> (source, key). Keys are the identifiers in each upstream library.
SOURCES: dict[str, tuple[str, str]] = {
    "OAc": ("2009IL", "ACE"),
    "BF4": ("2009IL", "BF4"),
    "BNZ": ("2009IL", "BNZ"),
    "Br": ("2009IL", "Br"),
    "Cl": ("2009IL", "Cl"),
    "ClO4": ("2009IL", "ClO4"),
    "DCA": ("2009IL", "DCA"),
    "HCOO": ("2009IL", "HCOO"),
    "MS": ("2009IL", "MS"),
    "NO3": ("2009IL", "NO3"),
    "TFSI": ("2009IL", "NTF2"),
    "PF6": ("2009IL", "PF6"),
    "PROP": ("2009IL", "PROP"),
    "SCN": ("2009IL", "SCN"),
    "TCM": ("2009IL", "TCM"),
    "OTf": ("2009IL", "TFO"),
    "FSI": ("oldrun", "an"),
    "TFA": ("fftool", "tfa"),
    "Beti": ("fftool", "beti"),
    "OTs": ("ligpargen", "[O-]S(=O)(=O)c1ccc(C)cc1"),
}

#: Aqvist cation parameters, as tabulated in CL&P il.ff.
#: name -> (element, mass, sigma_nm, epsilon_kJ)
CATIONS = {
    "Li": ("Li", 6.941, 0.2126, 0.07648),
    "Na": ("Na", 22.990, 0.3330, 0.01160),
}

SMILES_OF = {v: k for k, v in ANION_NAMES.items()}


# --------------------------------------------------------------- itp helpers
def parse_bonds(sections) -> list[tuple[int, int]]:
    out = []
    for name, lines in sections:
        if name != "bonds":
            continue
        for ln in lines:
            f = strip_comment(ln).split()
            if len(f) >= 2:
                try:
                    out.append((int(f[0]), int(f[1])))
                except ValueError:
                    pass
    return out


def n_atoms_of(sections) -> int:
    n = 0
    for name, lines in sections:
        if name != "atoms":
            continue
        for ln in lines:
            if strip_comment(ln).split():
                n += 1
    return n


def generate_pairs(bonds, natoms) -> list[tuple[int, int]]:
    """All atom pairs separated by exactly 3 bonds (OPLS 1-4 pairs)."""
    adj = collections.defaultdict(set)
    for i, j in bonds:
        adj[i].add(j)
        adj[j].add(i)
    pairs = set()
    for start in range(1, natoms + 1):
        # BFS to depth 3
        dist = {start: 0}
        queue = collections.deque([start])
        while queue:
            u = queue.popleft()
            if dist[u] == 3:
                continue
            for v in adj[u]:
                if v not in dist:
                    dist[v] = dist[u] + 1
                    queue.append(v)
        for v, d in dist.items():
            if d == 3:
                pairs.add((min(start, v), max(start, v)))
    return sorted(pairs)


def normalise_atomtypes(lines: list[str], prefix: str,
                        charge_scale: float = 1.0) -> list[str]:
    """Rewrite atomtype rows to `name at.nr mass charge ptype sigma epsilon`.

    Accepts the two layouts we encounter: the 6-column Acevedo form
    (name mass charge ptype sigma eps) and the 7-column form that already
    carries the atomic number.
    """
    out = []
    for ln in lines:
        f = strip_comment(ln).split()
        if not f:
            continue
        name = f[0]
        # locate ptype 'A' to disambiguate the layout
        try:
            p = f.index("A")
        except ValueError:
            continue
        sigma, eps = float(f[p + 1]), float(f[p + 2])
        charge = float(f[p - 1]) * charge_scale
        mass = float(f[p - 2])
        z = atomic_number(mass)
        out.append(f"  {prefix}{name:<12s} {z:3d} {mass:9.4f} {charge:9.4f}  A  "
                   f"{sigma:.5e} {eps:.5e}")
    return out


def rewrite_molecule(sections, prefix: str, molname: str, resname: str,
                     add_pairs: bool, charge_scale: float = 1.0) -> str:
    """Emit the bonded part of an .itp with namespaced types and fixed names."""
    natoms = n_atoms_of(sections)
    bonds = parse_bonds(sections)
    has_pairs = any(n == "pairs" for n, _ in sections)

    chunks: list[str] = []
    for name, lines in sections:
        if name in ("", "atomtypes", "defaults", "system", "molecules"):
            continue
        if name == "moleculetype":
            chunks.append("[ moleculetype ]\n; name  nrexcl\n"
                          f"{molname}   3\n")
            continue
        if name == "atoms":
            rows = ["[ atoms ]",
                    ";  nr    type  resnr  residue  atom  cgnr    charge      mass"]
            for ln in lines:
                f = strip_comment(ln).split()
                if len(f) < 7:
                    continue
                nr, typ, _resnr, _res, atom, cgnr, q = f[:7]
                mass = f[7] if len(f) > 7 else ""
                rows.append(f"{int(nr):5d}  {prefix}{typ:<10s} 1  {resname:<5s} "
                            f"{atom:<5s} {int(cgnr):4d} "
                            f"{float(q) * charge_scale:11.6f}  {mass}")
            chunks.append("\n".join(rows) + "\n")
            continue
        body = "\n".join(l for l in lines if strip_comment(l).strip() or l.strip().startswith(";"))
        chunks.append(f"[ {name} ]\n{body}\n")

    if add_pairs and not has_pairs:
        gp = generate_pairs(bonds, natoms)
        if gp:
            rows = ["[ pairs ]",
                    "; generated from the bond graph (3-bond separation);",
                    "; the OPLS-2009IL distribution omits this section"]
            rows += [f"{i:5d} {j:5d}     1" for i, j in gp]
            chunks.append("\n".join(rows) + "\n")
    return "\n".join(chunks)


# --------------------------------------------------------------- pdb helpers
def clean_pdb(src_text: str, resname: str) -> str:
    """Keep only ATOM/HETATM records and stamp a consistent residue name."""
    out = []
    n = 0
    for ln in src_text.splitlines():
        if ln.startswith(("ATOM", "HETATM")):
            n += 1
            ln = ln.ljust(80)
            ln = ln[:17] + f"{resname[:3]:>3s}" + ln[20:]
            out.append(ln.rstrip())
    out.append("END")
    return "\n".join(out) + "\n"


def xyz_to_pdb(xyz_text: str, resname: str,
               names: list[str] | None = None) -> str:
    """fftool .xyz -> .pdb, taking atom names from the topology when given.

    The xyz carries only element symbols, so naming the PDB atoms from it
    yields C, C, F, F ... while the .itp says C1, C2, F1, F2 ...  grompp then
    reports "non-matching atom names" for every atom of every such molecule
    and silently keeps the topology's. Harmless in itself, but it costs a
    grompp warning per stage, which is exactly the budget a low -maxwarn
    needs for real problems. Passing the topology names makes the two agree.
    """
    lines = xyz_text.splitlines()
    n = int(lines[0].split()[0])
    out = []
    for i, ln in enumerate(lines[2:2 + n], start=1):
        f = ln.split()
        el, x, y, z = f[0], float(f[1]), float(f[2]), float(f[3])
        nm = names[i - 1] if names and i <= len(names) else el
        out.append(f"HETATM{i:5d} {nm[:4]:<4s}{resname[:3]:>3s}"
                   f"     1    {x:8.3f}{y:8.3f}{z:8.3f}  1.00  0.00"
                   f"          {el[:2]:>2s}")
    out.append("END")
    return "\n".join(out) + "\n"


def atom_names(itp_text: str) -> list[str]:
    """Atom names from an .itp [atoms] block, in file order."""
    out = []
    for sname, lines in split_sections(itp_text):
        if sname != "atoms":
            continue
        for ln in lines:
            f = strip_comment(ln).split()
            if len(f) >= 5 and f[0].isdigit():
                out.append(f[4])
    return out


# --------------------------------------------------------------- builders
def build_from_2009il(name: str, key: str) -> tuple[str, str]:
    itp = (IL_ITP / f"{key}.itp").read_text()
    at_file = IL_ITP / f"{key}_atomtypes.itp"
    at_lines: list[str] = []
    for sname, lines in split_sections(at_file.read_text()):
        if sname == "atomtypes":
            at_lines = lines
    sections = split_sections(itp)
    types = normalise_atomtypes(at_lines, ANION_PREFIX)
    body = rewrite_molecule(sections, ANION_PREFIX, ANION_MOLNAME,
                            resname_of(name), add_pairs=True)
    pdb = clean_pdb((IL_PDB / f"{key}.pdb").read_text(), resname_of(name))
    return _assemble(name, types, body), pdb


def build_from_oldrun(name: str, key: str) -> tuple[str, str]:
    """FSI, lifted from a completed ce_solvation_md run.

    Those runs stored the anion ALREADY ECC-scaled (charges sum to -0.8), so we
    divide the archived charges back out to keep this library uniformly
    unscaled. setup_run.py re-applies whatever q_scale the campaign uses.
    """
    unscale = 1.0 / C.Q_SCALE_DEFAULT
    run = next(d for d in sorted(C.OLD_RUNS.glob("ce*"))
               if (d / "an.itp").exists() and (d / "system.top").exists())
    itp = (run / "an.itp").read_text()
    # its atomtypes live in that run's system.top, already `an_`-prefixed
    top = (run / "system.top").read_text()
    raw = []
    for sname, lines in split_sections(top):
        if sname == "atomtypes":
            raw = [l for l in lines
                   if strip_comment(l).split() and
                   strip_comment(l).split()[0].startswith(ANION_PREFIX) and
                   strip_comment(l).split()[0] != "an_Li"]
    types = normalise_atomtypes([l.replace(ANION_PREFIX, "", 1) for l in raw],
                                ANION_PREFIX, charge_scale=unscale)
    body = rewrite_molecule(split_sections(itp), "", ANION_MOLNAME,
                            resname_of(name), add_pairs=False,
                            charge_scale=unscale)
    pdb = clean_pdb((run / "an.pdb").read_text(), resname_of(name))
    return _assemble(name, types, body), pdb


def build_from_fftool(name: str, key: str) -> tuple[str, str]:
    """TFA / Beti, generated from the CL&P database through fftool."""
    with tempfile.TemporaryDirectory() as td:
        w = Path(td)
        for f in (f"{key}.zmat", "il.ff"):
            shutil.copy(CLANDP / f, w / f)
        run = lambda *a: subprocess.run(  # noqa: E731
            [str(C.PYTHON), str(FFTOOL), *a], cwd=w, check=True,
            capture_output=True, text=True)
        run("1", f"{key}.zmat", "--box", "30")
        subprocess.run([str(C.PACKMOL)], cwd=w, check=True,
                       stdin=(w / "pack.inp").open(), capture_output=True)
        run("1", f"{key}.zmat", "--box", "30", "--gmx")
        top = (w / "field.top").read_text()
        xyz = (w / f"{key}_pack.xyz").read_text()
    at_lines: list[str] = []
    for sname, lines in split_sections(top):
        if sname == "atomtypes":
            at_lines = lines
    types = normalise_atomtypes(at_lines, ANION_PREFIX)
    body = rewrite_molecule(split_sections(top), ANION_PREFIX, ANION_MOLNAME,
                            resname_of(name), add_pairs=False)
    assembled = _assemble(name, types, body)
    return assembled, xyz_to_pdb(xyz, resname_of(name),
                                 atom_names(assembled))


def build_from_ligpargen(name: str, smiles: str) -> tuple[str, str]:
    from ligpargen import fetch
    cache = FF_IONS / "_lpg"
    cache.mkdir(parents=True, exist_ok=True)
    prefix = cache / name
    fetch(smiles, prefix, net_charge=-1, verbose=True)
    itp = prefix.with_suffix(".itp").read_text()
    at_lines: list[str] = []
    for sname, lines in split_sections(itp):
        if sname == "atomtypes":
            at_lines = lines
    types = normalise_atomtypes(at_lines, ANION_PREFIX)
    body = rewrite_molecule(split_sections(itp), ANION_PREFIX, ANION_MOLNAME,
                            resname_of(name), add_pairs=False)
    pdb = clean_pdb(prefix.with_suffix(".pdb").read_text(), resname_of(name))
    return _assemble(name, types, body), pdb


def build_cation(name: str) -> tuple[str, str]:
    el, mass, sigma, eps = CATIONS[name]
    types = [f"  {el:<12s} {atomic_number(mass):3d} {mass:9.4f} {1.0:9.4f}  A  "
             f"{sigma:.5e} {eps:.5e}"]
    body = (f"[ moleculetype ]\n; name  nrexcl\n{el}+   3\n\n"
            "[ atoms ]\n"
            ";  nr    type  resnr  residue  atom  cgnr    charge      mass\n"
            f"    1  {el:<10s} 1  {resname_of(name):<5s} {el:<5s}    1    1.000000  {mass:.4f}\n")
    pdb = (f"HETATM    1 {el:<4s}{resname_of(name):>3s}     1       0.000   0.000"
           f"   0.000  1.00  0.00          {el:>2s}\nEND\n")
    return _assemble(name, types, body), pdb


def _assemble(name: str, types: list[str], body: str) -> str:
    head = (f"; {name} -- ion parameters for the ht1000 campaign\n"
            f"; assembled by build_ion_ff.py; charges UNSCALED "
            f"(ECC scaling applied at run-build time)\n\n"
            "[ atomtypes ]\n"
            ";  name       at.nr   mass    charge  ptype  sigma(nm)   epsilon(kJ/mol)\n")
    return head + "\n".join(types) + "\n\n" + body


BUILDERS = {
    "2009IL": build_from_2009il,
    "oldrun": build_from_oldrun,
    "fftool": build_from_fftool,
    "ligpargen": build_from_ligpargen,
}


def main() -> int:
    ap = argparse.ArgumentParser()
    ap.add_argument("--only", help="comma-separated ion names")
    ap.add_argument("--force", action="store_true")
    args = ap.parse_args()

    FF_IONS.mkdir(parents=True, exist_ok=True)
    wanted = set(args.only.split(",")) if args.only else None
    manifest: dict[str, dict] = {}
    mpath = FF_IONS / "ions_manifest.json"
    if mpath.exists():
        manifest = json.loads(mpath.read_text())

    targets: list[tuple[str, str, str]] = [
        *((n, *SOURCES[n]) for n in SOURCES),
        *((n, "builtin", n) for n in CATIONS),
    ]
    failures = []
    for name, source, key in targets:
        if wanted and name not in wanted:
            continue
        itp_path = FF_IONS / f"{name}.itp"
        if itp_path.exists() and not args.force:
            print(f"[skip] {name}")
            continue
        try:
            if source == "builtin":
                itp, pdb = build_cation(name)
            else:
                itp, pdb = BUILDERS[source](name, key)
        except Exception as e:                                # noqa: BLE001
            print(f"[FAIL] {name} ({source}): {type(e).__name__}: {e}")
            failures.append(name)
            continue

        itp_path.write_text(itp)
        (FF_IONS / f"{name}.pdb").write_text(pdb)
        q = C.total_charge(itp)
        expect = 1.0 if source == "builtin" else -1.0
        natoms = sum(1 for l in pdb.splitlines() if l.startswith(("ATOM", "HETATM")))
        nat_itp = n_atoms_of(split_sections(itp))
        ok = abs(q - expect) < 0.02 and natoms == nat_itp
        manifest[name] = {
            "source": source, "key": key,
            "smiles": SMILES_OF.get(name, {"Li": "[Li+]", "Na": "[Na+]"}.get(name, "")),
            "n_atoms": nat_itp, "charge_unscaled": round(q, 5),
            "pdb_atoms": natoms, "ok": ok,
        }
        flag = "" if ok else "   <-- CHECK"
        print(f"[ok]   {name:<6s} {source:<10s} atoms={nat_itp:<3d} q={q:+.4f}{flag}")
        if not ok:
            failures.append(name)

    mpath.write_text(json.dumps(manifest, indent=2, sort_keys=True))
    print(f"\nmanifest -> {mpath}")
    if failures:
        print(f"PROBLEMS: {sorted(set(failures))}")
        return 1
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
