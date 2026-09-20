#!/usr/bin/env python
"""Build one GROMACS run directory for a single electrolyte formulation.

Reconstructs the ce_solvation_md run layout (the original `elytefm` builder was
never backed up, so this is a re-implementation from its outputs):

    runs/el<hash>/
        sol1.itp  Li.itp|Na.itp  an.itp        force field, ECC-scaled ions
        sol1.pdb  Li.pdb|Na.pdb  an.pdb        packmol inputs
        packmol.inp  mixture.pdb  mixture.gro  packed starting box
        system.top                             merged topology
        em/nvt_heat/npt/prod .mdp              protocol (identical to the old runs)
        submit.sh                              SLURM job for Midway3
        summary.json                           formulation spec + counts

Usage:
    python setup_run.py --salt '[Li+].[N-](S(=O)(=O)F)S(=O)(=O)F' \\
                        --solvent 'COCCOC' --molarity 1.0
    python setup_run.py --from-csv index/pilot_rows.csv
"""
from __future__ import annotations

import argparse
import json
import shutil
import subprocess
import sys
from pathlib import Path

import pandas as pd

sys.path.insert(0, str(Path(__file__).resolve().parent))

import common as C  # noqa: E402
from common import (FF_IONS, FF_SOLV, RUNS, anion_name, canon,  # noqa: E402
                    cation_name, compute_counts, extract_atomtypes, fix_salt,
                    mol_weight, run_key, split_salt, split_sections,
                    strip_comment, total_charge)

SOLV_PREFIX = "sol1_"
SOLV_MOL = "sol1"
SOLV_RES = "SOL"


#: decimals used when writing charges into an .itp
Q_DECIMALS = 6


def neutralise(charges: list[float], target: float) -> list[float]:
    """Make the charges sum to `target` exactly *as written*.

    Two rounding errors compound here, and both scale with molecule count:

    1. LigParGen rounds CM1A charges to 4 decimals, leaving ~1e-4 e per
       molecule. Over a 924-solvent box that reached -0.092 e -- enough for
       grompp to warn about running Ewald on a charged system.
    2. Writing the corrected charges back out at 6 decimals reintroduces up to
       5e-7 e per atom, which over 695 molecules x 52 atoms came to -0.014 e.

    So: spread the bulk residual over all atoms (shifting each by ~1e-6 e, far
    below the precision of the charge model), round to the output precision,
    then put whatever is left on the largest-magnitude atom so the printed
    column sums exactly.
    """
    if not charges:
        return charges
    delta = (target - sum(charges)) / len(charges)
    out = [round(q + delta, Q_DECIMALS) for q in charges]
    residual = round(target - sum(out), Q_DECIMALS)
    if residual:
        i = max(range(len(out)), key=lambda k: abs(out[k]))
        out[i] = round(out[i] + residual, Q_DECIMALS)
    return out


# ------------------------------------------------------------------ topology
def prepare_solvent(itp_text: str, q_scale: float = 1.0) -> tuple[list[str], str]:
    """Namespace a LigParGen solvent .itp: UNK -> sol1/SOL, types -> sol1_*.

    The atomtype rows are left in LigParGen's own layout (name, bonded-type,
    mass, charge, ptype, sigma, eps) with only the name prefixed -- exactly what
    the old campaign's system.top contains. Solvent charges are NOT scaled;
    ECC scaling applies to ions only.
    """
    types, rest = extract_atomtypes(itp_text)
    out_types = []
    for ln in types:
        f = strip_comment(ln).split()
        if not f:
            continue
        f[0] = SOLV_PREFIX + f[0]
        out_types.append("  " + " ".join(f))

    chunks: list[str] = []
    for name, lines in split_sections(rest):
        if name in ("", "defaults", "system", "molecules"):
            continue
        if name == "moleculetype":
            chunks.append(f"[ moleculetype ]\n; name  nrexcl\n{SOLV_MOL}   3\n")
            continue
        if name == "atoms":
            parsed = [strip_comment(ln).split() for ln in lines]
            parsed = [f for f in parsed if len(f) >= 7]
            qs = neutralise([float(f[6]) * q_scale for f in parsed], 0.0)
            rows = ["[ atoms ]",
                    ";  nr    type  resnr  residue  atom  cgnr    charge      mass"]
            for f, q in zip(parsed, qs):
                nr, typ, atom, cgnr = f[0], f[1], f[4], f[5]
                mass = f[7] if len(f) > 7 else ""
                rows.append(f"{int(nr):5d}  {SOLV_PREFIX}{typ:<12s} 1  {SOLV_RES:<5s} "
                            f"{atom:<5s} {int(cgnr):4d} {q:11.6f}  {mass}")
            chunks.append("\n".join(rows) + "\n")
            continue
        body = "\n".join(lines)
        chunks.append(f"[ {name} ]\n{body}\n")
    return out_types, "\n".join(chunks)


def scale_ion(itp_text: str, q_scale: float, ion_charge: float) -> tuple[list[str], str]:
    """Split an ff_ions/*.itp and apply ECC charge scaling to [ atoms ].

    `ion_charge` is the formal charge (+1 / -1); charges are first forced to sum
    to it exactly, then scaled, so the ions contribute exactly +/-q_scale each.
    """
    types, rest = extract_atomtypes(itp_text)
    scaled_types = []
    for ln in types:
        f = strip_comment(ln).split()
        if not f:
            continue
        try:
            p = f.index("A")
            f[p - 1] = f"{float(f[p - 1]) * q_scale:.6f}"
        except (ValueError, IndexError):
            pass
        scaled_types.append("  " + " ".join(f))

    chunks: list[str] = []
    for name, lines in split_sections(rest):
        if name in ("", "defaults", "system", "molecules"):
            continue
        if name == "atoms":
            parsed = [strip_comment(ln).split() for ln in lines]
            parsed = [f for f in parsed if len(f) >= 7]
            # neutralise AFTER scaling, so the written column sums to exactly
            # +/-q_scale rather than to a rounded multiple of it
            qs = neutralise([float(f[6]) * q_scale for f in parsed],
                            ion_charge * q_scale)
            rows = ["[ atoms ]",
                    ";  nr    type  resnr  residue  atom  cgnr    charge      mass"]
            for f, q in zip(parsed, qs):
                nr, typ, res, atom, cgnr = f[0], f[1], f[3], f[4], f[5]
                mass = f[7] if len(f) > 7 else ""
                rows.append(f"{int(nr):5d}  {typ:<12s} 1  {res:<5s} {atom:<5s} "
                            f"{int(cgnr):4d} {q:11.6f}  {mass}")
            chunks.append("\n".join(rows) + "\n")
            continue
        chunks.append(f"[ {name} ]\n" + "\n".join(lines) + "\n")
    return scaled_types, "\n".join(chunks)


# ------------------------------------------------------------------ mdp files
def mdp_files(t_k: float, prod_ns: float) -> dict[str, str]:
    """The four .mdp stages, byte-for-byte the protocol of the old campaign."""
    return {
        "em.mdp": (
            "integrator       = steep\n"
            "emtol            = 1000.0\n"
            "emstep           = 0.01\n"
            "nsteps           = 5000\n"
            "nstlist          = 10\n"
            "cutoff-scheme    = Verlet\n"
            "coulombtype      = PME\n"
            "rcoulomb         = 1.2\n"
            "rvdw             = 1.2\n"
            "pbc              = xyz\n"),
        "nvt_heat.mdp": (
            "integrator           = md\n"
            "dt                   = 0.001\n"
            "nsteps               = 100000\n"
            "nstlist              = 20\n"
            "cutoff-scheme        = Verlet\n"
            "coulombtype          = PME\n"
            "rcoulomb             = 1.2\n"
            "rvdw                 = 1.2\n"
            "constraints          = h-bonds\n"
            "constraint_algorithm = LINCS\n"
            "gen_vel              = yes\n"
            "gen_temp             = 50.0\n"
            "gen_seed             = 73921\n"
            "tcoupl               = v-rescale\n"
            "tc-grps              = System\n"
            "tau_t                = 0.5\n"
            f"ref_t                = {t_k}\n"
            "nstenergy            = 5000\n"
            "nstlog               = 5000\n"
            "annealing            = single\n"
            "annealing_npoints    = 2\n"
            "annealing_time       = 0 100\n"
            f"annealing_temp       = 50 {t_k}\n"),
        "npt.mdp": (
            "integrator           = md\n"
            "dt                   = 0.001\n"
            "nsteps               = 2000000\n"
            "nstlist              = 20\n"
            "cutoff-scheme        = Verlet\n"
            "coulombtype          = PME\n"
            "rcoulomb             = 1.2\n"
            "rvdw                 = 1.2\n"
            "constraints          = h-bonds\n"
            "constraint_algorithm = LINCS\n"
            "gen_vel              = no\n"
            "tcoupl               = v-rescale\n"
            "tc-grps              = System\n"
            "tau_t                = 0.5\n"
            f"ref_t                = {t_k}\n"
            "pcoupl               = C-rescale\n"
            "pcoupltype           = isotropic\n"
            "tau_p                = 2.0\n"
            "ref_p                = 1.0\n"
            "compressibility      = 4.5e-5\n"
            "refcoord_scaling     = com\n"
            "nstenergy            = 1000\n"
            "nstlog               = 5000\n"),
        "prod.mdp": (
            f"; NVT production {prod_ns} ns -> solvation structure + ion MSD\n"
            "integrator           = md\n"
            "dt                   = 0.001\n"
            f"nsteps               = {int(prod_ns * 1_000_000)}\n"
            "nstlist              = 20\n"
            "cutoff-scheme        = Verlet\n"
            "coulombtype          = PME\n"
            "rcoulomb             = 1.2\n"
            "rvdw                 = 1.2\n"
            "constraints          = h-bonds\n"
            "constraint_algorithm = LINCS\n"
            "gen_vel              = no\n"
            "tcoupl               = v-rescale\n"
            "tc-grps              = System\n"
            "tau_t                = 0.5\n"
            f"ref_t                = {t_k}\n"
            "nstxout-compressed     = 5000\n"
            "compressed-x-precision = 10000\n"
            "nstenergy            = 1000\n"
            "nstlog               = 10000\n"),
    }


def submit_script(key: str, rundir: Path, hours: int, cation_mol: str) -> str:
    """Per-run SLURM wrapper. The protocol itself lives in scripts/run_md.sh so
    that single submissions and the job-array driver cannot drift apart."""
    scripts = Path(__file__).resolve().parent
    return f"""#!/bin/bash
#SBATCH --job-name=md-{key}
#SBATCH --account={C.SLURM_ACCOUNT}
#SBATCH --partition={C.SLURM_PARTITION}
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=8
#SBATCH --gres=gpu:1
#SBATCH --time={hours:02d}:00:00
#SBATCH --output={key}_%j.out
#SBATCH --error={key}_%j.err

bash {scripts}/run_md.sh {rundir}
"""


# ------------------------------------------------------------------ build
def build(salt_smi: str, solvent_smi: str, molarity: float, *,
          q_scale: float = C.Q_SCALE_DEFAULT, n_salt: int = C.N_SALT_DEFAULT,
          prod_ns: float = C.PROD_NS, t_k: float = C.T_K,
          hours: int = 8, force: bool = False,
          outroot: Path = RUNS) -> tuple[Path, dict]:
    salt_smi = fix_salt(salt_smi)
    solvent_smi = canon(solvent_smi)
    key = run_key(salt_smi, solvent_smi, molarity)
    d = outroot / key
    if (d / "mixture.gro").exists() and not force:
        return d, json.loads((d / "summary.json").read_text())
    d.mkdir(parents=True, exist_ok=True)

    cat_smi, an_smi = split_salt(salt_smi)
    cat, an = cation_name(cat_smi), anion_name(an_smi)
    skey = C.solvent_key(solvent_smi)

    sol_itp_src = FF_SOLV / f"{skey}.itp"
    if not sol_itp_src.exists():
        raise FileNotFoundError(
            f"no solvent FF for {solvent_smi} ({skey}); run fetch_solvent_ff.py")
    # Pre-flight: confirm the cached .itp really describes this SMILES. LigParGen
    # has been seen returning a different molecule's parameters, which produces
    # a perfectly valid topology for the wrong chemistry -- fail loudly here
    # rather than 4 ns later.
    ok, why = C.check_solvent_ff(sol_itp_src.read_text(), solvent_smi)
    if not ok:
        raise ValueError(f"solvent FF {skey} does not match {solvent_smi}: {why}; "
                         f"run verify_solvent_ff.py --fix")
    for p in (FF_IONS / f"{cat}.itp", FF_IONS / f"{an}.itp"):
        if not p.exists():
            raise FileNotFoundError(f"missing ion FF {p}; run build_ion_ff.py")

    # ---- counts and box
    counts = compute_counts(mol_weight(solvent_smi), mol_weight(salt_smi),
                            molarity, n_salt)
    n_sol, box_nm = counts["n_solvent"], counts["box_build_nm"]

    # ---- topology pieces
    sol_types, sol_body = prepare_solvent(sol_itp_src.read_text())
    cat_types, cat_body = scale_ion((FF_IONS / f"{cat}.itp").read_text(),
                                    q_scale, +1.0)
    an_types, an_body = scale_ion((FF_IONS / f"{an}.itp").read_text(),
                                  q_scale, -1.0)

    (d / "sol1.itp").write_text(sol_body)
    (d / f"{cat}.itp").write_text(cat_body)
    (d / "an.itp").write_text(an_body)
    shutil.copy(FF_SOLV / f"{skey}.pdb", d / "sol1.pdb")
    shutil.copy(FF_IONS / f"{cat}.pdb", d / f"{cat}.pdb")
    shutil.copy(FF_IONS / f"{an}.pdb", d / "an.pdb")

    # the solvent pdb is still LigParGen's UNK; make it match sol1.itp
    p = d / "sol1.pdb"
    p.write_text("\n".join(
        (ln.ljust(80)[:17] + f"{SOLV_RES:>3s}" + ln.ljust(80)[20:]).rstrip()
        if ln.startswith(("ATOM", "HETATM")) else ln
        for ln in p.read_text().splitlines()) + "\n")

    all_types = cat_types + an_types + sol_types
    (d / "system.top").write_text(
        f"; {key}  (q_scale={q_scale})\n"
        f"; salt={salt_smi}  solvent={solvent_smi}  M={molarity}\n\n"
        "[ defaults ]\n"
        "; nbfunc  comb-rule  gen-pairs  fudgeLJ  fudgeQQ\n"
        "1         3          yes        0.5      0.5\n\n"
        "[ atomtypes ]\n"
        ";  name    at.nr/btype   mass    charge  ptype  sigma(nm)  epsilon(kJ/mol)\n"
        + "\n".join(all_types) + "\n\n"
        f'#include "sol1.itp"\n#include "{cat}.itp"\n#include "an.itp"\n\n'
        f"[ system ]\n{key} electrolyte\n\n"
        "[ molecules ]\n"
        f"{SOLV_MOL:<12s} {n_sol}\n"
        f"{cat + '+':<12s} {n_salt}\n"
        f"{'an':<12s} {n_salt}\n")

    # ---- packmol
    box_a = box_nm * 10.0
    (d / "packmol.inp").write_text(
        f"tolerance 2.5\nfiletype pdb\noutput {d/'mixture.pdb'}\n"
        f"seed 12345\nnloop 20\nnloop0 30\n"
        f"structure {d/'sol1.pdb'}\n  number {n_sol}\n"
        f"  inside box 1.000 1.000 1.000 {box_a-1:.3f} {box_a-1:.3f} {box_a-1:.3f}\n"
        f"end structure\n"
        f"structure {d/(cat + '.pdb')}\n  number {n_salt}\n"
        f"  inside box 1.000 1.000 1.000 {box_a-1:.3f} {box_a-1:.3f} {box_a-1:.3f}\n"
        f"end structure\n"
        f"structure {d/'an.pdb'}\n  number {n_salt}\n"
        f"  inside box 1.000 1.000 1.000 {box_a-1:.3f} {box_a-1:.3f} {box_a-1:.3f}\n"
        f"end structure\n")

    # Packing is by far the slowest step (~1-3 min) and depends only on the
    # molecule counts and box, not on charges -- so a topology-only rebuild
    # reuses an existing box whose atom count still matches.
    n_expect = (n_sol * C.n_atoms_with_h(solvent_smi)
                + n_salt * (len(_atom_rows(cat_body)) + len(_atom_rows(an_body))))
    existing = (d / "mixture.pdb")
    reuse = (existing.exists() and sum(
        1 for ln in existing.read_text().splitlines()
        if ln.startswith(("ATOM", "HETATM"))) == n_expect)

    if not reuse:
        with (d / "packmol.inp").open() as fh:
            r = subprocess.run([str(C.PACKMOL)], cwd=d, stdin=fh,
                               capture_output=True, text=True)
        (d / "packmol.log").write_text(r.stdout + r.stderr)
    if not (d / "mixture.pdb").exists():
        raise RuntimeError(f"packmol produced no mixture.pdb (see {d/'packmol.log'})")
    # "ENDED WITHOUT PERFECT PACKING" is expected: we deliberately build the box
    # at 65% of target density and let NPT compress it. Only the atom count matters.
    n_pdb = sum(1 for ln in (d / "mixture.pdb").read_text().splitlines()
                if ln.startswith(("ATOM", "HETATM")))

    # ---- pdb -> gro with an explicit cubic box
    g = subprocess.run(
        ["bash", "-lc",
         f"module load {C.GMX_MODULE} >/dev/null 2>&1; cd {d} && "
         f"{C.GMX} editconf -f mixture.pdb -o mixture.gro "
         f"-box {box_nm:.4f} {box_nm:.4f} {box_nm:.4f}"],
        capture_output=True, text=True)
    if not (d / "mixture.gro").exists():
        raise RuntimeError(f"editconf failed: {g.stderr[-500:]}")

    for fn, txt in mdp_files(t_k, prod_ns).items():
        (d / fn).write_text(txt)
    sub = d / "submit.sh"
    sub.write_text(submit_script(key, d, hours, cat))
    sub.chmod(0o755)

    net_q = (n_sol * total_charge(sol_body)
             + n_salt * total_charge(cat_body)
             + n_salt * total_charge(an_body))

    summary = {
        "name": key,
        "T_K": t_k, "P_bar": C.P_BAR, "prod_ns": prod_ns, "q_scale": q_scale,
        "counts": {"sol1": n_sol, "cat": n_salt, "an": n_salt},
        "box_build_nm": round(box_nm, 4),
        "n_atoms": n_pdb,
        "n_atoms_expected": n_expect,
        "n_atomtypes": len(all_types),
        "net_charge": round(net_q, 4),
        "cation": cat, "anion": an,
        "spec": {
            "solvents": [{"smi": solvent_smi, "frac": 1.0}],
            "salt_smi": salt_smi, "salt_M": float(molarity),
            "T_K": t_k, "prod_ns": prod_ns, "n_salt": n_salt,
            "q_scale": q_scale, "name": key,
        },
        "comp_info": {"sol1": {"smi": solvent_smi, "frac": 1.0,
                               "mw": round(mol_weight(solvent_smi), 3),
                               "rho": C.RHO_SOLVENT,
                               "ff_key": skey}},
        "clamped_solvent_count": counts["clamped"],
    }
    (d / "summary.json").write_text(json.dumps(summary, indent=2))
    return d, summary


def _atom_rows(itp_body: str) -> list[str]:
    for name, lines in split_sections(itp_body):
        if name == "atoms":
            return [l for l in lines if strip_comment(l).split()]
    return []


def main() -> int:
    ap = argparse.ArgumentParser()
    ap.add_argument("--salt")
    ap.add_argument("--solvent")
    ap.add_argument("--molarity", type=float)
    ap.add_argument("--from-csv", type=Path)
    ap.add_argument("--q-scale", type=float, default=C.Q_SCALE_DEFAULT)
    ap.add_argument("--n-salt", type=int, default=C.N_SALT_DEFAULT)
    ap.add_argument("--prod-ns", type=float, default=C.PROD_NS)
    ap.add_argument("--hours", type=int, default=8)
    ap.add_argument("--force", action="store_true")
    args = ap.parse_args()

    if args.from_csv:
        rows = pd.read_csv(args.from_csv).to_dict("records")
    elif args.salt and args.solvent and args.molarity:
        rows = [{"salt": args.salt, "solvent": args.solvent,
                 "concentration": args.molarity}]
    else:
        ap.error("need --from-csv, or --salt/--solvent/--molarity")

    bad = 0
    for r in rows:
        try:
            d, s = build(r["salt"], r["solvent"], float(r["concentration"]),
                         q_scale=args.q_scale, n_salt=args.n_salt,
                         prod_ns=args.prod_ns, hours=args.hours,
                         force=args.force)
        except Exception as e:                                # noqa: BLE001
            print(f"[FAIL] {r.get('salt','?')[:24]} / {str(r.get('solvent'))[:28]}: "
                  f"{type(e).__name__}: {e}")
            bad += 1
            continue
        ok = s["n_atoms"] == s["n_atoms_expected"] and abs(s["net_charge"]) < 0.05
        print(f"[{'ok' if ok else '??'}] {s['name']}  {s['cation']}{s['anion']:<5s} "
              f"{s['spec']['salt_M']:>4.1f}M  nsol={s['counts']['sol1']:<5d} "
              f"box={s['box_build_nm']:.2f}nm  atoms={s['n_atoms']:<6d} "
              f"q={s['net_charge']:+.3f}"
              + ("" if ok else "   <-- CHECK"))
        if not ok:
            bad += 1
    return 1 if bad else 0


if __name__ == "__main__":
    raise SystemExit(main())
