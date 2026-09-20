#!/usr/bin/env python
"""Gather every run's summary/solvation/density into one results table.

Writes index/results.csv (one row per run) and prints a status breakdown. Safe
to run at any point during the campaign -- unfinished runs appear with their
stage and empty analysis columns.

Usage:  python collect_results.py [--out ../index/results.csv]
"""
from __future__ import annotations

import argparse
import json
import re
import sys
from pathlib import Path

import pandas as pd

sys.path.insert(0, str(Path(__file__).resolve().parent))

from common import DATASET, INDEX, N_AVOGADRO, RUNS, run_key  # noqa: E402

STAGES = ["em", "nvt_heat", "npt", "prod"]


def dataset_rows() -> dict[str, int]:
    """run key -> 0-based row number in the source CSV.

    Run directories are named by a content hash of (salt, solvent, molarity),
    which is stable but unreadable; this recovers the link back to the dataset
    so a result can be traced to the row that produced it.
    """
    if not DATASET.exists():
        print(f"warning: {DATASET} not found; row_index will be blank")
        return {}
    df = pd.read_csv(DATASET)
    out: dict[str, int] = {}
    for i, r in df.iterrows():
        out[run_key(r["salt"], r["solvent"], float(r["concentration"]))] = int(i)
    return out


def box_volume_nm3(gro: Path) -> float | None:
    """Box volume from the last line of a .gro file.

    The line holds v1x v2y v3z and, for triclinic cells, six more off-diagonal
    terms; the volume is the determinant, which for the upper-triangular form
    GROMACS writes is just the product of the diagonal.
    """
    if not gro.exists():
        return None
    last = gro.read_text().rstrip().rsplit("\n", 1)[-1].split()
    try:
        v = [float(x) for x in last[:3]]
    except (ValueError, IndexError):
        return None
    return v[0] * v[1] * v[2] if len(v) == 3 else None


def achieved_molarity(d: Path, n_salt: int) -> tuple[float | None, float | None]:
    """(molarity, box edge) actually sampled, from the equilibrated box.

    Production runs NVT, so its box is fixed at whatever NPT settled on -- that
    is the concentration the trajectory really represents. It differs from the
    requested molarity whenever the assumed component densities behind
    compute_counts() were off, and NPT cannot correct that: it relaxes the
    volume, not the molecule counts.
    """
    vol = box_volume_nm3(d / "prod.gro") or box_volume_nm3(d / "npt.gro")
    if not vol:
        return None, None
    return n_salt / (vol * N_AVOGADRO / 1e24), vol ** (1.0 / 3.0)


def mean_density(xvg: Path) -> float | None:
    """Average density in g/cm^3 from a gmx energy .xvg.

    GROMACS reports density in kg/m^3, so the raw column is ~1000x the value
    chemists expect; convert here rather than leaving a mislabelled number.
    """
    if not xvg.exists():
        return None
    vals = []
    for ln in xvg.read_text().splitlines():
        if ln.startswith(("#", "@")):
            continue
        f = ln.split()
        if len(f) >= 2:
            try:
                vals.append(float(f[1]))
            except ValueError:
                pass
    if not vals:
        return None
    # drop the first 20% as residual equilibration
    keep = vals[len(vals) // 5:] or vals
    return sum(keep) / len(keep) / 1000.0        # kg/m^3 -> g/cm^3


def perf_ns_day(log: Path) -> float | None:
    if not log.exists():
        return None
    m = re.findall(r"Performance:\s+([\d.]+)", log.read_text())
    return float(m[-1]) if m else None


def stage_reached(d: Path) -> str:
    last = "none"
    for s in STAGES:
        if (d / f"{s}.gro").exists():
            last = s
    return last


def main() -> int:
    ap = argparse.ArgumentParser()
    ap.add_argument("--out", type=Path, default=INDEX / "results.csv")
    args = ap.parse_args()

    row_of = dataset_rows()

    rows = []
    for d in sorted(RUNS.glob("el*")):
        if not (d / "summary.json").exists():
            continue
        s = json.loads((d / "summary.json").read_text())
        n_salt = s["counts"]["cat"]
        m_act, box_eq = achieved_molarity(d, n_salt)
        req = s["spec"]["salt_M"]
        row = {
            "key": s["name"],
            "row_index": row_of.get(s["name"]),
            "csv_line": (row_of[s["name"]] + 2) if s["name"] in row_of else None,
            "cation": s.get("cation"),
            "anion": s.get("anion"),
            "salt_smi": s["spec"]["salt_smi"],
            "solvent_smi": s["spec"]["solvents"][0]["smi"],
            "molarity": req,
            "molarity_actual": m_act,
            "molarity_err_pct": (100.0 * (m_act - req) / req
                                 if m_act is not None and req else None),
            "n_solvent": s["counts"]["sol1"],
            "n_salt": n_salt,
            "solv_per_ion_pair": s["counts"]["sol1"] / n_salt if n_salt else None,
            "n_atoms": s.get("n_atoms"),
            "box_build_nm": s.get("box_build_nm"),
            "box_equil_nm": box_eq,
            "net_charge": s.get("net_charge"),
            "q_scale": s.get("q_scale"),
            "stage": stage_reached(d),
            "has_xtc": (d / "prod.xtc").exists(),
            "density_g_cm3": mean_density(d / "density.xvg"),
            "prod_ns_day": perf_ns_day(d / "prod.log"),
        }
        sj = d / "solvation.json"
        if sj.exists():
            try:
                v = json.loads(sj.read_text())
                row.update({
                    "cn_solv_O": v.get("cn_solv_O"),
                    "cn_anion": v.get("cn_anion"),
                    "cn_total": v.get("cn_total"),
                    "ssip": v.get("ssip"),
                    "cip": v.get("cip"),
                    "agg": v.get("agg"),
                    "rdf_MO_peak_A": v.get("rdf_LiO_peak_A"),
                    "shell_cutoff_A": v.get("shell_cutoff_A"),
                    "n_frames_used": v.get("n_frames_used"),
                })
            except json.JSONDecodeError:
                pass

        cj = d / "clusters.json"
        if cj.exists():
            try:
                v = json.loads(cj.read_text())
                row.update({
                    "clust_cutoff_A": v.get("cutoff_A"),
                    "rdf_cat_an_peak_A": v.get("rdf_cation_anion_peak_A"),
                    # two conventions, differing by whether lone ions count:
                    # n_components is OVITO's, n_clusters is aggregates only
                    "n_components": v.get("n_components"),
                    "n_components_sd": v.get("n_components_sd"),
                    "n_clusters": v.get("n_clusters"),
                    "n_clusters_sd": v.get("n_clusters_sd"),
                    "free_cation_frac": v.get("free_cation_frac"),
                    "mean_cluster_size": v.get("mean_cluster_size"),
                    "max_cluster_size": v.get("max_cluster_size"),
                    "largest_cluster_frac": v.get("largest_cluster_frac"),
                    "percolating_frac": v.get("percolating_frac"),
                    "mean_cluster_charge": v.get("mean_cluster_charge"),
                })
            except json.JSONDecodeError:
                pass
        rows.append(row)

    if not rows:
        print("no runs found")
        return 0

    df = pd.DataFrame(rows)
    args.out.parent.mkdir(parents=True, exist_ok=True)
    df.to_csv(args.out, index=False)

    print(f"{len(df)} runs")
    print("\nstage reached:")
    print(df.stage.value_counts().to_string())
    done = df[df.get("cn_solv_O").notna()] if "cn_solv_O" in df else df.iloc[:0]
    print(f"\nwith solvation analysis: {len(done)}")
    if len(done):
        cols = ["key", "cation", "anion", "molarity", "molarity_actual",
                "solv_per_ion_pair", "density_g_cm3", "cn_solv_O", "cn_anion",
                "ssip", "cip", "agg", "clust_cutoff_A", "n_components",
                "n_clusters", "largest_cluster_frac"]
        print(done[[c for c in cols if c in done]].to_string(index=False))

    if "percolating_frac" in df and df.percolating_frac.notna().any():
        perc = df[df.percolating_frac > 0.5]
        print(f"\nclustering: {df.percolating_frac.notna().sum()} run(s) analysed, "
              f"{len(perc)} percolated (one network spanning >50% of the ions; "
              f"n_clusters is not meaningful for those)")

    drift = df["molarity_err_pct"].abs() if "molarity_err_pct" in df else None
    if drift is not None and drift.notna().any():
        bad = df[drift > 10]
        print(f"\nmolarity drift (requested -> equilibrated): "
              f"median {drift.median():.1f}%, max {drift.max():.1f}%")
        if len(bad):
            print(f"  {len(bad)} run(s) off by >10% -- the assumed component "
                  f"densities in compute_counts() are the usual cause")
    print(f"\n-> {args.out}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
