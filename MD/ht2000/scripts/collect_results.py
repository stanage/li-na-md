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

from common import INDEX, RUNS  # noqa: E402

STAGES = ["em", "nvt_heat", "npt", "prod"]


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

    rows = []
    for d in sorted(RUNS.glob("el*")):
        if not (d / "summary.json").exists():
            continue
        s = json.loads((d / "summary.json").read_text())
        row = {
            "key": s["name"],
            "cation": s.get("cation"),
            "anion": s.get("anion"),
            "salt_smi": s["spec"]["salt_smi"],
            "solvent_smi": s["spec"]["solvents"][0]["smi"],
            "molarity": s["spec"]["salt_M"],
            "n_solvent": s["counts"]["sol1"],
            "n_salt": s["counts"]["cat"],
            "n_atoms": s.get("n_atoms"),
            "box_build_nm": s.get("box_build_nm"),
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
        cols = ["key", "cation", "anion", "molarity", "density_g_cm3",
                "cn_solv_O", "cn_anion", "ssip", "cip", "agg"]
        print(done[[c for c in cols if c in done]].to_string(index=False))
    print(f"\n-> {args.out}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
