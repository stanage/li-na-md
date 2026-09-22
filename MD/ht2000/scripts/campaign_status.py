#!/usr/bin/env python
"""Progress, failures and ETA for the campaign. Safe to run at any time.

Reads the per-lane progress CSVs plus the run directories themselves, so it
reports what is actually on disk rather than trusting the logs alone.

Usage:  python campaign_status.py [--failures] [--lanes]
"""
from __future__ import annotations

import argparse
import glob
import sys
from datetime import timedelta
from pathlib import Path

import pandas as pd

sys.path.insert(0, str(Path(__file__).resolve().parent))

from common import HT, INDEX, RUNS  # noqa: E402

STAGES = ["em", "nvt_heat", "npt", "prod"]


def disk_state() -> pd.DataFrame:
    man = pd.read_csv(INDEX / "run_manifest.csv")
    rows = []
    for key in man.key:
        d = RUNS / key
        if not d.is_dir():
            rows.append("not_built"); continue
        if (d / "clusters.json").exists():
            rows.append("complete"); continue
        st = "built"
        for s in STAGES:
            if (d / f"{s}.gro").exists():
                st = s
        rows.append(st)
    man["state"] = rows
    return man


def main() -> int:
    ap = argparse.ArgumentParser()
    ap.add_argument("--failures", action="store_true", help="list failed runs")
    ap.add_argument("--lanes", action="store_true", help="per-lane breakdown")
    ap.add_argument("--remaining", action="store_true",
                    help="print only the number of runs left, for scripts")
    ap.add_argument("--lane", type=int, help="with --remaining: one lane only")
    ap.add_argument("--nlanes", type=int, default=3)
    args = ap.parse_args()

    man = disk_state()
    n = len(man)
    done = int((man.state == "complete").sum())
    runnable = int(man.solvent_ff.sum())

    # Machine-readable, for the self-chaining step in submit_campaign.sbatch:
    # a lane resubmits itself only while its own slice has work left.
    if args.remaining:
        m = man[man.solvent_ff.astype(bool)]
        if args.lane is not None:
            m = m[m.row_index % args.nlanes == args.lane]
        print(int((m.state != "complete").sum()))
        return 0

    print(f"campaign: {done}/{runnable} complete "
          f"({100*done/max(runnable,1):.1f}%)   "
          f"[{n - runnable} rows have no solvent force field]")
    print("\nstate on disk:")
    print(man.state.value_counts().to_string())

    prog = sorted(glob.glob(str(HT / "logs" / "campaign" / "lane_*.csv")))
    if not prog:
        print("\nno lane progress yet -- campaign has not run")
        return 0

    df = pd.concat([pd.read_csv(p) for p in prog], ignore_index=True)
    df = df[df.status.isin(["ok", "FAILED"])]
    if df.empty:
        print("\nlanes started but nothing has finished yet")
        return 0

    ok = df[df.status == "ok"]
    bad = df[df.status == "FAILED"]
    med = ok.seconds.median() if len(ok) else float("nan")
    print(f"\nfinished this campaign: {len(ok)} ok, {len(bad)} failed")
    if len(ok):
        print(f"median run time: {timedelta(seconds=int(med))}  "
              f"(p90 {timedelta(seconds=int(ok.seconds.quantile(0.9)))})")
        remaining = runnable - done
        lanes = max(len(prog), 1)
        eta = remaining * med / lanes
        print(f"remaining: {remaining} runs over {lanes} lanes "
              f"-> ~{timedelta(seconds=int(eta))} wall")

    if args.lanes:
        print("\nper lane:")
        for p in prog:
            d = pd.read_csv(p)
            d = d[d.status.isin(["ok", "FAILED"])]
            print(f"  {Path(p).stem}: {int((d.status=='ok').sum()):4d} ok, "
                  f"{int((d.status=='FAILED').sum()):3d} failed")

    if len(bad):
        print(f"\nfailures by last stage reached:")
        print(bad.stage.fillna("none").value_counts().to_string())
        if args.failures:
            print("\nfailed runs:")
            for _, r in bad.iterrows():
                print(f"  row {r.row_index:<5} {r.key}  stage={r.stage}  "
                      f"-> runs/{r.key}/pipeline.log")
        else:
            print("  (--failures to list them)")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
