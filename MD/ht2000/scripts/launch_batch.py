#!/usr/bin/env python
"""Build and submit a batch of formulations from the dataset.

Build and submission are separate concerns and are separately resumable:

    --build-only     construct run dirs (packmol is the slow part, ~1-3 min each,
                     so this parallelises across CPU workers)
    --submit-only    sbatch whatever is already built and unfinished
    (neither)        do both

Idempotent throughout: a run with solvation.json is considered finished and is
skipped, so the campaign survives preemption, cancelled chunks and partial
fetches. Nothing is submitted twice unless --force.

Usage:
    python launch_batch.py --from-csv ../index/pilot_rows.csv --build-only -j 8
    python launch_batch.py --from-csv ../index/pilot_rows.csv --submit-only
    python launch_batch.py --limit 100 --offset 0 -j 16      # slice of the 2000
"""
from __future__ import annotations

import argparse
import subprocess
import sys
from concurrent.futures import ProcessPoolExecutor, as_completed
from pathlib import Path

import pandas as pd

sys.path.insert(0, str(Path(__file__).resolve().parent))

import common as C  # noqa: E402
from common import (DATASET, FF_SOLV, INDEX, LOGS, RUNS, anion_name, canon,  # noqa: E402
                    cation_name, fix_salt, run_key, solvent_key, split_salt)
from ligpargen import have_all  # noqa: E402

import setup_run  # noqa: E402


def load_rows(args) -> pd.DataFrame:
    df = pd.read_csv(args.from_csv) if args.from_csv else pd.read_csv(DATASET)
    if "concentration" not in df.columns:
        raise SystemExit("input csv needs a 'concentration' column")
    df = df.iloc[args.offset:]
    if args.limit:
        df = df.iloc[:args.limit]
    df = df.copy()
    df["salt"] = df["salt"].map(fix_salt)
    df["solvent"] = df["solvent"].map(canon)
    df["key"] = [run_key(r.salt, r.solvent, r.concentration)
                 for r in df.itertuples()]
    return df


def status_of(key: str) -> str:
    d = RUNS / key
    if (d / "solvation.json").exists():
        return "done"
    if (d / "prod.gro").exists():
        return "prod-done"
    if (d / "mixture.gro").exists():
        return "built"
    return "new"


def _build_one(payload):
    salt, solvent, molarity, kw = payload
    try:
        d, s = setup_run.build(salt, solvent, molarity, **kw)
        return s["name"], "ok", (f"nsol={s['counts']['sol1']} atoms={s['n_atoms']} "
                                 f"q={s['net_charge']:+.3f}")
    except Exception as e:                                    # noqa: BLE001
        return run_key(salt, solvent, molarity), "fail", f"{type(e).__name__}: {e}"


def submit_build_array(df: pd.DataFrame, args) -> int:
    """Farm run-dir construction out to CPU nodes.

    Packing 2000 systems is ~160 CPU-hours, which has no business running on a
    login node and should not occupy a GPU. This splits the work into chunks and
    submits them to the caslake partition; each task builds its slice with the
    same setup_run.build() the pilot used.
    """
    scripts = Path(__file__).resolve().parent
    todo = df[df.status == "new"]
    if not len(todo):
        print("nothing to build")
        return 0

    rows_csv = INDEX / f"build_{args.tag}_rows.csv"
    todo[["salt", "solvent", "concentration"]].to_csv(rows_csv, index=False)
    n_chunks = max(1, -(-len(todo) // args.chunk))
    script = INDEX / f"build_{args.tag}.sh"
    script.write_text(f"""#!/bin/bash
#SBATCH --job-name=ht2000-build-{args.tag}
#SBATCH --account={C.SLURM_ACCOUNT}
#SBATCH --partition=caslake
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task={args.jobs}
#SBATCH --time=12:00:00
#SBATCH --array=0-{n_chunks - 1}%{args.build_throttle}
#SBATCH --output={LOGS}/build_{args.tag}_%A_%a.out
#SBATCH --error={LOGS}/build_{args.tag}_%A_%a.err

OFFSET=$(( SLURM_ARRAY_TASK_ID * {args.chunk} ))
{C.PYTHON} {scripts}/launch_batch.py --from-csv {rows_csv} \\
    --offset $OFFSET --limit {args.chunk} \\
    --build-only -j {args.jobs} \\
    --q-scale {args.q_scale} --n-salt {args.n_salt} --prod-ns {args.prod_ns}
""")
    script.chmod(0o755)
    print(f"  build array: {len(todo)} formulations in {n_chunks} chunks of "
          f"{args.chunk}, {args.jobs} workers each, {args.build_throttle} concurrent")
    print(f"  -> {script}")
    if args.dry_run:
        print("  (dry run, not submitted)")
        return 0
    p = subprocess.run(["sbatch", str(script)], capture_output=True, text=True)
    if p.returncode != 0:
        print(f"  sbatch FAILED: {p.stderr.strip()[:300]}")
        return 1
    print(f"  submitted as job array {p.stdout.strip().split()[-1]}")
    return 0


def submit_array(ready: pd.DataFrame, args) -> int:
    """Submit the whole batch as one throttled SLURM job array.

    Preferred over N individual jobs for the full campaign: the gpu QOS only
    lets 12 of this user's jobs run at once, so 2000 separate sbatch calls would
    just park 1988 pending jobs in a shared queue to no benefit.
    """
    scripts = Path(__file__).resolve().parent
    keys = list(ready.key)
    if not keys:
        print("nothing to submit")
        return 0

    keyfile = INDEX / f"array_{args.tag}_keys.txt"
    keyfile.write_text("\n".join(keys) + "\n")
    script = INDEX / f"array_{args.tag}.sh"
    script.write_text(f"""#!/bin/bash
#SBATCH --job-name=ht2000-{args.tag}
#SBATCH --account={C.SLURM_ACCOUNT}
#SBATCH --partition={C.SLURM_PARTITION}
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=8
#SBATCH --gres=gpu:1
#SBATCH --time={args.hours:02d}:00:00
#SBATCH --array=1-{len(keys)}%{args.throttle}
#SBATCH --output={LOGS}/{args.tag}_%A_%a.out
#SBATCH --error={LOGS}/{args.tag}_%A_%a.err

KEY=$(sed -n "${{SLURM_ARRAY_TASK_ID}}p" {keyfile})
if [[ -z "$KEY" ]]; then echo "no key for task $SLURM_ARRAY_TASK_ID"; exit 1; fi
if [[ -f "{RUNS}/$KEY/solvation.json" ]]; then
    echo "$KEY already finished, skipping"; exit 0
fi
bash {scripts}/run_md.sh {RUNS}/$KEY
""")
    script.chmod(0o755)
    print(f"  array of {len(keys)} tasks, throttled to {args.throttle} concurrent")
    print(f"  -> {script}")
    if args.dry_run:
        print("  (dry run, not submitted)")
        return 0
    p = subprocess.run(["sbatch", str(script)], capture_output=True, text=True)
    if p.returncode != 0:
        print(f"  sbatch FAILED: {p.stderr.strip()[:300]}")
        return 1
    jid = p.stdout.strip().split()[-1]
    print(f"  submitted as job array {jid}")
    pd.DataFrame({"key": keys, "jobid": f"{jid}"}).to_csv(
        INDEX / f"array_{args.tag}_submitted.csv", index=False)
    return 0


def main() -> int:
    ap = argparse.ArgumentParser()
    ap.add_argument("--from-csv", type=Path)
    ap.add_argument("--limit", type=int)
    ap.add_argument("--offset", type=int, default=0)
    ap.add_argument("-j", "--jobs", type=int, default=8,
                    help="parallel build workers (packmol is single-threaded)")
    ap.add_argument("--build-only", action="store_true")
    ap.add_argument("--submit-only", action="store_true")
    ap.add_argument("--q-scale", type=float, default=C.Q_SCALE_DEFAULT)
    ap.add_argument("--n-salt", type=int, default=C.N_SALT_DEFAULT)
    ap.add_argument("--prod-ns", type=float, default=C.PROD_NS)
    ap.add_argument("--hours", type=int, default=8)
    ap.add_argument("--max-submit", type=int,
                    help="cap how many jobs go to the queue in this call")
    ap.add_argument("--array", action="store_true",
                    help="submit as ONE throttled SLURM job array instead of N "
                         "individual jobs (use this for the full campaign)")
    ap.add_argument("--throttle", type=int, default=12,
                    help="max concurrently running array tasks; the gpu QOS "
                         "caps this user at 12 jobs / 16 GPUs")
    ap.add_argument("--tag", default="batch", help="name for the array files")
    ap.add_argument("--build-array", action="store_true",
                    help="submit run-dir construction to CPU nodes as a job "
                         "array instead of building here (use for the full 2000)")
    ap.add_argument("--chunk", type=int, default=100,
                    help="formulations per build-array task")
    ap.add_argument("--build-throttle", type=int, default=5,
                    help="max concurrent build-array tasks")
    ap.add_argument("--force", action="store_true")
    ap.add_argument("--dry-run", action="store_true")
    args = ap.parse_args()

    df = load_rows(args)
    print(f"{len(df)} formulations requested")

    # ---- skip anything whose solvent force field is missing
    missing = [s for s in df.solvent.unique() if not have_all(FF_SOLV / solvent_key(s))]
    if missing:
        print(f"WARNING: {len(missing)} solvents have no force field yet -- "
              f"run fetch_solvent_ff.py first. Affected rows are skipped.")
        (INDEX / "missing_ff.txt").write_text("\n".join(missing) + "\n")
        df = df[~df.solvent.isin(missing)]

    df["status"] = df.key.map(status_of)
    print("  " + "  ".join(f"{k}={v}" for k, v in
                           df.status.value_counts().to_dict().items()))

    # ---------------------------------------------------------------- build
    if args.build_array:
        return submit_build_array(df, args)

    if not args.submit_only:
        todo = df[df.status == "new"] if not args.force else df
        if len(todo):
            kw = dict(q_scale=args.q_scale, n_salt=args.n_salt,
                      prod_ns=args.prod_ns, hours=args.hours, force=args.force)
            payloads = [(r.salt, r.solvent, float(r.concentration), kw)
                        for r in todo.itertuples()]
            print(f"\nbuilding {len(payloads)} run dirs with {args.jobs} workers...")
            if args.dry_run:
                print("  (dry run, nothing built)")
            else:
                nfail = 0
                with ProcessPoolExecutor(max_workers=args.jobs) as ex:
                    futs = [ex.submit(_build_one, p) for p in payloads]
                    for i, f in enumerate(as_completed(futs), 1):
                        key, st, msg = f.result()
                        if st == "fail":
                            nfail += 1
                        print(f"  [{i}/{len(payloads)}] {st:4s} {key} {msg}")
                print(f"built {len(payloads)-nfail}, failed {nfail}")
            df["status"] = df.key.map(status_of)

    # ---------------------------------------------------------------- submit
    if args.build_only:
        return 0

    ready = df[df.status == "built"]
    if args.max_submit:
        ready = ready.iloc[:args.max_submit]
    print(f"\nsubmitting {len(ready)} jobs to {C.SLURM_PARTITION} "
          f"(account {C.SLURM_ACCOUNT})")

    if args.array:
        return submit_array(ready, args)

    submitted = []
    for r in ready.itertuples():
        sub = RUNS / r.key / "submit.sh"
        if not sub.exists():
            print(f"  skip {r.key}: no submit.sh")
            continue
        if args.dry_run:
            print(f"  would sbatch {sub}")
            continue
        p = subprocess.run(["sbatch", str(sub)], cwd=RUNS / r.key,
                           capture_output=True, text=True)
        if p.returncode != 0:
            print(f"  FAIL {r.key}: {p.stderr.strip()[:160]}")
            continue
        jid = p.stdout.strip().split()[-1]
        submitted.append({"key": r.key, "jobid": jid})
        print(f"  {r.key} -> job {jid}")

    if submitted and not args.dry_run:
        man = INDEX / "submitted.csv"
        old = pd.read_csv(man) if man.exists() else pd.DataFrame(columns=["key", "jobid"])
        pd.concat([old, pd.DataFrame(submitted)]).drop_duplicates(
            subset=["key", "jobid"]).to_csv(man, index=False)
        print(f"\n-> {man}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
