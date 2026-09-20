#!/usr/bin/env python
"""Fetch OPLS-AA/CM1A solvent parameters from LigParGen into ff_solvents/.

This is the answer to "where do the .itp/.gro files come from?" -- the previous
campaign pulled them, one SMILES at a time, from the Yale LigParGen web service.
None of the 490 solvents in this dataset overlap that old library, so all of them
have to be fetched fresh.

MUST RUN ON A LOGIN NODE: compute nodes have no outbound internet.

Files written per solvent, keyed by a hash of the canonical SMILES:

    ff_solvents/S<12hex>.itp  .gro  .pdb  .lmp  .smi

Idempotent -- re-running skips anything already complete, so it is safe to
interrupt and resume. Failures are recorded in index/ff_failures.csv along with
the server's error message.

Usage:
    python fetch_solvent_ff.py                 # every solvent in the dataset
    python fetch_solvent_ff.py --limit 10      # first N (pilot)
    python fetch_solvent_ff.py --smiles-file f # explicit list, one per line
    python fetch_solvent_ff.py --workers 3     # modest concurrency (default 2)
"""
from __future__ import annotations

import argparse
import csv
import sys
import threading
import time
from concurrent.futures import ThreadPoolExecutor, as_completed
from pathlib import Path

import pandas as pd

sys.path.insert(0, str(Path(__file__).resolve().parent))

import common as C  # noqa: E402
from common import (DATASET, FF_SOLV, INDEX, canon, check_solvent_ff,  # noqa: E402
                    solvent_key)
from ligpargen import LigParGenError, fetch, have_all  # noqa: E402

_print_lock = threading.Lock()


def log(msg: str):
    with _print_lock:
        print(msg, flush=True)


def targets_from_dataset(limit: int | None) -> list[str]:
    df = pd.read_csv(DATASET)
    smis: list[str] = []
    seen = set()
    for s in df["solvent"]:
        c = canon(s)
        if c and c not in seen:
            seen.add(c)
            smis.append(c)
    return smis[:limit] if limit else smis


def _rewrite_failures(fail_csv: Path, drop: set[str]):
    with fail_csv.open() as fh:
        rows = [r for r in csv.DictReader(fh) if r["smiles"] not in drop]
    with fail_csv.open("w", newline="") as fh:
        wr = csv.writer(fh)
        wr.writerow(["smiles", "key", "error"])
        for r in rows:
            wr.writerow([r["smiles"], r["key"], r["error"]])


def _discard(prefix: Path):
    for e in ("lmp", "pdb", "itp", "gro", "smi"):
        prefix.with_suffix(f".{e}").unlink(missing_ok=True)


def fetch_one(smi: str, pause: float, tries: int = 2) -> tuple[str, str, str]:
    """Returns (smiles, status, message). status in {ok, skip, fail}.

    Every download is verified against the requested SMILES before being kept.
    The LigParGen server has been observed returning a *different* molecule's
    files, so an unverified cache would quietly poison the whole campaign.
    """
    key = solvent_key(smi)
    prefix = FF_SOLV / key
    if have_all(prefix):
        return smi, "skip", ""

    last = ""
    for attempt in range(1, tries + 1):
        try:
            fetch(smi, prefix, net_charge=0, pause=pause, verbose=False)
        except LigParGenError as e:
            last = str(e)[:200]
            _discard(prefix)
            continue
        except Exception as e:                                # noqa: BLE001
            last = f"{type(e).__name__}: {e}"[:200]
            _discard(prefix)
            continue

        ok, why = check_solvent_ff(prefix.with_suffix(".itp").read_text(), smi)
        if ok:
            prefix.with_suffix(".smi").write_text(smi + "\n")
            return smi, "ok", ""
        last = f"{why} (attempt {attempt})"
        _discard(prefix)
        time.sleep(pause)
    return smi, "fail", last


def main() -> int:
    ap = argparse.ArgumentParser()
    ap.add_argument("--limit", type=int)
    ap.add_argument("--smiles-file", type=Path)
    # Concurrency corrupts results: the LigParGen backend runs BOSS in a shared
    # server-side working directory, so simultaneous submissions clobber one
    # another and come back as "error detected in your input data" even for
    # molecules that succeed perfectly well on their own. Measured: 7/10 spurious
    # failures at 2 workers, 0/10 at 1. Keep this at 1 unless that changes.
    ap.add_argument("--workers", type=int, default=1,
                    help="leave at 1; the server mis-handles concurrent requests")
    ap.add_argument("--pause", type=float, default=2.0)
    ap.add_argument("--retry-failed", action="store_true",
                    help="also retry SMILES recorded in ff_failures.csv")
    args = ap.parse_args()

    if args.workers > 1:
        log(f"WARNING: --workers {args.workers} -- LigParGen returns spurious "
            f"failures under concurrency; 1 is strongly recommended")

    FF_SOLV.mkdir(parents=True, exist_ok=True)
    INDEX.mkdir(parents=True, exist_ok=True)
    fail_csv = INDEX / "ff_failures.csv"

    if args.smiles_file:
        smis = [canon(l.strip()) for l in args.smiles_file.read_text().splitlines()
                if l.strip()]
        smis = [s for s in smis if s]
    else:
        smis = targets_from_dataset(args.limit)

    known_bad: set[str] = set()
    if fail_csv.exists() and not args.retry_failed:
        with fail_csv.open() as fh:
            known_bad = {r["smiles"] for r in csv.DictReader(fh)}
        # A SMILES that failed once (often only because of a transient server
        # error) may since have been fetched successfully. Drop those from the
        # blocklist so it reflects genuine failures rather than old history.
        healed = {s for s in known_bad if have_all(FF_SOLV / solvent_key(s))}
        if healed:
            known_bad -= healed
            _rewrite_failures(fail_csv, drop=healed)
            log(f"cleared {len(healed)} stale failure records (now cached)")
        if known_bad:
            log(f"skipping {len(known_bad)} previously-failed SMILES "
                f"(use --retry-failed to retry)")

    todo = [s for s in smis
            if not have_all(FF_SOLV / solvent_key(s)) and s not in known_bad]
    done_already = len(smis) - len(todo) - len(known_bad & set(smis))
    log(f"{len(smis)} unique solvents | {done_already} cached | {len(todo)} to fetch")
    if not todo:
        return 0

    t0 = time.time()
    results: list[tuple[str, str, str]] = []
    with ThreadPoolExecutor(max_workers=args.workers) as ex:
        futs = {ex.submit(fetch_one, s, args.pause): s for s in todo}
        for i, fut in enumerate(as_completed(futs), start=1):
            smi, status, msg = fut.result()
            results.append((smi, status, msg))
            rate = i / max(time.time() - t0, 1e-9)
            eta = (len(todo) - i) / rate if rate > 0 else 0
            mark = {"ok": "ok  ", "skip": "skip", "fail": "FAIL"}[status]
            log(f"[{i:4d}/{len(todo)}] {mark} {solvent_key(smi)} {smi[:52]:<52s}"
                f" eta {eta/60:5.1f}m" + (f"  {msg[:80]}" if msg else ""))

    failed = [(s, m) for s, st, m in results if st == "fail"]
    n_ok = sum(1 for _, st, _ in results if st == "ok")
    log(f"\nfetched {n_ok}, failed {len(failed)}, in {(time.time()-t0)/60:.1f} min")

    if failed:
        existing: dict[str, str] = {}
        if fail_csv.exists():
            with fail_csv.open() as fh:
                existing = {r["smiles"]: r["error"] for r in csv.DictReader(fh)}
        existing.update(dict(failed))
        with fail_csv.open("w", newline="") as fh:
            wr = csv.writer(fh)
            wr.writerow(["smiles", "key", "error"])
            for s, m in sorted(existing.items()):
                wr.writerow([s, solvent_key(s), m])
        log(f"failures -> {fail_csv}")

    # refresh the manifest of what the library now holds
    man = INDEX / "solvent_manifest.csv"
    with man.open("w", newline="") as fh:
        wr = csv.writer(fh)
        wr.writerow(["key", "smiles"])
        for p in sorted(FF_SOLV.glob("*.smi")):
            wr.writerow([p.stem, p.read_text().strip()])
    log(f"manifest -> {man}")
    return 1 if failed else 0


if __name__ == "__main__":
    raise SystemExit(main())
