#!/usr/bin/env python
"""Compute the solvation structure of one finished run -> solvation.json.

Re-implements `elytefm.analyze_solvation` (never backed up) against the schema
of the archived ce_solvation_md solvation.json files:

    n_li              cations analysed
    n_frames_used     production frames averaged
    cn_solv_O         cation first-shell coordination to solvent oxygen
    cn_anion          cation first-shell coordination to anion donor atoms
    ssip / cip / agg  fraction of cations with 0 / 1 / >=2 anions in the shell
    rdf_LiO_peak_A    first peak of the cation-O RDF, in Angstrom
    anion_resnames    e.g. ["fsi"]
    solvent_resnames  e.g. ["SOL"]

Fields added here (the old schema had no Na runs and one shell convention):
    cation            "Li" or "Na"
    shell_cutoff_A    first-minimum cutoff actually used
    shell_cutoff_from_rdf  False when it fell back to a hardcoded radius
    cn_solv_N, cn_total

The first-shell cutoff is taken from the first minimum of the cation-O RDF
rather than being hardcoded, because Na+ has a visibly larger first shell than
Li+ and a single fixed radius would bias the two halves of this dataset against
each other. If the minimum cannot be located the element default is used.

Usage:  python analyze_solvation.py --run-dir runs/el<hash>
"""
from __future__ import annotations

import argparse
import json
import sys
from pathlib import Path

import numpy as np

sys.path.insert(0, str(Path(__file__).resolve().parent))

#: fallback first-shell radii (Angstrom) if the RDF minimum is not resolvable
DEFAULT_CUTOFF = {"Li": 2.8, "Na": 3.3}

#: Which atoms actually coordinate the cation, per anion. Counting *every* anion
#: atom inside the shell would badly inflate cn_anion -- e.g. FSI's two sulfurs
#: and two fluorines sit close behind the coordinating oxygens without donating
#: to the cation, and no anion's carbon skeleton coordinates at all.
DONORS_BY_ANION = {
    "FSI": ("O", "N"), "TFSI": ("O", "N"), "Beti": ("O", "N"),
    "BF4": ("F",), "PF6": ("F",),
    "ClO4": ("O",), "NO3": ("O",), "OTf": ("O",), "OTs": ("O",),
    "MS": ("O",), "OAc": ("O",), "PROP": ("O",), "BNZ": ("O",),
    "HCOO": ("O",), "TFA": ("O",),
    "DCA": ("N",), "TCM": ("N",), "SCN": ("N", "S"),
    "Cl": ("CL",), "Br": ("BR",),
}
#: used when an anion is not in the table above
DEFAULT_DONORS = ("O", "N", "F", "CL", "BR", "S")

RDF_RMAX = 8.0
RDF_NBINS = 400


def _element_of(name: str) -> str:
    """Crude element guess from a GROMACS atom name (no element field in .tpr)."""
    n = name.upper().lstrip("0123456789")
    for two in ("CL", "BR", "SI"):
        if n.startswith(two):
            return two
    return n[:1]


#: a first peak must rise this far above the ideal-gas baseline to count as
#: coordination rather than noise. An unstructured liquid oscillates about
#: g(r) = 1; a real first solvation shell is several times that.
MIN_PEAK_HEIGHT = 1.8


def first_peak_and_min(r: np.ndarray, g: np.ndarray,
                       rmin_search: float = 1.2) -> tuple[float, float]:
    """Locate the FIRST RDF maximum and the following minimum.

    Two things this deliberately does not do.

    It does not take the global maximum. argmax over the whole range returns
    whichever peak is tallest, and when a solvent-separated shell outgrows the
    contact shell -- which happens for weakly coordinating solvents -- that is
    the wrong peak, and the minimum search then runs off the end. So we walk
    outward and stop at the first local maximum that clears MIN_PEAK_HEIGHT.

    It does not report a peak it cannot justify. On a structureless or
    unminimised configuration every candidate stays near g(r) = 1; returning
    nan there lets the caller fall back explicitly and say so, rather than
    quoting a cutoff read off noise.
    """
    m = r >= rmin_search
    rr, gg = r[m], g[m]
    if gg.size < 5 or gg.max() <= 0:
        return float("nan"), float("nan")
    # smooth lightly so shot noise does not create spurious extrema
    k = np.ones(5) / 5.0
    gs = np.convolve(gg, k, mode="same")

    ipk = None
    for i in range(1, len(gs) - 1):
        if gs[i] >= gs[i - 1] and gs[i] >= gs[i + 1] and gs[i] >= MIN_PEAK_HEIGHT:
            ipk = i
            break
    if ipk is None:
        return float("nan"), float("nan")
    peak = float(rr[ipk])

    # first minimum after that peak
    tail = gs[ipk:]
    imin = None
    for i in range(1, len(tail) - 1):
        if tail[i] <= tail[i - 1] and tail[i] <= tail[i + 1] and tail[i] < 0.6 * tail[0]:
            imin = ipk + i
            break
    rmin = float(rr[imin]) if imin is not None else float("nan")
    return peak, rmin


def compute(run_dir: Path, stride: int = 1, max_frames: int = 400) -> dict:
    import MDAnalysis as mda
    from MDAnalysis.analysis.rdf import InterRDF

    tpr, xtc = run_dir / "prod.tpr", run_dir / "prod.xtc"
    if not tpr.exists() or not xtc.exists():
        raise FileNotFoundError(f"need prod.tpr and prod.xtc in {run_dir}")

    summary = json.loads((run_dir / "summary.json").read_text())
    cation = summary.get("cation", "Li")
    anion = summary.get("anion", "")
    an_res = anion.lower()[:3]

    u = mda.Universe(str(tpr), str(xtc))

    cat = u.select_atoms(f"resname {cation.lower()[:3]}")
    if cat.n_atoms == 0:                       # fall back to atom name
        cat = u.select_atoms(f"name {cation.upper()} or name {cation}")
    sol = u.select_atoms("resname SOL")
    ani = u.select_atoms(f"resname {an_res}")
    if ani.n_atoms == 0:
        ani = u.select_atoms("resname an")

    sol_O = sol.select_atoms("name O*")
    sol_N = sol.select_atoms("name N*")
    donors = DONORS_BY_ANION.get(anion, DEFAULT_DONORS)
    an_donor = ani[[i for i, a in enumerate(ani.names)
                    if _element_of(a) in donors]]

    # ---- frames
    n_tot = len(u.trajectory)
    step = max(stride, n_tot // max_frames or 1)
    frames = list(range(0, n_tot, step))

    # ---- RDF cation-O(solvent) -> peak and shell cutoff
    peak = rmin = float("nan")
    if cat.n_atoms and sol_O.n_atoms:
        rdf = InterRDF(cat, sol_O, nbins=RDF_NBINS, range=(0.0, RDF_RMAX))
        rdf.run(start=0, stop=n_tot, step=step)
        peak, rmin = first_peak_and_min(rdf.results.bins, rdf.results.rdf)
    # 921 rows of this dataset have no solvent oxygen at all, so the RDF above
    # is never built and rmin is nan. The fallback is then a hardcoded radius,
    # which must not be indistinguishable from a measured one downstream.
    cutoff_from_rdf = bool(np.isfinite(rmin))
    cutoff = rmin if cutoff_from_rdf else DEFAULT_CUTOFF.get(cation, 3.0)

    # ---- coordination, per frame
    from MDAnalysis.lib.distances import capped_distance

    def count_within(group, cut):
        """Per-cation neighbour count within `cut` for this frame."""
        if group.n_atoms == 0 or cat.n_atoms == 0:
            return np.zeros(cat.n_atoms)
        pairs = capped_distance(cat.positions, group.positions, cut,
                                box=u.dimensions, return_distances=False)
        out = np.zeros(cat.n_atoms)
        if len(pairs):
            np.add.at(out, pairs[:, 0], 1.0)
        return out

    acc = {"O": [], "N": [], "an": [], "an_res": []}
    for fi in frames:
        u.trajectory[fi]
        acc["O"].append(count_within(sol_O, cutoff).mean() if sol_O.n_atoms else 0.0)
        acc["N"].append(count_within(sol_N, cutoff).mean() if sol_N.n_atoms else 0.0)
        c_an = count_within(an_donor, cutoff)
        acc["an"].append(c_an.mean() if an_donor.n_atoms else 0.0)

        # distinct anion RESIDUES in the shell -> SSIP / CIP / AGG
        if an_donor.n_atoms and cat.n_atoms:
            pairs = capped_distance(cat.positions, an_donor.positions, cutoff,
                                    box=u.dimensions, return_distances=False)
            per_cat = [set() for _ in range(cat.n_atoms)]
            if len(pairs):
                resids = an_donor.resids
                for ci, ai in pairs:
                    per_cat[ci].add(resids[ai])
            acc["an_res"].append([len(s) for s in per_cat])
        else:
            acc["an_res"].append([0] * cat.n_atoms)

    nres = np.array(acc["an_res"], dtype=float)
    ssip = float((nres == 0).mean()) if nres.size else float("nan")
    cip = float((nres == 1).mean()) if nres.size else float("nan")
    agg = float((nres >= 2).mean()) if nres.size else float("nan")

    return {
        "n_li": int(cat.n_atoms),
        "n_frames_used": len(frames),
        "cn_solv_O": round(float(np.mean(acc["O"])), 4),
        "cn_anion": round(float(np.mean(acc["an"])), 4),
        "ssip": round(ssip, 4),
        "cip": round(cip, 4),
        "agg": round(agg, 4),
        "rdf_LiO_peak_A": None if not np.isfinite(peak) else round(peak, 4),
        "anion_resnames": [an_res] if ani.n_atoms else [],
        "solvent_resnames": ["SOL"] if sol.n_atoms else [],
        # --- additions beyond the original schema
        "cation": cation,
        "anion": anion,
        "shell_cutoff_A": round(float(cutoff), 4),
        "shell_cutoff_from_rdf": cutoff_from_rdf,
        "anion_donor_elements": list(donors),
        "n_anion_donor_atoms": int(an_donor.n_atoms),
        "cn_solv_N": round(float(np.mean(acc["N"])), 4),
        "cn_total": round(float(np.mean(acc["O"]) + np.mean(acc["N"])
                                + np.mean(acc["an"])), 4),
    }


def main() -> int:
    ap = argparse.ArgumentParser()
    ap.add_argument("--run-dir", type=Path, required=True)
    ap.add_argument("--stride", type=int, default=1)
    ap.add_argument("--max-frames", type=int, default=400)
    ap.add_argument("--out", type=Path)
    args = ap.parse_args()

    res = compute(args.run_dir, args.stride, args.max_frames)
    out = args.out or (args.run_dir / "solvation.json")
    out.write_text(json.dumps(res, indent=2))
    print(json.dumps(res, indent=2))
    print(f"-> {out}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
