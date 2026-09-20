#!/usr/bin/env python
"""Ion-cluster analysis for one finished run -> clusters.json.

Answers "how many cation-anion clusters are there, and how big?" without
leaving the cluster: it reads prod.tpr/prod.xtc in place, so there is nothing
to download and no OVITO round-trip.

Three steps, in order:

  1. Cutoff from the data. The cation-anion RDF is built over the anion's
     *donor* atoms only (see DONORS_BY_ANION), its first peak is located, and
     the following minimum becomes the contact cutoff. That minimum is the
     edge of the first coordination shell, which is the physically meaningful
     definition of "in contact" -- not a round number picked by eye.

  2. Contact graph per frame. Nodes are ions (cations + whole anion
     residues); an edge joins a cation to an anion when any donor atom of that
     anion is within the cutoff. Distances use the minimum image convention,
     so clusters wrapped across the periodic boundary stay whole -- the usual
     failure mode of eyeballing this in a viewer.

  3. Connected components. Every maximal connected set of ions is one cluster.

Reported per frame and averaged over the trajectory:

    n_clusters        aggregates, i.e. components holding >= 2 ions
    n_components      ALL components, lone ions included -- this is the number
                      OVITO's Cluster Analysis reports, and it equals
                      n_clusters + n_free_cation + n_free_anion
    n_free_cation     cations contacting no anion
    n_free_anion      anions contacting no cation
    mean/max_cluster_size
    largest_cluster_frac   ions in the biggest cluster / all ions
    percolating_frac       frames whose largest cluster holds > 50% of ions
    size_histogram    cluster size -> mean count per frame
    mean_cluster_charge    <n_cat - n_an> over aggregates

Note on interpretation: once largest_cluster_frac approaches 1 the system has
percolated into one network, and "number of clusters" stops being a useful
descriptor -- read size_histogram and largest_cluster_frac instead.

Usage:  python cluster_analysis.py --run-dir runs/el<hash>
        python cluster_analysis.py --run-dir ... --cutoff 3.5   # override
"""
from __future__ import annotations

import argparse
import json
import sys
from collections import Counter
from pathlib import Path

import numpy as np

sys.path.insert(0, str(Path(__file__).resolve().parent))

from analyze_solvation import (DEFAULT_DONORS, DONORS_BY_ANION,  # noqa: E402
                               RDF_NBINS, RDF_RMAX, _element_of,
                               first_peak_and_min)

#: fallback cation-anion contact radii (Angstrom) if the RDF minimum is
#: unresolvable -- wider than the cation-O shell because anion donors sit
#: slightly further out.
DEFAULT_CUTOFF = {"Li": 3.0, "Na": 3.5}


def load(run_dir: Path):
    import MDAnalysis as mda

    tpr, xtc = run_dir / "prod.tpr", run_dir / "prod.xtc"
    if not tpr.exists() or not xtc.exists():
        raise FileNotFoundError(f"need prod.tpr and prod.xtc in {run_dir}")

    summary = json.loads((run_dir / "summary.json").read_text())
    cation = summary.get("cation", "Li")
    anion = summary.get("anion", "")

    u = mda.Universe(str(tpr), str(xtc))
    cat = u.select_atoms(f"resname {cation.lower()[:3]}")
    if cat.n_atoms == 0:
        cat = u.select_atoms(f"name {cation.upper()} or name {cation}")
    ani = u.select_atoms(f"resname {anion.lower()[:3]}")
    if ani.n_atoms == 0:
        ani = u.select_atoms("resname an")

    donors = DONORS_BY_ANION.get(anion, DEFAULT_DONORS)
    keep = [i for i, a in enumerate(ani.names) if _element_of(a) in donors]
    an_donor = ani[keep]
    return u, cat, ani, an_donor, cation, anion, donors


def contact_cutoff(u, cat, an_donor, frames, cation: str) -> tuple[float, float, bool]:
    """(cutoff, rdf first peak, whether the minimum was resolved)."""
    from MDAnalysis.analysis.rdf import InterRDF

    if not cat.n_atoms or not an_donor.n_atoms:
        return DEFAULT_CUTOFF.get(cation, 3.2), float("nan"), False
    rdf = InterRDF(cat, an_donor, nbins=RDF_NBINS, range=(0.0, RDF_RMAX))
    rdf.run(start=frames[0], stop=frames[-1] + 1, step=frames[1] - frames[0]
            if len(frames) > 1 else 1)
    peak, rmin = first_peak_and_min(rdf.results.bins, rdf.results.rdf)
    if np.isfinite(rmin):
        return float(rmin), float(peak), True
    return DEFAULT_CUTOFF.get(cation, 3.2), float(peak), False


def frame_clusters(cat, an_donor, an_resindex, n_an, cutoff, box):
    """Component labels and per-component (n_cat, n_an) for one frame."""
    from MDAnalysis.lib.distances import capped_distance
    from scipy.sparse import coo_matrix
    from scipy.sparse.csgraph import connected_components

    n_cat = cat.n_atoms
    n = n_cat + n_an
    rows: np.ndarray
    if n_cat and an_donor.n_atoms:
        pairs = capped_distance(cat.positions, an_donor.positions, cutoff,
                                box=box, return_distances=False)
    else:
        pairs = np.empty((0, 2), dtype=int)

    if len(pairs):
        rows = pairs[:, 0]
        cols = n_cat + an_resindex[pairs[:, 1]]
    else:
        rows = np.empty(0, dtype=int)
        cols = np.empty(0, dtype=int)

    adj = coo_matrix((np.ones(len(rows)), (rows, cols)), shape=(n, n))
    n_comp, labels = connected_components(adj, directed=False)

    is_cat = np.zeros(n, dtype=bool)
    is_cat[:n_cat] = True
    ncat = np.bincount(labels, weights=is_cat, minlength=n_comp)
    nan_ = np.bincount(labels, weights=~is_cat, minlength=n_comp)
    return ncat, nan_


def compute(run_dir: Path, stride: int = 1, max_frames: int = 200,
            cutoff_override: float | None = None) -> dict:
    u, cat, ani, an_donor, cation, anion, donors = load(run_dir)

    n_tot = len(u.trajectory)
    step = max(stride, n_tot // max_frames or 1)
    frames = list(range(0, n_tot, step))

    if cutoff_override is not None:
        cutoff, peak, resolved = float(cutoff_override), float("nan"), True
    else:
        cutoff, peak, resolved = contact_cutoff(u, cat, an_donor, frames, cation)

    # donor atom -> 0-based anion residue, so a multi-donor anion is one node
    uniq = {r: i for i, r in enumerate(np.unique(ani.resids))}
    an_resindex = np.array([uniq[r] for r in an_donor.resids], dtype=int)
    n_an = len(uniq)
    n_cat = cat.n_atoms
    n_ion = n_cat + n_an

    per_frame = []
    hist = Counter()
    for fi in frames:
        u.trajectory[fi]
        ncat, nan_ = frame_clusters(cat, an_donor, an_resindex, n_an,
                                    cutoff, u.dimensions)
        size = ncat + nan_
        agg = size >= 2                      # a cluster proper, not a lone ion
        sizes = size[agg]
        for s in sizes:
            hist[int(s)] += 1
        per_frame.append({
            "n_clusters": int(agg.sum()),
            # every connected component, monomers included -- OVITO's Cluster
            # Analysis counts this way, so quote it when comparing.
            "n_components": int(len(size)),
            "n_free_cation": int(((size == 1) & (ncat == 1)).sum()),
            "n_free_anion": int(((size == 1) & (nan_ == 1)).sum()),
            "mean_cluster_size": float(sizes.mean()) if sizes.size else 0.0,
            "max_cluster_size": int(sizes.max()) if sizes.size else 0,
            "largest_cluster_frac": float(sizes.max() / n_ion) if sizes.size else 0.0,
            "mean_cluster_charge": float((ncat[agg] - nan_[agg]).mean())
            if sizes.size else 0.0,
        })

    nf = len(per_frame)

    def avg(k):
        return float(np.mean([p[k] for p in per_frame])) if nf else float("nan")

    def sd(k):
        return float(np.std([p[k] for p in per_frame])) if nf else float("nan")

    return {
        "cation": cation,
        "anion": anion,
        "n_cation": int(n_cat),
        "n_anion": int(n_an),
        "n_frames_used": nf,
        "cutoff_A": round(float(cutoff), 4),
        "cutoff_from_rdf": bool(cutoff_override is None and resolved),
        "rdf_cation_anion_peak_A": None if not np.isfinite(peak) else round(peak, 4),
        "anion_donor_elements": list(donors),
        "n_clusters": round(avg("n_clusters"), 3),
        "n_clusters_sd": round(sd("n_clusters"), 3),
        "n_components": round(avg("n_components"), 3),
        "n_components_sd": round(sd("n_components"), 3),
        "n_free_cation": round(avg("n_free_cation"), 3),
        "n_free_anion": round(avg("n_free_anion"), 3),
        "free_cation_frac": round(avg("n_free_cation") / n_cat, 4) if n_cat else None,
        "mean_cluster_size": round(avg("mean_cluster_size"), 3),
        "max_cluster_size": round(avg("max_cluster_size"), 3),
        "largest_cluster_frac": round(avg("largest_cluster_frac"), 4),
        "percolating_frac": round(
            float(np.mean([p["largest_cluster_frac"] > 0.5 for p in per_frame])), 4)
        if nf else None,
        "mean_cluster_charge": round(avg("mean_cluster_charge"), 4),
        "size_histogram": {str(k): round(v / nf, 4) for k, v in sorted(hist.items())}
        if nf else {},
    }


def main() -> int:
    ap = argparse.ArgumentParser()
    ap.add_argument("--run-dir", type=Path, required=True)
    ap.add_argument("--stride", type=int, default=1)
    ap.add_argument("--max-frames", type=int, default=200)
    ap.add_argument("--cutoff", type=float,
                    help="contact cutoff in Angstrom; default is the "
                         "cation-anion RDF first minimum")
    ap.add_argument("--out", type=Path)
    args = ap.parse_args()

    res = compute(args.run_dir, args.stride, args.max_frames, args.cutoff)
    out = args.out or (args.run_dir / "clusters.json")
    out.write_text(json.dumps(res, indent=2))
    print(json.dumps(res, indent=2))
    print(f"-> {out}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
