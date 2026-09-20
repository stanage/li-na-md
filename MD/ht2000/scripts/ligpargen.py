"""LigParGen web-service client (OPLS-AA / 1.14*CM1A-LBCC).

Derived from ce_solvation_md's ligpargen_fetch.py, with two fixes needed here:

  1. **Net charge is honoured.** The original hardcoded the server's `dropcharge`
     field to 0, so every anion request failed with "error detected in your input
     data". We pass the real formal charge, which is what makes the 18 polyatomic
     anions in this dataset parameterisable.
  2. The charge/opt fallback ladder is retried per net-charge value, and the
     server's error text is extracted from the returned HTML for diagnostics.

Run from a login node: compute nodes have no outbound internet.
"""
from __future__ import annotations

import re
import time
from pathlib import Path

import requests
import urllib3

urllib3.disable_warnings(urllib3.exceptions.InsecureRequestWarning)

BASE = "https://zarbi.chem.yale.edu"
FORM_URL = f"{BASE}/cgi-bin/results_lpg.py"
DL_URL = f"{BASE}/cgi-bin/download_lpg.py"

EXTS = ("lmp", "pdb", "itp", "gro")

#: (chargetype, opt_iters) ladder, tried in order until the server yields files.
ATTEMPTS = (("cm1abcc", 0), ("cm1a", 1), ("cm1a", 2))


class LigParGenError(RuntimeError):
    pass


def _submit(smiles: str, opt_iters: int, chargetype: str, net_charge: int,
            timeout: int = 600) -> str:
    data = {
        "smiData": smiles,
        "checkopt": f" {opt_iters} ",
        "chargetype": chargetype,
        # The form sends this space-padded; the CGI is picky about the format.
        "dropcharge": f" {net_charge:+d} ".replace("+0", "0"),
    }
    files = {"molpdbfile": ("", b"", "application/octet-stream")}
    r = requests.post(FORM_URL, data=data, files=files, verify=False, timeout=timeout)
    r.raise_for_status()
    return r.text


def _extract_fileouts(html: str) -> dict[str, str]:
    out: dict[str, str] = {}
    for m in re.finditer(r'name="fileout"\s+value="(/tmp/[^"]+\.([A-Za-z]+))"', html):
        out.setdefault(m.group(2).lower(), m.group(1))
    return out


def server_error(html: str) -> str:
    """Pull the human-readable error out of a LigParGen failure page."""
    txt = re.sub(r"<[^>]*>", " ", html)
    txt = re.sub(r"\s+", " ", txt)
    m = re.search(r"(Sorry, an error[^.]*\.)", txt)
    if m:
        m2 = re.search(r"(Problem found in [^.]*\.)", txt)
        return m.group(1) + (" " + m2.group(1) if m2 else "")
    return txt[:300]


def _download(fileout: str, dst: Path, timeout: int = 120):
    r = requests.post(DL_URL, data={"fileout": fileout}, verify=False, timeout=timeout)
    r.raise_for_status()
    dst.write_bytes(r.content)


def have_all(prefix: Path) -> bool:
    return all(prefix.with_suffix(f".{e}").exists()
               and prefix.with_suffix(f".{e}").stat().st_size > 100 for e in EXTS)


def fetch(smiles: str, out_prefix: Path, net_charge: int = 0,
          pause: float = 2.0, verbose: bool = True) -> None:
    """Fetch OPLS-AA params for `smiles` into <out_prefix>.{lmp,pdb,itp,gro}.

    Idempotent: returns immediately if all four files already exist.
    Raises LigParGenError with the server's message if every attempt fails.
    """
    out_prefix = Path(out_prefix)
    out_prefix.parent.mkdir(parents=True, exist_ok=True)
    if have_all(out_prefix):
        if verbose:
            print(f"[skip] {out_prefix.name}")
        return

    last_html = ""
    for chargetype, opt in ATTEMPTS:
        if verbose:
            print(f"[submit] {smiles}  q={net_charge:+d}  {chargetype} opt={opt}")
        try:
            last_html = _submit(smiles, opt, chargetype, net_charge)
        except requests.RequestException as e:
            last_html = ""
            if verbose:
                print(f"  network error: {e}")
            time.sleep(pause * 3)
            continue
        files = _extract_fileouts(last_html)
        if all(k in files for k in EXTS):
            for ext in EXTS:
                _download(files[ext], out_prefix.with_suffix(f".{ext}"))
            if verbose:
                print(f"  saved {out_prefix.name}.{{{','.join(EXTS)}}}")
            return
        if verbose:
            print(f"  attempt failed (got {sorted(files)})")
        time.sleep(pause)

    out_prefix.with_suffix(".html").write_text(last_html)
    raise LigParGenError(server_error(last_html) if last_html else "no response")
