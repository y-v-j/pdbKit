"""Download structures from the RCSB PDB and the AlphaFold Protein Structure Database.

* **RCSB PDB** - ``https://files.rcsb.org/download/``: the deposited entry
  (``.cif`` / ``.pdb``) or a biological assembly
  (``<ID>-assembly<N>.cif`` / ``<ID>.pdb<N>``).
* **AlphaFold DB** - the model URLs are looked up through the AFDB API
  (``/api/prediction/<UniProt>``), so the current model version is used
  without hard-coding it; direct file URLs are tried as a fallback.

Large entries (for example ribosomes) have no PDB-format file at RCSB; asking
for ``fmt='pdb'`` then fails with a clear message suggesting mmCIF.
"""

from __future__ import annotations

import json
import os
import re
import sys
import tempfile
import time
import urllib.error
import urllib.request
from pathlib import Path
from typing import List, Optional

from .structure import PdbKitError

__all__ = ["fetch", "fetch_pdb", "fetch_alphafold", "alphafold_entries",
           "classify_identifier", "default_cache_dir"]

RCSB_URL = "https://files.rcsb.org/download/{name}"
AFDB_API = "https://alphafold.ebi.ac.uk/api/prediction/{acc}"
AFDB_FILE = "https://alphafold.ebi.ac.uk/files/AF-{acc}-F{frag}-model_v{ver}.{ext}"
AFDB_MISSING_HINT = ("AlphaFold DB does not cover viral proteins, and proteins longer than 2,700 residues "
                     "are available only as fragments in bulk proteome downloads; also check the accession.")
USER_AGENT = "pdbKit/1.0 (+https://github.com; structural-biology toolkit)"

PDB_ID_RE = re.compile(r"^[0-9][A-Za-z0-9]{3}$")
UNIPROT_RE = re.compile(r"^([OPQ][0-9][A-Z0-9]{3}[0-9]|[A-NR-Z][0-9]([A-Z][A-Z0-9]{2}[0-9]){1,2})$")
AF_ID_RE = re.compile(r"^AF-([A-Z0-9]+)(?:-F(\d+))?(?:-model_v(\d+))?$", re.IGNORECASE)


def default_cache_dir() -> str:
    return os.environ.get("PDBKIT_CACHE", os.path.join(Path.home(), ".cache", "pdbkit"))


def classify_identifier(ident: str) -> Optional[str]:
    """``'pdb'`` for a PDB ID, ``'alphafold'`` for an AF-/UniProt ID, else None."""
    s = ident.strip()
    if PDB_ID_RE.match(s):
        return "pdb"
    if AF_ID_RE.match(s) or UNIPROT_RE.match(s.upper()):
        return "alphafold"
    return None


def _http_get(url: str, timeout: float = 60, retries: int = 2) -> bytes:
    last: Optional[Exception] = None
    for attempt in range(retries + 1):
        try:
            req = urllib.request.Request(url, headers={"User-Agent": USER_AGENT})
            with urllib.request.urlopen(req, timeout=timeout) as resp:
                return resp.read()
        except urllib.error.HTTPError as exc:
            if exc.code in (400, 403, 404, 410):
                raise
            last = exc
        except (urllib.error.URLError, TimeoutError, ConnectionError, OSError) as exc:
            last = exc
        time.sleep(1.5 * (attempt + 1))
    raise PdbKitError(f"Download failed for {url}: {last}")


def _write_atomic(data: bytes, path: str) -> str:
    os.makedirs(os.path.dirname(os.path.abspath(path)), exist_ok=True)
    fd, tmp = tempfile.mkstemp(dir=os.path.dirname(os.path.abspath(path)), suffix=".part")
    with os.fdopen(fd, "wb") as fh:
        fh.write(data)
    os.replace(tmp, path)
    return path


def fetch_pdb(pdb_id: str, fmt: str = "cif", *, assembly: Optional[int] = None,
              outdir: str = ".", overwrite: bool = False, filename: Optional[str] = None) -> str:
    """Download an RCSB entry (or biological assembly ``assembly``). Returns the local path."""
    pid = pdb_id.strip().upper()
    if not PDB_ID_RE.match(pid):
        raise PdbKitError(f"{pdb_id!r} is not a 4-character PDB ID.")
    fmt = "cif" if fmt in ("cif", "mmcif") else "pdb"
    if assembly:
        remote = f"{pid}-assembly{assembly}.cif" if fmt == "cif" else f"{pid}.pdb{assembly}"
        local = filename or (f"{pid}-assembly{assembly}.cif" if fmt == "cif" else f"{pid}-assembly{assembly}.pdb")
    else:
        remote = f"{pid}.{fmt}"
        local = filename or remote
    path = os.path.join(outdir, local)
    if os.path.exists(path) and not overwrite:
        return path
    try:
        data = _http_get(RCSB_URL.format(name=remote))
    except urllib.error.HTTPError as exc:
        if exc.code == 404 and fmt == "pdb":
            raise PdbKitError(
                f"RCSB has no PDB-format file for {pid}{' assembly ' + str(assembly) if assembly else ''} "
                "(large or multi-character-chain entries are distributed as mmCIF only). "
                "Use --format cif.") from exc
        raise PdbKitError(f"RCSB returned HTTP {exc.code} for {remote} - check the ID"
                          f"{' and assembly number' if assembly else ''}.") from exc
    return _write_atomic(data, path)


def alphafold_entries(accession: str) -> List[dict]:
    """Model records for a UniProt accession from the AlphaFold DB API."""
    data = _http_get(AFDB_API.format(acc=accession.upper()))
    entries = json.loads(data.decode("utf-8"))
    return entries if isinstance(entries, list) else [entries]


def fetch_alphafold(identifier: str, fmt: str = "cif", *, fragment: int = 1,
                    version: Optional[int] = None, outdir: str = ".", overwrite: bool = False,
                    pae: bool = False) -> str:
    """Download an AlphaFold DB model by UniProt accession or ``AF-<acc>-F<n>`` ID.

    With ``pae=True`` the predicted-aligned-error JSON is saved next to the model.
    Returns the path of the structure file.
    """
    m = AF_ID_RE.match(identifier.strip())
    if m:
        acc = m.group(1).upper()
        fragment = int(m.group(2) or fragment)
        version = int(m.group(3)) if m.group(3) else version
    else:
        acc = identifier.strip().upper()
    fmt = "cif" if fmt in ("cif", "mmcif") else "pdb"
    model_id = f"AF-{acc}-F{fragment}"

    url = pae_url = None
    ver = version
    if version is None:
        try:
            for entry in alphafold_entries(acc):
                if entry.get("modelEntityId", "").upper() == model_id.upper() or \
                        entry.get("entryId", "").upper() == model_id.upper():
                    url = entry.get("cifUrl" if fmt == "cif" else "pdbUrl")
                    pae_url = entry.get("paeDocUrl")
                    ver = entry.get("latestVersion")
                    break
        except urllib.error.HTTPError as exc:
            if exc.code == 404:
                raise PdbKitError(f"AlphaFold DB has no model for {acc}. {AFDB_MISSING_HINT}") from exc
        except PdbKitError:
            pass  # API unreachable - fall back to direct URLs below
    local = os.path.join(outdir, f"{model_id}.{fmt}")
    if os.path.exists(local) and not overwrite:
        return local
    candidates = [url] if url else [AFDB_FILE.format(acc=acc, frag=fragment, ver=v, ext=fmt)
                                    for v in ([version] if version else range(6, 0, -1))]
    data = None
    for cand in candidates:
        try:
            data = _http_get(cand)
            url = cand
            break
        except urllib.error.HTTPError:
            continue
    if data is None:
        raise PdbKitError(f"AlphaFold DB has no {fmt} file for {model_id}"
                          f"{' version ' + str(version) if version else ''}. {AFDB_MISSING_HINT}")
    _write_atomic(data, local)
    if pae:
        if pae_url is None and url:
            pae_url = re.sub(r"model_v(\d+)\.\w+$", r"predicted_aligned_error_v\1.json", url)
        try:
            _write_atomic(_http_get(pae_url), os.path.join(outdir, f"{model_id}_pae.json"))
        except Exception as exc:  # noqa: BLE001
            print(f"pdbKit: could not download PAE for {model_id}: {exc}", file=sys.stderr)
    return local


def fetch(identifier: str, fmt: str = "cif", *, source: str = "auto", assembly: Optional[int] = None,
          outdir: str = ".", overwrite: bool = False, pae: bool = False) -> str:
    """Download a structure by PDB ID or AlphaFold/UniProt ID (``source='auto'`` decides)."""
    kind = source if source != "auto" else classify_identifier(identifier)
    if kind == "pdb":
        return fetch_pdb(identifier, fmt, assembly=assembly, outdir=outdir, overwrite=overwrite)
    if kind == "alphafold":
        if assembly:
            raise PdbKitError("AlphaFold models have no biological assemblies.")
        return fetch_alphafold(identifier, fmt, outdir=outdir, overwrite=overwrite, pae=pae)
    raise PdbKitError(f"Cannot tell whether {identifier!r} is a PDB ID or a UniProt accession; "
                      "pass --source pdb or --source alphafold.")
