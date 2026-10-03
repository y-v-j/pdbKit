"""Shared fixtures and oracle helpers for the pdbKit test suite.

Oracles used (all independent of pdbKit's own code paths):

* ``upstream/pdb-tools-2.0.0-rc1`` scripts, run as subprocesses;
* ``upstream/pdbtools-master`` (Harms lab) functions, imported from source;
* ``upstream/hetpdbnavi_parser.py`` and ``upstream/mmcif_domain_family.py``;
* PyMOL (``pymol`` module in the test interpreter) and gemmi (``py312`` env);
* a minimal fixed-column PDB parser defined here, so coordinate comparisons
  do not go through Biopython.
"""

from __future__ import annotations

import os
import shutil
import subprocess
import sys
import warnings
from pathlib import Path

import pytest

warnings.filterwarnings("ignore", message="You may be importing Biopython")

ROOT = Path(__file__).resolve().parents[1]
DATA = ROOT / "tests" / "data"
UPSTREAM = ROOT / "upstream"
PDBTOOLS = UPSTREAM / "pdb-tools-2.0.0-rc1" / "pdbtools"
HARMS = UPSTREAM / "pdbtools-master"
GEMMI_PY = os.environ.get("PDBKIT_GEMMI_PYTHON", "/home/kestrel/anaconda3/envs/py312/bin/python")
ENV = dict(os.environ, PYTHONWARNINGS="ignore")

ENTRIES = ["1AKE", "1BNA", "1BRS", "1CRN", "1D3Z", "1EHZ", "1K4C", "1L2Y", "1WHZ", "2GB1", "3NIR", "4HHB", "4MMK", "6LU7"]


@pytest.fixture
def data():
    return DATA


def run_pdbtool(tool: str, *args: str, stdin: str | None = None) -> str:
    """Run an upstream pdb-tools script and return stdout."""
    cmd = [sys.executable, str(PDBTOOLS / f"{tool}.py"), *args]
    res = subprocess.run(cmd, input=stdin, capture_output=True, text=True, env=ENV)
    if res.returncode not in (0,):
        raise RuntimeError(f"{tool} failed: {res.stderr}")
    return res.stdout


def harms():
    """Import the Harms-lab package from upstream/ (Python 3 compatible modules only)."""
    if str(HARMS) not in sys.path:
        sys.path.insert(0, str(HARMS))
    with warnings.catch_warnings():
        warnings.simplefilter("ignore")
        import pdbtools  # noqa: F401  (the Harms package)
    return sys.modules["pdbtools"]


def atom_records(text: str, all_models: bool = True):
    """Independent fixed-column parse of ATOM/HETATM lines.

    Returns tuples (model, record, name, altloc, resname, chain, resseq, icode,
    x, y, z, occ, b, element); serial numbers are ignored (writers renumber).
    """
    out, model = [], 1
    for ln in text.splitlines():
        if ln.startswith("MODEL"):
            model = int(ln[10:14])
        elif ln.startswith(("ATOM  ", "HETATM")):
            ln = ln.ljust(80)
            out.append((model if all_models else 1, ln[:6].strip(), ln[12:16].strip(), ln[16].strip(),
                        ln[17:20].strip(), ln[21], int(ln[22:26]), ln[26].strip(),
                        round(float(ln[30:38]), 3), round(float(ln[38:46]), 3), round(float(ln[46:54]), 3),
                        round(float(ln[54:60] or 0), 2), round(float(ln[60:66] or 0), 2), ln[76:78].strip()))
    return out


def core(records, with_record=False):
    """Comparable subset: identity + coordinates (+ occupancy/B)."""
    return [((r[1],) if with_record else ()) + r[2:13] for r in records]


def gemmi_run(code: str, *args: str) -> str:
    if not Path(GEMMI_PY).exists():
        pytest.skip("gemmi interpreter not available")
    res = subprocess.run([GEMMI_PY, "-c", code, *args], capture_output=True, text=True, env=ENV)
    if res.returncode:
        raise RuntimeError(res.stderr)
    return res.stdout


def have_pymol() -> bool:
    try:
        import pymol  # noqa: F401
        return True
    except ImportError:
        return False


def online() -> bool:
    if os.environ.get("PDBKIT_OFFLINE"):
        return False
    import urllib.request
    try:
        urllib.request.urlopen("https://files.rcsb.org/download/1CRN.pdb", timeout=10).read(10)
        return True
    except Exception:  # noqa: BLE001
        return False


_ONLINE = None


def pytest_collection_modifyitems(config, items):
    global _ONLINE
    if any("network" in it.keywords for it in items):
        _ONLINE = online()
        if not _ONLINE:
            skip = pytest.mark.skip(reason="no network")
            for it in items:
                if "network" in it.keywords:
                    it.add_marker(skip)
