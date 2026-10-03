"""The example notebooks in docs/examples must execute without errors (run in a temporary copy)."""

import shutil
from pathlib import Path

import pytest

from conftest import ROOT

nbformat = pytest.importorskip("nbformat")
nbclient = pytest.importorskip("nbclient")

EXAMPLES = ROOT / "docs" / "examples"


@pytest.mark.parametrize("name", ["pdbKit_protein_EGFR_erlotinib.ipynb", "pdbKit_dna_1BNA_Zif268.ipynb"])
def test_example_notebook_runs(name, tmp_path):
    shutil.copy(EXAMPLES / name, tmp_path / name)
    shutil.copytree(EXAMPLES / "data", tmp_path / "data")
    nb = nbformat.read(tmp_path / name, as_version=4)
    nbclient.NotebookClient(nb, timeout=600, kernel_name="python3",
                            resources={"metadata": {"path": str(tmp_path)}}).execute()
    errors = [o for c in nb.cells if c.cell_type == "code" for o in c.get("outputs", []) if o.output_type == "error"]
    assert not errors
    assert any((tmp_path / "output").rglob("*.pdb"))
