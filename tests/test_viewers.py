"""Requirement 9: pdbKit output must open correctly in molecular viewers.

PyMOL (headless, separate process) and gemmi (the parser behind many tools,
run from another interpreter) read the *original* wwPDB file and pdbKit's
PDB and mmCIF outputs. Atom counts, models, chains and PyMOL's own
polymer / organic / inorganic / solvent classification must agree.
"""

import json
import subprocess
import sys

import pytest

import pdbkit as pk
from conftest import DATA, ENV, gemmi_run, have_pymol

pytestmark = pytest.mark.oracle

PYMOL_SCRIPT = r"""
import json, sys
from pymol import cmd
out = {}
for path in sys.argv[1:]:
    cmd.reinitialize()
    cmd.load(path, "m")
    out[path] = {
        # 'present' excludes the coordinate-less CA placeholders PyMOL builds for
        # unmodelled residues from a wwPDB mmCIF's deposited sequence
        "atoms": cmd.count_atoms("m and present"), "states": cmd.count_states("m"),
        "chains": sorted(cmd.get_chains("m")),
        "polymer": cmd.count_atoms("m and present and polymer"),
        "protein": cmd.count_atoms("m and present and polymer.protein"),
        "nucleic": cmd.count_atoms("m and present and polymer.nucleic"),
        "organic": cmd.count_atoms("m and organic"),
        "inorganic": cmd.count_atoms("m and inorganic"), "solvent": cmd.count_atoms("m and solvent"),
        "ss": cmd.count_atoms("m and name CA and ss H+S"),
    }
print(json.dumps(out))
"""

GEMMI_SCRIPT = r"""
import json, sys, gemmi
out = {}
for path in sys.argv[1:]:
    st = gemmi.read_structure(path)
    st.setup_entities()
    out[path] = {"models": len(st), "atoms": sum(1 for m in [st[0]] for c in m for r in c for a in r),
                 "chains": sorted({c.name for c in st[0]}),
                 "entity_types": sorted(str(e.entity_type).split(".")[-1] for e in st.entities)}
    if path.endswith(".cif"):
        doc = gemmi.cif.read(path)             # strict CIF syntax check
        out[path]["blocks"] = len(doc)
print(json.dumps(out))
"""

CASES = ["4HHB", "1WHZ", "3NIR", "1BNA", "1D3Z", "6LU7", "4MMK", "1AKE", "1EHZ", "1K4C"]


def pymol(paths):
    if not have_pymol():
        pytest.skip("PyMOL not installed")
    res = subprocess.run([sys.executable, "-c", PYMOL_SCRIPT, *map(str, paths)],
                         capture_output=True, text=True, env=ENV)
    assert res.returncode == 0, res.stderr
    return json.loads(res.stdout.strip().splitlines()[-1])


def outputs(entry, tmp_path):
    s = pk.load(DATA / f"{entry}.cif")
    return {
        "orig_pdb": DATA / f"{entry}.pdb", "orig_cif": DATA / f"{entry}.cif",
        "ours_pdb": pk.save(s, tmp_path / f"{entry}_ours.pdb"),
        "ours_cif": pk.save(s, tmp_path / f"{entry}_ours.cif"),
        "ours_pdb_from_pdb": pk.save(pk.load(DATA / f"{entry}.pdb"), tmp_path / f"{entry}_p2p.pdb"),
        "ours_cif_from_pdb": pk.save(pk.load(DATA / f"{entry}.pdb"), tmp_path / f"{entry}_p2c.cif"),
    }


@pytest.mark.parametrize("entry", CASES)
def test_pymol_reads_outputs_like_the_originals(entry, tmp_path):
    files = outputs(entry, tmp_path)
    r = pymol(files.values())
    ref = r[str(files["orig_pdb"])]
    for key in ("ours_pdb", "ours_pdb_from_pdb"):
        got = r[str(files[key])]
        assert {k: got[k] for k in got if k != "ss"} == {k: ref[k] for k in ref if k != "ss"}, key
    ref_c = r[str(files["orig_cif"])]
    for key in ("ours_cif", "ours_cif_from_pdb"):
        got = r[str(files[key])]
        for k in ("atoms", "states", "chains", "polymer", "protein", "nucleic", "organic", "inorganic", "solvent"):
            assert got[k] == ref_c[k], (key, k, got[k], ref_c[k])


def test_pymol_reads_split_and_centred_files(tmp_path):
    s = pk.load(DATA / "4HHB.cif")
    paths = []
    for ext in ("pdb", "cif"):
        parts = pk.split_chains(s, atoms="ca", het="all", water=False)
        paths += pk.write_parts(parts, str(tmp_path), f"4HHB_{ext}", ext, "ca_het")
        paths.append(pk.save(pk.center(s, chains="B")[0], tmp_path / f"centred.{ext}"))
    r = pymol(paths)
    for p in paths:
        if "_A_" in p:
            assert r[p]["protein"] == 141 and r[p]["organic"] == 43 and r[p]["solvent"] == 0
        if "centred" in p:
            ref = pk.extract(s, chains="B")
            assert r[p]["chains"] == ["B"] and r[p]["atoms"] == len(list(pk.structure.iter_atoms_all_altlocs(ref)))


@pytest.mark.parametrize("entry", CASES)
def test_gemmi_reads_outputs(entry, tmp_path):
    files = outputs(entry, tmp_path)
    r = json.loads(gemmi_run(GEMMI_SCRIPT, *map(str, files.values())).strip().splitlines()[-1])
    ref = r[str(files["orig_cif"])]
    for key in ("ours_pdb", "ours_cif", "ours_pdb_from_pdb", "ours_cif_from_pdb"):
        got = r[str(files[key])]
        assert got["models"] == ref["models"] and got["atoms"] == ref["atoms"] and got["chains"] == ref["chains"], key
    for key in ("ours_cif", "ours_cif_from_pdb"):
        got = r[str(files[key])]
        assert got["blocks"] == 1
        assert set(got["entity_types"]) == set(ref["entity_types"]), key
