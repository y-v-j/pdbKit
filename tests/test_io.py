"""Reading and writing: format detection, parity between PDB and mmCIF, round trips."""

import collections
import gzip
import io

import numpy as np
import pytest

import pdbkit as pk
from pdbkit.cifio import cif_quote
from conftest import DATA, ENTRIES, atom_records, core


def _atoms(s):
    """(model, chain, resseq, icode, resname, atom, altloc) -> coord, over all altlocs."""
    out = {}
    for m in s:
        for r in pk.structure.iter_residues_all(m):      # includes microheterogeneity variants
            for a in r.get_unpacked_list():
                key = (m.serial_num, r.get_parent().id, r.id[1], r.id[2], r.get_resname(), a.get_id(), a.altloc)
                out[key] = np.round(a.coord.astype(float), 3)
    return out


@pytest.mark.parametrize("entry", ENTRIES)
def test_pdb_and_cif_inputs_give_the_same_atoms(entry):
    """The same entry read from .pdb and .cif yields identical atoms and coordinates."""
    a, b = _atoms(pk.load(DATA / f"{entry}.pdb")), _atoms(pk.load(DATA / f"{entry}.cif"))
    assert a.keys() == b.keys()
    assert all(np.allclose(a[k], b[k], atol=1e-3) for k in a)


@pytest.mark.parametrize("entry", ENTRIES)
@pytest.mark.parametrize("src", ["pdb", "cif"])
@pytest.mark.parametrize("dst", ["pdb", "cif"])
def test_round_trip_preserves_atoms_and_classification(entry, src, dst):
    s = pk.load(DATA / f"{entry}.{src}")
    s2 = pk.loads(pk.dumps(s, dst), dst)
    a, b = _atoms(s), _atoms(s2)
    assert a.keys() == b.keys()
    assert all(np.allclose(a[k], b[k], atol=1e-3) for k in a)
    kind = lambda st: collections.Counter(pk.residue_kind(r) for r in pk.structure.iter_residues_all(st))
    assert kind(s) == kind(s2)


def test_pdb_output_matches_input_records_exactly():
    """Independent fixed-column parse: all ATOM/HETATM fields survive PDB->PDB.

    Order may differ: Biopython groups each chain's HETATMs with that chain,
    where 4HHB.pdb lists them after all chains. Compared as sorted multisets.
    """
    for e in ENTRIES:
        text = (DATA / f"{e}.pdb").read_text()
        out = pk.dumps(pk.load(DATA / f"{e}.pdb"), "pdb")
        assert sorted(core(atom_records(text), True)) == sorted(core(atom_records(out), True)), e


def test_classification_parity_pdb_heuristic_vs_mmcif_labels():
    """Polymer/ligand/water from PDB heuristics equals mmCIF's label_seq_id-based truth."""
    for e in ENTRIES:
        ka = [(r.get_parent().id, r.id[1], pk.residue_kind(r)) for r in pk.load(DATA / f"{e}.pdb").get_residues()]
        kb = [(r.get_parent().id, r.id[1], pk.residue_kind(r)) for r in pk.load(DATA / f"{e}.cif").get_residues()]
        assert ka == kb, e


def test_modified_residue_mse_is_polymer():
    for ext in ("pdb", "cif"):
        s = pk.load(DATA / f"1WHZ.{ext}")
        mse = [r for r in s.get_residues() if r.get_resname() == "MSE"]
        assert mse and all(pk.residue_kind(r) == "polymer" for r in mse)
        assert all(r.id[0] == "H_MSE" for r in mse)       # still HETATM records


def test_capped_peptide_inhibitor_is_polymer():
    """6LU7 chain C (N3 inhibitor) incl. ester-linked C-terminal cap 010 is one polymer."""
    s = pk.load(DATA / "6LU7.pdb")
    assert {pk.residue_kind(r) for r in s[0]["C"]} == {"polymer"}


def test_format_detection_and_sources(tmp_path):
    text = (DATA / "1CRN.cif").read_text()
    assert pk.detect_format(text) == "cif"
    assert pk.detect_format((DATA / "1CRN.pdb").read_text()) == "pdb"
    gz = tmp_path / "x.cif.gz"
    with gzip.open(gz, "wt") as fh:
        fh.write(text)
    s = pk.load(gz)
    assert pk.source_info(s).fmt == "cif" and s.id == "1CRN"
    noext = tmp_path / "structure_without_extension"
    noext.write_text(text)
    assert pk.source_info(pk.load(noext)).fmt == "cif"
    s2 = pk.loads((DATA / "1CRN.pdb").read_text())
    assert len(list(s2.get_atoms())) == 327
    buf = io.StringIO()
    pk.save(s2, buf, "cif")
    assert buf.getvalue().startswith("data_1CRN")


def test_empty_input_is_an_error():
    with pytest.raises(pk.PdbKitError):
        pk.loads("REMARK nothing here\nEND\n", "pdb")


def test_pdb_writer_refuses_multichar_chains_and_remaps_on_request():
    s = pk.load(DATA / "1CRN.cif")
    s[0]["A"].id = "AB"
    with pytest.raises(pk.PdbKitError, match="chain IDs longer"):
        pk.dumps(s, "pdb")
    out = pk.dumps(s, "pdb", remap_chains=True)
    assert {r[5] for r in atom_records(out)} == {"B"} or len({r[5] for r in atom_records(out)}) == 1
    assert "AB" in pk.dumps(s, "cif")          # mmCIF keeps it


def test_pdb_writer_refuses_long_residue_names():
    s = pk.load(DATA / "1CRN.pdb")
    next(s.get_residues()).resname = "ABCDE"
    with pytest.raises(pk.PdbKitError, match="residue names"):
        pk.dumps(s, "pdb")


def test_cif_quoting():
    assert cif_quote("O5'") == '"O5\'"'
    assert cif_quote("A B") == "'A B'"
    assert cif_quote("_x") == "'_x'"
    assert cif_quote("data_x") == "'data_x'"
    assert cif_quote("") == "?"
    assert cif_quote("CA") == "CA"


def test_nucleic_acid_primes_survive_cif():
    s = pk.load(DATA / "1BNA.pdb")
    s2 = pk.loads(pk.dumps(s, "cif"), "cif")
    assert "O5'" in {a.get_id() for a in s2.get_atoms()}


def test_cif_writer_entities_and_label_seq_gaps():
    """Each ligand gets its own asym; label_seq_id keeps chain breaks visible."""
    s = pk.load(DATA / "4HHB.pdb")
    d = pk.metadata.as_mmcif_dict(pk.loads(pk.dumps(s, "cif"), "cif"))
    types = collections.Counter(d["_entity.type"])
    assert types == {"polymer": 2, "non-polymer": 2, "water": 1}     # alpha, beta, HEM, PO4, water
    assert len(d["_struct_asym.id"]) == 4 + 6 + 4                    # chains, ligand instances, waters
    # 1AKE chain A: preserved mmCIF label_seq_id equals the original file's values
    orig = pk.load(DATA / "1AKE.cif")
    d2 = pk.metadata.as_mmcif_dict(pk.loads(pk.dumps(orig, "cif"), "cif"))
    o = pk.metadata.as_mmcif_dict(orig)
    pairs = lambda dd: {(a, s_, l) for a, s_, l in zip(dd["_atom_site.auth_asym_id"], dd["_atom_site.auth_seq_id"],
                                                       dd["_atom_site.label_seq_id"]) if l != "."}
    assert pairs(o) == pairs(d2)


@pytest.mark.network
def test_load_by_id_fetches_and_caches(tmp_path, monkeypatch):
    monkeypatch.setenv("PDBKIT_CACHE", str(tmp_path))
    s = pk.load("1CRN")
    assert s.id == "1CRN" and (tmp_path / "1CRN.cif").exists()
