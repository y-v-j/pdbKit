"""Chain splitting / selection modes (features 3-4) and centre of mass (feature 5)."""

import itertools

import numpy as np
import pytest

import pdbkit as pk
from pdbkit.select import Selection, parse_ranges
from conftest import DATA, ENTRIES


def names(s):
    return [(r.get_parent().id, r.get_resname().strip(), a.get_id())
            for r in s.get_residues() for a in r.get_unpacked_list()]


@pytest.mark.parametrize("ext", ["pdb", "cif"])
class TestSplitModes:
    """The four representations asked for, on 4HHB (HEM, PO4, water) and 1WHZ (MSE)."""

    def test_ca_only(self, ext):
        parts = pk.split_chains(pk.load(DATA / f"4HHB.{ext}"), atoms="ca", het="none", water=False)
        assert sorted(parts) == ["A", "B", "C", "D"]
        for cid, p in parts.items():
            assert {a for _, _, a in names(p)} == {"CA"}
            assert {c for c, _, _ in names(p)} == {cid}
        assert len(names(parts["A"])) == 141 and len(names(parts["B"])) == 146

    def test_all_atom(self, ext):
        s = pk.load(DATA / f"4HHB.{ext}")
        p = pk.split_chains(s, atoms="all", het="none", water=False)["A"]
        ref = [(c, r, a) for c, r, a in names(s) if c == "A" and r not in ("HEM", "HOH", "PO4")]
        assert names(p) == ref

    def test_ca_with_all_het(self, ext):
        p = pk.split_chains(pk.load(DATA / f"4HHB.{ext}"), atoms="ca", het="all", water=False)["B"]
        res = {r for _, r, _ in names(p)}
        assert "HEM" in res and "PO4" in res and "HOH" not in res
        hem_atoms = [a for _, r, a in names(p) if r == "HEM"]
        assert len(hem_atoms) == 43                       # the whole ligand, not just a "CA"

    def test_all_atom_with_specific_het_and_water(self, ext):
        p = pk.split_chains(pk.load(DATA / f"4HHB.{ext}"), atoms="all", het=["HEM"], water=True)["A"]
        res = {r for _, r, _ in names(p)}
        assert "HEM" in res and "HOH" in res and "PO4" not in res

    def test_selected_chains_only(self, ext):
        parts = pk.split_chains(pk.load(DATA / f"4HHB.{ext}"), chains="A,C", atoms="ca", het="HEM")
        assert sorted(parts) == ["A", "C"]

    def test_unknown_chain_is_an_error(self, ext):
        with pytest.raises(pk.PdbKitError, match="not found"):
            pk.split_chains(pk.load(DATA / f"4HHB.{ext}"), chains="Z")

    def test_modified_residue_kept_in_ca_mode(self, ext):
        s = pk.load(DATA / f"1WHZ.{ext}")
        p = pk.extract(s, atoms="ca", het="none", water=False)
        n_poly = sum(1 for r in s.get_residues() if pk.residue_kind(r) == "polymer")
        assert len(names(p)) == n_poly                    # one CA per residue incl. the MSEs
        assert "MSE" in {r for _, r, _ in names(p)}

    def test_dna_trace_and_empty_ca(self, ext, capsys):
        s = pk.load(DATA / f"1BNA.{ext}")
        assert pk.split_chains(s, atoms="ca", het="none", water=False) == {}   # DNA has no C-alpha
        assert "skipped" in capsys.readouterr().err
        tr = pk.split_chains(s, atoms="trace", het="none", water=False)
        assert {a for _, _, a in names(tr["A"])} == {"P"} and len(names(tr["A"])) == 11   # 5' end has no P

    def test_nmr_split_keeps_models(self, ext):
        p = pk.split_chains(pk.load(DATA / f"1D3Z.{ext}"), atoms="ca")["A"]
        assert len(p) == 10


def test_split_outputs_are_identical_for_pdb_and_cif_input(tmp_path):
    """Every chain x representation: PDB and mmCIF input give byte-identical PDB output."""
    modes = [dict(atoms=a, het=h, water=w) for a, h, w in
             itertools.product(["ca", "all", "backbone"], ["none", "all", ["HEM"]], [False, True])]
    a, b = pk.load(DATA / "4HHB.pdb"), pk.load(DATA / "4HHB.cif")
    for m in modes:
        pa, pb = pk.split_chains(a, **m), pk.split_chains(b, **m)
        assert pa.keys() == pb.keys()
        for c in pa:
            assert pk.dumps(pa[c], "pdb") == pk.dumps(pb[c], "pdb"), (m, c)


def test_write_parts_names(tmp_path):
    s = pk.load(DATA / "4HHB.cif")
    sel = Selection(atoms="ca", het=["HEM", "PO4"], water=False)
    paths = pk.write_parts(pk.split_chains(s, selection=sel), str(tmp_path), "4HHB", "cif",
                           pk.select.describe(sel))
    assert [p.split("/")[-1] for p in paths] == [f"4HHB_{c}_ca_het-HEM-PO4.cif" for c in "ABCD"]


def test_parse_ranges():
    assert parse_ranges("1:10,20,30-40,50:,:5,1:100:5") == [
        (1, 10, 1), (20, 20, 1), (30, 40, 1), (50, None, 1), (None, 5, 1), (1, 100, 5)]
    assert parse_ranges("-5:-1") == [(-5, -1, 1)]
    with pytest.raises(pk.PdbKitError):
        parse_ranges("a:b")


# --------------------------------------------------------------------------- #
# centre of mass
# --------------------------------------------------------------------------- #

def independent_com(path, chains=None):
    """COM from a fixed-column parse and the IUPAC masses (no Biopython)."""
    from Bio.Data.IUPACData import atom_weights
    xyz, w = [], []
    for ln in open(path):
        if ln.startswith(("ATOM", "HETATM")) and (chains is None or ln[21] in chains) and ln[16] in " A":
            xyz.append([float(ln[30:38]), float(ln[38:46]), float(ln[46:54])])
            w.append(atom_weights[ln[76:78].strip().capitalize()])
    xyz, w = np.array(xyz), np.array(w)
    return (xyz * w[:, None]).sum(0) / w.sum()


@pytest.mark.parametrize("entry", ["1CRN", "4HHB", "4MMK", "1AKE"])
def test_com_against_independent_calculation(entry):
    ref = independent_com(DATA / f"{entry}.pdb")
    for ext in ("pdb", "cif"):
        assert np.allclose(pk.center_of_mass(pk.load(DATA / f"{entry}.{ext}")), ref, atol=2e-3)


@pytest.mark.parametrize("ext", ["pdb", "cif"])
def test_center_whole_structure_and_chain(ext, tmp_path):
    s = pk.load(DATA / f"4HHB.{ext}")
    new, com = pk.center(s)
    out = tmp_path / f"c.{ext}"
    pk.save(new, out)
    back = pk.load(out)
    assert np.allclose(pk.center_of_mass(back), 0, atol=1e-3)            # written file is centred
    assert np.allclose(com, independent_com(DATA / "4HHB.pdb"), atol=2e-3)

    # chain A only, written alone
    ca, com_a = pk.center(s, chains="A")
    assert {c.id for c in ca.get_chains()} == {"A"}
    assert np.allclose(pk.center_of_mass(ca), 0, atol=1e-4)
    assert np.allclose(com_a, independent_com(DATA / "4HHB.pdb", chains="A"), atol=2e-3)

    # chain A's COM to the origin, whole tetramer kept
    allc, _ = pk.center(s, chains="A", apply_to="all")
    assert len(list(allc.get_atoms())) == len(list(s.get_atoms()))
    assert np.allclose(pk.center_of_mass(allc, chains="A"), 0, atol=1e-4)
    # rigid translation: all interatomic distances unchanged
    x0, x1 = pk.geometry.coordinates(s), pk.geometry.coordinates(allc)
    assert np.allclose(x1 - x0, (x1 - x0)[0], atol=1e-4)


def test_center_ca_selection_and_geometric():
    s = pk.load(DATA / "1CRN.cif")
    new, com = pk.center(s, atoms="ca", weighted=False)
    xyz = pk.geometry.coordinates(new)
    assert np.allclose(xyz.mean(0), 0, atol=1e-4) and len(xyz) == 46


def test_center_moves_all_altlocs():
    s = pk.load(DATA / "3NIR.cif")
    new, com = pk.center(s)
    for a, b in zip(pk.structure.iter_atoms_all_altlocs(s), pk.structure.iter_atoms_all_altlocs(new)):
        assert np.allclose(b.coord, a.coord - com, atol=1e-4)


def test_center_nmr_ensemble_preserves_superposition():
    s = pk.load(DATA / "1D3Z.pdb")
    new, com = pk.center(s)
    first = np.array([a.coord for a in s[0].get_atoms()])
    for m_old, m_new in zip(s, new):
        d = np.array([a.coord for a in m_new.get_atoms()]) - np.array([a.coord for a in m_old.get_atoms()])
        assert np.allclose(d, -com, atol=1e-4)
    per, _ = pk.center(s, per_model=True)
    for m in per:
        xyz = np.array([a.coord for a in m.get_atoms()])
        w = pk.geometry.atom_masses(list(m.get_atoms()), warn=False)
        assert np.allclose((xyz * w[:, None]).sum(0) / w.sum(), 0, atol=1e-4)


def test_radius_of_gyration_translation_invariant():
    s = pk.load(DATA / "1CRN.pdb")
    assert pk.radius_of_gyration(s) == pytest.approx(pk.radius_of_gyration(pk.translate(s, [10, -20, 5])), abs=1e-6)


@pytest.mark.parametrize("entry", ["3NIR", "1AKE", "1BRS"])
def test_copies_keep_alternate_location_selection(entry):
    """Found while writing the protein notebook: after Biopython's copy(), a disordered atom's
    *selected* location still pointed at the original structure, so in-memory reads of a
    moved copy (get_atoms, centre of mass) returned the old coordinates for altloc atoms."""
    s = pk.load(DATA / f"{entry}.cif")
    moved, com = pk.center(s, apply_to="all")
    x0, x1 = pk.geometry.coordinates(s), pk.geometry.coordinates(moved)
    assert np.allclose(x1 - x0, -com, atol=1e-3)                       # every atom as read through wrappers
    assert np.linalg.norm(pk.center_of_mass(moved)) < 1e-3               # without writing / re-reading
    for a in moved.get_atoms():
        if a.is_disordered():
            assert any(a.selected_child is c for c in a.disordered_get_list())
    assert np.allclose(pk.geometry.coordinates(pk.load(DATA / f"{entry}.cif")), x0)   # input untouched
