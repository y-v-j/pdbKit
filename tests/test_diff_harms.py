"""Differential tests: pdbKit vs the Harms-lab ``pdbtools`` functions (imported from upstream/).

The upstream functions take a list of PDB lines; pdbKit runs on both the PDB
and the mmCIF version of the same entry. Where the original is wrong, the
test demonstrates the defect and checks pdbKit's corrected behaviour.
"""

import math
import re

import numpy as np
import pytest

import pdbkit as pk
from pdbkit.constants import ATOM_WEIGHTS
from pdbkit.select import Selection
from conftest import DATA, atom_records, harms

pytestmark = pytest.mark.oracle


def lines(entry):
    return (DATA / f"{entry}.pdb").read_text().splitlines(keepends=True)


def both(entry):
    return [pk.load(DATA / f"{entry}.pdb"), pk.load(DATA / f"{entry}.cif")]


# --------------------------------------------------------------------------- #
# centre of mass (feature 5)
# --------------------------------------------------------------------------- #

@pytest.mark.parametrize("entry", ["1CRN", "4HHB", "4MMK", "2GB1"])
@pytest.mark.parametrize("het", [False, True])
@pytest.mark.parametrize("mass", [True, False])
def test_centermass(entry, het, mass):
    h = harms()
    center, _ = h.centermass.pdbCentermass(lines(entry), include_hetatm=het, include_mass=mass)
    ref = np.array([float(c) for c in center])
    for s in both(entry):
        sel = None if het else Selection(only_het_records=False)
        com = pk.center_of_mass(s, sel, weighted=mass)
        assert np.allclose(com, ref, atol=2e-4), (com, ref)


def test_centermass_recentred_coordinates_match():
    h = harms()
    _, out = h.centermass.pdbCentermass(lines("4HHB"), write_coord=True, include_hetatm=True)
    ref = sorted((r[5], r[6], r[2], r[8], r[9], r[10]) for r in atom_records("\n".join(out)))
    for s in both("4HHB"):
        new, _ = pk.center(s)
        ours = sorted((r[5], r[6], r[2], r[8], r[9], r[10]) for r in atom_records(pk.dumps(new, "pdb")))
        assert len(ours) == len(ref)
        assert all(a[:3] == b[:3] and np.allclose(a[3:], b[3:], atol=1.1e-3) for a, b in zip(ref, ours))


def test_centermass_upstream_reads_segid_as_element():
    """Upstream takes the element from line[73:], which includes the segment-ID
    columns: with a segid present every atom is 'unknown' and weighted as
    carbon. pdbKit reads the element column."""
    h = harms()
    s = pk.set_segid(pk.load(DATA / "1CRN.pdb"), "PROT")
    text = pk.dumps(s, "pdb").splitlines(keepends=True)
    clean, _ = h.centermass.pdbCentermass(lines("1CRN"))
    segid, _ = h.centermass.pdbCentermass(text)
    geometric, _ = h.centermass.pdbCentermass(lines("1CRN"), include_mass=False)
    f = lambda c: np.array([float(x) for x in c])
    assert np.allclose(f(segid), f(geometric), atol=1e-3)       # all masses equal -> geometric centre
    assert not np.allclose(f(segid), f(clean), atol=1e-2)
    assert np.allclose(pk.center_of_mass(pk.loads("".join(text))), f(clean), atol=2e-4)


def test_tungsten_mass():
    h = harms()
    assert h.data.common.ATOM_WEIGHTS["W"] == 95.94             # upstream: molybdenum's mass
    assert ATOM_WEIGHTS["W"] == pytest.approx(183.84, abs=0.01)


def test_centerasu_matches_and_undoes_lattice_shift():
    h = harms()
    for entry in ("4HHB", "1AKE"):
        fvec, out = h.centerasu.pdbCenterasu(lines(entry))
        for s in both(entry):
            new, shift = pk.center_asu(s)
            assert np.allclose(shift, fvec, atol=1e-6)
            a = sorted(r[8:11] for r in atom_records("".join(out)))
            b = sorted(r[8:11] for r in atom_records(pk.dumps(new, "pdb")))
            assert np.allclose(np.array(a), np.array(b), atol=1.1e-3)
    # move 1AKE by two cell edges along a: centerasu must bring it back
    s = pk.load(DATA / "1AKE.cif")
    a_len = float(pk.metadata.first(pk.metadata.as_mmcif_dict(s), "_cell.length_a"))
    moved = pk.translate(s, [2 * a_len, 0, 0])
    back, shift = pk.center_asu(moved)
    ref, _ = pk.center_asu(s)
    assert np.allclose(pk.geometry.coordinates(back), pk.geometry.coordinates(ref), atol=5e-3)  # SCALE has 6 decimals


# --------------------------------------------------------------------------- #
# contacts
# --------------------------------------------------------------------------- #

@pytest.mark.parametrize("entry,target", [("4HHB", "HEM"), ("1CRN", "PHE"), ("1AKE", "AP5")])
def test_contacting_residues(entry, target):
    h = harms()
    out = h.contact.pdbContacting(lines(entry), target, 3.5)
    ref = set()
    for ln in out:
        tgt, *hits = ln.rstrip("\n").split("\t")
        t_res = tgt[5:].strip()
        for x in hits:
            if x != t_res:                                  # upstream also lists the target itself
                ref.add((re.sub(r"\s+", " ", t_res), re.sub(r"\s+", " ", x)))
    for s in both(entry):
        ours = set()
        for r in pk.contacts(s, target, 3.5, include_water=True):
            (tc, tn, ti), (rc, rn, ri) = r["target"].split(":"), r["residue"].split(":")
            ours.add((f"{tn} {tc} {ti}", f"{rn} {rc} {ri}"))
        assert ours == ref


# Inputs without hydrogens (upstream mis-slices 4-character names) and, for all-atom
# mode, without alternate locations (upstream pairs every altloc line; pdbKit uses
# one location per atom).
@pytest.mark.parametrize("entry,ca_only,cutoff", [("1CRN", False, 3.5), ("1CRN", True, 6.0), ("1AKE", True, 6.0)])
def test_close_contacts(entry, ca_only, cutoff):
    h = harms()
    ref = {(a.split()[0], int(a[-4:]), b.split()[0], int(b[-4:])): d
           for a, b, d in h.closecontacts.pdbCloseContacts(lines(entry), ca_only, cutoff)}
    for s in both(entry):
        ours = {}
        for r in pk.close_contacts(s, cutoff, ca_only):
            (c1, n1, r1, a1), (c2, n2, r2, a2) = (x.split(":") for x in (r["atom1"], r["atom2"]))
            ours[(a1, int(r1), a2, int(r2))] = r["distance"]
        assert ours.keys() == ref.keys()
        assert all(abs(ours[k] - ref[k]) < 1e-3 for k in ref)


def test_close_contacts_upstream_ignores_chain_breaks():
    """The upstream skip rule tests only residue numbers: in 1BNA, G12 (chain A) and
    C13 (chain B) count as 'consecutive', hiding their base-pair contacts."""
    h = harms()
    ref = {(a.split()[0], int(a[-4:]), b.split()[0], int(b[-4:]))
           for a, b, d in h.closecontacts.pdbCloseContacts(lines("1BNA"), False, 3.5)}
    for s in both("1BNA"):
        ours = set()
        for r in pk.close_contacts(s, 3.5):
            (c1, n1, r1, a1), (c2, n2, r2, a2) = (x.split(":") for x in (r["atom1"], r["atom2"]))
            ours.add((a1, int(r1), a2, int(r2)))
        assert ours - ref == {("O6", 12, "N3", 13), ("O6", 12, "N4", 13)}
        assert ref <= ours


def test_disulfide_upstream_cannot_run_on_python3():
    """pdb_disulfide calls clean.stripACS, which uses Python 2's string.letters."""
    h = harms()
    with pytest.raises(AttributeError, match="letters"):
        h.disulfide.pdbDisulfide(lines("1CRN"), 4.0)


def test_disulfide_distance_rule_matches_independent_calculation():
    """Harms' rule (all SG pairs closer than the cutoff), recomputed from raw lines."""
    for entry in ("1CRN", "4MMK", "1BRS"):
        sg = [(ln[21], int(ln[22:26]), np.array([float(ln[30:38]), float(ln[38:46]), float(ln[46:54])]))
              for ln in lines(entry) if ln.startswith("ATOM") and ln[12:16] == " SG " and ln[16] in " A"]
        for cutoff in (2.5, 3.0, 4.0):
            ref = {tuple(sorted([(a[0], a[1]), (b[0], b[1])])) for i, a in enumerate(sg) for b in sg[i + 1:]
                   if np.linalg.norm(a[2] - b[2]) < cutoff}
            for s in both(entry):
                ours = {tuple(sorted((x.split(":")[0], int(x.split(":")[2])) for x in (p["cys1"], p["cys2"])))
                        for p in pk.disulfides(s, cutoff)["bonds"]}
                assert ours == ref


def test_disulfide_cutoff_agrees_with_ssbond_records():
    """With the 3.0 A default, detected bonds equal the file's SSBOND/_struct_conn disulf records."""
    for entry in ("1CRN", "4MMK", "1BRS", "1K4C"):
        for s in both(entry):
            d = pk.metadata.as_mmcif_dict(s)
            recorded = {tuple(sorted([(a, int(b)), (c, int(e))])) for t, a, b, c, e in zip(
                pk.metadata.column(d, "_struct_conn.conn_type_id"),
                pk.metadata.column(d, "_struct_conn.ptnr1_auth_asym_id"),
                pk.metadata.column(d, "_struct_conn.ptnr1_auth_seq_id"),
                pk.metadata.column(d, "_struct_conn.ptnr2_auth_asym_id"),
                pk.metadata.column(d, "_struct_conn.ptnr2_auth_seq_id")) if t == "disulf"}
            found = set()
            for p in pk.disulfides(s)["bonds"]:
                ends = [(x.split(":")[0], int(x.split(":")[2])) for x in (p["cys1"], p["cys2"])]
                found.add(tuple(sorted(ends)))
            assert found == recorded, entry


@pytest.mark.parametrize("entry", ["1CRN"])   # upstream is O(n^2) in pure Python and fails on H atoms
def test_water_contacts(entry):
    h = harms()
    ref = h.watercontact.pdbWaterContact(lines(entry))
    ref = {re.sub(r"\s+", "", k): {re.sub(r"\s+", "", v) for v in vs} - {re.sub(r"\s+", "", k)}
           for k, vs in ref.items()}                        # upstream lists a residue as its own neighbour
    for s in both(entry):
        ours = {}
        for k, vs in pk.water_contacts(s).items():
            key = lambda x: x.split(":")[0] + x.split(":")[2]
            ours[key(k)] = {key(v) for v in vs}
        assert ours == ref


def test_water_contacts_upstream_fails_on_hydrogens():
    h = harms()
    with pytest.raises(KeyError):
        h.watercontact.pdbWaterContact(lines("2GB1"))
    assert pk.water_contacts(pk.load(DATA / "2GB1.cif"))           # works


@pytest.mark.parametrize("entry", ["1CRN", "2GB1", "1AKE"])
def test_neighbors(entry):
    h = harms()
    names, counts = h.neighbors.pdbNeighbors(lines(entry))
    ref = [(int(n[10:14]), c) for n, c in zip(names, counts)]
    for s in both(entry):
        ours = [(int(r["residue"].split(":")[2]), r["neighbors"]) for r in pk.neighbors(s)]
        assert ours == ref


def test_neighbors_upstream_misparses_large_negative_coordinates():
    """pdb_neighbors reads CB coordinates from column 32 (not 31), dropping the
    minus sign of values <= -100: counts change under a pure translation."""
    h = harms()
    s = pk.load(DATA / "1CRN.pdb")
    moved = pk.dumps(pk.translate(s, [-300.0, 0, 0]), "pdb").splitlines(keepends=True)
    _, before = h.neighbors.pdbNeighbors(lines("1CRN"))
    _, after = h.neighbors.pdbNeighbors(moved)
    assert before != after
    ours = lambda st: [r["neighbors"] for r in pk.neighbors(st)]
    assert ours(s) == ours(pk.loads("".join(moved))) == before


@pytest.mark.parametrize("entry", ["1CRN", "2GB1", "1AKE"])
def test_torsions(entry):
    h = harms()
    dihedrals, labels = h.torsion.pdbTorsion(lines(entry))
    ref = {(lab[4], int(lab[5:9])): d for lab, d in zip(labels, dihedrals)}
    for s in both(entry):
        ours = {(r["residue"].split(":")[0], int(r["residue"].split(":")[2])): (r["phi"], r["psi"])
                for r in pk.torsions(s)}
        n = 0
        for key, (phi, psi) in ref.items():
            if key in ours and ours[key][0] is not None and ours[key][1] is not None:
                assert ours[key][0] == pytest.approx(phi, abs=0.02)
                assert ours[key][1] == pytest.approx(psi, abs=0.02)
                n += 1
        assert n >= len(ref) - 2


def test_torsions_upstream_crosses_chain_boundaries():
    """4HHB: upstream reports phi/psi for the last residue of chain A using chain B's first N."""
    h = harms()
    _, labels = h.torsion.pdbTorsion(lines("4HHB"))
    assert any(lab[4] == "A" and int(lab[5:9]) == 141 for lab in labels)
    last_a = [r for r in pk.torsions(pk.load(DATA / "4HHB.pdb")) if r["residue"] == "A:ARG:141"][0]
    assert last_a["psi"] is None


# --------------------------------------------------------------------------- #
# sequence parameters and electrostatics
# --------------------------------------------------------------------------- #

@pytest.mark.parametrize("entry", ["1CRN", "2GB1", "1AKE"])
def test_param_mw_and_pi(entry):
    h = harms()
    chains, _ = h.seq.pdbSeq(lines(entry), use_atoms=True)
    names = [n for c in sorted(chains) for n in chains[c]]
    mw, pi = h.param.calcMW(names), h.param.pdbPi(names)
    for s in both(entry):
        p = pk.param(s)
        assert p["molecular_weight"] == pytest.approx(mw, abs=0.01)
        assert p["pI"] == pytest.approx(pi, abs=0.001)


def test_param_upstream_cannot_run():
    h = harms()
    with pytest.raises(NameError, match="pdb_seq"):
        h.param.pdbParam(lines("1CRN"))


@pytest.mark.parametrize("entry", ["1CRN", "2GB1", "4HHB"])
def test_coulomb(entry):
    h = harms()
    coord, pka, charge = h.coulomb.readPDB(str(DATA / f"{entry}.pdb"))
    ref = h.coulomb.pdbCoulomb(coord, pka, charge, 40.0, 100.0, 7.0, 298.0)
    for s in both(entry):
        assert pk.coulomb(s) == pytest.approx(ref, rel=1e-5)


@pytest.mark.parametrize("entry", ["1CRN", "4HHB"])
def test_iondist(entry):
    h = harms()
    hist = h.iondist.pdbIonDist(lines(entry), 0.5)
    for s in both(entry):
        ours = pk.ion_distance_histogram(s, 0.5)
        assert [r["acid_acid"] for r in ours] == hist[0]
        assert [r["acid_base"] for r in ours] == hist[1]
        assert [r["base_base"] for r in ours] == hist[2]


# --------------------------------------------------------------------------- #
# header-derived reports
# --------------------------------------------------------------------------- #

@pytest.mark.parametrize("entry", ["4HHB", "1AKE", "1BRS", "4MMK"])
def test_ligands(entry):
    h = harms()
    ref = set(filter(None, h.ligand.pdbLigand(lines(entry)).split(";")))
    for s in both(entry):
        ours = {r["resname"] for r in pk.ligands(s)} - {"HOH", "CA", "SO4", "IOD", "NA", "CL", "GOL", "PO4"}
        assert ours == ref


@pytest.mark.parametrize("entry", ["4HHB", "1BRS", "1AKE", "4MMK"])
def test_oligomer_wysiwyg(entry):
    h = harms()
    out = h.oligomer.pdbOligomer(lines(entry))
    ref = [block.split("\t")[2] == "True" for block in out.split("|") if block.strip()]
    for s in both(entry):
        ours = [a["in_asymmetric_unit"] for a in pk.assemblies(s)]
        assert ours == ref


@pytest.mark.parametrize("entry", ["1AKE", "4MMK", "1BRS"])
def test_exper(entry):
    h = harms()
    out = h.exper.pdbExperiment(lines(entry))
    num = lambda x: None if x.strip() == "na" else float(x)
    res, rval, rfree = num(out[-30:-20]), num(out[-20:-10]), num(out[-10:])
    for s in both(entry):
        e = pk.experiment(s)
        assert float(e["resolution"]) == pytest.approx(res, abs=0.01)
        if rval is not None:                 # upstream cannot parse every REMARK 3 layout
            assert 100 * float(e["r_work"]) == pytest.approx(rval, abs=0.06)
        if rfree is not None:
            assert 100 * float(e["r_free"]) == pytest.approx(rfree, abs=0.06)


@pytest.mark.parametrize("entry", ["4HHB", "1AKE", "4MMK", "1BNA"])
def test_seqres_sequences(entry):
    h = harms()
    fasta = h.seq.pdbSeq2Fasta(lines(entry), chain="all")
    ref = {}
    for block in fasta.split(">")[1:]:
        head, *seq = block.splitlines()
        ref[head.split("_")[0][-1]] = "".join(seq)
    for s in both(entry):
        ours = pk.sequences(s, "seqres")
        if entry == "1BNA":       # Harms maps nucleotides to X; compare lengths
            assert {c: len(v) for c, v in ours.items()} == {c: len(v) for c, v in ref.items()}
        else:
            assert ours == ref


# --------------------------------------------------------------------------- #
# editing
# --------------------------------------------------------------------------- #

def test_bfactor_from_data(tmp_path):
    h = harms()
    data = tmp_path / "vals.txt"
    data.write_text("# residue value\n" + "".join(f"{i} {i * 0.5:.2f}\n" for i in range(1, 47, 2)))
    ref = h.bfactor.pdbBfactor(lines("1CRN"), h.bfactor.loadDataFile(str(data)))
    ref_b = {(ln[22:26].strip(), ln[12:16].strip()): ln[60:66].strip() for ln in ref if ln.startswith("ATOM")}
    for s in both("1CRN"):
        out = pk.bfactor_from_data(s, str(data))
        for r in atom_records(pk.dumps(out, "pdb")):
            v = ref_b[(str(r[6]), r[2])]
            assert (r[12] == 0.0) if v == "NA" else r[12] == pytest.approx(float(v))


@pytest.mark.parametrize("remove_gaps", [False, True])
def test_residue_renumber(remove_gaps):
    h = harms()
    out, _ = h.residue_renumber.pdbResidueRenumber(lines("4HHB"), start_res=101, remove_gaps=remove_gaps)
    ref = sorted((r[5], r[4], r[6], r[2]) for r in atom_records("".join(out)))
    for s in both("4HHB"):
        new = pk.renumber_residues(s, start=101, keep_gaps=not remove_gaps)
        ours = sorted((r[5], r[4], r[6], r[2]) for r in atom_records(pk.dumps(new, "pdb")))
        assert ours == ref


def test_offset_and_atom_renumber(tmp_path):
    h = harms()
    ref = h.offset.pdbOffset(str(DATA / "1CRN.pdb"), 7)
    for s in both("1CRN"):
        ours = pk.renumber_residues(s, shift=7)
        assert sorted(r[2:7] for r in atom_records("".join(ref))) == \
            sorted(r[2:7] for r in atom_records(pk.dumps(ours, "pdb")))
    ref = h.atom_renumber.pdbAtomRenumber(lines("1CRN"))
    serials_ref = [int(ln[6:11]) for ln in ref if ln.startswith(("ATOM", "HETATM"))]
    out = pk.dumps(pk.renumber_atoms(pk.load(DATA / "1CRN.cif"), 1), "pdb")
    assert [int(ln[6:11]) for ln in out.splitlines() if ln.startswith(("ATOM", "HETATM"))] == serials_ref
