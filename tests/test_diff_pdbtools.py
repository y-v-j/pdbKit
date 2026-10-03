"""Differential tests: pdbKit vs the original HADDOCK pdb-tools 2.0.0-rc1 scripts.

Each case runs the upstream script on the PDB file and pdbKit on *both* the
PDB and the mmCIF file of the same entry; the atom records (fixed-column
parse, serial numbers ignored) must be identical as multisets.
"""

import pytest

import pdbkit as pk
from pdbkit.select import Selection
from conftest import DATA, atom_records, core, run_pdbtool

pytestmark = pytest.mark.oracle


def same(upstream_text, structure, with_models=False):
    ours = pk.dumps(structure, "pdb")
    a = sorted(core(atom_records(upstream_text, with_models), True))
    b = sorted(core(atom_records(ours, with_models), True))
    assert a == b, f"{len(a)} upstream atoms vs {len(b)} pdbKit atoms"


def both(entry):
    return [pk.load(DATA / f"{entry}.pdb"), pk.load(DATA / f"{entry}.cif")]


SELECTIONS = [
    # (tool, args, entry, pdbKit selection)
    ("pdb_selchain", ["-A,C"], "4HHB", dict(chains="A,C")),
    ("pdb_delchain", ["-A"], "4HHB", dict(exclude_chains="A")),
    ("pdb_selresname", ["-HEM,HOH"], "4HHB", dict(resnames="HEM,HOH")),
    ("pdb_delresname", ["-HOH"], "4HHB", dict(exclude_resnames="HOH")),
    ("pdb_selelem", ["-C,N"], "4HHB", dict(elements="C,N")),
    ("pdb_delelem", ["-C"], "4HHB", dict(exclude_elements="C")),
    ("pdb_delhetatm", [], "4HHB", dict(only_het_records=False)),
    ("pdb_delhetatm", [], "1WHZ", dict(only_het_records=False)),     # also drops MSE, like upstream
    # pdb_selatom matches names everywhere: water "O" atoms are kept too
    ("pdb_selatom", ["-CA,C,N,O"], "4HHB", dict(atoms="CA,C,N,O", het="none", water=True)),
    ("pdb_selatom", ["-CA"], "1CRN", dict(atoms="ca")),
    ("pdb_selres", ["-10:20"], "1CRN", dict(residues="10:20")),
    ("pdb_selres", ["-::5"], "1CRN", dict(residues="::5")),
    ("pdb_selres", ["-1:30:4"], "1CRN", dict(residues="1:30:4")),
    ("pdb_selres", ["-40:"], "1CRN", dict(residues="40:")),
    ("pdb_delres", ["-1:10"], "1CRN", dict(exclude_residues="1:10")),
    ("pdb_selaltloc", [], "1AKE", dict(altloc="best")),
]


def test_selaltloc_on_file_with_anisou_and_microheterogeneity():
    """3NIR has ANISOU records and two microheterogeneity sites (PRO/SER22, LEU/ILE25).

    Upstream pdb_selaltloc crashes on the ANISOU records (it reads U values as
    occupancy); after stripping them with pdb_keepcoord it keeps *both*
    residue identities at 22 and 25 as overlapping atoms, and with -A/-B it
    crashes on residues whose atoms carry only the other label. pdbKit keeps
    one variant. Everything outside residues 22/25 must agree exactly.
    """
    with pytest.raises(RuntimeError, match="could not convert"):
        run_pdbtool("pdb_selaltloc", str(DATA / "3NIR.pdb"))
    stripped = run_pdbtool("pdb_keepcoord", str(DATA / "3NIR.pdb"))
    for opt in ("-A", "-B"):
        with pytest.raises(RuntimeError, match="index out of range"):
            run_pdbtool("pdb_selaltloc", opt, stdin=stripped)
    up = atom_records(run_pdbtool("pdb_selaltloc", stdin=stripped))
    assert {r[4] for r in up if r[6] == 22} == {"PRO", "SER"}        # both identities kept upstream
    for s in both("3NIR"):
        ours = atom_records(pk.dumps(pk.extract(s, altloc="best"), "pdb"))
        outside = lambda recs: sorted(core([r for r in recs if r[6] not in (22, 25)], True))
        assert outside(up) == outside(ours)
        assert {r[4] for r in ours if r[6] == 22} == {"PRO"}           # occupancy 0.63 > 0.37
        assert all(r[3] == "" for r in ours)                          # no altloc labels left
        for label in ("A", "B"):
            sel = atom_records(pk.dumps(pk.extract(s, altloc=label), "pdb"))
            assert len({r[4] for r in sel if r[6] == 22}) == 1


@pytest.mark.parametrize("tool,args,entry,sel", SELECTIONS, ids=[f"{t}{''.join(a)}-{e}" for t, a, e, _ in SELECTIONS])
def test_selection_tools(tool, args, entry, sel):
    up = run_pdbtool(tool, *args, str(DATA / f"{entry}.pdb"))
    for s in both(entry):
        same(up, pk.extract(s, Selection(**sel)))


EDITS = [
    ("pdb_shiftres", ["-10"], "4HHB", lambda s: pk.renumber_residues(s, shift=10)),
    ("pdb_shiftres", ["--5"], "1CRN", lambda s: pk.renumber_residues(s, shift=-5)),
    ("pdb_reres", ["-5"], "1CRN", lambda s: pk.renumber_residues(s, start=5, continuous=True)),
    ("pdb_reres", ["-1"], "2GB1", lambda s: pk.renumber_residues(s, start=1, continuous=True)),
    ("pdb_rplchain", ["-A:Z"], "4HHB", lambda s: pk.rename_chains(s, "A:Z")),
    ("pdb_rplresname", ["-HEM:HEE"], "4HHB", lambda s: pk.rename_residues(s, "HEM:HEE")),
    ("pdb_chain", ["-X"], "1CRN", lambda s: pk.set_chain(s, "X")),
    ("pdb_b", ["-5.5"], "4HHB", lambda s: pk.set_bfactor(s, 5.5)),
    ("pdb_occ", ["-0.5"], "4HHB", lambda s: pk.set_occupancy(s, 0.5)),
]


@pytest.mark.parametrize("tool,args,entry,fn", EDITS, ids=[f"{t}{''.join(a)}-{e}" for t, a, e, _ in EDITS])
def test_edit_tools(tool, args, entry, fn):
    up = run_pdbtool(tool, *args, str(DATA / f"{entry}.pdb"))
    for s in both(entry):
        same(up, fn(s))


def test_splitchain(tmp_path):
    import os
    import subprocess
    import sys
    from conftest import PDBTOOLS, ENV
    subprocess.run([sys.executable, str(PDBTOOLS / "pdb_splitchain.py"), str(DATA / "4HHB.pdb")],
                   cwd=tmp_path, check=True, env=ENV)
    for s in both("4HHB"):
        parts = pk.split_chains(s)                   # defaults keep everything
        for cid, part in parts.items():
            same((tmp_path / f"4HHB_{cid}.pdb").read_text(), part)
        assert sorted(parts) == sorted(p.split("_")[1][0] for p in os.listdir(tmp_path))


def test_splitmodel(tmp_path):
    import subprocess
    import sys
    from conftest import PDBTOOLS, ENV
    subprocess.run([sys.executable, str(PDBTOOLS / "pdb_splitmodel.py"), str(DATA / "1D3Z.pdb")],
                   cwd=tmp_path, check=True, env=ENV)
    for s in both("1D3Z"):
        parts = pk.split_models(s)
        assert len(parts) == 10
        for n, part in parts.items():
            same((tmp_path / f"1D3Z_{n}.pdb").read_text(), part)


@pytest.mark.parametrize("entry", ["1CRN", "2GB1", "1D3Z", "1L2Y"])
def test_tofasta(entry):
    up = run_pdbtool("pdb_tofasta", "-multi", str(DATA / f"{entry}.pdb"))
    up_seqs = {}
    for block in up.strip().split(">PDB|")[1:]:
        lines = block.splitlines()
        up_seqs[lines[0]] = "".join(lines[1:])
    for s in both(entry):
        assert pk.sequences(s, "atom") == up_seqs


@pytest.mark.parametrize("entry", ["1CRN", "4HHB", "1D3Z", "3NIR"])
def test_wc_counts(entry):
    up = run_pdbtool("pdb_wc", str(DATA / f"{entry}.pdb"))
    vals = dict(line.split(":", 1) for line in up.strip().splitlines())
    for s in both(entry):
        info = pk.summary(s)
        assert int(vals["No. models"].strip()) == info["models"]
        assert int(vals["No. HETATM"].strip()) == \
            sum(1 for m in s for a in pk.structure.iter_atoms_all_altlocs(m) if a.get_parent().id[0] != " ")


def test_reres_upstream_numbers_in_record_order():
    """pdb_reres follows raw record order: 4HHB's HETATM block (after all chains)
    is numbered after chain D. pdbKit numbers per chain; polymer residues agree."""
    up = atom_records(run_pdbtool("pdb_reres", "-1", str(DATA / "4HHB.pdb")))
    ours = atom_records(pk.dumps(pk.renumber_residues(pk.load(DATA / "4HHB.pdb"), start=1, continuous=True,
                                                      include_het=False), "pdb"))
    poly = lambda recs: sorted(r[2:] for r in recs if r[1] == "ATOM")
    assert poly(up) == poly(ours)
    hem_up = {r[6] for r in up if r[4] == "HEM"}
    assert min(hem_up) > max(r[6] for r in up if r[1] == "ATOM")   # after every polymer residue


def test_gap_upstream_reports_altlocs_as_gaps():
    up = run_pdbtool("pdb_gap", str(DATA / "3NIR.pdb"))
    assert "A:THR1 < Seq. Gap > A:THR1" in up           # same residue, two altlocs
    assert len([ln for ln in up.splitlines() if "<" in ln]) == 18
    for s in both("3NIR"):
        assert pk.gaps(s) == []                          # no real break: every C-N bond is intact


def test_gap_upstream_false_positives_on_stretched_ca_ca():
    """4HHB (1984): pdb_gap reports 36 'gaps' with CA-CA 4.00-4.26 A across intact
    peptide bonds (C-N < 2 A). pdbKit uses the peptide bond and reports none."""
    import math
    up = [ln for ln in run_pdbtool("pdb_gap", str(DATA / "4HHB.pdb")).splitlines() if "<" in ln]
    assert len(up) == 36
    d = [math.sqrt(float(ln.split("<")[1].split("A")[0])) for ln in up]
    assert 3.99 <= min(d) and max(d) < 4.3
    for s in both("4HHB"):
        assert pk.gaps(s) == []


@pytest.mark.parametrize("entry", ["1AKE", "1BRS", "6LU7", "4MMK", "1CRN"])
def test_gap_detects_the_same_breaks_and_upstream_prints_squared_distances(entry):
    up = run_pdbtool("pdb_gap", str(DATA / f"{entry}.pdb"))
    up_lines = [ln for ln in up.splitlines() if "<" in ln]
    for s in both(entry):
        ours = pk.gaps(s)
        assert len(ours) == len(up_lines)
        for ln, g in zip(up_lines, ours):
            if g["type"] == "distance":
                printed = float(ln.split("<")[1].split("A")[0])
                assert printed == pytest.approx(g["distance"] ** 2, abs=0.05), "upstream prints d^2"


def test_merge_intersect_mkensemble(tmp_path):
    a = tmp_path / "a.pdb"
    b = tmp_path / "b.pdb"
    s = pk.load(DATA / "1CRN.pdb")
    pk.save(pk.extract(s, residues="1:20"), a)
    pk.save(pk.rename_chains(pk.extract(s, residues="15:46"), "A:B"), b)
    same(run_pdbtool("pdb_merge", str(a), str(b)), pk.merge([a, b]))
    same(run_pdbtool("pdb_intersect", str(DATA / "1CRN.pdb"), str(a)), pk.intersect([DATA / "1CRN.pdb", a]))
    up = run_pdbtool("pdb_mkensemble", str(DATA / "1CRN.pdb"), str(DATA / "1CRN.pdb"))
    same(up, pk.make_ensemble([DATA / "1CRN.pdb", DATA / "1CRN.cif"]), with_models=True)
