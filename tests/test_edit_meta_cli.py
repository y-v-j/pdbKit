"""Editing, metadata parity, FASTA, CLI and fetching."""

import json
import os
import subprocess
import sys

import numpy as np
import pytest

import pdbkit as pk
from pdbkit.cli import main
from conftest import DATA, ENTRIES, ENV, atom_records


def variants(s, num):
    return sorted({r.get_resname() for r in pk.structure.iter_residues_all(s) if r.id[1] == num})


# --------------------------------------------------------------------------- #
# editing
# --------------------------------------------------------------------------- #

def test_edits_reach_microheterogeneity_variants():
    s = pk.load(DATA / "3NIR.cif")
    assert variants(s, 22) == ["PRO", "SER"]
    shifted = pk.loads(pk.dumps(pk.renumber_residues(s, shift=100), "pdb"), "pdb")
    assert variants(shifted, 122) == ["PRO", "SER"] and variants(shifted, 22) == []
    b = pk.loads(pk.dumps(pk.set_bfactor(s, 7.0), "cif"), "cif")
    assert {a.bfactor for a in pk.structure.iter_atoms_all_altlocs(b)} == {7.0}
    ca = pk.extract(s, atoms="ca", het="none", water=False)
    assert {a.get_id() for a in pk.structure.iter_atoms_all_altlocs(ca)} == {"CA"}


def test_rename_chain_swap_and_merge_error():
    s = pk.load(DATA / "4HHB.pdb")
    sw = pk.rename_chains(s, "A:B,B:A")
    assert len(sw[0]["A"]) == len(s[0]["B"]) and len(sw[0]["B"]) == len(s[0]["A"])
    with pytest.raises(pk.PdbKitError, match="already exists and residue"):
        pk.set_chain(s, "X")                    # chains share residue numbers
    assert [c.id for c in pk.set_chain(pk.load(DATA / "1CRN.pdb"), "Z")[0]] == ["Z"]


def test_renumber_modes():
    s = pk.load(DATA / "1AKE.pdb")
    first = lambda st, c: [r.id[1] for r in st[0][c]][:3]
    assert first(pk.renumber_residues(s, start=10), "B") == [10, 11, 12]            # per chain
    kept = pk.renumber_residues(s, start=101, keep_gaps=True)
    assert first(kept, "A") == [101, 102, 103]
    poly_only = pk.renumber_residues(s, start=1000, include_het=False)
    ligs = [r.id[1] for r in poly_only[0]["A"] if pk.residue_kind(r) != "polymer"]
    assert ligs == [r.id[1] for r in s[0]["A"] if pk.residue_kind(r) != "polymer"]


def test_merge_ensemble_intersect(tmp_path):
    a = pk.load(DATA / "1CRN.pdb")
    with pytest.raises(pk.PdbKitError, match="more than one input"):
        pk.merge([a, a])
    m = pk.merge([a, a], rename_duplicates=True)
    assert [c.id for c in m[0]] == ["A", "B"]
    ens = pk.make_ensemble([DATA / "1CRN.pdb", DATA / "1CRN.cif"])
    assert len(ens) == 2 and pk.check_ensemble(ens) == []
    bad = pk.make_ensemble([a, pk.extract(a, residues="1:40")])
    assert pk.check_ensemble(bad)
    assert pk.check_ensemble(pk.load(DATA / "1D3Z.pdb")) == []


def test_bfactor_from_chain_qualified_file(tmp_path):
    f = tmp_path / "v.txt"
    f.write_text("A:5 9.5\nB:5 1.5\n")
    s = pk.bfactor_from_data(pk.load(DATA / "4HHB.pdb"), str(f), default=-1)
    assert s[0]["A"][5]["CA"].bfactor == 9.5 and s[0]["B"][5]["CA"].bfactor == 1.5
    assert s[0]["A"][6]["CA"].bfactor == -1


def test_segid_round_trip():
    s = pk.chain_to_segid(pk.load(DATA / "4HHB.pdb"))
    text = pk.dumps(s, "pdb")
    assert {ln[72:76].strip() for ln in text.splitlines() if ln.startswith("ATOM")} == set("ABCD")
    back = pk.segid_to_chain(pk.set_chain(pk.extract(s, chains="A"), "Q"))
    assert [c.id for c in back[0]] == ["A"]
    assert sorted(pk.split_segments(s)) == list("ABCD")


# --------------------------------------------------------------------------- #
# metadata parity (PDB header adapter vs mmCIF)
# --------------------------------------------------------------------------- #

@pytest.mark.parametrize("entry", ["4HHB", "1AKE", "4MMK", "1BRS", "6LU7", "1BNA", "3NIR"])
def test_metadata_parity(entry):
    a, b = pk.load(DATA / f"{entry}.pdb"), pk.load(DATA / f"{entry}.cif")
    assert pk.sequences(a, "seqres") == pk.sequences(b, "seqres")
    assert pk.sequences(a, "atom") == pk.sequences(b, "atom")
    ea, eb = pk.experiment(a), pk.experiment(b)
    assert ea["method"] == eb["method"] and ea["title"].upper() == eb["title"].upper()
    assert float(ea["resolution"]) == pytest.approx(float(eb["resolution"]), abs=0.006)
    for k in ("r_work", "r_free"):
        if eb[k] != "NA":
            assert float(ea[k]) == pytest.approx(float(eb[k]), abs=0.0006), k
    strip = lambda rows: [(r["id"], r["chains"], r["n_operators"], r["in_asymmetric_unit"]) for r in rows]
    assert strip(pk.assemblies(a)) == strip(pk.assemblies(b))


def test_assembly_is_asymmetric_unit():
    assert [x["is_asymmetric_unit"] for x in pk.assemblies(pk.load(DATA / "4HHB.cif"))] == [True]
    brs = pk.assemblies(pk.load(DATA / "1BRS.cif"))
    assert [x["chains"] for x in brs] == [["A", "D"], ["B", "E"], ["C", "F"]]
    assert all(x["in_asymmetric_unit"] and not x["is_asymmetric_unit"] for x in brs)


# --------------------------------------------------------------------------- #
# FASTA (feature added on request: 4MMK.cif -> 4MMK.fasta)
# --------------------------------------------------------------------------- #

def read_fasta(path):
    seqs, name = {}, None
    for ln in open(path):
        ln = ln.strip()
        if ln.startswith(">"):
            name = ln[1:].split()[0]
            seqs[name] = ""
        elif ln:
            seqs[name] += ln
    return seqs


@pytest.mark.parametrize("ext", ["cif", "pdb"])
def test_fasta_cli_writes_stem_fasta(ext, tmp_path):
    src = tmp_path / f"4MMK.{ext}"
    src.write_text((DATA / f"4MMK.{ext}").read_text())
    assert main(["fasta", str(src)]) == 0
    seqs = read_fasta(tmp_path / "4MMK.fasta")
    assert sorted(seqs) == [f"4MMK_{c}" for c in "ABCDEFGHIJKL"]
    assert all(len(s) == 82 and s.startswith("MSKGHSLADPYLNTLRKERVPVSIYLVNG") for s in seqs.values())
    header = open(tmp_path / "4MMK.fasta").readline()
    assert "mol:protein" in header and "length:82" in header


def test_fasta_nucleic_and_mixed(tmp_path):
    out = tmp_path / "x.fasta"
    pk.write_fasta(DATA / "1BNA.cif", out)
    assert set(read_fasta(out).values()) == {"CGCGAATTCGCG"}
    assert "mol:dna" in out.read_text()
    atom = pk.fasta(pk.load(DATA / "4MMK.pdb"), source="atom")
    assert "source:atom" in atom and "length:69" in atom


def test_fasta_without_seqres_falls_back_to_atoms(tmp_path):
    p = tmp_path / "noseq.pdb"
    pk.save(pk.load(DATA / "1CRN.pdb"), p)          # our writer emits no SEQRES
    text = pk.fasta(pk.load(p))
    assert "source:atom" in text and "TTCCPSIVARSNFNVCRLPGTPEAICATYTGCIIIPGATCPGDYAN" in text.replace("\n", "")


@pytest.mark.network
@pytest.mark.parametrize("entry", ["4MMK", "4HHB", "1BNA", "6LU7"])
def test_fasta_matches_rcsb_fasta_service(entry):
    """Oracle: RCSB's own FASTA (https://www.rcsb.org/fasta/entry/<ID>)."""
    import urllib.request
    text = urllib.request.urlopen(f"https://www.rcsb.org/fasta/entry/{entry}", timeout=30).read().decode()
    rcsb = {}
    for block in text.strip().split(">")[1:]:
        head, *seq = block.splitlines()
        chains = head.split("|")[1].replace("Chains ", "").replace("Chain ", "")
        for c in chains.split(","):
            c = c.strip().split("[auth ")[-1].rstrip("]")
            rcsb[c] = "".join(seq)
    for ext in ("pdb", "cif"):
        ours = pk.sequences(pk.load(DATA / f"{entry}.{ext}"), "seqres")
        for c, s in ours.items():
            assert s.replace("X", "") == rcsb[c].replace("X", "") or len(s) == len(rcsb[c]), (entry, c)


# --------------------------------------------------------------------------- #
# CLI
# --------------------------------------------------------------------------- #

def run_cli(*args, stdin=None):
    res = subprocess.run([sys.executable, "-m", "pdbkit.cli", *map(str, args)], input=stdin,
                         capture_output=True, text=True, env=ENV)
    return res


def test_cli_pipeline_extract_center(tmp_path):
    a = run_cli("extract", DATA / "4HHB.cif", "--chains", "A", "--atoms", "ca", "--het", "HEM", "--no-water", "-f", "pdb")
    assert a.returncode == 0, a.stderr
    b = run_cli("center", "-", "-o", "-", stdin=a.stdout)
    assert b.returncode == 0, b.stderr
    s = pk.loads(b.stdout)
    assert np.allclose(pk.center_of_mass(s), 0, atol=1e-3)
    assert {r.get_resname() for r in s.get_residues() if pk.residue_kind(r) == "ligand"} == {"HEM"}


def test_cli_split_and_center_default_names(tmp_path):
    cwd = os.getcwd()
    os.chdir(tmp_path)
    try:
        assert main(["split", str(DATA / "4HHB.pdb"), "--atoms", "ca", "--het", "all", "-d", "out"]) == 0
        assert sorted(os.listdir("out")) == [f"4HHB_{c}_ca_het.pdb" for c in "ABCD"]
        assert main(["split", str(DATA / "4HHB.cif"), "-c", "B", "--het", "HEM,PO4", "--water", "-d", "o2"]) == 0
        assert os.listdir("o2") == ["4HHB_B_all_het-HEM-PO4_wat.cif"]
        assert main(["center", str(DATA / "1CRN.cif")]) == 0
        assert os.path.exists("1CRN_centered.cif")
        assert main(["center", str(DATA / "4HHB.pdb"), "-c", "A", "--keep-all", "-o", "c.pdb"]) == 0
        assert len(list(pk.load("c.pdb").get_chains())) == 4
    finally:
        os.chdir(cwd)


@pytest.mark.parametrize("cmd", [
    ["info"], ["validate"], ["gaps"], ["ligands"], ["disulfide"], ["contacts", "-t", "HEM"], ["closecontacts"],
    ["watercontacts"], ["neighbors"], ["torsion"], ["param"], ["coulomb"], ["iondist"], ["sasa"], ["exper"],
    ["assembly"], ["com", "-c", "A"], ["distances", "-r", "A:87:NE2", "--atom", "FE"], ["convert", "-f", "cif"],
    ["sort"], ["renumber", "--start", "1"], ["rename", "--chain-map", "A:X"], ["setfield", "--bfactor", "1"],
    ["centerasu"],
])
@pytest.mark.parametrize("ext", ["pdb", "cif"])
def test_cli_commands_run_on_both_formats(cmd, ext, capsys):
    rc = main([cmd[0], str(DATA / f"4HHB.{ext}"), *cmd[1:]])
    out = capsys.readouterr().out
    assert rc == 0 and out.strip()


def test_cli_table_formats(capsys):
    main(["ligands", str(DATA / "4HHB.cif"), "--table", "json"])
    rows = json.loads(capsys.readouterr().out)
    assert {r["resname"] for r in rows} == {"HEM", "PO4"}
    main(["sasa", str(DATA / "1CRN.pdb"), "--table", "csv"])
    assert capsys.readouterr().out.startswith("residue,sasa,relative")


def test_cli_errors_are_clean():
    r = run_cli("extract", DATA / "1CRN.pdb", "--chains", "Q")
    assert r.returncode == 2 and "pdbkit: error: The selection matched no atoms" in r.stderr
    r = run_cli("extract", "no_such_file.pdb")
    assert r.returncode == 2 and "not an existing file" in r.stderr


def test_sasa_sanity():
    rows = pk.sasa(pk.load(DATA / "1CRN.cif"))
    total = sum(r["sasa"] for r in rows)
    assert 2500 < total < 3500                                       # crambin ~2.9-3.0e3 A^2
    assert all(0 <= (r["relative"] or 0) < 1.5 for r in rows)


def test_validate():
    assert pk.validate_file(str(DATA / "4HHB.pdb"))["ok"]
    assert pk.validate_file(str(DATA / "4HHB.cif"))["ok"]
    bad = pk.validate_text("ATOM      1  N   ALA A   1      xx.xxx  12.000  13.000  1.00  0.00           N\n")
    assert not bad["ok"]
    ours = pk.validate_text(pk.dumps(pk.load(DATA / "3NIR.cif"), "cif"))
    assert ours["ok"] and not any("_entity" in w for w in ours["warnings"])


# --------------------------------------------------------------------------- #
# fetching
# --------------------------------------------------------------------------- #

@pytest.mark.network
def test_fetch_rcsb(tmp_path):
    p = pk.fetch("1crn", "pdb", outdir=str(tmp_path))
    assert p.endswith("1CRN.pdb") and pk.load(p).id == "1CRN"
    c = pk.fetch("1BRS", "cif", assembly=1, outdir=str(tmp_path))
    assert sorted(ch.id for ch in pk.load(c).get_chains()) == ["A", "D"]
    p1 = pk.fetch("1BRS", "pdb", assembly=1, outdir=str(tmp_path))
    assert len(list(pk.load(p1).get_chains())) == 2
    with pytest.raises(pk.PdbKitError, match="mmCIF only|--format cif"):
        pk.fetch("4V6X", "pdb", outdir=str(tmp_path))
    with pytest.raises(pk.PdbKitError, match="HTTP 404"):
        pk.fetch("0ZZZ", "cif", outdir=str(tmp_path))


@pytest.mark.network
def test_fetch_alphafold(tmp_path):
    p = pk.fetch("P69905", "pdb", outdir=str(tmp_path), pae=True)
    s = pk.load(p)
    assert os.path.basename(p) == "AF-P69905-F1.pdb"
    assert pk.sequences(s)["A"].startswith("MVLSPADKTNVKAAWGKVGAHAGEYGAEALERMF")
    assert (tmp_path / "AF-P69905-F1_pae.json").exists()
    c = pk.fetch("AF-P69905-F1", "cif", outdir=str(tmp_path))
    assert np.allclose(pk.center_of_mass(pk.load(c)), pk.center_of_mass(s), atol=1e-3)
    with pytest.raises(pk.PdbKitError):
        pk.fetch("AF-P69905-F9", "cif", outdir=str(tmp_path))      # no such fragment


@pytest.mark.network
def test_cli_fetch_to_stdout_and_file(tmp_path):
    r = run_cli("fetch", "1CRN", "-f", "pdb", "-o", "-")
    assert r.returncode == 0 and "ATOM" in r.stdout
    r = run_cli("fetch", "AF-P69905-F1", "-o", tmp_path / "hba.cif")
    assert r.returncode == 0 and (tmp_path / "hba.cif").exists()


# --------------------------------------------------------------------------- #
# regressions found by the case-study run (testcases/)
# --------------------------------------------------------------------------- #

def test_casestudy_close_contacts_ca_only_on_nucleic_acids():
    """--ca-only on DNA/RNA used to crash inside NeighborSearch (empty atom list)."""
    for entry in ("1BNA", "1EHZ"):
        assert pk.close_contacts(pk.load(DATA / f"{entry}.cif"), 3.5, ca_only=True) == []
    assert main(["closecontacts", str(DATA / "1EHZ.cif"), "--ca-only"]) == 0


@pytest.mark.parametrize("ext", ["pdb", "cif"])
def test_casestudy_long_deposited_disulfides(ext):
    """1K4C annotates SG-SG bonds of 2.70-2.72 A; the 3.0 A default finds all three."""
    s = pk.load(DATA / f"1K4C.{ext}")
    bonds = {(p["cys1"], p["cys2"]) for p in pk.disulfides(s)["bonds"]}
    assert bonds == {("A:CYS:22", "A:CYS:96"), ("B:CYS:23", "B:CYS:88"), ("B:CYS:134", "B:CYS:194")}


def test_casestudy_fasta_modified_nucleotides():
    """1EHZ (tRNA-Phe) has PSU, 2MG, M2G, OMC, OMG, 5MC, 7MG, 5MU, 1MA, H2U, YYG:
    SEQRES output equals the file's canonical sequence; modelled residues map to parents."""
    s = pk.load(DATA / "1EHZ.cif")
    can = "".join(pk.metadata.first(pk.metadata.as_mmcif_dict(s), "_entity_poly.pdbx_seq_one_letter_code_can").split())
    assert len(can) == 76 and "X" not in can
    p = pk.load(DATA / "1EHZ.pdb")
    for st in (s, p):
        for source in ("seqres", "atom"):
            assert pk.sequences(st, source)["A"] == can, source


def test_casestudy_pdb_ter_records_and_80_columns():
    """TER follows the last polymer residue of each chain (not the waters); records are 80 columns."""
    text = pk.dumps(pk.load(DATA / "4HHB.cif"), "pdb")
    lines = text.splitlines()
    assert all(len(ln) == 80 for ln in lines)
    ters = [i for i, ln in enumerate(lines) if ln.startswith("TER")]
    assert len(ters) == 4
    for i in ters:
        prev = lines[i - 1]
        assert prev.startswith("ATOM") and lines[i][17:27] == prev[17:27]
    serials = [int(ln[6:11]) for ln in lines if ln[:6] in ("ATOM  ", "HETATM", "TER   ")]
    assert serials == list(range(1, len(serials) + 1))


def test_casestudy_centerasu_refuses_non_crystal():
    with pytest.raises(pk.PdbKitError, match="placeholder cell"):
        pk.center_asu(pk.load(DATA / "1D3Z.cif"))


def test_casestudy_split_copies_only_selected_chains():
    """Chain splitting must not deep-copy the whole structure per chain (was quadratic)."""
    import time
    s = pk.load(DATA / "4MMK.cif")
    t = time.time()
    parts = pk.split_chains(s)
    assert len(parts) == 12 and time.time() - t < 5
    assert sum(len(list(pk.structure.iter_atoms_all_altlocs(p))) for p in parts.values()) == \
        len(list(pk.structure.iter_atoms_all_altlocs(s)))


# regressions found by the AlphaFold case study (testcases/alphafold/)

@pytest.mark.network
def test_casestudy_fetch_output_across_filesystems():
    """`fetch -o FILE` downloaded into /tmp and os.replace()d the file, which fails when the
    destination is on another filesystem (here: the home directory vs tmpfs /tmp)."""
    import tempfile
    dest_dir = tempfile.mkdtemp(dir=os.path.dirname(os.path.abspath(__file__)))
    try:
        r = run_cli("fetch", "AF-P69905-F1", "-o", os.path.join(dest_dir, "hba.cif"))
        assert r.returncode == 0, r.stderr
        assert pk.sequences(pk.load(os.path.join(dest_dir, "hba.cif")))["A"].startswith("MVLSPADKTN")
    finally:
        import shutil
        shutil.rmtree(dest_dir)


def test_casestudy_hetnavi_on_predicted_model_is_a_clean_error(tmp_path):
    s = pk.load(DATA / "1CRN.cif")
    pk.source_info(s).entry_id = "AF-P69905-F1"
    p = tmp_path / "AF-P69905-F1.cif"
    p.write_text(pk.dumps(s, "cif").replace("data_1CRN", "data_AF-P69905-F1").replace("_entry.id 1CRN", "_entry.id AF-P69905-F1"))
    r = run_cli("hetnavi", p)
    assert r.returncode == 2 and "Het-PDB Navi covers PDB entries only" in r.stderr and "Traceback" not in r.stderr


@pytest.mark.network
@pytest.mark.parametrize("ident", ["P00720", "Q8WZ42", "AF-P69905-F2"])
def test_casestudy_alphafold_missing_models_explained(ident, tmp_path):
    """Viral (T4 lysozyme), >2,700-residue (titin) and non-existent fragments give a clear reason."""
    with pytest.raises(pk.PdbKitError, match="viral proteins"):
        pk.fetch(ident, "cif", outdir=str(tmp_path))
