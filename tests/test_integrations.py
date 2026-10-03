"""Het-PDB Navi parser and domain/family annotation vs the original scripts."""

import csv
import filecmp
import os
import subprocess
import sys

import pytest

import pdbkit as pk
from conftest import DATA, ENV, UPSTREAM


def run_script(name, *args, cwd=None):
    res = subprocess.run([sys.executable, str(UPSTREAM / name), *args], capture_output=True, text=True,
                         env=ENV, cwd=cwd)
    assert res.returncode == 0, res.stderr
    return res.stdout


def read_csv(path):
    with open(path, newline="") as fh:
        return list(csv.DictReader(fh))


# --------------------------------------------------------------------------- #
# Het-PDB Navi
# --------------------------------------------------------------------------- #

@pytest.mark.oracle
@pytest.mark.parametrize("entry", ["4HHB", "1AKE", "1M17"])
def test_hetnavi_outputs_match_original(entry, tmp_path):
    """Every file the original writes is byte-identical in pdbKit's output. For data
    whose #HET lines lack the H-bond column (1AKE, 1M17) the original writes empty
    tables; pdbKit additionally keeps those contacts (one CSV row per #HET line)."""
    dat = str(DATA / f"{entry}_assembly1.dat")
    run_script("hetpdbnavi_parser.py", entry, "--file", dat, "--outdir", str(tmp_path / "orig"))
    records, written = pk.hetnavi.interactions(entry, dat_file=dat, outdir=str(tmp_path / "ours"))
    n_het = sum(1 for ln in open(dat) if ln.startswith("FT\t#HET"))
    assert len(read_csv(tmp_path / "ours" / f"{entry}_all_het_interactions.csv")) == n_het
    full_layout = all(len(ln.rstrip("\n").split("\t")) == 13 for ln in open(dat) if ln.startswith("FT\t#HET"))
    orig = sorted(os.listdir(tmp_path / "orig"))
    assert set(orig) <= {os.path.basename(p) for p in written}
    for name in orig:
        same = filecmp.cmp(tmp_path / "orig" / name, tmp_path / "ours" / name, shallow=False)
        if full_layout or "het" not in name.lower():
            assert same, name


@pytest.mark.parametrize("ext", ["pdb", "cif"])
def test_hetnavi_takes_the_id_from_a_structure_file(ext, tmp_path):
    dat = str(DATA / "4HHB_assembly1.dat")
    records, written = pk.hetnavi.interactions(DATA / f"4HHB.{ext}", dat_file=dat, outdir=str(tmp_path))
    assert any(p.endswith("4HHB_HEM_142_A_interaction.csv") for p in written)
    rows = pk.hetnavi.het_contacts_table(records)
    assert {r["heteroatom"] for r in rows} >= {"HEM:142:A", "HEM:148:B"}


@pytest.mark.network
def test_hetnavi_download(tmp_path):
    records, written = pk.hetnavi.interactions("4HHB", outdir=str(tmp_path), save_raw=True)
    assert records and any(p.endswith(".dat") for p in written)


# --------------------------------------------------------------------------- #
# Domain / family annotation
# --------------------------------------------------------------------------- #

DOMAIN_ENTRIES = ["4HHB", "1AKE", "4MMK", "6LU7", "1WHZ", "3NIR", "1BNA", "1CRN"]


@pytest.mark.oracle
@pytest.mark.parametrize("entry", DOMAIN_ENTRIES)
def test_domains_mmcif_identical_to_original_offline(entry, tmp_path):
    run_script("mmcif_domain_family.py", str(DATA / f"{entry}.cif"), "-o", str(tmp_path / "orig"), "-q")
    ours = pk.domains.run([str(DATA / f"{entry}.cif")], str(tmp_path / "ours"), skip_lookups=True, log=lambda m: None)
    assert filecmp.cmp(tmp_path / "orig" / f"{entry}_domain_ss_info.csv", ours[0], shallow=False)


PARITY_COLUMNS = ["Chain_ID", "EC_Number", "EC_Type", "Enzyme Type", "UniProt_Accession",
                  "Secondary_Structure_Order", "Secondary_Structure_Residue_Ranges", "Heteroatoms_In_Chain"]


@pytest.mark.parametrize("entry", DOMAIN_ENTRIES)
def test_domains_pdb_input_matches_mmcif(entry, tmp_path):
    """The original reported N/A for most columns of .pdb input; pdbKit reads the PDB
    header records and gets the same values as from the mmCIF file."""
    a = pk.domains.annotate(DATA / f"{entry}.pdb", skip_lookups=True)
    b = pk.domains.annotate(DATA / f"{entry}.cif", skip_lookups=True)
    assert len(a) == len(b)
    for ra, rb in zip(a, b):
        for col in PARITY_COLUMNS:
            assert ra[col] == rb[col], (entry, ra["Chain_ID"], col, ra[col], rb[col])
        assert ra["Title"].upper() == rb["Title"].upper()
        assert float(ra["Resolution"]) == pytest.approx(float(rb["Resolution"]), abs=0.006)  # PDB keeps 2 decimals
        assert ra["Experimental_Method"] == rb["Experimental_Method"]


def test_domains_original_pdb_input_reports_na(tmp_path):
    run_script("mmcif_domain_family.py", str(DATA / "4HHB.pdb"), "-o", str(tmp_path), "-q")
    rows = read_csv(tmp_path / "4HHB_domain_ss_info.csv")
    assert all(r["Secondary_Structure_Order"] == "N/A" and r["UniProt_Accession"] == "N/A" for r in rows)
    ours = pk.domains.annotate(DATA / "4HHB.pdb", skip_lookups=True)
    assert ours[0]["UniProt_Accession"] == "P69905" and ours[0]["Secondary_Structure_Order"].startswith("Coil-Helix")


def test_domains_bonding_partners_pdb_vs_cif():
    """Covalent/metal/disulfide partners agree; mmCIF-only hydrogen-bond records are the only extra."""
    import re
    for entry in ("4HHB", "1AKE", "6LU7", "1BNA", "4MMK"):
        a = pk.domains.annotate(DATA / f"{entry}.pdb", skip_lookups=True)
        b = pk.domains.annotate(DATA / f"{entry}.cif", skip_lookups=True)
        for ra, rb in zip(a, b):
            strip_h = lambda s: re.sub(r"[^,>]*\(chain \w+, hydrog\),? ?", "", s)
            assert ra["Heteroatom_Bonding_Partners"] == strip_h(rb["Heteroatom_Bonding_Partners"]).replace(", ;", ";"), \
                (entry, ra["Chain_ID"])


@pytest.mark.network
def test_domains_with_live_lookups_identical_to_original(tmp_path):
    run_script("mmcif_domain_family.py", str(DATA / "1AKE.cif"), "-o", str(tmp_path / "orig"))
    ours = pk.domains.run([str(DATA / "1AKE.cif")], str(tmp_path / "ours"), log=lambda m: None)
    assert filecmp.cmp(tmp_path / "orig" / "1AKE_domain_ss_info.csv", ours[0], shallow=False)
    rows = read_csv(ours[0])
    assert "Pfam" in rows[0]["Family_Domain_DB_Xrefs"] and rows[0]["EC_Number"] == "2.7.4.3"


@pytest.mark.oracle
def test_hetnavi_lines_without_hbond_count(tmp_path):
    """Case-study finding: 1M17 (EGFR + erlotinib) lists its 15 AQ4 contacts on #HET lines
    without the H-bond column. The original parser drops them all; pdbKit keeps them."""
    dat = DATA / "1M17_assembly1.dat"
    n_lines = sum(1 for ln in open(dat) if ln.startswith("FT\t#HET"))
    run_script("hetpdbnavi_parser.py", "1M17", "--file", str(dat), "--outdir", str(tmp_path / "orig"))
    assert len(read_csv(tmp_path / "orig" / "1M17_all_het_interactions.csv")) == 0          # upstream: nothing
    records, written = pk.hetnavi.interactions("1M17", dat_file=str(dat), outdir=str(tmp_path / "ours"))
    rows = read_csv(tmp_path / "ours" / "1M17_all_het_interactions.csv")
    assert len(rows) == n_lines == 15
    assert {r["heteroatom"] for r in rows} == {"AQ4:999:A"}
    assert (tmp_path / "ours" / "1M17_AQ4_999_A_interaction.csv").exists()
    met769 = [r for r in rows if r["res_end"] == "769"][0]           # hinge Met769, the key erlotinib contact
    assert met769["res_name"] == "MET" and met769["num_contacts"] == "13" and met769["num_hbonds"] == ""
