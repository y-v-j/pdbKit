---
title: pdbKit User Manual
version: 1.0.0
created: 2026-10-01
tags: [pdbkit, structural-biology, pdb, mmcif, alphafold, manual]
aliases: [pdbKit manual, pdbkit]
---

# pdbKit User Manual

pdbKit is a single Python package (and `pdbkit` command) for three-dimensional macromolecular structure files. It works the same way on **PDB** (`.pdb`, `.ent`, `.pdbN`) and **PDBx/mmCIF** (`.cif`, `.mmcif`) files, gzipped or not. It merges:

| Source | Licence | What came from it |
|---|---|---|
| HADDOCK **pdb-tools** 2.0.0-rc1 | Apache-2.0 | selection, deletion, renumbering, renaming, splitting, merging, gaps, FASTA, validation, fetch |
| Harms-lab **pdbtools** | GPL-3.0 | centre of mass, ASU centring, contacts, disulfides, neighbours, torsions, pI/MW, Coulomb energy, ion distances, ligands, experiment/oligomer reports, B-factor loading |
| `hetpdbnavi_parser.py` | project | Het-PDB Navi protein–ligand / subunit interaction tables |
| `mmcif_domain_family.py` | project | domain, family, EC, secondary-structure and bonding annotation |

New in pdbKit: download from RCSB **and AlphaFold DB**, chain splitting and selection in four representations, centre-of-mass centring, FASTA for all chains, a standards-conforming mmCIF writer, a PDB-header → mmCIF translator (so every header-based feature works for both formats), and a test suite that checks everything against independent oracles (chapter [[#11. Testing and results]]).

> [!info] Contents
> 1. [[#1. Installation]]
> 2. [[#2. Concepts you need once]]
> 3. [[#3. Quick start]]
> 4. [[#4. Fetching structures (RCSB, AlphaFold)]]
> 5. [[#5. Splitting chains and choosing a representation]]
> 6. [[#6. Selecting chains and atoms]]
> 7. [[#7. Centre of mass and centring]]
> 8. [[#8. FASTA sequences for all chains]]
> 9. [[#9. Het-PDB Navi interactions]]
> 10. [[#10. Domain, family and secondary-structure annotation]]
> 11. [[#11. Testing and results]]
> 12. [[#12. Editing tools (pdb-tools equivalents)]]
> 13. [[#13. Analysis tools (Harms-lab equivalents)]]
> 14. [[#14. Using pdbKit from Python and Jupyter]]
> 15. [[#15. Viewer compatibility (PyMOL, Mol*, ChimeraX)]]
> 16. [[#16. Differences from the original tools]]
> 17. [[#17. Troubleshooting]]
> 18. [[#18. Reference: all commands]]
> 19. [[#19. Worked terminal examples (case-study structures)]]

---

## 1. Installation

pdbKit needs **Python ≥ 3.10**, **Biopython ≥ 1.80** and **NumPy**. pandas is optional (for DataFrames).

```bash
conda activate biopyc                      # the development environment
cd /path/to/github-new/pdbKit
pip install -e .                           # editable install; adds the `pdbkit` command
pdbkit --version                           # pdbKit 1.0.0
```

In any other environment:

```bash
pip install /path/to/pdbKit                # or: pip install -e /path/to/pdbKit
pip install "pdbKit[dataframes,test]"      # + pandas and pytest (from the package folder)
```

Uninstall with `pip uninstall pdbKit`.

> [!warning] A warning you may see in `biopyc`
> Every Biopython import in `biopyc` prints *"BiopythonWarning: You may be importing Biopython from inside the source tree … you have a pyproject.toml file in an unexpected directory: …/site-packages"*. This does not come from pdbKit. The `dash-bootstrap-components` package left its `pyproject.toml` directly in `site-packages`. It is harmless; to remove it, delete (or rename) `…/envs/biopyc/lib/python3.14/site-packages/pyproject.toml`, or reinstall `dash-bootstrap-components`. To silence it for one shell: `export PYTHONWARNINGS="ignore:You may be importing Biopython"`.

**Download cache.** When you pass a PDB or AlphaFold ID where a file is expected (for example `pdbkit info 4HHB`), the file is downloaded once into `~/.cache/pdbkit`. Set `PDBKIT_CACHE=/some/dir` to change this.

---

## 2. Concepts you need once

### 2.1 Inputs and outputs

Every structure command accepts:

- a file: `.pdb .ent .pdb1 … .cif .mmcif`, optionally `.gz`; files without an extension are sniffed;
- a **PDB ID** (`4HHB`) or **AlphaFold ID** (`AF-P69905-F1`, or a UniProt accession `P69905`), fetched and cached;
- `-` for **stdin**, so commands can be piped.

The output format is decided in this order: `-f/--format pdb|cif` → the extension of `-o` → **the input format**. Without `-o`, structures go to stdout (except `center`, `split` and `fasta`, which write files).

```bash
pdbkit fetch 1brs -f pdb -o - | pdbkit extract - -c A,D --no-water | pdbkit center - -o 1brs_AD_centred.pdb
```

### 2.2 How residues are classified

Every residue is classified once, when the file is read. This is what "C-alpha", "all-atom" and "heteroatoms" mean everywhere in pdbKit.

| Class | Rule |
|---|---|
| **polymer** | mmCIF: the file gives the residue a `label_seq_id`. PDB: standard `ATOM` residues, plus `HETATM` residues that are **covalently chained** to a neighbour (previous residue's C bonded to them, or bonded to the next residue's N, or an O3′–P link). Modified residues such as selenomethionine `MSE` and capped/modified peptide residues are therefore *polymer*. |
| **ligand** ("heteroatoms") | every other non-water residue: cofactors, ions, sugars, inhibitors. |
| **water** | `HOH WAT DOD H2O TIP TIP3 SOL …` |

> [!tip] PDB and mmCIF agree
> On all 14 test entries and 28 case-study entries, classification from the PDB file (heuristic) is identical, residue by residue, to the mmCIF file's own `label_seq_id` statement ([[#11. Testing and results]]).

### 2.3 Representations and selection options

These options are shared by `split`, `extract`, `com` and `center`:

| Option | Values | Meaning |
|---|---|---|
| `-a/--atoms` | `all` · `ca` · `backbone` · `trace` · `N,CA,C,O…` | atoms kept from **polymer** residues. `trace` = CA for protein, P for nucleic acids. |
| `--het` | `all` · `none` · `HEM,ZN,…` | ligands kept (whole residues). |
| `--water` / `--no-water` | | keep or drop water. |
| `-c/--chains`, `--exclude-chains` | `A` · `A,C` | chains. |
| `--residues`, `--exclude-residues` | `1:10` `5:` `:20` `1:100:5` `42` `30-40` | author residue numbers; `start:end:step` as in `pdb_selres`. |
| `--resnames`, `--exclude-resnames` | `ALA,GLY` | residue names. |
| `--elements`, `--exclude-elements`, `--no-hydrogens` | `C,N` | elements. |
| `--models` | `1` · `1,3` | NMR models (1-based). |
| `--altloc` | `best` · `A` | resolve alternate locations (`best` = highest occupancy); also picks one residue identity at microheterogeneity sites. |
| `--hetatm-only` / `--no-hetatm` | | by record type, as `pdb_delhetatm`. |
| `--segids` | | PDB segment IDs. |

> [!note] Defaults differ on purpose
> `split` defaults to **`--het none --no-water`** (a clean chain). `extract`, `com` and `center` keep **everything** unless you restrict it.

---

## 3. Quick start

```bash
pdbkit fetch 4HHB                             # -> ./4HHB.cif
pdbkit info 4HHB.cif                          # counts, chains, ligands, method, resolution
pdbkit split 4HHB.cif --atoms ca --het all -d chains
pdbkit center 4HHB.cif                        # -> 4HHB_centered.cif
pdbkit fasta 4HHB.cif                         # -> 4HHB.fasta
pdbkit hetnavi 4HHB -d hetnavi_4HHB           # ligand-contact CSVs
pdbkit domains 4HHB.cif -d annotation         # 4HHB_domain_ss_info.csv
```

More real terminal sessions, on 28 structures from ribosomes to tRNA, are in [[#19. Worked terminal examples (case-study structures)]].

```python
import pdbkit as pk
s = pk.load("4HHB")                            # fetch + parse
chains = pk.split_chains(s, atoms="ca", het="all")
pk.save(chains["A"], "4HHB_A_ca_het.pdb")
```

---

## 4. Fetching structures (RCSB, AlphaFold)

```bash
pdbkit fetch 4HHB                       # mmCIF (default)          -> 4HHB.cif
pdbkit fetch 4HHB 1CRN -f pdb -d pdbs   # several, PDB format      -> pdbs/4HHB.pdb, pdbs/1CRN.pdb
pdbkit fetch 1BRS --assembly 1          # biological assembly 1    -> 1BRS-assembly1.cif
pdbkit fetch 1BRS --assembly 1 -f pdb   #                          -> 1BRS-assembly1.pdb
pdbkit fetch P69905                     # AlphaFold model by UniProt accession -> AF-P69905-F1.cif
pdbkit fetch AF-P69905-F1 -f pdb --pae  # model + predicted aligned error JSON
pdbkit fetch 1CRN -f pdb -o -           # to stdout, for pipes
pdbkit fetch 1CRN -o crambin.cif        # choose the file name
```

- IDs are recognised automatically. Use `--source pdb|alphafold` if needed.
- AlphaFold file URLs come from the AFDB API, so the **current model version** is used (v6 at the time of writing). Fragments of long proteins: `AF-Q8WZ42-F2`.
- Large entries (ribosomes, multi-character chain IDs) have **no PDB-format file** at RCSB. pdbKit then says so and suggests `-f cif`.
- Existing files are not downloaded again unless `--overwrite` is given.

Python: `pk.fetch("4HHB", "cif", outdir=".")`, `pk.fetch("P69905", "pdb", pae=True)`, `pk.fetch_pdb("1BRS", "cif", assembly=1)`.

---

## 5. Splitting chains and choosing a representation

`pdbkit split` writes **one file per chain**. The representation is chosen with `--atoms`, `--het` and `--water`. File names are `<stem>_<chain>_<tag>.<ext>`; the tag records the representation (disable it with `--no-tag`).

| You want | Command | File for chain A |
|---|---|---|
| C-alpha only | `pdbkit split 4HHB.cif --atoms ca` | `4HHB_A_ca.cif` |
| all-atom | `pdbkit split 4HHB.cif` | `4HHB_A_all.cif` |
| C-alpha + all heteroatoms | `pdbkit split 4HHB.cif --atoms ca --het all` | `4HHB_A_ca_het.cif` |
| C-alpha + specific heteroatoms | `pdbkit split 4HHB.cif --atoms ca --het HEM` | `4HHB_A_ca_het-HEM.cif` |
| all-atom + all heteroatoms | `pdbkit split 4HHB.cif --het all` | `4HHB_A_all_het.cif` |
| all-atom + specific heteroatoms + water | `pdbkit split 4HHB.pdb --het HEM,PO4 --water` | `4HHB_A_all_het-HEM-PO4_wat.pdb` |
| backbone / trace | `--atoms backbone` · `--atoms trace` | `…_backbone…` · `…_trace…` |
| PDB output from mmCIF input | add `-f pdb` | `….pdb` |
| one file per NMR model / segment | `--by model` · `--by segment` | `1D3Z_model1.pdb` … |

Real session (4HHB, haemoglobin):

```text
$ pdbkit split 4HHB.cif --atoms ca --het all -d chains
chains/4HHB_A_ca_het.cif
chains/4HHB_B_ca_het.cif
chains/4HHB_C_ca_het.cif
chains/4HHB_D_ca_het.cif

$ pdbkit split 4HHB.cif -c B --het ZN -d chains
pdbKit: ligand(s) ZN not present in the structure (ligands here: HEM, PO4)
chains/4HHB_B_all_het-ZN.cif
```

> [!note] Things split handles for you
> - Ligands are kept **whole**: in C-alpha mode chain A gets 141 CA atoms plus the 43 atoms of its haem, not a haem "CA".
> - Modified residues (MSE, …) are kept in C-alpha mode because they are polymer residues ([[#2.2 How residues are classified]]).
> - A DNA/RNA chain has no C-alpha. In `--atoms ca` mode it is **skipped with a message**; use `--atoms trace` (phosphates).
> - NMR ensembles keep all models in each chain file (`--models 1` for one).
> - Chain IDs longer than one character cannot go into PDB format: write mmCIF, or use `--remap-chains`.

Python:

```python
parts = pk.split_chains(s, atoms="ca", het=["HEM"], water=False)       # {"A": Structure, ...}
paths = pk.write_parts(parts, "chains", "4HHB", "pdb", "ca_het-HEM")
```

---

## 6. Selecting chains and atoms

`pdbkit split -c` selects chains **and** writes them separately. `pdbkit extract` writes the selection as **one** file. All representation options of [[#5. Splitting chains and choosing a representation]] apply.

```bash
pdbkit split   4HHB.pdb -c A,C --atoms ca --het HEM -d out   # out/4HHB_A_ca_het-HEM.pdb, out/4HHB_C_...
pdbkit extract 4HHB.pdb -c A,C --atoms ca --het HEM -o 4HHB_AC_ca_hem.pdb    # one file, two chains
pdbkit extract 4HHB.cif --exclude-chains B,D --no-water -o alpha_chains.cif
pdbkit extract 1CRN.pdb --residues 1:10 --atoms backbone -o 1crn_bb.pdb
pdbkit extract 3NIR.cif --altloc best -o 3NIR_single_conformer.cif
pdbkit extract 4HHB.cif --hetatm-only -o ligands_and_water.cif
pdbkit extract 4HHB.cif --resnames HEM -f pdb                # stdout
```

An unknown chain is an error that lists the chains present. An empty selection is an error, not an empty file.

Python: `pk.extract(s, chains="A,C", atoms="ca", het=["HEM"], water=False)` or `pk.extract(s, pk.Selection(...))`.

---

## 7. Centre of mass and centring

`pdbkit com` **reports** the centre of mass (and radius of gyration). `pdbkit center` **moves** it to the origin and writes a **new file**.

```text
$ pdbkit com 4HHB.cif
x	10.3794
y	13.2648
z	39.9448
weighting	mass
radius_of_gyration	23.5844

$ pdbkit com 4HHB.pdb -c A --atoms ca --geometric     # C-alpha geometric centre of chain A
x	23.752
y	18.428
z	29.065

$ pdbkit center 4HHB.cif                              # whole structure
centre of mass before: 10.3794 13.2648 39.9448
wrote 4HHB_centered.cif

$ pdbkit center 4HHB.pdb -c A --keep-all -o 4HHB_A0.pdb   # chain A's COM -> origin, keep all chains
$ pdbkit com 4HHB_A0.pdb -c A
x	-0.0004
y	-0.0004
z	-0.0002
```

| Goal | Command |
|---|---|
| centre the whole structure | `pdbkit center in.cif` → `in_centered.cif` |
| centre one chain, write only that chain | `pdbkit center in.cif -c A -o A_centred.cif` |
| put chain A's COM at the origin, keep everything | `pdbkit center in.cif -c A --keep-all -o out.cif` |
| protein only (ignore ligands and water) | add `--het none --no-water` |
| C-alpha centroid | `--atoms ca --geometric` |
| NMR: centre each model separately | `--per-model` (default: every model moves by model 1's COM, preserving the superposition) |
| output format | `-o name.pdb` / `-f pdb` |

> [!note] Definition and precision
> **One conformer per atom:** where an atom has alternate locations, the highest-occupancy location is used, with the full element mass (occupancies are not used as weights). This is pdbKit's fixed definition. The case-study check recomputes it independently with gemmi. For comparison, weighting every alternate location by occupancy moves the centre of 1A6M (myoglobin, 897 alternate-location atoms) by 0.06 Å.
> Masses come from the element column (IUPAC atomic weights). The residual after centring (≤ 0.0005 Å per axis above) is the 3-decimal rounding of PDB/mmCIF coordinates. All alternate locations are moved.

Python:

```python
com = pk.center_of_mass(s, chains="A")                       # numpy array
new, com_before = pk.center(s, chains="A", apply_to="all")   # apply_to="selection" (default) or "all"
pk.save(new, "4HHB_A_origin.cif")
```

`pdbkit centerasu` (Harms `pdb_centerasu`) moves the whole asymmetric unit into the unit cell. It uses SCALE records, `_atom_sites`, or the cell. NMR and cryo-EM files carry a placeholder 1 × 1 × 1 Å cell; for those `centerasu` refuses with an explanation.

---

## 8. FASTA sequences for all chains

`pdbkit fasta` writes one FASTA record per polymer chain, protein **or nucleic acid**. `4MMK.cif` produces `4MMK.fasta` **next to the input file**. For a fetched ID the file goes to the current directory; use `-o`/`-d` to choose.

```text
$ pdbkit fasta 4MMK.cif
./4MMK.fasta
$ head -3 4MMK.fasta
>4MMK_A mol:protein length:82 source:seqres Protein hfq
MSKGHSLADPYLNTLRKERVPVSIYLVNGIKLQGQIESFDQFVILLKNTVSQMVYKHAIS
TVVPSRPVRLPSGDQPAEPGNA
$ pdbkit fasta 1BNA.cif -o -
>1BNA_A mol:dna length:12 source:seqres DNA (5'-D(*CP*GP*CP*GP*AP*AP*TP*TP*CP*GP*CP*G)-3')
CGCGAATTCGCG
>1BNA_B mol:dna length:12 source:seqres DNA (5'-D(*CP*GP*CP*GP*AP*AP*TP*TP*CP*GP*CP*G)-3')
CGCGAATTCGCG
```

| `--source` | Sequence written |
|---|---|
| `auto` (default) | the deposited sequence, including unmodelled residues, when the file has it; otherwise the modelled residues. For mmCIF this is the file's own canonical sequence (`_entity_poly.pdbx_seq_one_letter_code_can`), identical to RCSB's FASTA. For PDB files it is built from SEQRES. |
| `seqres` | deposited sequence only (error if absent) |
| `atom` | residues present in the coordinates (4MMK chain A: 69 of 82) |

Other options: `-c A,B` (only some chains), `--width 80`. The header shows the molecule type (`protein`, `dna`, `rna`), the length, which source was used, and the entity description. Modified residues are written as their parent residue (MSE → M, PSU → U, 2MG → G, …), using MODRES, `_pdbx_struct_mod_residue`, `_chem_comp.mon_nstd_parent_comp_id` and a table of common modified nucleotides; residues without a known parent are written as `X`.

Python: `pk.write_fasta("4MMK.cif")`, `text = pk.fasta(s, source="atom")`, `pk.sequences(s, "seqres")` → `{"A": "MSKG…"}`.

---

## 9. Het-PDB Navi interactions

Integrated from `hetpdbnavi_parser.py`. It downloads the interaction file of an entry from [Het-PDB Navi](https://hetpdbnavi.nagahama-i-bio.ac.jp/) (or reads a saved `.dat`) and writes CSV tables:

| File | Content |
|---|---|
| `<ID>_<CODE>_<RESNUM>_<CHAIN>_interaction.csv` | one per ligand instance: the protein residues it contacts (the site's "Select hetatm" menu) |
| `<ID>_all_het_interactions.csv` | all protein–ligand contacts, with a `heteroatom` column |
| `<ID>_subunit_interactions.csv` | chain–chain interface contacts |
| `<ID>_disorder.csv` | unresolved residue ranges |

```bash
pdbkit hetnavi 4HHB -d hetnavi_4HHB                          # download + parse
pdbkit hetnavi 4HHB.cif -d out                               # ID taken from the file's header
pdbkit hetnavi 4HHB --file 4HHB_assembly1.dat -d out         # offline, from a saved file
pdbkit hetnavi 4HHB --assembly assembly2 --save-raw -d out   # other assembly; keep the .dat
```

```python
records, paths = pk.hetnavi.interactions("4HHB", outdir="out")
rows = pk.hetnavi.het_contacts_table(records)      # list of dicts -> pk.to_dataframe(rows)
```

> [!warning] Fixed: contacts silently dropped by the original parser
> Het-PDB Navi writes protein–ligand contact lines (`#HET`) in two layouts: with and **without** a trailing hydrogen-bond count. The original `hetpdbnavi_parser.py` required the count and silently skipped every other line. In the case study that emptied the ligand tables of **22 of the 24 entries** that have ligand contacts (e.g. all 15 erlotinib contacts of 1M17). pdbKit keeps these lines and leaves `num_hbonds` empty. Where the data do include the count, the CSV files are byte-identical to the original script's.

Column names are kept from the original: `res_start` is the residue's **position in the sequence** and `res_end` its **author (PDB) residue number** (1M17: `29, 694, LEU` is Leu694).

---

## 10. Domain, family and secondary-structure annotation

Integrated from `mmcif_domain_family.py`. It writes `<ID>_domain_ss_info.csv` with one row per chain: title, method, resolution, EC number/type/KEGG class, enzyme type and name, UniProt accession, family, UniProt domains mapped to PDB numbering, Pfam/InterPro/CATH/SCOP cross-references, secondary-structure order and ranges, heteroatoms and their bonding partners.

```bash
pdbkit domains 1AKE                          # download mmCIF, live UniProt + KEGG lookups
pdbkit domains 1AKE.cif 4HHB.pdb -d ann      # local files, both formats
pdbkit domains ./my_structures/ -d ann       # every structure file in a folder
pdbkit domains 1AKE.cif -q                   # offline: skip UniProt/KEGG
```

> [!tip] PDB files now get the full annotation
> The original script reported most columns as N/A for `.pdb` input, stating the data were mmCIF-only. wwPDB PDB files carry the same information in `HELIX`/`SHEET`, `SSBOND`/`LINK`, `DBREF` and `COMPND EC:` records. pdbKit translates them, and the tests check that PDB and mmCIF give the same values. The one inherent difference: hydrogen-bond connections (`hydrog`) exist only in mmCIF.

Python: `rows = pk.domains.annotate("1AKE.cif")` (list of dicts), `pk.domains.run([...], outdir)`.

---

## 11. Testing and results

This chapter records how pdbKit was tested and what the tests found (requirement 11). All numbers below are from real runs on 2026-10-01.

### 11.1 Test environment

| Item | Value |
|---|---|
| OS | Linux 7.2.7 (Fedora-based Bazzite), x86-64, 15 GB RAM |
| Primary interpreter | `biopyc` conda env: Python 3.14.6, Biopython 1.87, NumPy 2.5.1, pytest 9.1.1 |
| Second interpreter | `mayavipy` env: Python 3.12.13, Biopython 1.88 (package not installed there, run from `src/`) |
| Viewer oracle | PyMOL 3.1.0 (headless, in `biopyc`) |
| Parser oracle | gemmi 0.7.5 (`py312` env) |
| Test data | 14 wwPDB entries in **both** formats (`tests/data`), 3 Het-PDB Navi files |

Test entries, chosen to exercise specific cases:

| Entry | Why |
|---|---|
| 1CRN | small protein, 3 disulfides |
| 4HHB | 4 chains, HEM + PO4 ligands, water; low-precision 1984 model |
| 1AKE | 2 chains, bound inhibitor AP5, alternate locations |
| 1BRS | 6 chains, 3 biological assemblies, chain breaks |
| 1D3Z · 1L2Y | NMR ensembles (10 and 38 models), hydrogens |
| 1BNA | DNA double helix |
| 1WHZ | selenomethionine (MSE) modified residues |
| 3NIR | 478 alternate-location atoms, ANISOU records, **microheterogeneity** (PRO/SER 22, LEU/ILE 25) |
| 4MMK | 12 chains, metal ions (Zn, K, Na) |
| 6LU7 | covalent peptide inhibitor with non-standard residues and an ester-linked cap |
| 2GB1 | NMR structure with hydrogens |
| 1EHZ | tRNA-Phe with eleven kinds of modified nucleotides (added after the case study) |
| 1K4C | KcsA–Fab complex with long deposited disulfides (added after the case study) |

### 11.2 How to run the tests

```bash
cd pdbKit
export PYTHONWARNINGS="ignore:You may be importing Biopython"
python -m pytest tests                    # everything (network tests skip themselves when offline)
python -m pytest tests -m "not network"   # offline only
python -m pytest tests -m oracle          # only comparisons against external programs
PDBKIT_OFFLINE=1 python -m pytest tests   # force offline
```

Runs were memory-capped (`ulimit -v 6000000`).

### 11.3 Test design: independent oracles

Tests that compare pdbKit with *itself* prove little, so almost every test compares it with something independent:

| Oracle | Used for |
|---|---|
| The original **pdb-tools** scripts, run as subprocesses | every selection, deletion, renumbering, renaming, B/occupancy, split, merge, intersect, ensemble, FASTA, count and gap tool |
| The original **Harms-lab** functions, imported from `upstream/` | centre of mass (4 weightings), ASU centring, contacts, close contacts, water contacts, neighbours, torsions, MW/pI, Coulomb, ion distances, ligands, oligomer, experiment, SEQRES sequences, B-factor loading, renumbering |
| The original **hetpdbnavi_parser.py** and **mmcif_domain_family.py** | byte-for-byte comparison of every CSV |
| **PDB vs mmCIF** of the same entry | every feature runs on both; results must be identical |
| **PyMOL** and **gemmi** | reading pdbKit output exactly as they read the original wwPDB files |
| **RCSB FASTA service** | chain sequences |
| A fixed-column PDB parser in `tests/conftest.py` (no Biopython) | field-by-field comparison of written files |
| Physics / geometry | translation invariance, rigid-body distances, SSBOND records, NumPy recomputation |

### 11.4 Results

**Final run (Python 3.14.6, Biopython 1.87, network available): 369 passed, 0 failed, 0 skipped.**
Wall time 3 min 46 s, peak memory 280 MB.
**Python 3.12.13 / Biopython 1.88, offline: 350 passed, 14 network tests deselected, 0 failed** (3 min 17 s, 143 MB).
**Python 3.10 grammar:** every source and test file parses with `ast.parse(..., feature_version=(3, 10))`.
**Notebook:** `docs/examples/pdbKit_demo.ipynb` executes end-to-end under `jupyter nbconvert --execute` with 0 errors.

| Test file | What it checks | Tests | Passed | Time (s) |
|---|---|---|---|---|
| `test_diff_harms.py` | vs original Harms-lab functions | 68 | 68 | 28 |
| `test_diff_pdbtools.py` | vs original pdb-tools scripts | 45 | 45 | 26 |
| `test_edit_meta_cli.py` | editing, metadata parity, FASTA, CLI, fetching, case-study regressions | 91 | 91 | 46 |
| `test_integrations.py` | Het-PDB Navi & domain annotation vs original scripts | 26 | 26 | 23 |
| `test_io.py` | Reading/writing, PDB↔mmCIF parity, round trips, writer limits | 82 | 82 | 46 |
| `test_notebooks.py` | the example notebooks execute without errors | 2 | 2 | 13 |
| `test_select_geometry.py` | chain split modes, selection, centre of mass, copies with alternate locations | 34 | 34 | 17 |
| `test_viewers.py` | PyMOL and gemmi read outputs like the originals | 21 | 21 | 20 |
| **total** | | **369** | **369** | |

What the numbers cover, by requirement:

| Requirement | Evidence |
|---|---|
| 1 · one package, CLI + Python + notebooks | 24 commands run through the CLI on 4HHB in both formats, plus CLI tests of `fetch`, `fasta`, `split`, `center` and an `extract` → `center` stdin/stdout pipe; `merge`, `ensemble`, `intersect`, `hetnavi` and `domains` are tested through the functions their commands call; executed notebook |
| 2 · fetch from PDB and AlphaFold | RCSB entry (PDB + mmCIF), assemblies (both formats), missing-PDB-format message, 404 message; AlphaFold by accession and AF-ID, PAE, missing fragment |
| 3 · split chains in 4 representations | CA, all-atom, CA+all/specific het, all-atom+all/specific het, backbone, trace; **PDB and mmCIF input give byte-identical output** for 18 representation × 4 chain combinations; matches `pdb_splitchain`, `pdb_splitmodel` |
| 4 · select chains | matches `pdb_selchain`, `pdb_delchain`; error on unknown chain |
| 5 · centre of mass, set to zero | matches Harms `pdb_centermass` in all 4 weighting modes on 4 entries (≤ 2×10⁻⁴ Å); matches an independent NumPy calculation; re-read file has COM within 10⁻³ Å of 0; chain-only and keep-all modes; altlocs and NMR models moved correctly |
| 6–8 · Het-PDB Navi and domain annotation integrated | domain CSVs byte-identical to the original script (8 entries offline, 1 live here; all 28 case-study entries offline and live); Het-PDB Navi CSVs byte-identical wherever the original parser reads the data, plus the contact lines it silently dropped; PDB input now matches mmCIF |
| 9 · viewer compatibility | PyMOL: identical atom/state/chain counts and identical polymer/protein/nucleic/organic/inorganic/solvent classification for 10 entries × 4 output variants; gemmi: identical models/atoms/chains and entity types; strict CIF syntax parse |
| FASTA (added on request) | 12 chains of 4MMK; DNA; SEQRES vs ATOM; fallback; identical to RCSB's FASTA for 4 entries |

### 11.5 Bugs the tests found in pdbKit (all fixed)

Recorded honestly, because they show what the oracles catch:

1. **Microheterogeneity lost in mmCIF output.** In 3NIR, residues 22 and 25 each have two identities (PRO/SER, LEU/ILE). Biopython groups them in a `DisorderedResidue` whose atoms list shows only one variant, so the first mmCIF writer silently dropped 33 atoms. The PyMOL comparison exposed it; the first round-trip test had missed it because it iterated residues the same way. Now every writer, selection and editing function expands all variants, and renumbering re-keys each variant.
2. **Ester-capped peptide misclassified.** In 6LU7 the inhibitor's C-terminal cap (`010`) is attached by an ester bond, not a peptide bond, so it came back as a ligand after a PDB round trip. The linkage rule now accepts any atom bonded to the previous residue's carbonyl carbon.
3. **False chain break at a modified residue** (6LU7 LEU4→PJE5, CA–CA 4.17 Å across a 1.42 Å peptide bond). Gap detection now uses the C–N bond.
4. **R-work taken from the wrong REMARK 3 line** (4MMK: 0.229 "working + test" instead of 0.216 "working set"). Both values are now mapped to their own mmCIF items.
5. Ambiguous residue labels (`A:AP5215`) → now `A:AP5:215`. Model numbers differed between formats for single-model files (0 vs 1) → now always 1-based. `--het NAME` for an absent ligand produced a silently ligand-free file → now warns.

6. The case-study run on 28 PDB entries and 8 AlphaFold models found 13 more problems (performance on large assemblies, nucleic-acid edge cases, disulfide cutoff, FASTA of modified residues, TER placement, network robustness, `fetch -o` across filesystems, a bug inherited from the original Het-PDB Navi parser, …). All are fixed and listed in [[#11.8 Case-study run on 28 structures and 8 AlphaFold models]].
7. Writing the example notebooks exposed a Biopython pitfall: after `copy()`, an atom with alternate locations still *read* its coordinates from the original structure, so in-memory checks on a moved copy (e.g. centre of mass without saving) saw old positions for those atoms. Files were always written correctly. Every copy in pdbKit now repairs the selection (`structure.deep_copy`), with a regression test.

### 11.6 Defects confirmed in the original tools

Each is demonstrated by a test (`tests/test_diff_*.py`) that runs the original code.

| Tool | Defect | pdbKit |
|---|---|---|
| pdb-tools `pdb_gap` | prints the **squared** distance labelled as Å | Å |
| pdb-tools `pdb_gap` | reports each alternate location of a CA as a gap (18 false gaps in 3NIR, which has none) | altloc-aware |
| pdb-tools `pdb_gap` | 36 false gaps in 4HHB (CA–CA 4.00–4.26 Å across intact peptide bonds) | uses the C–N bond: 0 |
| pdb-tools `pdb_selaltloc` | crashes on files with ANISOU records (reads U values as occupancy) | handled |
| pdb-tools `pdb_selaltloc` | at microheterogeneity keeps **both** residue identities as overlapping atoms; `-A`/`-B` crash | keeps one variant |
| pdb-tools `pdb_reres` | numbers in raw record order: ligands listed after all chains are numbered after chain D; numbering continues across NMR models | per chain, per model (`--continuous` for the old behaviour) |
| Harms `centermass` | reads the element from `line[73:]`, which includes the segment-ID columns: with a segid every atom is weighted as carbon | element column |
| Harms `data/common.py` | tungsten given molybdenum's mass (95.94; W = 183.84) | IUPAC masses |
| Harms `param` | cannot run at all (`pdb_seq` never imported → NameError) | ported formulas |
| Harms `disulfide` | cannot run on Python 3 (`string.letters`) | ported rule |
| Harms `dist_filter` | raises the misspelt `IOErro` on its error path | – |
| Harms `neighbors` | reads CB coordinates from column 32, dropping the minus sign of values ≤ −100: counts change under pure translation | invariant |
| Harms `torsion` | computes φ/ψ across chain boundaries (last residue of chain A uses chain B's N) | only across real peptide bonds |
| Harms `closecontacts` | "consecutive residue" test ignores chain IDs, hiding G12(A)–C13(B) base-pair contacts in 1BNA; mis-slices 4-character atom names | chain-aware |
| Harms `watercontact` | KeyError on any atom missing from its radius table (all hydrogens); lists each residue as its own neighbour | element fallback radii |
| Harms `bfactor` | writes the text `NA` into the numeric B-factor column; keys use only 3 digits of the residue number | numeric default, full numbers |
| `mmcif_domain_family.py` | states that secondary structure, bonding, UniProt and EC data are mmCIF-only; returns N/A for PDB input | read from PDB records |
| `hetpdbnavi_parser.py` | silently drops every `#HET` contact line that has no hydrogen-bond column (22 of 24 case-study entries with contacts got empty tables) | all lines kept |

### 11.7 Limits of the testing

- Python 3.10 and 3.11 were checked **statically** (grammar) only; no 3.10/3.11 interpreter was available. Runtime tests ran on 3.12 and 3.14.
- Mol* and ChimeraX were not run. PyMOL and gemmi were, and pdbKit's mmCIF follows the wwPDB conventions those viewers rely on ([[#15. Viewer compatibility (PyMOL, Mol*, ChimeraX)]]).
- Network tests depend on live services (RCSB, AlphaFold DB, Het-PDB Navi, UniProt, KEGG); their data can change.
- Close contacts use one location per atom; the upstream comparison therefore used altloc-free inputs.


### 11.8 Case-study run on 28 structures and 8 AlphaFold models

In addition to the unit tests, every feature was run through the `pdbkit` command line on **28 PDB entries** (the mmCIF files of `mayaviBiomolecules/casestudy/`, from a 566-atom DNA dodecamer to the 149,338-atom cryo-EM ribosome 7K00) and on **8 AlphaFold DB models** downloaded with `pdbkit fetch`. Every output was verified independently: gemmi and PyMOL read every written file with equal atom counts; centre of mass recomputed with gemmi; FASTA compared with RCSB and UniProt; disulfides compared with the files' own records; domain and Het-PDB Navi tables compared with the original scripts; RCSB's PDB-format files of the same entries give identical results.

**Result: 1197/1197 commands behaved as expected; 587/587 independent checks passed.** 2168 files, 272 MB.

All outputs, one folder per structure with a terminal transcript (`commands.log`) and the verification record (`checks.json`), are in `testcases/`; `testcases/README.md` has the per-entry table, and chapter [[#19. Worked terminal examples (case-study structures)]] shows excerpts.

The case study found 13 problems, all fixed and covered by new unit tests:

| # | Problem | Fix |
|---|---|---|
| 1 | chain splitting deep-copied the whole structure per chain (7K00: 213 s) | only selected chains copied (3.7 s) |
| 2 | `neighbors`, `coulomb`, `iondist` were pure-Python O(n²) (up to 35 s) | KD-tree / NumPy, identical results (< 1 s) |
| 3 | `closecontacts --ca-only` crashed on DNA/RNA | handled (no pairs) |
| 4 | 2.5 Å disulfide default missed deposited 2.70–2.72 Å bonds (1K4C) | default 3.0 Å |
| 5 | modified residues written as `X` in FASTA (7K00, 1EHZ) | file's canonical sequence; parents of modified residues |
| 6 | PDB output: TER after waters, 81-column TER lines | TER after the polymer, 80-column records |
| 7 | `centerasu` accepted the 1 × 1 × 1 Å placeholder cell of NMR/EM files | refuses with explanation |
| 8 | Het-PDB Navi contacts without an H-bond column silently dropped (**bug of the original parser**, 22 of 24 entries affected) | all contact lines kept |
| 9 | `rename` error named an internal temporary chain | names the chains involved |
| 10 | Het-PDB Navi timeout crashed with a traceback | retries, clean error |
| 11 | `fetch -o FILE` failed across filesystems ("Invalid cross-device link") | `shutil.move` |
| 12 | `hetnavi` on an AlphaFold model crashed with a traceback | clean error explaining why |
| 13 | AlphaFold "no model" message gave no reason | explains viral proteins / >2,700-residue proteins |

### 11.9 Output files of one case-study entry (1M17)

Each folder in `testcases/` contains about 50 files. The table below explains them, using
`testcases/1M17/` as the example: the EGFR tyrosine kinase domain with the drug erlotinib (ligand
code `AQ4`). It is an X-ray structure at 2.60 Å with one chain (A), 312 modelled residues, one
ligand and 20 waters.

In the **Command** column, `1m17.cif` stands for the input file
`mayaviBiomolecules/casestudy/1m17.cif`. Files written with `-d` take their stem from the input
file name, so they are lower case (`1m17_A_ca.cif`). Files named with `-o` are upper case
(`1M17_centered.cif`). Commands without `-o` print to the terminal. For those, the case-study
script saved the printed output in the file listed (marked *stdout*).

#### Run record

| File | Command | Contents |
|---|---|---|
| `commands.log` | – | Terminal transcript: each `$ pdbkit …` command, its output (shortened if long), exit code and run time. |
| `checks.json` | – | Machine-readable record: `id`, every command with its exit code and time (`commands`), and every independent verification with result and detail (`checks`). |

#### Summary and validation

| File | Command | Contents |
|---|---|---|
| `info.json` | `info 1m17.cif --table json` | Overview: models, chain IDs, number of polymer residues, ligands and waters, atom and HETATM counts, flags for alternate locations, insertion codes and sequence gaps, ligand codes, title, method, resolution, R-work/R-free. |
| `validate.txt` | `validate 1m17.cif` (*stdout*) | Format check result: `OK: cif file, 1 model(s), 2546 atoms`. |
| `exper.tsv` | `exper 1m17.cif` (*stdout*) | Experiment summary: ID, title, classification, organism, method, resolution, R-work, R-free. |
| `assembly.tsv` | `assembly 1m17.cif` | Biological assemblies: assembly ID, oligomeric state and count, chains, number of symmetry operators, and whether the file already *is* the assembly (1M17: monomer, chain A). |

#### Sequences

| File | Command | Contents |
|---|---|---|
| `1m17.fasta` | `fasta 1m17.cif -d .` | **Deposited** sequence of every polymer chain (from SEQRES / `_entity_poly`): 333 residues, including the residues not seen in the density. |
| `1M17_modelled.fasta` | `fasta … --source atom` | Sequence of the residues that have coordinates: 312 residues. The comparison with `1m17.fasta` shows the disordered stretches. |

#### Centre of mass and centring

| File | Command | Contents |
|---|---|---|
| `com_whole.tsv` | `com 1m17.cif` (*stdout*) | Mass-weighted centre of mass (x, y, z) of all atoms and the radius of gyration. |
| `com_chainA_ca.tsv` | `com … -c A --atoms ca --geometric` | Geometric centre (unweighted) and radius of gyration of the C-alpha atoms of chain A. |
| `1M17_centered.cif` | `center 1m17.cif` | Whole structure translated so that its centre of mass is at the origin. |
| `1M17_chainA_com0.cif` | `center … -c A --keep-all` | Whole structure (protein, drug, water) translated so that **chain A's** centre of mass is at the origin. |
| `1M17_chainA_protein_centered.cif` | `center … -c A --het none --no-water` | Chain A protein only, without ligand or water, centred at the origin. |

#### Chain splitting (one file per chain)

The file-name tag records the representation: `ca` = C-alpha only, `all` = all atoms,
`het` = all heteroatoms, `het-AQ4` = only ligand AQ4, `wat` = water included.

| File | Command | Contents |
|---|---|---|
| `split_ca/1m17_A_ca.cif` | `split … --atoms ca` | Chain A, C-alpha atoms only, no ligand or water. |
| `split_ca_het/1m17_A_ca_het.cif` | `split … --atoms ca --het all` | Chain A, C-alpha atoms plus all ligands of the chain (the whole erlotinib molecule). |
| `split_all_atom_het_water/1m17_A_all_het_wat.cif` | `split … --het all --water` | Chain A, all atoms, with all ligands and water. |
| `split_selected/1m17_A_all_het-AQ4.cif` | `split … -c A --het AQ4` | Chain A, all atoms, with erlotinib (AQ4) only. |
| `split_selected/1m17_A_ca_het-AQ4.pdb` | `split … -c A --atoms ca --het AQ4 -f pdb` | Chain A, C-alpha atoms plus erlotinib, in **PDB** format. |
| `split_ca_pdb/1m17_A_ca.pdb` | `split … --atoms ca -f pdb` | Chain A, C-alpha only, in PDB format. |

#### Selections

| File | Command | Contents |
|---|---|---|
| `1M17_backbone.cif` | `extract … --atoms backbone --het none --no-water` | Protein backbone atoms (N, CA, C, O) only. |
| `1M17_hetatm_only.cif` | `extract … --hetatm-only` | HETATM records only: the 29 erlotinib atoms and 20 waters. |
| `1M17_single_conformer.cif` | `extract … --altloc best --no-hydrogens` | One conformation per atom (the alternate location with highest occupancy), without hydrogens. |

#### Format conversion

| File | Command | Contents |
|---|---|---|
| `1M17.pdb` | `convert 1m17.cif -o 1M17.pdb` | The structure in PDB format. |
| `1M17_roundtrip.cif` | `convert 1M17.pdb -o 1M17_roundtrip.cif` | Converted back to mmCIF. Used to check that a mmCIF → PDB → mmCIF round trip loses nothing. |

#### Editing tools

| File | Command | Contents |
|---|---|---|
| `1M17_renumbered.cif` | `renumber … --start 1` | Residues renumbered from 1 (the protein was originally numbered 672–995). |
| `1M17_chainA_as_z.cif` | `rename … --chain-map A:z` | Chain A renamed to `z`. |
| `1M17_b0_occ1.cif` | `setfield … --bfactor 0 --occupancy 1` | All B-factors set to 0 and occupancies to 1. |
| `1M17_sorted.cif` | `sort 1m17.cif` | Chains sorted by ID; residues by number, polymer first, then ligands, then water. |
| `1M17_asu_in_cell.cif` | `centerasu 1m17.cif` | Asymmetric unit moved into the crystallographic unit cell (Harms `pdb_centerasu`). |

#### Analyses (tab-separated tables)

| File | Command | Contents |
|---|---|---|
| `gaps.tsv` | `gaps 1m17.cif` | Chain breaks: model, chain, residues on either side, and the distance. 1M17 has one: HIS964 → LEU977, 7.19 Å. |
| `ligands.tsv` | `ligands 1m17.cif` | Ligands: chain, residue name, number, insertion code, atom count (AQ4 999, 29 atoms). |
| `disulfide.tsv` | `disulfide 1m17.cif` | Every cysteine pair within the cutoff (3.0 Å S–S) as `bonded`, and every other cysteine as `free`. All six cysteines of 1M17 are free. |
| `contacts.tsv` | `contacts … -t AQ4 --cutoff 4.0` | Residues with any atom within 4.0 Å of erlotinib, with their minimum distance: the binding site (15 residues, e.g. Leu694, Ala719, Met769). |
| `closecontacts.tsv` | `closecontacts … --ca-only` | Pairs of C-alpha atoms from non-adjacent residues closer than 3.5 Å. **Empty** for 1M17, because no such pairs exist. |
| `watercontacts.tsv` | `watercontacts 1m17.cif` | For each residue: the number and list of neighbouring residues, found from van der Waals radii plus a water probe (Harms `watercontact`). |
| `neighbors.tsv` | `neighbors 1m17.cif` | For each residue: the number of CB atoms within 10 Å of its CB (pseudo-CB for Gly), excluding residues within 4 positions in sequence. This is a burial measure. |
| `torsion.tsv` | `torsion 1m17.cif` | Backbone φ and ψ for each residue, in degrees; empty where undefined (chain ends, breaks). Use it for a Ramachandran plot. |
| `sasa.tsv` | `sasa 1m17.cif` | Solvent-accessible surface area per residue (Å², Shrake–Rupley) and relative SASA (fraction of the maximum for that residue type). |
| `distances.tsv` | `distances … -r A:672:CA --atom CA -c A` | Distance from the CA of residue A672 to every CA in chain A. |
| `iondist.tsv` | `iondist … --step 0.5` | Histogram of distances between titratable groups, in 0.5 Å bins, with counts for acid–acid, acid–base and base–base pairs. |
| `param.tsv` | `param 1m17.cif` (*stdout*) | Sequence parameters: residue count, molecular weight, pI, fractions of titratable, acidic and basic residues, net charge, amino-acid composition, non-standard residues. |
| `coulomb.tsv` | `coulomb 1m17.cif` (*stdout*) | Screened Coulomb energy of the titratable groups in kcal/mol (Harms `pdb_coulomb` defaults). |

#### Annotation and interactions

| File | Command | Contents |
|---|---|---|
| `domains/1M17_domain_ss_info.csv` | `domains 1m17.cif -d domains` | One row per chain: title, method, resolution, entity, EC number and KEGG classification, enzyme type and name, UniProt accession and family, UniProt domains in PDB numbering, Pfam/InterPro/CATH/SCOP cross-references, secondary-structure order and ranges, heteroatoms and their bonding partners. Uses live UniProt and KEGG lookups. |
| `domains_offline/1M17_domain_ss_info.csv` | `domains … -q` | The same table built offline from the file alone. Columns that need UniProt/KEGG read `N/A (skipped …)`. |
| `hetnavi/1M17_assembly1.dat` | `hetnavi 1M17 --save-raw` | The raw Het-PDB Navi interaction file, as downloaded. |
| `hetnavi/1M17_AQ4_999_A_interaction.csv` | ″ | Protein residues in contact with ligand AQ4 999 of chain A: sequence position (`res_start`), author residue number (`res_end`), residue name, contact type, number of contacts and hydrogen bonds. One such file is written per ligand instance. |
| `hetnavi/1M17_all_het_interactions.csv` | ″ | All protein–ligand contacts of the entry in one table, with a `heteroatom` column (`AQ4:999:A`). |
| `hetnavi/1M17_disorder.csv` | ″ | Unresolved (disordered) residue ranges per chain, by sequence position: 1–6, 300–311, 331–333. |

> [!note] Files that depend on the entry
> Some outputs appear only when the structure has the matching content. For example:
> - `hetnavi/<ID>_subunit_interactions.csv` appears only for entries with chain–chain interfaces (e.g. 1AON, 1AAY). 1M17 has a single chain.
> - `split_models/` appears only for NMR ensembles (1BCT).
> - Entries with more chains get one file per chain in each `split_*` folder.

---

## 12. Editing tools (pdb-tools equivalents)

| pdb-tools | pdbKit command | Python |
|---|---|---|
| `pdb_fetch` | `pdbkit fetch` | `pk.fetch` |
| `pdb_selchain` / `pdb_delchain` | `extract -c` / `--exclude-chains` | `pk.extract(chains=…)` |
| `pdb_selatom` | `extract -a CA,C,N,O` | `atoms=[…]` |
| `pdb_selelem` / `pdb_delelem` | `--elements` / `--exclude-elements` | |
| `pdb_selres` / `pdb_delres` | `--residues` / `--exclude-residues` | |
| `pdb_selresname` / `pdb_delresname` | `--resnames` / `--exclude-resnames` | |
| `pdb_delhetatm` | `--no-hetatm` | `only_het_records=False` |
| `pdb_selaltloc` | `--altloc best` / `--altloc A` | `altloc=` |
| `pdb_selseg` | `--segids` | |
| `pdb_splitchain` / `pdb_splitmodel` / `pdb_splitseg` | `split --by chain/model/segment` | `pk.split_chains/…` |
| `pdb_reres` / `pdb_shiftres` | `renumber --start N [--continuous]` / `--shift K` | `pk.renumber_residues` |
| `pdb_reatom` | `renumber --atoms-from N` | `pk.renumber_atoms` |
| `pdb_chain` / `pdb_rplchain` | `rename --set-chain X` / `--chain-map A:B` | `pk.set_chain` / `pk.rename_chains` |
| `pdb_rplresname` | `rename --residue-map HSD:HIS` | `pk.rename_residues` |
| `pdb_seg` / `pdb_chainxseg` / `pdb_segxchain` | `rename --set-segid / --chain-to-segid / --segid-to-chain` | |
| `pdb_b` / `pdb_occ` | `setfield --bfactor 10 / --occupancy 1` | `pk.set_bfactor/…` |
| `pdb_sort` | `sort` | `pk.sort_structure` |
| `pdb_merge` | `merge a.pdb b.pdb [--rename-duplicates]` | `pk.merge` |
| `pdb_mkensemble` / `pdb_chkensemble` | `ensemble a b c` / `ensemble --check x.pdb` | `pk.make_ensemble` / `pk.check_ensemble` |
| `pdb_intersect` | `intersect a b` | `pk.intersect` |
| `pdb_tofasta` | `fasta` | `pk.fasta` |
| `pdb_tocif` / `pdb_fromcif` | `convert` (to the other format by default) | `pk.save(s, "x.cif")` |
| `pdb_tidy` / `pdb_keepcoord` / `pdb_element` | `convert -f pdb` (clean coordinates, TER/END, elements) | |
| `pdb_wc` | `info` | `pk.summary` |
| `pdb_gap` | `gaps` | `pk.gaps` |
| `pdb_validate` | `validate` | `pk.validate_file` |
| `pdb_head` | use `head` | |

Examples:

```bash
pdbkit renumber in.cif --start 1 -o r.cif                  # 1..n in every chain
pdbkit renumber in.pdb --start 101 --keep-gaps -c A        # offset chain A, keep gaps
pdbkit rename in.pdb --chain-map A:H,B:L -o hl.pdb         # swaps are allowed: A:B,B:A
pdbkit setfield in.pdb --bfactor-file scores.txt -o colored.pdb   # "42 0.87" or "A:42 0.87" per line
pdbkit merge receptor.pdb ligand.pdb --rename-duplicates -o complex.pdb
pdbkit validate in.pdb
```

`setfield --bfactor-file` puts per-residue values (conservation, predictions, …) into the B-factor column for colouring in PyMOL (`spectrum b`) or Mol*.

---

## 13. Analysis tools (Harms-lab equivalents)

| Harms script | pdbKit | Output |
|---|---|---|
| `pdb_centermass` | `com`, `center` | COM, recentred file |
| `pdb_centerasu` | `centerasu` | file |
| `pdb_contact` | `contacts -t HEM` (name, `42` or `A:42`) | target, residue, min distance |
| `pdb_closecontacts` | `closecontacts [--ca-only] --cutoff 3.5` | atom pairs |
| `pdb_watercontact` | `watercontacts --probe 1.4` | neighbour lists |
| `pdb_neighbors` | `neighbors --cutoff 10 --seq-cutoff 4` | CB neighbour counts |
| `pdb_disulfide` | `disulfide [--cutoff 3.0]` | bonded pairs, free Cys |
| `pdb_torsion` | `torsion` | φ, ψ |
| `pdb_param` | `param [-c A] [--source seqres]` | MW, pI, composition, charge |
| `pdb_coulomb` | `coulomb --dielectric 40 --ionic-strength 100 --ph 7` | kcal/mol |
| `pdb_iondist` | `iondist --step 0.1` | histogram |
| `pdb_dist_filter` | `distances -r A:87:NE2 --atom FE` | distances |
| `pdb_ligand` | `ligands [--skip-boring]` | ligand instances |
| `pdb_exper` | `exper` | title, organism, method, resolution, R-work, R-free |
| `pdb_oligomer` | `assembly` | assemblies; `in_asymmetric_unit` (= Harms "wysiwyg"), `is_asymmetric_unit` (exact) |
| `pdb_seq` | `fasta --source seqres` | FASTA |
| `pdb_sasa` (needed NACCESS) | `sasa [--level atom] [--het]` | Shrake–Rupley SASA and relative SASA (Tien 2013) |
| `pdb_bfactor`, `pdb_offset`, `pdb_residue_renumber`, `pdb_atom_renumber`, `pdb_subset`, `pdb_splitnmr`, `pdb_download` | `setfield`, `renumber`, `extract`, `split --by model`, `fetch` | |

Table commands print tab-separated text; add `--table csv` or `--table json`, and `-o file` to save.

```text
$ pdbkit disulfide 1CRN.pdb
cys1	cys2	distance	state
A:CYS:3	A:CYS:40	2.004	bonded
A:CYS:4	A:CYS:32	2.035	bonded
A:CYS:16	A:CYS:26	2.051	bonded

$ pdbkit assembly 1BRS.cif
id	details	oligomeric_details	oligomeric_count	chains	n_operators	in_asymmetric_unit	is_asymmetric_unit
1	author_defined_assembly	dimeric	2	A,D	1	True	False
2	author_defined_assembly	dimeric	2	B,E	1	True	False
3	author_defined_assembly	dimeric	2	C,F	1	True	False

$ pdbkit gaps 1BRS.cif
model	chain	from	to	type	distance
1	D	D:THR:63	D:GLY:66	distance	4.88
1	E	E:THR:63	E:GLY:66	distance	4.95
```

> [!warning] Not ported
> `pdb_addH`, `pdb_clean` and `pdb_mutator` drive the commercial **CHARMM** program. `pdb_satk` needs compiled Fortran and NACCESS. `pdb_moment` needs the UHBD-format hydrogens that `pdb_addH` writes. They cannot run without those programs, so they are not part of pdbKit. The originals remain in `upstream/pdbtools-master`.

---

## 14. Using pdbKit from Python and Jupyter

Everything the CLI does is a function. Structures are ordinary Biopython `Structure` objects, so pdbKit mixes freely with Biopython code. Functions **return new structures** and never modify their input.

```python
import pdbkit as pk

s = pk.load("1AKE.cif")                    # file, "1AKE", "AF-P69905-F1", or "-" (stdin)
pk.summary(s); pk.experiment(s)            # dicts
a = pk.extract(s, chains="A", het="none", water=False)
parts = pk.split_chains(s, atoms="ca", het="all")
new, com = pk.center(s, chains="A", apply_to="all")
pk.save(new, "1AKE_A0.pdb")                # format from the extension
text = pk.dumps(a, "cif")                  # string

df = pk.to_dataframe(pk.contacts(s, "AP5", 4.0))    # pandas DataFrame
pk.to_dataframe(pk.sasa(s))
pk.to_dataframe(pk.domains.annotate(s, skip_lookups=True))
records, paths = pk.hetnavi.interactions("1AKE", outdir="hn")
```

Main functions:

| Area | Functions |
|---|---|
| I/O | `load`, `loads`, `save`, `dumps`, `fetch`, `fetch_pdb`, `fetch_alphafold`, `write_fasta` |
| selection | `Selection`, `extract`, `split_chains`, `split_models`, `split_segments`, `write_parts`, `chain_ids`, `residue_kind` |
| geometry | `center_of_mass`, `center`, `translate`, `radius_of_gyration`, `center_asu` |
| editing | `renumber_residues`, `renumber_atoms`, `rename_chains`, `set_chain`, `rename_residues`, `set_bfactor`, `set_occupancy`, `bfactor_from_data`, `set_segid`, `chain_to_segid`, `segid_to_chain`, `sort_structure`, `merge`, `make_ensemble`, `check_ensemble`, `intersect` |
| analysis | `summary`, `sequences`, `fasta`, `gaps`, `ligands`, `disulfides`, `contacts`, `close_contacts`, `water_contacts`, `neighbors`, `torsions`, `param`, `coulomb`, `ion_distance_histogram`, `distances`, `sasa`, `experiment`, `assemblies` |
| metadata | `pk.metadata.as_mmcif_dict(s)` – mmCIF-style dict for **either** format; `seqres`, `chain_descriptions` |
| integrations | `pk.hetnavi.*`, `pk.domains.*` |
| validation | `validate_file`, `validate_text` |

Worked notebooks are in `docs/examples/`; each explains the code and how to read every output:

| Notebook | Content |
|---|---|
| `pdbKit_protein_EGFR_erlotinib.ipynb` | protein–drug complex EGFR/erlotinib (1M17): representations, binding site (contacts, buried surface, Het-PDB Navi), Ramachandran plot, drug centred for docking, annotation, crystal vs AlphaFold (pLDDT) |
| `pdbKit_dna_1BNA_Zif268.ipynb` | DNA: Watson–Crick pairs, minor-groove width, spine of hydration (1BNA); zinc sites and protein–DNA recognition contacts of Zif268 (1AAY) |
| `pdbKit_demo.ipynb` | short tour on haemoglobin (4HHB) |

They run offline from `docs/examples/data/` and are executed by the test suite.

---

## 15. Viewer compatibility (PyMOL, Mol*, ChimeraX)

**PDB output** is written with Biopython's `PDBIO` and then tidied: every record is exactly 80 columns, **TER follows each chain's last polymer residue** (not its waters, as `PDBIO` does), atom serials count the TER records, and element columns, MODEL/ENDMDL and END are written. If a structure cannot be represented in PDB format (chain IDs longer than one character, residue names longer than three characters, residue numbers outside −999…9999, more than 99 999 atoms), pdbKit **refuses with an explanation** rather than writing a corrupt file.

**mmCIF output** uses pdbKit's own writer, which follows the wwPDB conventions viewers use to classify molecules:

- `_entity` (polymer / non-polymer / water), `_entity_poly`, `_struct_asym`, `_chem_comp` (peptide/DNA/RNA linking or non-polymer), `_atom_type`;
- a separate `label_asym_id` for every polymer chain, every ligand instance, and each chain's waters;
- `label_seq_id` preserved from the source mmCIF; otherwise derived so that **chain breaks stay visible** (viewers do not draw cartoon across missing residues);
- all author fields (`auth_asym_id`, `auth_seq_id`, …), alternate locations and residue variants preserved; values quoted where CIF requires it (`"O5'"`).

Tests load every output in PyMOL and gemmi and compare with the original wwPDB file ([[#11. Testing and results]]).

> [!note] Why PyMOL shows extra atoms for a wwPDB mmCIF
> For original wwPDB mmCIF files, PyMOL adds coordinate-less CA placeholders for **unmodelled** residues from the deposited sequence (158 in 4MMK). pdbKit output contains only modelled atoms, so its count matches the number of real atoms (`count_atoms present` in PyMOL).

---

## 16. Differences from the original tools

- **One structure model for both formats.** pdb-tools and Harms tools edit PDB text lines; pdbKit parses into a structure and writes a clean file. Atom serial numbers are renumbered on output (use `renumber --atoms-from` to control them). Within a chain, ligands and waters are written after that chain's polymer.
- **Semantics that differ deliberately:** see the pdbKit column of [[#11.6 Defects confirmed in the original tools]]. Also: `disulfide` defaults to 3.0 Å. The ideal S–S bond is ~2.05 Å, but deposited disulfides in poorly restrained regions reach 2.7 Å (1K4C), while non-bonded SG–SG contacts start near 3.4 Å. With this default it reproduces the SSBOND / `_struct_conn` records of every test and case-study entry. Harms used 4.0 Å (`--cutoff 4`). `pdb_selatom` kept any atom with the name everywhere (including a calcium ion named CA); pdbKit's `--atoms` applies to polymer residues, and ligands follow `--het`.
- **Header records** (REMARKs, CRYST1, SEQRES …) are not copied into edited coordinate files. Metadata stays available through `pk.metadata` while the structure is loaded.

---

## 17. Troubleshooting

| Message | Meaning / fix |
|---|---|
| `cannot be written in PDB format (chain IDs longer than 1 character …)` | write mmCIF (`-f cif`), or add `--remap-chains` (mapping is printed) |
| `RCSB has no PDB-format file for …` | large entry: fetch with `-f cif` |
| `The selection matched no atoms` | check `pdbkit info file` for chain IDs and ligand names |
| `chain 'X' has no atoms matching the selection - skipped` | e.g. a DNA chain in `--atoms ca` mode: use `--atoms trace` |
| `ligand(s) ZN not present …` | the message lists the ligands that are present |
| `No SEQRES / _pdbx_poly_seq_scheme records` | use `--source atom` (the default `auto` falls back by itself) |
| `unknown element(s) given the mass of carbon` | the element column is missing or unusual; run `pdbkit convert` or fix the file |
| BiopythonWarning about `pyproject.toml` | environment issue, see [[#1. Installation]] |

---

## 18. Reference: all commands

```text
fetch  convert  info  validate  fasta                     input / output
extract  split                                            selection (features 3, 4)
com  center  centerasu                                    geometry (feature 5)
renumber  rename  setfield  sort  merge  ensemble  intersect      editing
gaps  ligands  disulfide  contacts  closecontacts  watercontacts  neighbors
torsion  param  coulomb  iondist  distances  sasa  exper  assembly   analysis
hetnavi  domains                                          integrations (features 6-8)
```

`pdbkit <command> -h` shows every option. Project layout:

```text
pdbKit/
├── pyproject.toml, README.md, LICENSE (GPL-3.0), NOTICE
├── src/pdbkit/        structure.py (I/O, classification) · cifio.py (mmCIF writer) · fetch.py
│                      select.py · geometry.py · edit.py · analysis.py · metadata.py
│                      hetnavi.py · domains.py · validate.py · cli.py
├── tests/             test_*.py, conftest.py (oracles), data/ (14 entries × 2 formats)
├── docs/              pdbKit_User_Manual.md (this file), examples/pdbKit_demo.ipynb
└── upstream/          the original packages and scripts, unchanged (used by the tests)
```

**Licence.** pdbKit includes GPL-3.0 code from Harms-lab pdbtools, so the package is distributed under GPL-3.0-or-later. pdb-tools' Apache-2.0 code may be combined into it (see `NOTICE`).

---

## 19. Worked terminal examples (case-study structures)

These sessions were run in a folder containing structures from the case study (`mayaviBiomolecules/casestudy/`), plus AlphaFold models downloaded with `pdbkit fetch` ([[#19.9 AlphaFold models]]). The output shown is exactly what pdbKit printed. The complete transcripts of all 28 structures, about 43 commands each, are in `testcases/<ID>/commands.log`, and every output file is next to its log ([[#11.8 Case-study run on 28 structures and 8 AlphaFold models]]).

### 19.1 A 21-chain complex: GroEL/GroES (1AON)

```text
$ pdbkit info 1aon.cif
id	1AON
models	1
chains	["A", "B", "C", "D", "E", "F", "G", "H", "I", "J", "K", "L", "M", "N", "O", "P", "Q", "R", "S", "T", "U"]
polymer_residues	8015
ligand_residues	14
water_molecules	0
atoms	58870
hetatm_atoms	196
...
ligands	{"MG": 7, "ADP": 7}
exp_title	CRYSTAL STRUCTURE OF THE ASYMMETRIC CHAPERONIN COMPLEX GROEL/GROES/(ADP)7
exp_method	X-RAY DIFFRACTION
exp_resolution	3.0

$ pdbkit assembly 1aon.cif
id	details	oligomeric_details	oligomeric_count	chains	n_operators	in_asymmetric_unit	is_asymmetric_unit
1	author_and_software_defined_assembly	21-meric	21	A,B,C,D,E,F,G,H,I,J,K,L,M,N,O,P,Q,R,S,T,U	1	True	True

$ pdbkit split 1aon.cif -c A,H,O --atoms ca --het ADP,MG -d groel
groel/1aon_A_ca_het-ADP-MG.cif
groel/1aon_H_ca_het-ADP-MG.cif
groel/1aon_O_ca_het-ADP-MG.cif
```

Three chains were selected (one from each GroEL ring and one GroES), each as a C-alpha trace with its ADP and Mg²⁺ kept as whole ligands.

### 19.2 Protein–DNA complex: Zif268 zinc fingers (1AAY)

```text
$ pdbkit split 1aay.cif --atoms trace -d zif268       # CA for the protein, P for the DNA strands
zif268/1aay_B_trace.cif
zif268/1aay_C_trace.cif
zif268/1aay_A_trace.cif

$ pdbkit fasta 1aay.cif                               # writes 1aay.fasta next to the input
./1aay.fasta
$ head -2 1aay.fasta
>1AAY_B mol:dna length:11 source:seqres DNA (5'-D(*AP*GP*CP*GP*TP*GP*GP*GP*CP*GP*T)-3')
AGCGTGGGCGT

$ pdbkit center 1aay.cif -c A --het none --no-water -o zif268_protein_centered.pdb
centre of mass before: 4.5500 10.3275 48.1012
wrote zif268_protein_centered.pdb
```

`--atoms ca` would skip the DNA chains (they have no C-alpha); `--atoms trace` keeps one atom per nucleotide.

### 19.3 Centring a high-resolution structure with alternate locations (1A6M, myoglobin)

```text
$ pdbkit com 1a6m.cif
x	14.0222
y	21.0561
z	8.3082
weighting	mass
radius_of_gyration	15.6306

$ pdbkit center 1a6m.cif
centre of mass before: 14.0222 21.0561 8.3082
wrote 1a6m_centered.cif

$ pdbkit com 1a6m_centered.cif
x	0.0002
y	0.0001
z	0.0002
weighting	mass
radius_of_gyration	15.6306
```

The residual is coordinate rounding (3 decimals). The radius of gyration is unchanged because the move is a pure translation.

### 19.4 Modified nucleotides: tRNA-Phe (1EHZ)

```text
$ pdbkit fasta 1ehz.cif -o -
>1EHZ_A mol:rna length:76 source:seqres TRANSFER RNA (PHE)
GCGGAUUUAGCUCAGUUGGGAGAGCGCCAGACUGAAGAUCUGGAGGUCCUGUGUUCGAUC
CACAGAAUUCGCACCA
```

The eleven kinds of modified nucleotides (PSU, 2MG, M2G, OMC, OMG, 7MG, 5MU, 5MC, 1MA, H2U, YYG) appear as their parent bases. The sequence is identical to RCSB's FASTA.

### 19.5 Disulfides, ligand contacts and ligand lists

```text
$ pdbkit disulfide 1k4c.cif                 # KcsA-Fab complex
cys1	cys2	distance	state
A:CYS:22	A:CYS:96	2.724	bonded
B:CYS:23	B:CYS:88	2.705	bonded
B:CYS:134	B:CYS:194	2.71	bonded
A:CYS:145			free
A:CYS:200			free
C:CYS:90			free

$ pdbkit contacts 4nxa.cif -t HEM --cutoff 3.5
target	residue	min_distance
A:HEM:201	A:LYS:42	3.253
A:HEM:201	A:PHE:43	3.429
A:HEM:201	A:ARG:45	3.01
A:HEM:201	A:SER:92	2.64
A:HEM:201	A:HIS:93	2.108
A:HEM:201	A:HIS:97	2.618
A:HEM:201	A:ILE:99	3.429

$ pdbkit ligands 2cb8.cif --skip-boring
chain	resname	resseq	icode	atoms
A	MYA	1088		27
A	ZN	1089		1
...
B	MXE	1092		5
```

HIS93 at 2.1 Å from the haem is the proximal iron ligand. The three disulfides of 1K4C are exactly those in the file's `_struct_conn` records.

### 19.6 Ligand interactions from Het-PDB Navi: EGFR + erlotinib (1M17)

```text
$ pdbkit hetnavi 1m17.cif -d hetnavi_1m17      # the PDB ID is read from the file
Parsed 1 chain record(s).
hetnavi_1m17/1M17_AQ4_999_A_interaction.csv
hetnavi_1m17/1M17_all_het_interactions.csv
hetnavi_1m17/1M17_disorder.csv

$ head -6 hetnavi_1m17/1M17_AQ4_999_A_interaction.csv
protein_chain,res_start,res_end,res_name,contact_type,num_contacts,num_hbonds,het_code,het_resnum,het_chain,het_index
A,29,694,LEU,A,7,,AQ4,999,A,1
A,54,719,ALA,S,3,,AQ4,999,A,1
A,56,721,LYS,S,7,,AQ4,999,A,1
A,73,738,GLU,S,1,,AQ4,999,A,1
A,99,764,LEU,A,3,,AQ4,999,A,1
```

The original parser wrote no rows for this entry ([[#9. Het-PDB Navi interactions]]).

### 19.7 NMR ensemble (1BCT) and a cryo-EM ribosome (7K00)

```text
$ pdbkit split 1bct.cif --by model -d nmr_models
nmr_models/1bct_model1.cif
...
nmr_models/1bct_model14.cif

$ pdbkit centerasu 1bct.cif -o x.cif
pdbkit: error: No crystallographic unit cell (SOLUTION NMR files carry a 1 x 1 x 1 A placeholder cell); centerasu only applies to crystal structures.

$ pdbkit convert 7k00.cif -f pdb -o 7k00.pdb
pdbkit: error: This structure cannot be written in PDB format (329 residue number(s) outside -999..9999; 149338 atoms per model (PDB format limit is 99999)). Write mmCIF instead (--format cif / a .cif file name).

$ pdbkit split 7k00.cif -c A,a --atoms trace -f pdb -d ribosome     # single chains do fit in PDB format
ribosome/7k00_A_trace.pdb
ribosome/7k00_a_trace.pdb
```

When an operation does not apply, pdbKit exits with status 2 and an explanation; it does not write a wrong file.

### 19.8 Sequence parameters and domain annotation (1CA2, carbonic anhydrase II)

```text
$ pdbkit param 1ca2.cif
residues	256
molecular_weight	28766.13
pI	7.118
fraction_titratable	0.3203
fraction_acidic	0.125
fraction_basic	0.1602
charge_fully_ionised	9
...

$ pdbkit domains 1ca2.cif -d ann
  1ca2.cif -> ann/1CA2_domain_ss_info.csv
ann/1CA2_domain_ss_info.csv
```

### 19.9 AlphaFold models

```text
$ pdbkit fetch P69905 -d .                           # UniProt accession -> current model, mmCIF
./AF-P69905-F1.cif
$ pdbkit fetch AF-P69905-F1 -f pdb -d .              # AlphaFold ID, PDB format
./AF-P69905-F1.pdb
$ pdbkit fetch P69905 --pae -f pdb --overwrite -d with_pae   # + predicted aligned error (with_pae/AF-P69905-F1_pae.json)
with_pae/AF-P69905-F1.pdb

$ pdbkit fetch Q8WZ42 -d .                           # titin
pdbkit: error: AlphaFold DB has no model for Q8WZ42. AlphaFold DB does not cover viral proteins, and proteins longer than 2,700 residues are available only as fragments in bulk proteome downloads; also check the accession.

$ pdbkit fasta AF-P00533-F1.cif -d .                 # EGFR, 1210 residues
./AF-P00533-F1.fasta
$ head -2 AF-P00533-F1.fasta
>AF-P00533-F1_A mol:protein length:1210 source:seqres Epidermal growth factor receptor
MRPSGTAGAALLALLAALCPASRALEEKKVCQGTSNKLTQLGTFEDHFLSLQRMFNNCEV

$ pdbkit com AF-P00533-F1.cif
x	-6.7514
y	6.6987
z	3.4504
weighting	mass
radius_of_gyration	40.0615

$ pdbkit hetnavi AF-P00533-F1.cif -d hetnavi
pdbkit: error: AF-P00533-F1.cif has no 4-character PDB ID (entry ID 'AF-P00533-F1'). Het-PDB Navi covers PDB entries only - predicted models (AlphaFold) have no interaction data there. Pass a PDB ID if the file is a copy of one.
```

AlphaFold models keep their per-residue confidence (pLDDT) in the B-factor column. Every pdbKit output preserves it, so `spectrum b` in PyMOL colours split, centred or converted files by confidence, exactly as it does the original model. The case study checked this on 8 models ([[#11.8 Case-study run on 28 structures and 8 AlphaFold models]]).

