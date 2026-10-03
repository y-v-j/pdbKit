# pdbKit

pdbKit is a toolkit for three-dimensional macromolecular structure files. It can be used as a command-line program (`pdbkit`), as a Python library, or in IPython/Jupyter notebooks. It treats **PDB** (`.pdb`, `.ent`) and **PDBx/mmCIF** (`.cif`, `.mmcif`) files the same way- compressed and non-compressed.

With pdbKit you can:

- **fetch** experimental structures from RCSB PDB and predicted models from AlphaFold DB;
- **select** and **split** chains, in C-alpha, backbone, trace or all-atom form, with or without
  all or specific heteroatoms and water;
- compute the **centre of mass** of a whole structure or of selected chains, and move it to the
  origin;
- **edit** structures: renumber, rename, set B-factor/occupancy, sort, merge, build ensembles, and
  convert between PDB and mmCIF;
- **analyse** structures: FASTA, gaps, ligands, disulfides, contacts, torsions, SASA, pI/MW, and
  experiment and assembly reports;
- **annotate** structures with Het-PDB Navi protein–ligand and subunit interactions, and with
  domain, family, EC and secondary-structure information.

Output files open in PyMOL, Mol*, ChimeraX and other molecular viewers.

## Acknowledgements and origin

pdbKit merges two existing open-source toolkits into one package and builds on their code and
algorithms. Many thanks to both labs:

| Project | Lab | Licence | Used in pdbKit for |
|---|---|---|---|
| [**pdb-tools** 2.0.0-rc1](https://github.com/haddocking/pdb-tools) | **Bonvin lab** (HADDOCK team; João Rodrigues et al.), Utrecht University | Apache-2.0 | selection, deletion, renumbering, renaming, splitting, merging, gaps, FASTA, validation, fetching |
| [**pdbtools**](https://github.com/harmslab/pdbtools) | **Harms lab** (Michael J. Harms), University of Oregon | GPL-3.0 | centre of mass, ASU centring, contacts, disulfides, neighbours, torsions, pI/MW, Coulomb energy, ion distances, ligands, experiment/oligomer reports, B-factor loading |

The project scripts `hetpdbnavi_parser.py` (Het-PDB Navi parser) and `mmcif_domain_family.py`
(domain/family annotator) are also part of the package. The original sources are kept, with their
licences, in [upstream/](upstream/). See [NOTICE](NOTICE) for details.

If you use pdbKit in published work, please also cite the original tools, for example
*Rodrigues JPGLM, Teixeira JMC, Trellet M, Bonvin AMJJ. pdb-tools: a swiss army knife for
molecular structures. F1000Research 2018, 7:1961*.

## What pdbKit changes and adds

pdbKit is not a copy of the original tools. pdb-tools and the Harms-lab pdbtools are separate
command-line scripts that edit PDB text line by line. pdbKit puts them into a single package built
around a common structure model, and adds the following.

### New features

- **PDB and mmCIF support throughout.** Every command reads both formats, gzipped or not, and can
  write either one. The original tools work only on PDB-format text.
- **Standards-conforming mmCIF writer.** It writes `_entity`, `_entity_poly`, `_struct_asym`,
  `_chem_comp` and the other wwPDB categories that viewers use, keeps author numbering and
  alternate locations, and leaves chain breaks visible. Output loads correctly in PyMOL, Mol* and
  ChimeraX.
- **PDB-header → mmCIF translator** (`pk.metadata`). It reads `HELIX`/`SHEET`, `SSBOND`/`LINK`,
  `DBREF`, `SEQRES`, `COMPND EC:` and `REMARK` records into mmCIF-style data. Features that need
  header information therefore work the same on PDB and mmCIF input.
- **Fetching** from RCSB PDB (entries and biological assemblies) and from **AlphaFold DB** (by
  AlphaFold ID or UniProt accession, with optional PAE). The current model version is found
  through the AFDB API, and downloads are cached.
- **Chain splitting and selection in several representations:** C-alpha, backbone, trace
  (C-alpha/P) or all-atom. Each can be combined with all, specific or no heteroatoms, and with or
  without water.
- **Centre of mass** of a whole structure or of selected chains, and **centring** it at the origin
  in a new file.
- **FASTA for all polymer chains** (protein, DNA, RNA), taken from SEQRES when available.
- **Consistent residue classification** (polymer / ligand / water). Modified residues such as MSE
  count as polymer in both formats.
- **Python, IPython and Jupyter API:** every command is also a function, results can be turned
  into pandas DataFrames with `to_dataframe`, and the package includes example notebooks.
- **SASA by Shrake–Rupley** (through Biopython), replacing the Harms-lab `sasa`, which needed the
  external NACCESS program.

### Integrated project scripts

- **Het-PDB Navi parser** (`hetpdbnavi_parser.py` → `pdbkit hetnavi`, `pk.hetnavi`). It
  downloads or reads Het-PDB Navi interaction files and writes per-ligand, all-ligand,
  subunit-interface and disorder CSV tables. The original parser **silently dropped** every
  protein–ligand contact line that had no hydrogen-bond count. In a 28-entry case study, this
  emptied the ligand tables of 22 of the 24 entries that had ligand contacts. pdbKit keeps these
  lines.
- **Domain and family annotation** (`mmcif_domain_family.py` → `pdbkit domains`,
  `pk.domains`). It writes one row per chain with EC number and KEGG class, UniProt accession,
  family, UniProt domains mapped to PDB numbering, Pfam/InterPro/CATH/SCOP cross-references,
  secondary structure, and heteroatoms with their bonding partners. The original script returned
  N/A for most columns on PDB input. pdbKit takes the same information from the PDB header
  records, so **PDB and mmCIF input give the same annotation**. There is one exception:
  hydrogen-bond records exist only in mmCIF. A `-q` option skips the online UniProt/KEGG lookups.

### Fixes to defects in the original tools

Each defect below is confirmed by a test in `tests/test_diff_*.py` that runs the original code.

| Original tool | Defect | In pdbKit |
|---|---|---|
| pdb-tools `pdb_gap` | reports the squared distance as Å; reports false gaps at alternate locations and across intact peptide bonds (36 in 4HHB) | gaps found from the C–N bond, distance in Å, alternate locations handled |
| pdb-tools `pdb_selaltloc` | crashes on files with ANISOU records; keeps both residue identities at microheterogeneity sites | handled; keeps one variant |
| pdb-tools `pdb_reres` | numbers residues in raw record order, continuing across ligands and NMR models | numbers per chain and per model (`--continuous` gives the old behaviour) |
| Harms `centermass` | weights every atom as carbon when a segment ID is present | reads the element column |
| Harms `data/common.py` | tungsten has molybdenum's mass | IUPAC masses |
| Harms `param`, `disulfide` | do not run (missing import; Python 2 code) | ported and working |
| Harms `neighbors` | loses the minus sign of coordinates ≤ −100, so counts change when the structure is translated | translation-invariant |
| Harms `torsion` | computes φ/ψ across chain boundaries | only across real peptide bonds |
| Harms `closecontacts` | ignores chain IDs when skipping consecutive residues; mis-reads 4-character atom names | chain-aware |
| Harms `watercontact` | crashes on hydrogens; lists each residue as its own neighbour | fallback radii by element |
| Harms `bfactor` | writes the text `NA` into the numeric B-factor column; uses only 3 digits of the residue number | numeric default, full numbers |

### Behaviour that differs on purpose

- Files are parsed into a structure and written out cleanly, so atom serial numbers are
  renumbered. Header records (REMARK, CRYST1, SEQRES, …) are not copied into edited coordinate
  files, but stay available through `pk.metadata` while the structure is loaded.
- PDB output always has 80-column records and a TER after each chain's last polymer residue. pdbKit
  refuses, with an explanation, to write PDB format for structures that cannot be represented in
  it, such as multi-character chain IDs or more than 99,999 atoms.
- `--atoms` applies only to polymer residues and ligands follow `--het`. The original
  `pdb_selatom -CA` also kept calcium ions, whose atom name is `CA`.
- `disulfide` uses a 3.0 Å cutoff by default; Harms used 4.0 Å. With 3.0 Å it reproduces the
  deposited SSBOND records of all test and case-study entries.

The [User Manual](docs/pdbKit_User_Manual.md) has the full list (sections 11.5–11.8 and 16),
including the problems that testing found in pdbKit itself and how they were fixed.

## Requirements

| Dependency | Version | Required? |
|---|---|---|
| Python | ≥ 3.10 | yes |
| [Biopython](https://biopython.org) | ≥ 1.80 | yes |
| [NumPy](https://numpy.org) | ≥ 1.22 | yes |
| pandas | any recent | optional, for `to_dataframe` |
| pytest | any recent | optional, to run the tests |

You need internet access only to fetch structures, for Het-PDB Navi, and for the UniProt/KEGG
lookups of `domains`. Everything else works offline.

## Installation

### 1. Install the dependencies

With conda (recommended):

```bash
conda create -n pdbkit -c conda-forge python=3.12 biopython numpy pandas pytest
conda activate pdbkit
```

Or with pip in a virtual environment:

```bash
python -m venv .venv
source .venv/bin/activate          # Windows: .venv\Scripts\activate
pip install "biopython>=1.80" "numpy>=1.22" pandas pytest
```

### 2. Install pdbKit

From the pdbKit folder (the one that contains `pyproject.toml`):

```bash
pip install .                      # normal install
# or
pip install -e .                   # editable install, for development
pip install -e ".[dataframes,test]"   # also installs pandas and pytest
```

pip also installs the required dependencies if they are missing, so step 1 is optional.

Check the installation:

```bash
pdbkit --version                   # pdbKit 1.0.0
pdbkit --help                      # list of all commands
python -c "import pdbkit; print(pdbkit.__version__)"
```

To uninstall: `pip uninstall pdbKit`.

## Usage

### Inputs and outputs

Any command that reads a structure accepts:

- a local file (`.pdb`, `.ent`, `.cif`, `.mmcif`, optionally `.gz`);
- a **PDB ID** such as `4HHB`, or an **AlphaFold ID** such as `AF-P69905-F1` or a UniProt
  accession such as `P69905`. The file is downloaded once and cached in `~/.cache/pdbkit`; set
  `PDBKIT_CACHE` to change the location;
- `-` for standard input, so commands can be piped.

pdbKit chooses the output format from `-f pdb|cif`, then from the extension of `-o`, and
otherwise uses the input format.

### Command line

**Fetch structures**

```bash
pdbkit fetch 4HHB                        # -> 4HHB.cif (mmCIF is the default)
pdbkit fetch 4HHB 1CRN -f pdb -d pdbs    # several entries, PDB format, into pdbs/
pdbkit fetch 1BRS --assembly 1           # biological assembly 1
pdbkit fetch P69905 --pae                # AlphaFold model + predicted aligned error
```

**Inspect**

```bash
pdbkit info 4HHB.cif                     # chains, residue/atom counts, ligands, method, resolution
pdbkit validate 4HHB.pdb                 # format check
pdbkit fasta 4HHB.cif                    # -> 4HHB.fasta, all protein and nucleic-acid chains
```

**Split chains into separate files**

```bash
pdbkit split 4HHB.cif                         # all-atom, one file per chain
pdbkit split 4HHB.cif --atoms ca              # C-alpha only
pdbkit split 4HHB.cif --atoms ca --het all    # C-alpha + all heteroatoms
pdbkit split 4HHB.cif --atoms ca --het HEM    # C-alpha + specific heteroatoms
pdbkit split 4HHB.cif --het HEM,PO4 --water -f pdb -d chains   # all-atom + HEM/PO4 + water, PDB output
```

```text
$ pdbkit split 4HHB.cif --atoms ca --het HEM -d chains -f pdb
chains/4HHB_A_ca_het-HEM.pdb
chains/4HHB_B_ca_het-HEM.pdb
chains/4HHB_C_ca_het-HEM.pdb
chains/4HHB_D_ca_het-HEM.pdb
```

**Select chains, residues and atoms**

```bash
pdbkit extract 1BRS.cif -c A,D --no-water -o 1BRS_AD.cif
pdbkit extract 4HHB.cif -c A --atoms ca --het HEM --no-water -o 4HHB_A_ca_HEM.pdb
pdbkit extract 1AKE.cif -c A --residues 1:100 --no-hydrogens -o 1AKE_A_1-100.pdb
```

**Centre of mass**

```bash
pdbkit com 4HHB.cif                      # centre of mass + radius of gyration
pdbkit com 4HHB.cif -c A                 # one chain only
pdbkit center 1CRN.pdb                   # move the centre of mass to the origin -> 1CRN_centered.pdb
```

**Pipes**

```bash
pdbkit fetch 1BRS -f pdb -o - | pdbkit extract - -c A,D --no-water | pdbkit center - -o 1BRS_AD_centred.pdb
```

**Editing and conversion**

```bash
pdbkit convert 4HHB.pdb -o 4HHB.cif      # PDB -> mmCIF (and back)
pdbkit renumber 1CRN.pdb --start 1 -o 1CRN_renum.pdb
pdbkit rename 4HHB.pdb --chain-map A:X,B:Y -o 4HHB_XY.pdb
pdbkit merge a.pdb b.pdb -o complex.pdb
```

**Analysis**

```bash
pdbkit ligands 4HHB.cif
pdbkit disulfide 1CRN.pdb
pdbkit contacts 1AKE.cif -t AP5 --cutoff 4.0   # residues within 4 Å of ligand AP5
pdbkit torsion 1CRN.pdb                  # phi/psi
pdbkit sasa 1CRN.pdb
pdbkit param 1CRN.pdb                    # MW, pI, composition, charge
```

**Interactions and annotation**

```bash
pdbkit hetnavi 4HHB -d hetnavi_4HHB      # Het-PDB Navi ligand/subunit contacts -> CSV
pdbkit domains 1AKE -d annotation        # domain/family/EC/secondary structure -> CSV
pdbkit domains 1AKE.cif -q               # offline: skip UniProt/KEGG lookups
```

Run `pdbkit <command> --help` for every option of a command.

### Python, IPython and Jupyter

All CLI functions are also available from Python. Structures are ordinary Biopython `Structure`
objects. pdbKit functions return new structures and never modify their input.

```python
import pdbkit as pk

s = pk.load("4HHB")                      # a file path, a PDB/AlphaFold ID, or "-"
pk.summary(s)                            # dict: chains, counts, ligands, ...

# split chains: C-alpha + haem, then save one chain as PDB
parts = pk.split_chains(s, atoms="ca", het="HEM")
pk.save(parts["A"], "4HHB_A_ca_HEM.pdb")

# select chain A without ligands or water
a = pk.extract(s, chains="A", het="none", water=False)

# centre the whole structure on the centre of mass of chain A
centred, com = pk.center(s, chains="A", apply_to="all")
pk.save(centred, "4HHB_centred.cif")

print(pk.fasta(s))                       # FASTA text for all polymer chains

# tables as pandas DataFrames (needs pandas)
pk.to_dataframe(pk.ligands(s))
pk.to_dataframe(pk.contacts(pk.load("1AKE"), "AP5", 4.0))
pk.to_dataframe(pk.domains.annotate("1AKE.cif", skip_lookups=True))

# Het-PDB Navi interactions (needs internet)
records, paths = pk.hetnavi.interactions("4HHB", outdir="hetnavi_4HHB")
```

Example notebooks that run offline are in [docs/examples/](docs/examples/):
a protein–drug complex (EGFR/erlotinib, 1M17), DNA and a protein–DNA complex (1BNA, Zif268), and
a short tour on haemoglobin (4HHB).

## Documentation

The full [User Manual](docs/pdbKit_User_Manual.md) covers every command and option, how residues
are classified, viewer compatibility, the differences from the original tools, troubleshooting,
the test results, and worked examples on 28 structures.

## Running the tests

```bash
pip install -e ".[test]"
pytest -m "not network"                  # offline tests
pytest                                   # all tests, including those that need internet
```

## Licence

pdbKit is distributed under the **GNU General Public License v3.0 or later** (see
[LICENSE](LICENSE)), because it includes GPL-3.0 code from the Harms lab's pdbtools.
The Apache-2.0 code from the Bonvin lab's pdb-tools can be combined into a GPL-3.0 work, and its
licence and attribution are preserved in [NOTICE](NOTICE) and [upstream/](upstream/).
