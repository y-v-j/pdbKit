# pdbKit case-study test cases

Every feature of pdbKit was run through the `pdbkit` command line on the 28 mmCIF files in `computational_geometry_projects/mayaviBiomolecules/casestudy/`, and every output was checked independently (gemmi, PyMOL, RCSB, and the original scripts). Only outputs are stored here; the input files were not copied.

**Result: 1197/1197 commands behaved as expected; 587/587 independent checks passed.** 2168 files, 272 MB.

**Result: PDB entries - 1197/1197 commands behaved as expected, 587/587 independent checks passed. AlphaFold models - 244/244 commands, 118/118 checks.** 2423 files, 296 MB.

| File / folder | Content |
|---|---|
| `commands.log` | terminal transcript: each `$ pdbkit ...` command, its output (long outputs shortened) and exit code |
| `checks.json` | every command (exit code, time) and every verification with its result and detail |
| `info.json`, `exper.tsv`, `assembly.tsv`, `validate.txt` | summary, experiment, assemblies, format check |
| `<id>.fasta`, `<ID>_modelled.fasta` | FASTA of all chains: deposited sequence / modelled residues |
| `com_*.tsv`, `<ID>_centered.cif`, `<ID>_chainX_com0.cif`, `<ID>_chainX_protein_centered.cif` | centre of mass; centred structure; chain X's COM at the origin with all chains kept; chain X alone (no ligands/water) centred |
| `split_ca/`, `split_ca_het/`, `split_all_atom_het_water/`, `split_selected/`, `split_ca_pdb/`, `split_models/` | chain splitting in the different representations (file-name tag = representation) |
| `<ID>_backbone.cif`, `<ID>_hetatm_only.cif`, `<ID>_single_conformer.cif` | selections |
| `<ID>.pdb`, `<ID>_roundtrip.cif` | format conversion and round trip |
| `<ID>_renumbered.cif`, `<ID>_chainX_as_z.cif`, `<ID>_b0_occ1.cif`, `<ID>_sorted.cif`, `<ID>_asu_in_cell.cif` | editing tools |
| `gaps.tsv`, `ligands.tsv`, `disulfide.tsv`, `contacts.tsv`, `closecontacts.tsv`, `watercontacts.tsv`, `neighbors.tsv`, `torsion.tsv`, `iondist.tsv`, `sasa.tsv`, `distances.tsv`, `param.tsv`, `coulomb.tsv` | analyses |
| `domains/`, `domains_offline/` | domain / family / EC / secondary-structure CSV (live lookups / offline) |
| `hetnavi/` | Het-PDB Navi interaction CSVs and the raw `.dat` |

## Results per entry

| ID | Method | Models | Chains | Atoms | Ligands | Commands OK | Checks passed | Time (s) | Files | MB |
|---|---|---|---|---|---|---|---|---|---|---|
| [1A6M](1A6M/commands.log) | X-RAY DIFFRACTION | 1 | 1 | 1584 | SO4;HEM;OXY | 43/43 | 21/21 | 29.2 | 50 | 1.7 |
| [1AAY](1AAY/commands.log) | X-RAY DIFFRACTION | 1 | 3 | 1330 | ZN | 43/43 | 21/21 | 25.1 | 52 | 1.3 |
| [1AON](1AON/commands.log) | X-RAY DIFFRACTION | 1 | 21 | 58870 | MG;ADP | 43/43 | 21/21 | 198.8 | 144 | 60.5 |
| [1BCT](1BCT/commands.log) | SOLUTION NMR | 14 | 1 | 1103 | - | 43/43 | 21/21 | 54.7 | 57 | 15.9 |
| [1BNA](1BNA/commands.log) | X-RAY DIFFRACTION | 1 | 2 | 566 | - | 42/42 | 21/21 | 21.8 | 44 | 0.6 |
| [1BRS](1BRS/commands.log) | X-RAY DIFFRACTION | 1 | 6 | 5153 | - | 42/42 | 21/21 | 35.5 | 68 | 5.2 |
| [1CA2](1CA2/commands.log) | X-RAY DIFFRACTION | 1 | 1 | 2207 | ZN | 43/43 | 21/21 | 32.6 | 48 | 2.5 |
| [1EHZ](1EHZ/commands.log) | X-RAY DIFFRACTION | 1 | 1 | 1821 | MG;MN | 43/43 | 21/21 | 25.4 | 51 | 1.9 |
| [1HVR](1HVR/commands.log) | X-RAY DIFFRACTION | 1 | 2 | 1890 | XK2 | 43/43 | 21/21 | 26.7 | 54 | 2.1 |
| [1IFC](1IFC/commands.log) | X-RAY DIFFRACTION | 1 | 1 | 2590 | - | 42/42 | 21/21 | 27.7 | 46 | 2.7 |
| [1K4C](1K4C/commands.log) | X-RAY DIFFRACTION | 1 | 3 | 4534 | F09;K;DGA | 43/43 | 21/21 | 35.0 | 64 | 5.0 |
| [1L87](1L87/commands.log) | X-RAY DIFFRACTION | 1 | 1 | 1434 | CL;BME | 43/43 | 21/21 | 27.5 | 49 | 1.6 |
| [1M17](1M17/commands.log) | X-RAY DIFFRACTION | 1 | 1 | 2560 | AQ4 | 43/43 | 21/21 | 31.7 | 48 | 3.0 |
| [1MJ5](1MJ5/commands.log) | X-RAY DIFFRACTION | 1 | 1 | 5802 | MG;CL | 43/43 | 21/21 | 44.5 | 49 | 6.0 |
| [1VFB](1VFB/commands.log) | X-RAY DIFFRACTION | 1 | 3 | 3456 | - | 42/42 | 21/21 | 30.7 | 55 | 3.5 |
| [2BFN](2BFN/commands.log) | X-RAY DIFFRACTION | 1 | 1 | 2743 | D2P;CL;CA | 43/43 | 21/21 | 31.8 | 51 | 3.0 |
| [2CB8](2CB8/commands.log) | X-RAY DIFFRACTION | 1 | 2 | 1918 | MYA;ZN;SO4;MXE | 43/43 | 21/21 | 29.2 | 64 | 2.1 |
| [2F9Q](2F9Q/commands.log) | X-RAY DIFFRACTION | 1 | 4 | 14421 | HEM;SO4 | 43/43 | 21/21 | 64.5 | 68 | 15.3 |
| [2YXP](2YXP/commands.log) | X-RAY DIFFRACTION | 1 | 1 | 2878 | - | 42/42 | 21/21 | 31.3 | 45 | 3.2 |
| [3A2M](3A2M/commands.log) | X-RAY DIFFRACTION | 1 | 5 | 5200 | CL;GLC;FRU | 43/43 | 21/21 | 38.2 | 68 | 5.7 |
| [3ERT](3ERT/commands.log) | X-RAY DIFFRACTION | 1 | 1 | 2070 | OHT | 43/43 | 21/21 | 27.6 | 48 | 2.4 |
| [3U1T](3U1T/commands.log) | X-RAY DIFFRACTION | 1 | 2 | 5164 | CL;MLI | 43/43 | 21/21 | 43.7 | 60 | 5.7 |
| [4F4L](4F4L/commands.log) | X-RAY DIFFRACTION | 1 | 4 | 2759 | - | 42/42 | 21/21 | 28.5 | 60 | 2.9 |
| [4F4X](4F4X/commands.log) | X-RAY DIFFRACTION | 1 | 3 | 3648 | CA;DCP | 43/43 | 21/21 | 47.2 | 57 | 4.1 |
| [4HZG](4HZG/commands.log) | X-RAY DIFFRACTION | 1 | 1 | 2642 | CL | 43/43 | 21/21 | 32.6 | 48 | 3.0 |
| [4MMK](4MMK/commands.log) | X-RAY DIFFRACTION | 1 | 12 | 6928 | ZN;K;NA;SO4 | 43/43 | 21/21 | 41.4 | 128 | 7.0 |
| [4NXA](4NXA/commands.log) | X-RAY DIFFRACTION | 1 | 1 | 1353 | HEM;ZN;SO4;XE | 43/43 | 21/21 | 25.5 | 58 | 1.5 |
| [7K00](7K00/commands.log) | ELECTRON MICROSCOPY | 1 | 56 | 149338 | PAR;MG;SPD;SPM;ZN | 42/42 | 20/20 | 529.5 | 534 | 103.1 |

## Failed checks

None.

## Commands that correctly refused (exit 2 with a message)

These operations do not apply to the entry; pdbKit refuses with an explanation instead of writing a wrong file.

| ID | Command | Message |
|---|---|---|
| 1AAY | distances | pdbkit: error: Reference atom 'B:1:CA' not found. |
| 1BCT | extract_hetatm | pdbkit: error: The selection matched no atoms. |
| 1BCT | centerasu | pdbkit: error: No crystallographic unit cell (SOLUTION NMR files carry a 1 x 1 x 1 A placeholder cell); centerasu only applies to crystal structures. |
| 1BNA | distances | pdbkit: error: Reference atom 'A:1:CA' not found. |
| 1BNA | param | pdbkit: error: No protein residues selected. |
| 1EHZ | distances | pdbkit: error: Reference atom 'A:1:CA' not found. |
| 1EHZ | param | pdbkit: error: No protein residues selected. |
| 7K00 | convert_pdb | pdbkit: error: This structure cannot be written in PDB format (329 residue number(s) outside -999..9999; 149338 atoms per model (PDB format limit is 99999)). Wr |
| 7K00 | distances | pdbkit: error: Reference atom 'A:2:CA' not found. |
| 7K00 | centerasu | pdbkit: error: No crystallographic unit cell (ELECTRON MICROSCOPY files carry a 1 x 1 x 1 A placeholder cell); centerasu only applies to crystal structures. |
| 7K00 | rename | pdbkit: error: Cannot rename chain 'A' to 'z': chain 'z' already exists and residue 2 is in both. Rename or renumber one of them first. |

## AlphaFold DB models (`alphafold/`)

Models were downloaded with `pdbkit fetch` and run through the same features (29 commands each). Download modes are tested in `alphafold/fetch_tests.log`: accession and AF-ID, mmCIF and PDB, `--pae`, `-o` file, `-o -` (stdout), several IDs, lower-case input, cache loading by ID, and the refusals for a missing fragment, titin (fragments only), T4 lysozyme (viral, not in AlphaFold DB) and `--assembly`.

**Download modes: 12/12 commands as expected; 6/6 checks passed.**

| Model | Protein | Residues | Atoms | Mean pLDDT | Commands OK | Checks passed |
|---|---|---|---|---|---|---|
| [AF-P02185-F1](alphafold/AF-P02185-F1/commands.log) | myoglobin, sperm whale (1A6M) | 154 | 1225 | 97.5 | 29/29 | 14/14 |
| [AF-P00918-F1](alphafold/AF-P00918-F1/commands.log) | carbonic anhydrase II, human (1CA2) | 260 | 2073 | 97.4 | 29/29 | 14/14 |
| [AF-P0DP23-F1](alphafold/AF-P0DP23-F1/commands.log) | calmodulin-1, human | 149 | 1174 | 85.2 | 29/29 | 14/14 |
| [AF-P0A6F5-F1](alphafold/AF-P0A6F5-F1/commands.log) | GroEL, E. coli (1AON) | 548 | 4001 | 90.8 | 29/29 | 14/14 |
| [AF-P03372-F1](alphafold/AF-P03372-F1/commands.log) | estrogen receptor alpha, human (3ERT) | 595 | 4637 | 66.4 | 29/29 | 14/14 |
| [AF-P00533-F1](alphafold/AF-P00533-F1/commands.log) | EGFR, human (1M17) | 1210 | 9392 | 75.9 | 29/29 | 14/14 |
| [AF-P0A334-F1](alphafold/AF-P0A334-F1/commands.log) | KcsA potassium channel (1K4C) | 160 | 1250 | 89.6 | 29/29 | 14/14 |
| [AF-P69905-F1](alphafold/AF-P69905-F1/commands.log) | haemoglobin alpha, human (4HHB) | 142 | 1077 | 98.1 | 29/29 | 14/14 |

AlphaFold checks: model sequence = UniProt sequence = AlphaFold DB API sequence = FASTA output; PAE matrix is L x L; pLDDT (B-factor column) preserved through conversion, centring, renumbering and splitting; mean pLDDT equals AlphaFold DB's global value; PDB and mmCIF downloads give identical atoms and centre of mass; centring and C-alpha centroid ~0 (gemmi); no gaps or ligands; domain annotation finds the UniProt entry; every output readable by gemmi and PyMOL with equal atom counts and valid under `pdbkit validate`; `centerasu` and `hetnavi` refuse with an explanation.

## What each check verifies

| Check | Meaning |
|---|---|
| `info_atoms_equal_gemmi` | atom count reported by `pdbkit info` equals gemmi's count of the input |
| `info_models_equal_gemmi` | model count equals gemmi's |
| `all_outputs_parse_in_gemmi` | every structure file written (.cif/.pdb) is read by gemmi |
| `all_outputs_load_in_pymol` | every structure file written loads in PyMOL (headless) |
| `gemmi_and_pymol_agree_on_atom_counts` | gemmi and PyMOL count the same atoms in every output file |
| `split_all_is_a_partition` | all-atom + all-het + water split: chain files together hold every input atom exactly once |
| `split_ca_counts_equal_gemmi` | C-alpha split: one file per protein chain, CA count per chain equals gemmi's |
| `split_ca_pdb_equals_cif` | C-alpha split written as PDB has the same atoms as the mmCIF split |
| `centered_com_zero_gemmi` | centred file: centre of mass recomputed with gemmi (pdbKit definition) is ~0 |
| `centered_com_zero_pdbkit` | centred file re-read by pdbKit: |COM| < 0.001 A |
| `keepall_chain_com_zero` | `center -c X --keep-all`: chain X's COM ~0 (gemmi) and no atom lost |
| `keepall_is_rigid_translation` | `--keep-all` output is a pure translation of the input (all atoms shifted equally) |
| `fasta_equals_rcsb` | every FASTA record equals RCSB's FASTA service for that chain |
| `disulfides_equal_struct_conn` | disulfides found equal the file's own _struct_conn disulf records (read by gemmi) |
| `cif_pdb_cif_round_trip_identical` | mmCIF -> PDB -> mmCIF keeps every atom, altloc and coordinate |
| `all_outputs_validate` | every structure file written passes `pdbkit validate` |
| `domains_offline_identical_to_original` | domain/family CSV (-q) byte-identical to the original mmcif_domain_family.py |
| `domains_live_identical_to_original` | domain/family CSV with live UniProt+KEGG byte-identical to the original script |
| `hetnavi_matches_original` | Het-PDB Navi: one CSV row per contact line of the raw data; CSVs byte-identical to the original hetpdbnavi_parser.py where its parser reads the data (see finding 8) |
| `pdb_input_parity` | RCSB's PDB-format file of the entry gives the same FASTA, sequences, split files, COM and domain columns |
| `all_commands_exit_as_expected` | every command exited as expected (0, or 2 with a clear message where the operation does not apply) |

## Errors found by this case study and fixed in pdbKit

1. **Chain splitting was quadratic**: each chain extraction deep-copied the whole structure (7K00: 213 s). Now only the selected chains are copied (3.7 s).
2. **Slow pure-Python pair loops**: `neighbors` (35 s on 1AON), `coulomb` (12 s), `iondist` (7 s) now use a KD-tree / NumPy (< 1 s, identical results).
3. **`closecontacts --ca-only` crashed on nucleic acids** (1BNA, 1EHZ: no C-alpha atoms).
4. **Disulfides deposited at 2.70-2.72 A were missed** (1K4C) by the 2.5 A default; default is now 3.0 A (non-bonded SG-SG contacts start near 3.4 A).
5. **FASTA of modified residues**: modified nucleotides and amino acids (7K00 ribosome, 1EHZ tRNA) were written as X. The deposited sequence now uses the file's own canonical sequence (identical to RCSB), and modelled sequences map modifications to their parent residue.
6. **PDB output**: TER records were placed after a chain's waters and were 81 columns long; TER now follows the last polymer residue and every record is 80 columns.
7. **`centerasu` on NMR / cryo-EM entries** ran on the 1 x 1 x 1 A placeholder cell and produced a meaningless shift; it now refuses with an explanation.
8. **Het-PDB Navi contacts silently dropped (bug in the original parser)**: most entries list protein-ligand contacts on `#HET` lines without the hydrogen-bond column; the original `hetpdbnavi_parser.py` required that column and wrote empty tables for 22 of the 24 entries with contacts (e.g. 1M17: all 15 erlotinib contacts lost). pdbKit now keeps them.
9. **Misleading error text**: `rename` refusing a clash named an internal temporary chain; it now names the chains involved.
10. **Het-PDB Navi download not robust**: a server timeout (seen for 1A6M and 1AAY during this run) crashed `pdbkit hetnavi` with a Python traceback. Downloads now retry transient errors and fail with a one-line message; both entries were re-run successfully.
11. **`pdbkit fetch -o FILE` failed across filesystems** (AlphaFold study): the file was downloaded to /tmp and moved with `os.replace`, which cannot cross filesystems ("Invalid cross-device link"). Now moved with `shutil.move`.
12. **`pdbkit hetnavi` on a predicted model crashed** with a Python traceback; it now explains that Het-PDB Navi covers PDB entries only.
13. **Unhelpful AlphaFold 'no model' message**: now states the usual reasons (viral proteins are not covered; proteins over 2,700 residues exist only as fragments in bulk downloads).

Centre-of-mass definition (as decided): one conformer per atom, the highest-occupancy alternate location, full element mass. The gemmi check recomputes exactly that definition independently.
