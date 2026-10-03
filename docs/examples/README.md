# pdbKit example notebooks

| Notebook | What it covers |
|---|---|
| `pdbKit_protein_EGFR_erlotinib.ipynb` | **Protein**: EGFR kinase with the drug erlotinib (1M17). Structure overview, chain-splitting representations, the drug's binding site (contacts, buried surface, Het-PDB Navi), Ramachandran plot, moving the drug to the origin, domain annotation and FASTA, and a comparison with the AlphaFold model using pLDDT. |
| `pdbKit_dna_1BNA_Zif268.ipynb` | **DNA**: the Dickerson dodecamer (1BNA) — strands, phosphate trace, Watson–Crick pairs, minor-groove width, spine of hydration, centring — and the Zif268 zinc-finger–DNA complex (1AAY) — zinc sites and the protein–DNA recognition contacts. |
| `pdbKit_demo.ipynb` | A short tour of the main functions on haemoglobin (4HHB). |

Each section explains the code, shows the output, and explains how to read it; the equivalent `pdbkit` terminal command is given alongside.

**Running them.** Open in JupyterLab/Notebook from this folder with the environment where pdbKit is installed (`conda activate biopyc`). The input files are in `data/` (no internet needed); output files are written to `output/protein/` and `output/dna/`. The notebooks are executed by the test suite (`tests/test_notebooks.py`), so their outputs reflect the current code.
