"""Chemical and naming constants shared by all pdbKit modules.

Most tables are carried over from the Harms-lab ``pdbtools`` (``data/common.py``)
so that ported calculations (pI, molecular weight, Coulomb energy, ion
distances, water contacts) reproduce the original numbers. Two deliberate
changes are documented inline: element masses come from Biopython's complete
IUPAC table, and tungsten is no longer given the mass of molybdenum.
"""

from __future__ import annotations

from Bio.Data.IUPACData import atom_weights as _IUPAC_WEIGHTS

# --------------------------------------------------------------------------- #
# Element masses (g/mol), keyed by upper-case element symbol
# --------------------------------------------------------------------------- #
# The Harms table mapped "W" to 95.94 (molybdenum); tungsten is 183.84.
ATOM_WEIGHTS = {sym.upper(): mass for sym, mass in _IUPAC_WEIGHTS.items()}
ATOM_WEIGHTS.setdefault("D", 2.0141)
CARBON_MASS = ATOM_WEIGHTS["C"]

# --------------------------------------------------------------------------- #
# Residue classes
# --------------------------------------------------------------------------- #
WATER_NAMES = frozenset({"HOH", "WAT", "DOD", "H2O", "TIP", "TIP3", "TIP4", "SOL", "OH2", "SPC"})

AMINO_ACIDS = frozenset({
    "ALA", "ARG", "ASN", "ASP", "CYS", "GLN", "GLU", "GLY", "HIS", "ILE",
    "LEU", "LYS", "MET", "PHE", "PRO", "SER", "THR", "TRP", "TYR", "VAL",
    "SEC", "PYL", "ASX", "GLX", "UNK",
    # common force-field protonation-state names
    "HSD", "HSE", "HSP", "HID", "HIE", "HIP", "CYX", "CYM", "ASH", "GLH", "LYN",
})
NUCLEOTIDES = frozenset({"A", "C", "G", "U", "I", "N", "DA", "DC", "DG", "DT", "DU", "DI", "DN"})
DNA_NUCLEOTIDES = frozenset({"DA", "DC", "DG", "DT", "DU", "DI", "DN"})

BACKBONE_ATOMS = ("N", "CA", "C", "O")
# "trace" atoms: one representative atom per polymer residue
TRACE_ATOMS = {"protein": "CA", "nucleic": "P"}

# Ligands that pdb_ligand (Harms) treats as uninteresting crystallisation additives
BORING_LIGANDS = frozenset({"HOH", "CA", "SO4", "IOD", "NA", "CL", "GOL", "PO4"})

# --------------------------------------------------------------------------- #
# One-letter codes
# --------------------------------------------------------------------------- #
AA3_TO_AA1 = {
    "ALA": "A", "CYS": "C", "ASP": "D", "GLU": "E", "PHE": "F", "GLY": "G",
    "HIS": "H", "ILE": "I", "LYS": "K", "LEU": "L", "MET": "M", "ASN": "N",
    "PRO": "P", "GLN": "Q", "ARG": "R", "SER": "S", "THR": "T", "VAL": "V",
    "TRP": "W", "TYR": "Y", "SEC": "U", "PYL": "O", "ASX": "B", "GLX": "Z",
    "UNK": "X",
    # aliases used by the Harms tools and by MD force fields
    "HSE": "H", "HSD": "H", "HSP": "H", "HID": "H", "HIE": "H", "HIP": "H",
    "CYX": "C", "CYM": "C", "ASH": "D", "GLH": "E", "LYN": "K", "MSE": "M",
}
NUC_TO_1 = {
    "A": "A", "C": "C", "G": "G", "U": "U", "I": "I", "N": "N",
    "DA": "A", "DC": "C", "DG": "G", "DT": "T", "DU": "U", "DI": "I", "DN": "N",
}

# --------------------------------------------------------------------------- #
# Electrostatics / titration (Harms-lab values, unchanged)
# --------------------------------------------------------------------------- #
PKA_DICT = {"ASP": 4.0, "CYS": 8.5, "GLU": 4.4, "TYR": 10.0, "CTERM": 3.1,
            "ARG": 12.0, "HIS": 6.5, "LYS": 10.4, "NTERM": 8.0}
CHARGE_DICT = {"ASP": -1.0, "CYS": -1.0, "GLU": -1.0, "TYR": -1.0, "CTERM": -1.0,
               "ARG": 1.0, "HIS": 1.0, "LYS": 1.0, "NTERM": 1.0}
# Atom on which the titratable charge is placed
TITR_ATOM = {"ASP": "CG", "GLU": "CD", "TYR": "OH", "ARG": "CZ", "HIS": "NE2", "LYS": "NZ"}

# Free amino-acid molecular weights; one water is removed per peptide bond.
MW_DICT = {"ALA": 89.09, "ARG": 174.20, "ASN": 132.12, "ASP": 133.10, "ASX": 132.61,
           "CYS": 121.15, "GLN": 146.15, "GLU": 147.13, "GLX": 146.64, "GLY": 75.07,
           "HIS": 155.16, "ILE": 131.17, "LEU": 131.17, "LYS": 146.19, "MET": 149.21,
           "PHE": 165.19, "PRO": 115.13, "SER": 105.09, "THR": 119.12, "TRP": 204.23,
           "TYR": 181.19, "VAL": 117.15, "IVA": 85.12, "STA": 157.15, "ACE": 43.30}
MW_H2O = 18.0

# van der Waals radii by atom name (Harms-lab table, used by water contacts)
VDW_DICT = {
    "N": 1.650, "CA": 1.870, "C": 1.760, "O": 1.400, "AD1": 1.700, "AD2": 1.700,
    "AE1": 1.700, "AE2": 1.700, "CB": 1.870, "CD": 1.870, "CD1": 1.760, "CD2": 1.760,
    "CG": 1.870, "CG1": 1.870, "CG2": 1.870, "CE": 1.870, "CE1": 1.760, "CE2": 1.760,
    "CE3": 1.760, "CH2": 1.760, "CZ": 1.870, "CZ2": 1.760, "CZ3": 1.760, "ND": 1.650,
    "ND1": 1.650, "ND2": 1.650, "NE": 1.650, "NE1": 1.650, "NE2": 1.650, "NH1": 1.650,
    "NH2": 1.650, "NZ": 1.650, "OD": 1.400, "OD1": 1.400, "OD2": 1.400, "OT2": 1.400,
    "OE": 1.400, "OE1": 1.400, "OE2": 1.400, "OH2": 1.400, "OG": 1.400, "OG1": 1.400,
    "OT1": 1.400, "OXT": 1.400, "OH": 1.400, "SD": 1.850, "SG": 1.850, "P": 1.900,
}
