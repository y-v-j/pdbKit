"""Structure analyses (ported from the Harms-lab ``pdbtools`` and HADDOCK ``pdb-tools``).

All functions take a Biopython Structure (so PDB and mmCIF input behave the
same) and return plain Python data - usually a list of dicts, one per row -
which :func:`to_dataframe` turns into a pandas DataFrame when pandas is
installed.

Ported calculations keep the original algorithms and default parameters so
their numbers can be compared; the few places where the original was wrong
are fixed and marked "differs from the original".
"""

from __future__ import annotations

import math
from collections import OrderedDict
from typing import Dict, List, Optional, Sequence, Tuple

import numpy as np
from Bio.PDB.NeighborSearch import NeighborSearch
from Bio.PDB.Structure import Structure

from .constants import (AA3_TO_AA1, BORING_LIGANDS, CHARGE_DICT, MW_DICT, MW_H2O, NUC_TO_1,
                        PKA_DICT, TITR_ATOM, VDW_DICT)
from .metadata import experiment as _experiment, seqres as _seqres, assemblies as _assemblies
from .select import _as_list
from .structure import PdbKitError, iter_atoms_all_altlocs, residue_kind, structure_id

__all__ = ["summary", "sequences", "to_fasta", "fasta", "gaps", "ligands", "disulfides", "contacts",
           "close_contacts", "water_contacts", "neighbors", "torsions", "param", "coulomb",
           "ion_distance_histogram", "distances", "sasa", "experiment", "assemblies", "to_dataframe"]


def to_dataframe(rows):
    """List of dicts -> pandas DataFrame (pandas is optional)."""
    import pandas as pd
    return pd.DataFrame(list(rows))


def _first_model(structure: Structure):
    return next(iter(structure))


def _res_label(res) -> str:
    chain = res.get_parent().id
    return f"{chain}:{res.get_resname().strip()}:{res.id[1]}{res.id[2].strip()}"


# --------------------------------------------------------------------------- #
# Counts (pdb_wc)
# --------------------------------------------------------------------------- #

def summary(structure: Structure) -> Dict[str, object]:
    """Counts of models, chains, residues (by kind), atoms, altlocs, insertion codes, gaps."""
    models = list(structure)
    m0 = models[0]
    kinds = {"polymer": 0, "ligand": 0, "water": 0}
    for r in m0.get_residues():
        kinds[residue_kind(r)] += 1
    atoms = list(iter_atoms_all_altlocs(m0))
    ligand_names = OrderedDict()
    for r in m0.get_residues():
        if residue_kind(r) == "ligand":
            ligand_names[r.get_resname().strip()] = ligand_names.get(r.get_resname().strip(), 0) + 1
    return {
        "id": structure_id(structure),
        "models": len(models),
        "chains": [c.id for c in m0],
        "polymer_residues": kinds["polymer"],
        "ligand_residues": kinds["ligand"],
        "water_molecules": kinds["water"],
        "atoms": len(atoms),
        "hetatm_atoms": sum(1 for a in atoms if a.get_parent().id[0] != " "),
        "has_altlocs": any(a.altloc.strip() for a in atoms),
        "has_insertion_codes": any(r.id[2].strip() for r in m0.get_residues()),
        "has_sequence_gaps": bool(gaps(structure)),
        "ligands": dict(ligand_names),
    }


# --------------------------------------------------------------------------- #
# Sequences (pdb_tofasta, Harms pdb_seq)
# --------------------------------------------------------------------------- #

def _one_letter(name: str) -> str:
    name = name.strip()
    return AA3_TO_AA1.get(name) or NUC_TO_1.get(name) or "X"


# Common modified nucleotides -> parent (used when the file itself gives no parent)
MODIFIED_NUCLEOTIDE_PARENT = {
    "PSU": "U", "5MU": "U", "H2U": "U", "4SU": "U", "OMU": "U", "5BU": "U", "UR3": "U", "2MU": "U",
    "5MC": "C", "OMC": "C", "4OC": "C", "CBR": "C", "5CM": "DC",
    "2MG": "G", "M2G": "G", "7MG": "G", "OMG": "G", "1MG": "G", "YG": "G", "QUO": "G",
    "1MA": "A", "6MA": "A", "MA6": "A", "2MA": "A", "6IA": "A", "T6A": "A", "MIA": "A",
}


def _parent_map(structure: Structure) -> Dict[str, str]:
    """Modified residue -> parent residue, from the file (MODRES, _pdbx_struct_mod_residue,
    _chem_comp.mon_nstd_parent_comp_id), then a built-in table of modified nucleotides."""
    from .metadata import as_mmcif_dict, column
    d = as_mmcif_dict(structure)
    parent = dict(MODIFIED_NUCLEOTIDE_PARENT)
    for cid, par in zip(column(d, "_chem_comp.id"), column(d, "_chem_comp.mon_nstd_parent_comp_id")):
        if par not in ("?", ".", "") and "," not in par:
            parent[cid] = par
    parent.update(zip(column(d, "_pdbx_struct_mod_residue.auth_comp_id") or column(d, "_pdbx_struct_mod_residue.label_comp_id"),
                      column(d, "_pdbx_struct_mod_residue.parent_comp_id")))
    return parent


def _canonical_seqres(structure: Structure) -> Dict[str, str]:
    """mmCIF's own canonical one-letter sequence per author chain (what RCSB's FASTA shows)."""
    from .metadata import column, source_info
    d = source_info(structure).mmcif or {}
    out: Dict[str, str] = {}
    for strands, seq in zip(column(d, "_entity_poly.pdbx_strand_id"),
                            column(d, "_entity_poly.pdbx_seq_one_letter_code_can")):
        if seq in ("?", ".", "") or strands in ("?", "."):
            continue
        clean = "".join(seq.split())
        for c in strands.split(","):
            if c.strip():
                out[c.strip()] = clean
    return out


def sequences(structure: Structure, source: str = "atom", chains: Optional[Sequence[str]] = None) -> Dict[str, str]:
    """One-letter sequence per chain.

    ``source="atom"`` - observed polymer residues (modified residues such as
    MSE included, ligands and water excluded); ``"seqres"`` - the full
    deposited sequence, including unobserved residues. For mmCIF this is the
    file's own canonical sequence (``_entity_poly.pdbx_seq_one_letter_code_can``,
    identical to RCSB's FASTA); for PDB files, SEQRES. Modified residues are
    written as their parent (MODRES, ``_pdbx_struct_mod_residue``,
    ``_chem_comp.mon_nstd_parent_comp_id``, common modified nucleotides),
    otherwise ``X``.
    """
    parent = _parent_map(structure)
    wanted = _as_list(chains)
    out: Dict[str, str] = OrderedDict()
    if source == "seqres":
        canonical = _canonical_seqres(structure)
        for c, names in _seqres(structure).items():
            if wanted is None or c in wanted:
                out[c] = canonical.get(c) or "".join(_one_letter(parent.get(n, n)) for n in names)
        if not out:
            raise PdbKitError("No SEQRES / _pdbx_poly_seq_scheme records in this file; use source='atom'.")
        return out
    for chain in _first_model(structure):
        if wanted is not None and chain.id not in wanted:
            continue
        seq = [_one_letter(parent.get(r.get_resname().strip(), r.get_resname().strip()))
               for r in chain if residue_kind(r) == "polymer"]
        if seq:
            out[chain.id] = "".join(seq)
    return out


def to_fasta(seqs: Dict[str, str], name: str = "", width: int = 60,
             descriptions: Optional[Dict[str, str]] = None) -> str:
    """FASTA text; headers are ``>NAME_CHAIN description`` (or ``>CHAIN`` without a name)."""
    lines = []
    for c, s in seqs.items():
        header = f">{name}_{c}" if name else f">{c}"
        if descriptions and descriptions.get(c):
            header += f" {descriptions[c]}"
        lines.append(header)
        lines.extend(s[i:i + width] for i in range(0, len(s), width))
    return "\n".join(lines) + "\n"


def _chain_moltype(structure: Structure) -> Dict[str, str]:
    out = {}
    for chain in _first_model(structure):
        kinds = {r.xtra.get("pdbkit_poly") for r in chain if residue_kind(r) == "polymer"}
        kinds.discard(None)
        if kinds == {"nucleic"}:
            names = [r.get_resname().strip() for r in chain if residue_kind(r) == "polymer"]
            out[chain.id] = "dna" if all(n.startswith("D") for n in names) else "rna"
        elif kinds:
            out[chain.id] = "protein" if kinds == {"protein"} else "mixed"
    return out


def fasta(structure: Structure, source: str = "auto", chains: Optional[Sequence[str]] = None,
          width: int = 60) -> str:
    """FASTA text for every polymer chain (protein, DNA, RNA).

    ``source``: ``"seqres"`` - full deposited sequence; ``"atom"`` - residues
    present in the coordinates; ``"auto"`` (default) - SEQRES when the file
    has it, otherwise the coordinates. Headers look like
    ``>4MMK_A mol:protein length:292 <entity description>``.
    """
    from .metadata import chain_descriptions
    if source == "auto":
        try:
            seqs = sequences(structure, "seqres", chains)
            used = "seqres"
        except PdbKitError:
            seqs, used = sequences(structure, "atom", chains), "atom"
    else:
        seqs, used = sequences(structure, source, chains), source
    moltype = _chain_moltype(structure)
    names = chain_descriptions(structure)
    desc = {}
    for c, s in seqs.items():
        parts = [f"mol:{moltype.get(c, 'unknown')}", f"length:{len(s)}", f"source:{used}"]
        if names.get(c):
            parts.append(names[c])
        desc[c] = " ".join(parts)
    return to_fasta(seqs, structure_id(structure), width, desc)


# --------------------------------------------------------------------------- #
# Gaps (pdb_gap)
# --------------------------------------------------------------------------- #

def gaps(structure: Structure, cutoff: float = 4.0) -> List[Dict[str, object]]:
    """Chain breaks between consecutive protein residues, and numbering jumps.

    A break is a missing peptide bond: C(i)-N(i+1) > 2.0 A when both atoms
    exist, else CA-CA > ``cutoff`` (the pdb_gap criterion). ``distance`` is
    the CA-CA distance in A.

    Differs from the original pdb_gap, which (a) printed the *squared*
    distance labelled as A, (b) reported every alternate location of a CA as
    a gap, and (c) flagged modified residues whose "CA" is not an alpha
    carbon; this version is altloc-aware and uses the peptide bond.
    """
    out = []
    for model in structure:
        mnum = getattr(model, "serial_num", None) or model.id + 1
        for chain in model:
            prev = None
            for res in chain:
                if residue_kind(res) != "polymer" or "CA" not in res:
                    continue
                if prev is not None:
                    d = float(np.linalg.norm(res["CA"].coord - prev["CA"].coord))
                    if "C" in prev and "N" in res:
                        broken = float(np.linalg.norm(prev["C"].coord - res["N"].coord)) > 2.0
                    else:
                        broken = d > cutoff
                    if broken:
                        out.append({"model": mnum, "chain": chain.id, "from": _res_label(prev),
                                    "to": _res_label(res), "type": "distance", "distance": round(d, 2)})
                    elif prev.id[1] + 1 != res.id[1] and not (prev.id[1] == res.id[1] and res.id[2] != prev.id[2]) \
                            and not res.id[2].strip():
                        out.append({"model": mnum, "chain": chain.id, "from": _res_label(prev),
                                    "to": _res_label(res), "type": "numbering", "distance": round(d, 2)})
                prev = res
    return out


# --------------------------------------------------------------------------- #
# Ligands (Harms pdb_ligand)
# --------------------------------------------------------------------------- #

def ligands(structure: Structure, skip_boring: bool = False) -> List[Dict[str, object]]:
    """Ligand instances in the first model (from coordinates, so mmCIF works too)."""
    out = []
    for chain in _first_model(structure):
        for res in chain:
            if residue_kind(res) != "ligand":
                continue
            name = res.get_resname().strip()
            if skip_boring and name in BORING_LIGANDS:
                continue
            out.append({"chain": chain.id, "resname": name, "resseq": res.id[1], "icode": res.id[2].strip(),
                        "atoms": len(res.get_unpacked_list())})
    return out


# --------------------------------------------------------------------------- #
# Disulfides (Harms pdb_disulfide)
# --------------------------------------------------------------------------- #

def disulfides(structure: Structure, cutoff: float = 3.0) -> Dict[str, list]:
    """Cys-Cys pairs with SG-SG distance below ``cutoff`` (A) and the free cysteines.

    An ideal S-S bond is ~2.05 A, but deposited disulfides in poorly
    restrained regions reach ~2.7-2.9 A (e.g. 1K4C), while non-bonded SG-SG
    contacts start near 3.4-3.6 A; 3.0 A separates the two. Harms
    ``pdb_disulfide`` used 4.0 A (pass ``cutoff=4.0`` for the original).
    """
    sgs = [r["SG"] for r in _first_model(structure).get_residues()
           if r.get_resname().strip() in ("CYS", "CYX", "CYM") and "SG" in r]
    bonded, pairs = set(), []
    for i in range(len(sgs)):
        for j in range(i + 1, len(sgs)):
            d = float(np.linalg.norm(sgs[i].coord - sgs[j].coord))
            if d < cutoff:
                a, b = sgs[i].get_parent(), sgs[j].get_parent()
                pairs.append({"cys1": _res_label(a), "cys2": _res_label(b), "distance": round(d, 3)})
                bonded.update((id(a), id(b)))
    free = [_res_label(s.get_parent()) for s in sgs if id(s.get_parent()) not in bonded]
    return {"bonds": pairs, "free": free}


# --------------------------------------------------------------------------- #
# Contacts (Harms pdb_contact, pdb_closecontacts, pdb_watercontact, pdb_neighbors)
# --------------------------------------------------------------------------- #

def contacts(structure: Structure, target: str = "TRP", cutoff: float = 3.5,
             include_water: bool = False) -> List[Dict[str, object]]:
    """Residues within ``cutoff`` A (any atom pair, any alternate location) of each target residue.

    ``target`` is a residue name (``"TRP"``, ``"HEM"``), a residue number
    (``"42"``) or chain-qualified number (``"A:42"``).
    """
    model = _first_model(structure)
    atoms = [a for a in iter_atoms_all_altlocs(model) if include_water or residue_kind(a.get_parent()) != "water"]
    ns = NeighborSearch(atoms) if atoms else None
    t = target.strip()
    if ":" in t:
        tc, tn = t.split(":", 1)
        targets = [r for r in model[tc] if str(r.id[1]) == tn] if tc in model else []
    elif t.lstrip("-").isdigit():
        targets = [r for r in model.get_residues() if r.id[1] == int(t)]
    else:
        targets = [r for r in model.get_residues() if r.get_resname().strip() == t.upper()]
    if not targets:
        raise PdbKitError(f"No residue matches target {target!r}.")
    out = []
    for tres in targets:
        best: Dict[int, Tuple[object, float]] = {}
        for a in (tres.get_unpacked_list() if ns is not None else []):
            for b in ns.search(a.coord, cutoff):
                r = b.get_parent()
                if r is tres:
                    continue
                d = float(np.linalg.norm(a.coord - b.coord))
                if id(r) not in best or d < best[id(r)][1]:
                    best[id(r)] = (r, d)
        for r, d in sorted(best.values(), key=lambda x: (x[0].get_parent().id, x[0].id[1])):
            out.append({"target": _res_label(tres), "residue": _res_label(r), "min_distance": round(d, 3)})
    return out


def close_contacts(structure: Structure, cutoff: float = 3.5, ca_only: bool = False) -> List[Dict[str, object]]:
    """Atom pairs from different polymer residues closer than ``cutoff`` A.

    Pairs across consecutive residues (i, i+1) are skipped when atom i's name
    starts with C and atom i+1's with N, O with N, or C with CA - Harms
    ``pdb_closecontacts``' rule for the bonded/geometrically fixed backbone
    neighbours.
    """
    model = _first_model(structure)
    atoms = [a for r in model.get_residues() if residue_kind(r) == "polymer" for a in r.get_atoms()
             if not ca_only or a.get_id() == "CA"]
    order = {id(a): i for i, a in enumerate(atoms)}
    out = []
    if len(atoms) < 2:          # e.g. --ca-only on a nucleic acid: nothing to pair
        return out
    for a, b in NeighborSearch(atoms).search_all(cutoff):
        if order[id(a)] > order[id(b)]:
            a, b = b, a
        ra, rb = a.get_parent(), b.get_parent()
        if ra is rb:
            continue
        na, nb = a.get_id(), b.get_id()
        if ra.get_parent() is rb.get_parent() and rb.id[1] == ra.id[1] + 1 and (
                (na[0] == "C" and nb[0] == "N") or (na[0] == "O" and nb[0] == "N") or
                (na[0] == "C" and nb[:2] == "CA")):
            continue
        out.append({"atom1": f"{_res_label(ra)}:{na}", "atom2": f"{_res_label(rb)}:{nb}",
                    "distance": round(float(np.linalg.norm(a.coord - b.coord)), 3)})
    out.sort(key=lambda r: (r["atom1"], r["atom2"]))
    return out


_ELEMENT_VDW = {"C": 1.76, "N": 1.65, "O": 1.40, "S": 1.85, "P": 1.90, "H": 1.10}


def water_contacts(structure: Structure, probe: float = 1.4) -> Dict[str, List[str]]:
    """Residue neighbour lists: residues touch if any atoms are closer than vdW1 + vdW2 + 2*probe.

    Differs from the original: atoms missing from the Harms vdW table
    (hydrogens, ligand atoms) raised KeyError there; here they fall back to an
    element radius.
    """
    model = _first_model(structure)
    atoms = [a for r in model.get_residues() if residue_kind(r) == "polymer" for a in r.get_atoms()]
    radius = {id(a): VDW_DICT.get(a.get_id(), _ELEMENT_VDW.get((a.element or "C").upper(), 1.8)) for a in atoms}
    max_r = max(radius.values()) if radius else 2.0
    neighbours: Dict[str, List[str]] = OrderedDict((_res_label(r), []) for r in model.get_residues()
                                                   if residue_kind(r) == "polymer")
    for a, b in (NeighborSearch(atoms).search_all(2 * max_r + 2 * probe) if len(atoms) > 1 else []):
        ra, rb = a.get_parent(), b.get_parent()
        if ra is rb:
            continue
        if np.linalg.norm(a.coord - b.coord) <= radius[id(a)] + radius[id(b)] + 2 * probe:
            la, lb = _res_label(ra), _res_label(rb)
            if lb not in neighbours[la]:
                neighbours[la].append(lb)
                neighbours[lb].append(la)
    return neighbours


def _gly_cbeta(n, ca, c) -> np.ndarray:
    """Harms-lab pseudo-CB for glycine: CA->C rotated 240 deg about CA->N."""
    v = c - ca
    axis = (n - ca) / np.linalg.norm(n - ca)
    theta = math.radians(240.0)
    x, y, z = axis
    ct, st, oc = math.cos(theta), math.sin(theta), 1 - math.cos(theta)
    mat = np.array([[ct + x * x * oc, x * y * oc - z * st, x * z * oc + y * st],
                    [x * y * oc + z * st, ct + y * y * oc, y * z * oc - x * st],
                    [x * z * oc - y * st, y * z * oc + x * st, ct + z * z * oc]])
    return v @ mat + ca     # row-vector convention, as in Harms geometry.arbRot


def neighbors(structure: Structure, cutoff: float = 10.0, seq_cutoff: int = 4) -> List[Dict[str, object]]:
    """Number of CB atoms within ``cutoff`` A of each residue's CB (pseudo-CB for Gly).

    Residues in the same chain closer than ``seq_cutoff`` in sequence are not
    counted, as in Harms ``pdb_neighbors``.
    """
    cbs = []
    for chain in _first_model(structure):
        for r in chain:
            if residue_kind(r) != "polymer":
                continue
            if r.get_resname().strip() == "GLY":
                if all(x in r for x in ("N", "CA", "C")):
                    cbs.append((r, _gly_cbeta(r["N"].coord.astype(float), r["CA"].coord.astype(float),
                                              r["C"].coord.astype(float))))
            elif "CB" in r:
                cbs.append((r, r["CB"].coord.astype(float)))
    xyz = np.array([c for _, c in cbs]) if cbs else np.zeros((0, 3))
    out = []
    if not cbs:
        return out
    from Bio.PDB.kdtrees import KDTree
    tree = KDTree(np.ascontiguousarray(xyz, dtype="d"), 10)
    chains = [r.get_parent().id for r, _ in cbs]
    nums = [r.id[1] for r, _ in cbs]
    for i, (r, c) in enumerate(cbs):
        n = 0
        for p in tree.search(np.ascontiguousarray(c, dtype="d"), float(cutoff) + 1e-6):
            j = p.index
            if chains[j] == chains[i] and abs(nums[j] - nums[i]) < seq_cutoff:
                continue
            if np.linalg.norm(xyz[j] - c) <= cutoff:
                n += 1
        out.append({"residue": _res_label(r), "neighbors": n})
    return out


# --------------------------------------------------------------------------- #
# Backbone torsions (Harms pdb_torsion)
# --------------------------------------------------------------------------- #

def torsions(structure: Structure) -> List[Dict[str, object]]:
    """phi/psi (degrees) per protein residue, IUPAC sign convention.

    Differs from the original: phi/psi are only computed across genuine
    peptide bonds (C-N <= 2.0 A); pdb_torsion also paired residues across
    chain breaks and between chains.
    """
    from Bio.PDB.vectors import calc_dihedral
    out = []
    for chain in _first_model(structure):
        res = [r for r in chain if residue_kind(r) == "polymer" and all(x in r for x in ("N", "CA", "C"))]
        for i, r in enumerate(res):
            prev = res[i - 1] if i > 0 else None
            nxt = res[i + 1] if i + 1 < len(res) else None
            phi = psi = None
            if prev is not None and np.linalg.norm(prev["C"].coord - r["N"].coord) <= 2.0:
                phi = math.degrees(calc_dihedral(prev["C"].get_vector(), r["N"].get_vector(),
                                                 r["CA"].get_vector(), r["C"].get_vector()))
            if nxt is not None and np.linalg.norm(r["C"].coord - nxt["N"].coord) <= 2.0:
                psi = math.degrees(calc_dihedral(r["N"].get_vector(), r["CA"].get_vector(),
                                                 r["C"].get_vector(), nxt["N"].get_vector()))
            out.append({"residue": _res_label(r), "phi": None if phi is None else round(phi, 2),
                        "psi": None if psi is None else round(psi, 2)})
    return out


# --------------------------------------------------------------------------- #
# Sequence parameters (Harms pdb_param)
# --------------------------------------------------------------------------- #

def _charge_state(counts: Dict[str, int], ph: float) -> float:
    q = 0.0
    for aa, n in counts.items():
        pka, charge = PKA_DICT[aa], CHARGE_DICT[aa]
        q += n * charge / (1 + 10 ** (charge * (ph - pka)))
    return q


def _pi(names: List[str], step: float = 2.0, cutoff: float = 0.001) -> float:
    """Isoelectric point by the Harms bisection (model-compound pKa values)."""
    counts = {k: sum(1 for n in names if n == k) for k in PKA_DICT}
    counts["NTERM"] = 1
    counts["CTERM"] = 1
    ph = 0.0
    while step > cutoff:
        if _charge_state(counts, ph) > 0:
            ph += step
        else:
            step *= 0.5
            ph -= step
    return ph


def param(structure: Structure, chains: Optional[Sequence[str]] = None, source: str = "atom") -> Dict[str, object]:
    """Molecular weight, pI, composition and fully-ionised charge (Harms ``pdb_param``).

    The original ``pdb_param`` could not run (it called an unimported module);
    this port uses its formulas: MW = sum of free amino-acid masses minus one
    water per peptide bond; pI by bisection on Henderson-Hasselbalch charges.
    """
    from .metadata import seqres as _sr
    if source == "seqres":
        per_chain = _sr(structure)
    else:
        per_chain = {c.id: [r.get_resname().strip() for r in c if residue_kind(r) == "polymer"
                            and r.xtra.get("pdbkit_poly") != "nucleic"]
                     for c in _first_model(structure)}
    wanted = _as_list(chains)
    names = [n for c, seq in per_chain.items() if wanted is None or c in wanted for n in seq]
    names = ["MET" if n == "MSE" else n for n in names]
    if not names:
        raise PdbKitError("No protein residues selected.")
    unknown = sorted({n for n in names if n not in MW_DICT})
    mw = sum(MW_DICT.get(n, 0.0) for n in names) - (len(names) - 1) * MW_H2O
    comp = {k: names.count(k) for k in sorted(set(names))}
    titr = [n for n in names if n in ("ASP", "GLU", "HIS", "LYS", "ARG", "TYR", "CYS")]
    acid = [n for n in names if n in ("ASP", "GLU")]
    base = [n for n in names if n in ("HIS", "LYS", "ARG")]
    return {
        "residues": len(names),
        "molecular_weight": round(mw, 2),
        "pI": round(_pi(names), 3),
        "fraction_titratable": round(len(titr) / len(names), 4),
        "fraction_acidic": round(len(acid) / len(names), 4),
        "fraction_basic": round(len(base) / len(names), 4),
        "charge_fully_ionised": len(base) - len(acid),
        "composition": comp,
        "non_standard_residues": unknown,
    }


# --------------------------------------------------------------------------- #
# Electrostatics (Harms pdb_coulomb, pdb_iondist)
# --------------------------------------------------------------------------- #

def _titratable(structure: Structure, skip: Sequence[str] = ()):
    coords, pkas, charges, labels = [], [], [], []
    for r in _first_model(structure).get_residues():
        name = r.get_resname().strip()
        if residue_kind(r) != "polymer" or name not in TITR_ATOM or name in skip:
            continue
        atom = TITR_ATOM[name]
        if atom in r:
            coords.append(r[atom].coord.astype(float))
            pkas.append(PKA_DICT[name])
            charges.append(CHARGE_DICT[name])
            labels.append(_res_label(r))
    return np.array(coords).reshape(-1, 3), pkas, charges, labels


def coulomb(structure: Structure, dielectric: float = 40.0, ionic_strength: float = 100.0,
            ph: float = 7.0, temperature: float = 298.0) -> float:
    """Screened Coulomb energy (kcal/mol) of titratable groups at model-compound pKa values.

    ``ionic_strength`` is in mM. Same formula and defaults as Harms ``pdb_coulomb``.
    """
    xyz, pkas, charges, _ = _titratable(structure)
    kappa = 50.29 * math.sqrt((ionic_strength / 1000.0) / (dielectric * temperature))
    q = np.array([c / (1 + 10 ** (c * (ph - p))) for p, c in zip(pkas, charges)])
    e = 0.0
    for start in range(0, len(q), 1000):            # row blocks keep memory bounded
        block = xyz[start:start + 1000]
        r = np.linalg.norm(block[:, None, :] - xyz[None, :, :], axis=2)
        i = np.arange(start, start + len(block))[:, None]
        j = np.arange(len(q))[None, :]
        mask = j > i
        with np.errstate(divide="ignore", invalid="ignore"):
            terms = 332 * q[start:start + len(block), None] * q[None, :] / (r * dielectric) * np.exp(-kappa * r)
        e += float(terms[mask].sum())
    return e


def ion_distance_histogram(structure: Structure, step: float = 0.1, skip: Sequence[str] = ("TYR",)) -> List[Dict[str, object]]:
    """Histogram of distances between titratable atoms: acid/acid, acid/base, base/base (Harms ``pdb_iondist``)."""
    xyz, _, charges, _ = _titratable(structure, skip)
    n = len(charges)
    if n < 2:
        return []
    ch = np.array(charges)
    dmax, kinds_bins = 0.0, []
    for start in range(0, n, 1000):
        block = xyz[start:start + 1000]
        d = np.linalg.norm(block[:, None, :] - xyz[None, :, :], axis=2)
        i = np.arange(start, start + len(block))[:, None]
        mask = np.arange(n)[None, :] > i
        dmax = max(dmax, float(d.max()))
        kind = ((ch[start:start + len(block), None] + ch[None, :]) / 2 + 1).astype(int)
        kinds_bins.append((kind[mask], np.rint(d[mask] / step).astype(int)))
    nbins = int(round(dmax / step)) + 1
    hist = np.zeros((3, nbins), dtype=int)
    for kind, b in kinds_bins:
        np.add.at(hist, (kind, b), 1)
    return [{"distance": round(k * step, 6), "acid_acid": int(hist[0, k]), "acid_base": int(hist[1, k]),
             "base_base": int(hist[2, k])} for k in range(nbins)]


def distances(structure: Structure, reference: str, atom_name: Optional[str] = None,
              chains: Optional[Sequence[str]] = None) -> List[Dict[str, object]]:
    """Distances from a reference atom ``"A:42:NZ"`` (or ``"42:N"``) to every atom named ``atom_name``.

    Generalises Harms ``pdb_dist_filter`` (distance from one residue's atom to
    the same atom type elsewhere).
    """
    parts = reference.split(":")
    if len(parts) == 3:
        rc, rn, ra = parts
    elif len(parts) == 2:
        rc, (rn, ra) = None, parts
    else:
        raise PdbKitError("Reference must be CHAIN:RESNUM:ATOM or RESNUM:ATOM")
    model = _first_model(structure)
    ref = None
    for r in model.get_residues():
        if str(r.id[1]) == rn and (rc is None or r.get_parent().id == rc) and ra in r:
            ref = r[ra]
            break
    if ref is None:
        raise PdbKitError(f"Reference atom {reference!r} not found.")
    name = atom_name or ra
    wanted = _as_list(chains)
    out = []
    for r in model.get_residues():
        if wanted is not None and r.get_parent().id not in wanted:
            continue
        if name in r and r[name] is not ref:
            out.append({"atom": f"{_res_label(r)}:{name}",
                        "distance": round(float(np.linalg.norm(r[name].coord - ref.coord)), 3)})
    return out


# --------------------------------------------------------------------------- #
# Solvent accessibility (replaces Harms pdb_sasa, which needed NACCESS)
# --------------------------------------------------------------------------- #

# Tien et al. 2013 (theoretical) maximum residue SASA, A^2
MAX_ASA = {"ALA": 129.0, "ARG": 274.0, "ASN": 195.0, "ASP": 193.0, "CYS": 167.0, "GLN": 225.0,
           "GLU": 223.0, "GLY": 104.0, "HIS": 224.0, "ILE": 197.0, "LEU": 201.0, "LYS": 236.0,
           "MET": 224.0, "PHE": 240.0, "PRO": 159.0, "SER": 155.0, "THR": 172.0, "TRP": 285.0,
           "TYR": 263.0, "VAL": 174.0}


def sasa(structure: Structure, level: str = "residue", probe: float = 1.4, n_points: int = 100,
         include_het: bool = False) -> List[Dict[str, object]]:
    """Shrake-Rupley solvent-accessible surface area (A^2), per residue or per atom.

    Waters are always removed first; ligands only count when ``include_het``.
    """
    from Bio.PDB.SASA import ShrakeRupley
    from .select import extract
    s = extract(structure, het="all" if include_het else "none", water=False, models=[
        getattr(_first_model(structure), "serial_num", None) or 1])
    model = _first_model(s)
    ShrakeRupley(probe_radius=probe, n_points=n_points).compute(model, level="A")
    out = []
    for r in model.get_residues():
        if level == "atom":
            for a in r.get_atoms():
                out.append({"residue": _res_label(r), "atom": a.get_id(), "sasa": round(float(a.sasa), 2)})
        else:
            total = float(sum(a.sasa for a in r.get_atoms()))
            name = r.get_resname().strip()
            rel = round(total / MAX_ASA[name], 3) if name in MAX_ASA else None
            out.append({"residue": _res_label(r), "sasa": round(total, 2), "relative": rel})
    return out


# --------------------------------------------------------------------------- #
# Header-derived reports (Harms pdb_exper, pdb_oligomer)
# --------------------------------------------------------------------------- #

def experiment(structure: Structure) -> Dict[str, str]:
    return _experiment(structure)


def assemblies(structure: Structure) -> List[Dict[str, object]]:
    return _assemblies(structure)
