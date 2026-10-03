"""Centre of mass, recentring and other rigid-body geometry.

``center_of_mass`` / ``center`` implement pdbKit feature 5 and replace
Harms-lab ``pdb_centermass``: the reference can be the whole structure or any
selection (chains, C-alpha only, with or without ligands/water), and the
translation can be applied to the selection alone or to the whole structure.

Masses come from the element column (Biopython infers the element from the
atom name when the column is blank). Unknown elements get the mass of carbon,
as in the original, and are reported.
"""

from __future__ import annotations

import sys
from typing import Dict, List, Optional, Tuple

import numpy as np
from Bio.PDB.Structure import Structure

from .constants import ATOM_WEIGHTS, CARBON_MASS
from .select import Selection, extract
from .structure import PdbKitError, deep_copy, iter_atoms_all_altlocs

__all__ = ["atom_masses", "center_of_mass", "center", "translate", "radius_of_gyration",
           "center_asu", "coordinates"]


def coordinates(entity) -> np.ndarray:
    """(N, 3) array of coordinates (one location per atom)."""
    return np.array([a.coord for a in entity.get_atoms()], dtype=float).reshape(-1, 3)


def atom_masses(atoms, warn: bool = True) -> np.ndarray:
    masses, unknown = [], {}
    for a in atoms:
        el = (a.element or "").strip().upper()
        m = ATOM_WEIGHTS.get(el)
        if m is None:
            unknown[el or "?"] = unknown.get(el or "?", 0) + 1
            m = CARBON_MASS
        masses.append(m)
    if unknown and warn:
        print("pdbKit: unknown element(s) given the mass of carbon: " +
              ", ".join(f"{k} x{v}" for k, v in unknown.items()), file=sys.stderr)
    return np.array(masses, dtype=float)


def _first_model(structure: Structure):
    return next(iter(structure))


def center_of_mass(structure: Structure, selection: Optional[Selection] = None, *,
                   weighted: bool = True, model: Optional[int] = None, **kwargs) -> np.ndarray:
    """Centre of mass (``weighted=True``) or geometric centre of a selection.

    By default the whole first model is used. Selection keywords are those of
    :class:`pdbkit.select.Selection`, e.g. ``chains="A"``, ``atoms="ca"``,
    ``het="none"``, ``water=False``.
    """
    sub = extract(structure, selection, **kwargs) if (selection is not None or kwargs) else structure
    if model is not None:
        target = next((m for m in sub if (getattr(m, "serial_num", None) or m.id + 1) == model), None)
        if target is None:
            raise PdbKitError(f"Model {model} not found.")
    else:
        target = _first_model(sub)
    atoms = list(target.get_atoms())
    if not atoms:
        raise PdbKitError("No atoms in the selection.")
    xyz = np.array([a.coord for a in atoms], dtype=float)
    if not weighted:
        return xyz.mean(axis=0)
    w = atom_masses(atoms)
    return (xyz * w[:, None]).sum(axis=0) / w.sum()


def translate(structure: Structure, vector, *, inplace: bool = False) -> Structure:
    """Shift every atom (all models, all alternate locations) by ``vector``."""
    s = structure if inplace else deep_copy(structure)
    v = np.asarray(vector, dtype=float)
    for atom in iter_atoms_all_altlocs(s):
        atom.coord = atom.coord.astype(float) + v
    return s


def center(structure: Structure, selection: Optional[Selection] = None, *, weighted: bool = True,
           apply_to: str = "selection", per_model: bool = False, **kwargs) -> Tuple[Structure, np.ndarray]:
    """Move the centre of mass of a selection to the origin.

    Parameters
    ----------
    apply_to:
        ``"selection"`` - output only the selected atoms, recentred (default);
        ``"all"`` - output the whole structure, translated by the selection's
        centre of mass (e.g. put chain A's COM at the origin, keep everything).
    per_model:
        Centre each model separately. By default every model is shifted by the
        first model's centre, which preserves an NMR ensemble's superposition.

    Returns ``(new_structure, com_before)``.
    """
    has_sel = selection is not None or bool(kwargs)
    sub = extract(structure, selection, **kwargs) if has_sel else deep_copy(structure)
    base = sub if apply_to == "selection" else deep_copy(structure)
    if apply_to not in ("selection", "all"):
        raise PdbKitError("apply_to must be 'selection' or 'all'")
    com = center_of_mass(sub, weighted=weighted)
    if not per_model:
        return translate(base, -com, inplace=True), com
    for model in base:
        serial = getattr(model, "serial_num", None) or model.id + 1
        ref = next((m for m in sub if (getattr(m, "serial_num", None) or m.id + 1) == serial), None)
        if ref is None:
            continue
        atoms = list(ref.get_atoms())
        xyz = np.array([a.coord for a in atoms], dtype=float)
        w = atom_masses(atoms, warn=False) if weighted else np.ones(len(atoms))
        c = (xyz * w[:, None]).sum(0) / w.sum()
        for atom in iter_atoms_all_altlocs(model):
            atom.coord = atom.coord.astype(float) - c
    return base, com


def radius_of_gyration(structure: Structure, selection: Optional[Selection] = None,
                       weighted: bool = True, **kwargs) -> float:
    sub = extract(structure, selection, **kwargs) if (selection is not None or kwargs) else structure
    atoms = list(_first_model(sub).get_atoms())
    xyz = np.array([a.coord for a in atoms], dtype=float)
    w = atom_masses(atoms, warn=False) if weighted else np.ones(len(atoms))
    c = (xyz * w[:, None]).sum(0) / w.sum()
    return float(np.sqrt((w * ((xyz - c) ** 2).sum(1)).sum() / w.sum()))


def center_asu(structure: Structure) -> Tuple[Structure, np.ndarray]:
    """Move the asymmetric unit into the unit cell (Harms-lab ``pdb_centerasu``).

    Coordinates are fractionalised with SCALEn (PDB) or ``_atom_sites``
    (mmCIF) - or, failing both, the cell - and the whole ASU is shifted by the
    integer lattice vector that brings its fractional centroid into [0, 1).
    Returns ``(new_structure, fractional_shift)``.
    """
    from .metadata import as_mmcif_dict, column, first, fractional_matrix
    d = as_mmcif_dict(structure)
    try:
        cell = [float(first(d, f"_cell.length_{x}")) for x in "abc"]
    except (TypeError, ValueError):
        cell = None
    if cell is not None and all(abs(c - 1.0) < 1e-6 for c in cell):
        method = "; ".join(column(d, "_exptl.method")) or "this entry"
        raise PdbKitError(f"No crystallographic unit cell ({method} files carry a 1 x 1 x 1 A placeholder cell); "
                          "centerasu only applies to crystal structures.")
    fm = fractional_matrix(structure)
    if fm is None:
        raise PdbKitError("No SCALE records, _atom_sites matrix or unit cell: cannot fractionalise.")
    fmx, fvec = fm
    xyz = np.array([a.coord for a in iter_atoms_all_altlocs(structure)], dtype=float)
    frac = xyz @ fmx.T + fvec
    fc = frac.mean(axis=0)
    shift = (fc % 1.0) - fc
    cart = np.linalg.inv(fmx) @ shift
    return translate(structure, cart), shift
