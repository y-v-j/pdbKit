"""Atom selection, chain extraction and splitting.

One :class:`Selection` describes what to keep. It covers everything the
``pdb_sel*`` / ``pdb_del*`` tools of pdb-tools did, plus the representation
modes asked for in chain splitting:

========================  ====================================================
``atoms``                 ``"all"`` (all-atom), ``"ca"`` (C-alpha only),
                          ``"backbone"`` (N, CA, C, O / nucleic backbone),
                          ``"trace"`` (CA for protein, P for nucleic acids),
                          or a list of atom names.  Applies to *polymer*
                          residues.
``het``                   ``"all"``, ``"none"`` or a list of residue names
                          (e.g. ``["HEM", "ZN"]``).  Applies to ligands
                          (non-polymer, non-water "heteroatom" groups).
``water``                 keep waters (HOH/WAT/DOD...)?
========================  ====================================================

Modified residues that are part of the chain (MSE, SEP, ...) are *polymer*
residues, not heteroatoms - see :func:`pdbkit.structure.annotate`.
"""

from __future__ import annotations

import os
import re
import sys
from dataclasses import dataclass, field, replace
from typing import Dict, Iterable, List, Optional, Sequence, Tuple, Union

from Bio.PDB.Model import Model
from Bio.PDB.Structure import Structure

from .constants import BACKBONE_ATOMS
from .structure import PdbKitError, deep_copy, residue_kind, save, source_info

__all__ = ["Selection", "extract", "split_chains", "split_models", "split_segments",
           "parse_ranges", "chain_ids", "write_parts", "describe"]

NUCLEIC_BACKBONE = ("P", "OP1", "OP2", "O5'", "C5'", "C4'", "C3'", "O3'")
AtomSpec = Union[str, Sequence[str]]


def _as_list(value) -> Optional[List[str]]:
    if value is None:
        return None
    if isinstance(value, str):
        return [v.strip() for v in value.split(",") if v.strip()]
    return [str(v).strip() for v in value]


# --------------------------------------------------------------------------- #
# Residue ranges: "1:10", "1-10", "5:", ":20", "1:100:5", "3,7,9" (comma lists)
# --------------------------------------------------------------------------- #

def parse_ranges(spec: Optional[str]) -> Optional[List[Tuple[Optional[int], Optional[int], int]]]:
    """Parse residue ranges into ``(start, end, step)`` triples (``None`` = open end).

    ``start:end:step`` follows pdb-tools' ``pdb_selres``; ``a-b`` is accepted
    for non-negative numbers. The step counts residues, not numbers.
    """
    if spec is None or str(spec).strip() == "":
        return None
    out = []
    for token in str(spec).split(","):
        t = token.strip()
        if not t:
            continue
        if ":" in t:
            parts = t.split(":")
            if len(parts) > 3 or not all(re.fullmatch(r"\s*(-?\d+)?\s*", p) for p in parts):
                raise PdbKitError(f"Bad residue range {t!r}; use e.g. 1:10, 5:, :20, 1:100:5 or 42")
            start = int(parts[0]) if parts[0].strip() else None
            end = int(parts[1]) if len(parts) > 1 and parts[1].strip() else None
            step = int(parts[2]) if len(parts) > 2 and parts[2].strip() else 1
            if step < 1:
                raise PdbKitError(f"Range step must be >= 1 in {t!r}")
        elif re.fullmatch(r"\d+-\d*", t):
            a, b = t.split("-")
            start, end, step = int(a), (int(b) if b else None), 1
        elif re.fullmatch(r"-?\d+", t):
            start = end = int(t)
            step = 1
        else:
            raise PdbKitError(f"Bad residue range {t!r}; use e.g. 1:10, 5:, :20, 1:100:5 or 42")
        out.append((start, end, step))
    return out


@dataclass
class Selection:
    """What to keep. The defaults keep everything."""

    chains: Optional[Sequence[str]] = None
    exclude_chains: Optional[Sequence[str]] = None
    atoms: AtomSpec = "all"
    exclude_atoms: Optional[Sequence[str]] = None
    het: AtomSpec = "all"
    water: bool = True
    hydrogens: bool = True
    models: Optional[Sequence[int]] = None          # 1-based model numbers
    residues: Optional[str] = None                  # range spec, see parse_ranges
    exclude_residues: Optional[str] = None
    resnames: Optional[Sequence[str]] = None
    exclude_resnames: Optional[Sequence[str]] = None
    elements: Optional[Sequence[str]] = None
    exclude_elements: Optional[Sequence[str]] = None
    segids: Optional[Sequence[str]] = None
    altloc: Optional[str] = None                    # None=keep all, 'best'=highest occupancy, or a label
    only_het_records: Optional[bool] = None         # True: HETATM only, False: drop HETATM records

    def normalized(self) -> "Selection":
        s = replace(self)
        for name in ("chains", "exclude_chains", "exclude_atoms", "resnames", "exclude_resnames", "segids"):
            setattr(s, name, _as_list(getattr(s, name)))
        s.elements = [e.upper() for e in _as_list(s.elements)] if s.elements is not None else None
        s.exclude_elements = [e.upper() for e in _as_list(s.exclude_elements)] if s.exclude_elements is not None else None
        if isinstance(s.atoms, str) and s.atoms.lower() in ("all", "ca", "backbone", "trace"):
            s.atoms = s.atoms.lower()
        else:
            s.atoms = _as_list(s.atoms)
        if isinstance(s.het, str) and s.het.lower() in ("all", "none"):
            s.het = s.het.lower()
        else:
            s.het = [h.upper() for h in _as_list(s.het)]
        if s.models is not None:
            s.models = [int(m) for m in (_as_list(s.models) if isinstance(s.models, str) else s.models)]
        return s


def _range_filter(residues: List, ranges) -> set:
    """ids of residues selected by a range list (order-aware for steps)."""
    keep = set()
    for start, end, step in ranges:
        span = [r for r in residues
                if (start is None or r.id[1] >= start) and (end is None or r.id[1] <= end)]
        keep.update(id(r) for r in span[::step])
    return keep


def _atom_ok(atom, residue, kind: str, sel: Selection) -> bool:
    name = atom.get_id()
    element = (atom.element or "").strip().upper()
    if kind == "polymer":
        mode = sel.atoms
        if mode == "ca":
            if name != "CA":
                return False
        elif mode == "backbone":
            if name not in BACKBONE_ATOMS and name not in NUCLEIC_BACKBONE:
                return False
        elif mode == "trace":
            ptype = residue.xtra.get("pdbkit_poly")
            if name != ("P" if ptype == "nucleic" else "CA"):
                return False
        elif mode != "all" and name not in mode:
            return False
    if sel.exclude_atoms and name in sel.exclude_atoms:
        return False
    if not sel.hydrogens and element in ("H", "D"):
        return False
    if sel.elements is not None and element not in sel.elements:
        return False
    if sel.exclude_elements and element in sel.exclude_elements:
        return False
    return True


def _residue_ok(residue, kind: str, sel: Selection) -> bool:
    name = residue.get_resname().strip()
    if kind == "water" and not sel.water:
        return False
    if kind == "ligand":
        if sel.het == "none":
            return False
        if isinstance(sel.het, list) and name.upper() not in sel.het:
            return False
    if sel.resnames is not None and name not in sel.resnames:
        return False
    if sel.exclude_resnames and name in sel.exclude_resnames:
        return False
    if sel.segids is not None and residue.segid.strip() not in sel.segids:
        return False
    if sel.only_het_records is True and residue.id[0] == " ":
        return False
    if sel.only_het_records is False and residue.id[0] != " ":
        return False
    return True


def _pick_residue_variant(chain, dres, mode: str):
    """Microheterogeneity (two residue identities at one position, e.g. PRO/SER):
    keep one variant - the highest mean occupancy, or the one carrying label ``mode``."""
    variants = dres.disordered_get_list()

    def mean_occ(r):
        occ = [a.occupancy or 0.0 for a in r.get_unpacked_list()]
        return sum(occ) / len(occ) if occ else 0.0

    chosen = None
    if mode != "best":
        chosen = next((v for v in variants if any(a.altloc == mode for a in v.get_unpacked_list())), None)
    if chosen is None:
        chosen = max(variants, key=mean_occ)
    pos = chain.child_list.index(dres)
    chain.detach_child(dres.id)
    chosen.detach_parent()
    chain.insert(pos, chosen)
    return chosen


def _resolve_altlocs(residue, mode: str) -> None:
    """Replace disordered atoms by a single location (label cleared, like pdb_selaltloc)."""
    for atom in list(residue.child_list):
        if not atom.is_disordered():
            if atom.altloc.strip():
                if mode != "best" and atom.altloc != mode:
                    residue.detach_child(atom.get_id())
                else:
                    atom.altloc = " "
            continue
        children = atom.disordered_get_list()
        if mode == "best":
            chosen = max(children, key=lambda a: (a.occupancy or 0.0, -ord(a.altloc or " ")))
        else:
            matches = [a for a in children if a.altloc == mode]
            if not matches:
                residue.detach_child(atom.get_id())
                continue
            chosen = matches[0]
        chosen.altloc = " "
        chosen.disordered_flag = 0
        residue.detach_child(atom.get_id())
        chosen.detach_parent()
        residue.add(chosen)
    if not any(a.is_disordered() for a in residue.child_list):
        residue.disordered = 0
    # keep original atom order
    residue.child_list.sort(key=lambda a: getattr(a, "serial_number", 0) or 0)


def _warn_missing_het(structure: Structure, sel: Selection) -> None:
    if not isinstance(sel.het, list):
        return
    present = {r.get_resname().strip().upper() for r in structure.get_residues() if residue_kind(r) == "ligand"}
    absent = [h for h in sel.het if h not in present]
    if absent:
        print(f"pdbKit: ligand(s) {', '.join(absent)} not present in the structure"
              f"{' (ligands here: ' + ', '.join(sorted(present)) + ')' if present else ' (it has no ligands)'}",
              file=sys.stderr)


def extract(structure: Structure, selection: Optional[Selection] = None, *, _check_het: bool = True,
            **kwargs) -> Structure:
    """Return a new structure containing only the selected atoms.

    Accepts a :class:`Selection` or its fields as keyword arguments::

        extract(s, chains="A,B", atoms="ca", het=["HEM"])
    """
    sel = (selection or Selection(**kwargs)).normalized()
    if selection is not None and kwargs:
        sel = replace(sel, **kwargs).normalized()
    ranges = parse_ranges(sel.residues)
    exclude_ranges = parse_ranges(sel.exclude_residues)
    if _check_het:
        _warn_missing_het(structure, sel)
    # copy only the models and chains that can contribute (a full deep copy of a
    # large assembly for every chain made splitting quadratic)
    new = Structure(structure.id)
    new.header = getattr(structure, "header", {})
    new.xtra = dict(structure.xtra)
    for model in structure:
        serial = getattr(model, "serial_num", None) or model.id + 1
        if sel.models is not None and serial not in sel.models:
            continue
        m = Model(model.id, serial)
        m.xtra = dict(model.xtra)
        for chain in model:
            if (sel.chains is not None and chain.id not in sel.chains) or \
                    (sel.exclude_chains and chain.id in sel.exclude_chains):
                continue
            m.add(deep_copy(chain))
        new.add(m)
    for model in list(new):
        for chain in list(model):
            residues = list(chain)
            in_range = _range_filter(residues, ranges) if ranges else None
            out_range = _range_filter(residues, exclude_ranges) if exclude_ranges else set()
            for res in residues:
                kind = residue_kind(res)
                if (in_range is not None and id(res) not in in_range) or id(res) in out_range \
                        or not _residue_ok(res, kind, sel):
                    chain.detach_child(res.id)
                    continue
                if sel.altloc:
                    if res.is_disordered() == 2:
                        res = _pick_residue_variant(chain, res, sel.altloc)
                    _resolve_altlocs(res, sel.altloc)
                variants = res.disordered_get_list() if res.is_disordered() == 2 else [res]
                for v in variants:                  # every microheterogeneity variant
                    for atom in list(v.child_list):
                        if not _atom_ok(atom, v, residue_kind(v), sel):
                            v.detach_child(atom.get_id())
                if all(len(v) == 0 for v in variants):
                    chain.detach_child(res.id)
            if len(chain) == 0:
                model.detach_child(chain.id)
        if len(model) == 0:
            new.detach_child(model.id)
    if len(new) == 0:
        raise PdbKitError("The selection matched no atoms.")
    return new


def chain_ids(structure: Structure) -> List[str]:
    """Chain IDs in file order (union over models)."""
    seen: List[str] = []
    for model in structure:
        for chain in model:
            if chain.id not in seen:
                seen.append(chain.id)
    return seen


def split_chains(structure: Structure, chains: Optional[Sequence[str]] = None,
                 selection: Optional[Selection] = None, **kwargs) -> Dict[str, Structure]:
    """One structure per chain, each filtered by the selection.

    ``chains`` limits which chains are produced (feature: *select chains by
    ID*). Chains that end up empty - e.g. a DNA chain in C-alpha mode - are
    skipped with a message on stderr.
    """
    sel = selection or Selection(**kwargs)
    wanted = _as_list(chains) if chains is not None else chain_ids(structure)
    present = chain_ids(structure)
    missing = [c for c in wanted if c not in present]
    if missing:
        raise PdbKitError(f"Chain(s) {', '.join(missing)} not found; available: {', '.join(present)}")
    _warn_missing_het(structure, sel.normalized())
    out: Dict[str, Structure] = {}
    for cid in wanted:
        try:
            out[cid] = extract(structure, replace(sel, chains=[cid]), _check_het=False)
        except PdbKitError:
            print(f"pdbKit: chain {cid!r} has no atoms matching the selection - skipped", file=sys.stderr)
    return out


def split_models(structure: Structure) -> Dict[int, Structure]:
    """One single-model structure per MODEL (keyed by model number)."""
    out = {}
    for model in structure:
        serial = getattr(model, "serial_num", None) or model.id + 1
        out[serial] = extract(structure, Selection(models=[serial]))
    return out


def split_segments(structure: Structure) -> Dict[str, Structure]:
    """One structure per segment identifier (PDB columns 73-76)."""
    segs: List[str] = []
    for res in structure.get_residues():
        s = res.segid.strip()
        if s not in segs:
            segs.append(s)
    return {s: extract(structure, Selection(segids=[s])) for s in segs}


def describe(selection: Selection) -> str:
    """Short file-name tag for a selection, e.g. ``ca_het-HEM`` or ``all_het_wat``."""
    sel = selection.normalized()
    atoms = sel.atoms if isinstance(sel.atoms, str) else "-".join(sel.atoms)
    parts = [atoms]
    if sel.het == "all":
        parts.append("het")
    elif isinstance(sel.het, list):
        parts.append("het-" + "-".join(sel.het))
    if sel.water:
        parts.append("wat")
    return "_".join(parts)


def write_parts(parts: Dict, outdir: str, stem: str, fmt: str, tag: str = "",
                remap_chains: bool = False) -> List[str]:
    """Write ``{key: structure}`` as ``<outdir>/<stem>_<key>[_<tag>].<ext>``."""
    os.makedirs(outdir, exist_ok=True)
    ext = "cif" if fmt == "cif" else "pdb"
    written = []
    for key, struct in parts.items():
        safe_key = re.sub(r"[^A-Za-z0-9_.-]", "_", str(key)) or "blank"
        name = f"{stem}_{safe_key}" + (f"_{tag}" if tag else "") + f".{ext}"
        written.append(save(struct, os.path.join(outdir, name), fmt, remap_chains=remap_chains))
    return written
