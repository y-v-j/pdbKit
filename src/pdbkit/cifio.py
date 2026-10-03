"""A small, spec-conforming PDBx/mmCIF coordinate writer.

Biopython's ``MMCIFIO`` writes only ``_atom_site`` and gives every residue of
an author chain (polymer, ligands and waters) the same ``label_asym_id``.
Mol*, PyMOL, ChimeraX and gemmi all read that, but they then have to guess
which atoms are polymer, ligand or water. This writer follows the wwPDB
conventions instead:

* one ``label_asym_id`` per polymer chain, per ligand instance and per
  chain's waters;
* ``_entity`` (polymer / non-polymer / water), ``_entity_poly``,
  ``_struct_asym``, ``_chem_comp`` and ``_atom_type`` categories;
* ``label_seq_id`` taken from the source mmCIF when available, otherwise
  derived from author numbering so that chain breaks remain visible as gaps
  (viewers do not draw cartoon across them);
* author fields (``auth_asym_id``, ``auth_seq_id``, ...) preserved exactly.
"""

from __future__ import annotations

import re
from typing import Dict, List, Optional, Tuple

from Bio.PDB.Structure import Structure

from .structure import iter_residues_all, residue_kind, source_info

_RESERVED = re.compile(r"^(data_|loop_|save_|global_|stop_)", re.IGNORECASE)


def cif_quote(value) -> str:
    """Quote a value for mmCIF output when it needs it."""
    s = "" if value is None else str(value)
    if s == "":
        return "?"
    needs = any(ch.isspace() or ch in "'\"" for ch in s) or s[0] in "_#$;[]" or _RESERVED.match(s)
    if not needs:
        return s
    if "'" not in s:
        return f"'{s}'"
    if '"' not in s:
        return f'"{s}"'
    return f";{s}\n;"  # text field; both quote kinds present


def _asym_names():
    """A, B, ..., Z, AA, BA, CA, ... (wwPDB style)."""
    letters = [chr(c) for c in range(ord("A"), ord("Z") + 1)]
    for c in letters:
        yield c
    n = 2
    while True:
        # base-26 little-endian names of length n
        for i in range(26 ** n):
            name, k = "", i
            for _ in range(n):
                name += letters[k % 26]
                k //= 26
            yield name
        n += 1


def _chem_comp_type(residue) -> str:
    kind = residue_kind(residue)
    if kind != "polymer":
        return "non-polymer"
    ptype = residue.xtra.get("pdbkit_poly")
    name = residue.get_resname().strip()
    if ptype == "nucleic":
        return "DNA linking" if name.startswith("D") and len(name) == 2 else "RNA linking"
    return "peptide linking" if name == "GLY" else "L-peptide linking"


def _entity_poly_type(residues) -> str:
    ptypes = {r.xtra.get("pdbkit_poly") for r in residues}
    if ptypes == {"nucleic"}:
        dna = all(r.get_resname().strip().startswith("D") for r in residues)
        return "polydeoxyribonucleotide" if dna else "polyribonucleotide"
    if "nucleic" in ptypes and "protein" in ptypes:
        return "other"
    return "polypeptide(L)"


def _derived_label_seq(residues) -> List[int]:
    """label_seq_id per residue: the source mmCIF's values when every residue has
    one, else sequential numbering that keeps author-numbering gaps as gaps.
    Microheterogeneity variants (same number and insertion code) share a value."""
    stored = [r.xtra.get("label_seq_id") for r in residues]
    if all(s is not None for s in stored):
        return [int(s) for s in stored]
    out: List[int] = []
    prev = None
    for r in residues:
        if prev is None:
            out.append(1)
        elif r.id[1:] == prev.id[1:]:
            out.append(out[-1])
        else:
            step = r.id[1] - prev.id[1]
            out.append(out[-1] + (step if step > 0 else 1))
        prev = r
    return out


def _descriptions(structure: Structure) -> Dict[str, str]:
    """Chain -> polymer description from the source header, when there is one."""
    try:
        from .metadata import chain_descriptions
        return chain_descriptions(structure)
    except Exception:  # noqa: BLE001 - descriptions are cosmetic
        return {}


def dump_cif(structure: Structure, data_name: Optional[str] = None) -> str:
    info = source_info(structure)
    block = re.sub(r"[^A-Za-z0-9_.-]", "_", data_name or info.entry_id or str(structure.id) or "pdbkit")
    descriptions = _descriptions(structure)

    # ---- asym units and entities (keyed so that all models agree) -------- #
    asym_of: Dict[Tuple, str] = {}
    asym_entity: Dict[str, str] = {}
    entity_rows: List[Dict[str, str]] = []
    entity_key_to_id: Dict[Tuple, str] = {}
    entity_poly: Dict[str, Dict[str, object]] = {}
    chem_comp: Dict[str, str] = {}
    names = _asym_names()
    label_seq_of: Dict[int, int] = {}   # id(residue) -> label_seq_id

    def entity_for(key: Tuple, etype: str, description: str) -> str:
        if key not in entity_key_to_id:
            eid = str(len(entity_key_to_id) + 1)
            entity_key_to_id[key] = eid
            entity_rows.append({"id": eid, "type": etype, "description": description})
        return entity_key_to_id[key]

    for model in structure:
        for chain in model:
            residues = chain.get_unpacked_list()          # every microheterogeneity variant
            poly = [r for r in residues if residue_kind(r) == "polymer"]
            if poly:
                seq_numbers = _derived_label_seq(poly)
                for r, n in zip(poly, seq_numbers):
                    label_seq_of[id(r)] = n
                key = (chain.id, "poly")
                if key not in asym_of:
                    asym_of[key] = next(names)
                    ptype = _entity_poly_type(poly)
                    ekey = ("poly", ptype, tuple(r.get_resname() for r in poly))
                    eid = entity_for(ekey, "polymer", descriptions.get(chain.id, "?"))
                    asym_entity[asym_of[key]] = eid
                    ep = entity_poly.setdefault(eid, {"type": ptype, "chains": []})
                    if chain.id not in ep["chains"]:
                        ep["chains"].append(chain.id)
            for r in residues:
                kind = residue_kind(r)
                name = r.get_resname().strip()
                chem_comp.setdefault(name, _chem_comp_type(r))
                if kind == "ligand":
                    key = (chain.id, "lig", r.id[1], r.id[2], name)
                    if key not in asym_of:
                        asym_of[key] = next(names)
                        asym_entity[asym_of[key]] = entity_for(("lig", name), "non-polymer", name)
                elif kind == "water":
                    key = (chain.id, "wat")
                    if key not in asym_of:
                        asym_of[key] = next(names)
                        asym_entity[asym_of[key]] = entity_for(("wat",), "water", "water")

    out: List[str] = [f"data_{block}", "#", f"_entry.id {cif_quote(block)}", "#"]

    out += ["loop_", "_entity.id", "_entity.type", "_entity.pdbx_description"]
    out += [f"{e['id']} {e['type']} {cif_quote(e['description'])}" for e in entity_rows]
    out.append("#")

    if entity_poly:
        out += ["loop_", "_entity_poly.entity_id", "_entity_poly.type",
                "_entity_poly.nstd_linkage", "_entity_poly.pdbx_strand_id"]
        for eid, ep in entity_poly.items():
            out.append(f"{eid} {cif_quote(ep['type'])} no {cif_quote(','.join(ep['chains']))}")
        out.append("#")

    out += ["loop_", "_struct_asym.id", "_struct_asym.pdbx_blank_PDB_chainid_flag",
            "_struct_asym.pdbx_modified", "_struct_asym.entity_id", "_struct_asym.details"]
    for key, aid in asym_of.items():
        blank = "Y" if key[0].strip() == "" else "N"
        out.append(f"{aid} {blank} N {asym_entity[aid]} ?")
    out.append("#")

    out += ["loop_", "_chem_comp.id", "_chem_comp.type"]
    out += [f"{cif_quote(n)} {cif_quote(t)}" for n, t in sorted(chem_comp.items())]
    out.append("#")

    elements = sorted({(a.element or "X").strip().upper() or "X"
                       for r in iter_residues_all(structure) for a in r.get_unpacked_list()})
    out += ["loop_", "_atom_type.symbol"] + [cif_quote(e) for e in elements] + ["#"]

    cols = ["group_PDB", "id", "type_symbol", "label_atom_id", "label_alt_id", "label_comp_id",
            "label_asym_id", "label_entity_id", "label_seq_id", "pdbx_PDB_ins_code",
            "Cartn_x", "Cartn_y", "Cartn_z", "occupancy", "B_iso_or_equiv", "pdbx_formal_charge",
            "auth_seq_id", "auth_comp_id", "auth_asym_id", "auth_atom_id", "pdbx_PDB_model_num"]
    out.append("loop_")
    out += [f"_atom_site.{c}" for c in cols]
    serial = 0
    for model in structure:
        mnum = getattr(model, "serial_num", None) or model.id + 1
        for chain in model:
            auth_chain = cif_quote(chain.id if chain.id.strip() else " ")
            for r in chain.get_unpacked_list():
                kind = residue_kind(r)
                name = r.get_resname().strip()
                if kind == "polymer":
                    aid = asym_of[(chain.id, "poly")]
                    lseq = str(label_seq_of[id(r)])
                elif kind == "ligand":
                    aid = asym_of[(chain.id, "lig", r.id[1], r.id[2], name)]
                    lseq = "."
                else:
                    aid = asym_of[(chain.id, "wat")]
                    lseq = "."
                group = "ATOM" if r.id[0] == " " else "HETATM"
                icode = r.id[2].strip() or "?"
                cname = cif_quote(name)
                for a in r.get_unpacked_list():
                    serial += 1
                    x, y, z = (float(v) for v in a.coord)
                    occ = "?" if a.occupancy is None else f"{a.occupancy:.2f}"
                    b = "?" if a.bfactor is None else f"{a.bfactor:.2f}"
                    elem = (a.element or "X").strip().upper() or "X"
                    aname = cif_quote(a.get_id())
                    alt = a.altloc.strip() or "."
                    out.append(
                        f"{group} {serial} {elem} {aname} {alt} {cname} {aid} {asym_entity[aid]} {lseq} "
                        f"{icode} {x:.3f} {y:.3f} {z:.3f} {occ} {b} ? {r.id[1]} {cname} {auth_chain} "
                        f"{aname} {mnum}")
    out.append("#")
    return "\n".join(out) + "\n"
