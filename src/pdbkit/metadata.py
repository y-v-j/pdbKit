"""Format-independent access to header metadata.

mmCIF files carry their metadata as categories (``_struct``, ``_entity``,
``_struct_conf``, ...). Legacy PDB files carry the same information in header
records (``TITLE``, ``COMPND``, ``DBREF``, ``HELIX``, ``SHEET``, ``SSBOND``,
``LINK``, ``CRYST1``, ``SCALE``, ``REMARK 2/3/350``, ``SEQRES``). The
:func:`as_mmcif_dict` adapter translates the PDB records into the mmCIF
category names, so every metadata consumer in pdbKit (domain/family
annotation, experiment summary, assemblies, SEQRES sequences, the mmCIF
writer's entity names) has exactly one code path for both formats.
"""

from __future__ import annotations

import re
from collections import OrderedDict
from typing import Dict, List, Optional, Tuple

import numpy as np
from Bio.PDB.Structure import Structure

from .structure import source_info

__all__ = ["as_mmcif_dict", "pdb_header_to_mmcif", "chain_descriptions", "experiment",
           "seqres", "assemblies", "fractional_matrix", "first", "column"]

METAL_ELEMENTS = {
    "LI", "NA", "K", "RB", "CS", "BE", "MG", "CA", "SR", "BA", "AL", "GA", "IN", "TL", "SN", "PB",
    "SC", "TI", "V", "CR", "MN", "FE", "CO", "NI", "CU", "ZN", "Y", "ZR", "NB", "MO", "TC", "RU",
    "RH", "PD", "AG", "CD", "HF", "TA", "W", "RE", "OS", "IR", "PT", "AU", "HG", "LA", "CE", "PR",
    "ND", "SM", "EU", "GD", "TB", "DY", "HO", "ER", "TM", "YB", "LU", "U",
}


def _is_null(v) -> bool:
    return v is None or v in ("?", ".", "")


def first(d: dict, key: str, default=None):
    v = d.get(key)
    if isinstance(v, list):
        v = v[0] if v else None
    return default if _is_null(v) else v


def column(d: dict, key: str) -> List[str]:
    v = d.get(key)
    if v is None:
        return []
    return v if isinstance(v, list) else [v]


# --------------------------------------------------------------------------- #
# PDB header -> mmCIF-style dict
# --------------------------------------------------------------------------- #

def _continued_text(lines: List[str], record: str) -> str:
    """Join a continued text record (TITLE, COMPND, SOURCE, KEYWDS...)."""
    parts = [ln[10:80].rstrip() for ln in lines if ln[:6].rstrip() == record]
    text = ""
    for p in parts:
        p = p.strip()
        if not text:
            text = p
        elif text.endswith("-") and not text.endswith(" -"):
            text += p          # hyphenated word split across lines
        else:
            text += " " + p
    return text


def _spec_list(text: str) -> List[Tuple[str, str]]:
    """Parse ``KEY: value; KEY: value;`` specification lists (COMPND/SOURCE)."""
    out = []
    for item in re.split(r";\s*(?=[A-Z_]+:)|;\s*$", text):
        if ":" not in item:
            if out and item.strip():
                k, v = out[-1]
                out[-1] = (k, (v + " " + item.strip()).strip())
            continue
        k, v = item.split(":", 1)
        out.append((k.strip(), v.strip().rstrip(";").strip()))
    return out


def _int_field(s: str) -> Optional[str]:
    s = s.strip()
    return s if re.fullmatch(r"-?\d+", s) else None


def _add_rows(d: dict, category: str, rows: List[Dict[str, str]]) -> None:
    if not rows:
        return
    keys = list(OrderedDict((k, None) for r in rows for k in r))
    for k in keys:
        d[f"{category}.{k}"] = [r.get(k, "?") if not _is_null(r.get(k)) else "?" for r in rows]


def pdb_header_to_mmcif(lines: List[str], structure: Optional[Structure] = None) -> Dict[str, List[str]]:
    """Translate PDB header records into mmCIF category/item names."""
    d: Dict[str, List[str]] = {}
    lines = [ln.rstrip("\n") for ln in lines]

    header = next((ln.ljust(80) for ln in lines if ln.startswith("HEADER")), None)
    if header:
        if header[62:66].strip():
            d["_entry.id"] = [header[62:66].strip()]
        if header[10:50].strip():
            d["_struct_keywords.pdbx_keywords"] = [header[10:50].strip()]
    title = _continued_text(lines, "TITLE")
    if title:
        d["_struct.title"] = [title]
    expdta = _continued_text(lines, "EXPDTA")
    if expdta:
        d["_exptl.method"] = [m.strip() for m in expdta.split(";") if m.strip()]

    remark2 = " ".join(ln[10:].strip() for ln in lines if ln.startswith("REMARK   2"))
    m = re.search(r"RESOLUTION\.\s*([0-9.]+)\s*ANGSTROM", remark2)
    if m:
        d["_refine.ls_d_res_high"] = [m.group(1)]
    remark3 = [ln[10:].strip() for ln in lines if ln.startswith("REMARK   3")]
    for ln in remark3:
        for key, rx in (("_refine.ls_R_factor_R_work", r"R VALUE\s+\(WORKING SET\)\s*:\s*([0-9.]+)"),
                        ("_refine.ls_R_factor_obs", r"R VALUE\s+\(WORKING \+ TEST SET\)\s*:\s*([0-9.]+)"),
                        ("_refine.ls_R_factor_R_free", r"FREE R VALUE\s*:\s*([0-9.]+)")):
            m = re.match(rx, ln)
            if m and key not in d:
                d[key] = [m.group(1)]

    # ---- entities from COMPND ------------------------------------------ #
    entities: List[Dict[str, str]] = []
    chain_to_entity: Dict[str, str] = {}
    current: Optional[Dict[str, str]] = None
    for key, value in _spec_list(_continued_text(lines, "COMPND")):
        if key == "MOL_ID":
            current = {"id": value, "type": "polymer", "pdbx_description": "?", "pdbx_ec": "?", "_chains": ""}
            entities.append(current)
        elif current is not None:
            if key == "MOLECULE":
                current["pdbx_description"] = value
            elif key == "EC":
                current["pdbx_ec"] = ",".join(e.strip() for e in value.split(",") if e.strip())
            elif key == "CHAIN":
                current["_chains"] = ",".join(c.strip() for c in value.split(",") if c.strip())
    for e in entities:
        for c in e["_chains"].split(","):
            if c:
                chain_to_entity[c] = e["id"]
    if not entities and structure is not None:
        # minimal files: one entity per polymer chain
        from .structure import residue_kind
        model = next(iter(structure), None)
        for i, chain in enumerate(model or [], start=1):
            if any(residue_kind(r) == "polymer" for r in chain):
                entities.append({"id": str(i), "type": "polymer", "pdbx_description": "?",
                                 "pdbx_ec": "?", "_chains": chain.id})
                chain_to_entity[chain.id] = str(i)
    _add_rows(d, "_entity", [{k: v for k, v in e.items() if not k.startswith("_")} for e in entities])
    _add_rows(d, "_entity_poly", [{"entity_id": e["id"], "pdbx_strand_id": e["_chains"] or "?"} for e in entities])

    organisms = {}
    mol = None
    for key, value in _spec_list(_continued_text(lines, "SOURCE")):
        if key == "MOL_ID":
            mol = value
        elif key == "ORGANISM_SCIENTIFIC" and mol is not None:
            organisms[mol] = value
    if organisms:
        _add_rows(d, "_entity_src_gen", [{"entity_id": k, "pdbx_gene_src_scientific_name": v}
                                         for k, v in organisms.items()])

    # ---- DBREF / DBREF1+2 -> _struct_ref / _struct_ref_seq ------------------ #
    refs: List[Dict[str, str]] = []
    ref_seqs: List[Dict[str, str]] = []
    ref_key_to_id: Dict[Tuple[str, str], str] = {}

    def add_ref(chain, db, acc, code, seq_beg, seq_end, db_beg, db_end):
        eid = chain_to_entity.get(chain, "?")
        key = (eid, acc)
        if key not in ref_key_to_id:
            ref_key_to_id[key] = str(len(refs) + 1)
            refs.append({"id": ref_key_to_id[key], "entity_id": eid, "db_name": db, "db_code": code,
                         "pdbx_db_accession": acc})
        ref_seqs.append({"align_id": str(len(ref_seqs) + 1), "ref_id": ref_key_to_id[key],
                         "pdbx_strand_id": chain, "db_align_beg": db_beg, "db_align_end": db_end,
                         "pdbx_auth_seq_align_beg": seq_beg, "pdbx_auth_seq_align_end": seq_end})

    dbref1: Dict[str, List[str]] = {}
    for ln in lines:
        rec = ln[:6]
        ln = ln.ljust(80)
        if rec == "DBREF ":
            add_ref(ln[12], ln[26:32].strip(), ln[33:41].strip(), ln[42:54].strip(),
                    ln[14:18].strip(), ln[20:24].strip(), ln[55:60].strip(), ln[62:67].strip())
        elif rec == "DBREF1":
            dbref1[ln[12]] = [ln[14:18].strip(), ln[20:24].strip(), ln[26:32].strip(), ln[47:67].strip()]
        elif rec == "DBREF2" and ln[12] in dbref1:
            sb, se, db, code = dbref1.pop(ln[12])
            add_ref(ln[12], db, ln[18:40].strip(), code, sb, se, ln[45:55].strip(), ln[57:67].strip())
    _add_rows(d, "_struct_ref", refs)
    _add_rows(d, "_struct_ref_seq", ref_seqs)

    # ---- secondary structure ------------------------------------------------ #
    helices = []
    for ln in lines:
        if ln.startswith("HELIX "):
            ln = ln.ljust(80)
            helices.append({"id": ln[11:14].strip(), "conf_type_id": "HELX_P",
                            "beg_auth_asym_id": ln[19], "beg_auth_seq_id": ln[21:25].strip(),
                            "end_auth_asym_id": ln[31], "end_auth_seq_id": ln[33:37].strip()})
        elif ln.startswith("TURN  "):
            ln = ln.ljust(80)
            helices.append({"id": ln[11:14].strip(), "conf_type_id": "TURN_P",
                            "beg_auth_asym_id": ln[19], "beg_auth_seq_id": ln[20:24].strip(),
                            "end_auth_asym_id": ln[30], "end_auth_seq_id": ln[31:35].strip()})
    _add_rows(d, "_struct_conf", helices)
    sheets = []
    for ln in lines:
        if ln.startswith("SHEET "):
            ln = ln.ljust(80)
            sheets.append({"sheet_id": ln[11:14].strip(), "id": ln[7:10].strip(),
                           "beg_auth_asym_id": ln[21], "beg_auth_seq_id": ln[22:26].strip(),
                           "end_auth_asym_id": ln[32], "end_auth_seq_id": ln[33:37].strip()})
    _add_rows(d, "_struct_sheet_range", sheets)

    # ---- connectivity: SSBOND + LINK -> _struct_conn ------------------------ #
    conns = []
    elements = {}
    if structure is not None:
        for atom in structure.get_atoms():
            res = atom.get_parent()
            elements[(res.get_parent().id, res.id[1], atom.get_id())] = (atom.element or "").upper()
    for ln in lines:
        if ln.startswith("SSBOND"):
            ln = ln.ljust(80)
            conns.append({"id": f"disulf{len(conns) + 1}", "conn_type_id": "disulf",
                          "ptnr1_auth_asym_id": ln[15], "ptnr1_auth_comp_id": ln[11:14].strip(),
                          "ptnr1_auth_seq_id": ln[17:21].strip(), "ptnr1_label_atom_id": "SG",
                          "ptnr2_auth_asym_id": ln[29], "ptnr2_auth_comp_id": ln[25:28].strip(),
                          "ptnr2_auth_seq_id": ln[31:35].strip(), "ptnr2_label_atom_id": "SG"})
        elif ln.startswith("LINK  ") or ln.startswith("LINKR "):
            ln = ln.ljust(80)
            c1, s1, a1 = ln[21], ln[22:26].strip(), ln[12:16].strip()
            c2, s2, a2 = ln[51], ln[52:56].strip(), ln[42:46].strip()
            e1 = elements.get((c1, int(s1) if _int_field(s1) else None, a1), "")
            e2 = elements.get((c2, int(s2) if _int_field(s2) else None, a2), "")
            ctype = "metalc" if (e1 in METAL_ELEMENTS or e2 in METAL_ELEMENTS) else "covale"
            conns.append({"id": f"{ctype}{len(conns) + 1}", "conn_type_id": ctype,
                          "ptnr1_auth_asym_id": c1, "ptnr1_auth_comp_id": ln[17:20].strip(),
                          "ptnr1_auth_seq_id": s1, "ptnr1_label_atom_id": a1,
                          "ptnr2_auth_asym_id": c2, "ptnr2_auth_comp_id": ln[47:50].strip(),
                          "ptnr2_auth_seq_id": s2, "ptnr2_label_atom_id": a2})
    _add_rows(d, "_struct_conn", conns)

    # ---- crystal ------------------------------------------------------------ #
    cryst = next((ln.ljust(80) for ln in lines if ln.startswith("CRYST1")), None)
    if cryst:
        try:
            for key, sl in (("length_a", slice(6, 15)), ("length_b", slice(15, 24)), ("length_c", slice(24, 33)),
                            ("angle_alpha", slice(33, 40)), ("angle_beta", slice(40, 47)),
                            ("angle_gamma", slice(47, 54))):
                d[f"_cell.{key}"] = [str(float(cryst[sl]))]
            d["_symmetry.space_group_name_H-M"] = [cryst[55:66].strip()]
        except ValueError:
            pass
    for ln in lines:
        if ln[:5] == "SCALE" and ln[5:6] in "123":
            i = ln[5]
            ln = ln.ljust(80)
            try:
                vals = [float(ln[10:20]), float(ln[20:30]), float(ln[30:40]), float(ln[45:55])]
            except ValueError:
                continue
            for j in range(3):
                d[f"_atom_sites.fract_transf_matrix[{i}][{j + 1}]"] = [str(vals[j])]
            d[f"_atom_sites.fract_transf_vector[{i}]"] = [str(vals[3])]

    # ---- SEQRES / MODRES -> _pdbx_poly_seq_scheme / _pdbx_struct_mod_residue -- #
    seq_rows = []
    counters: Dict[str, int] = {}
    for ln in lines:
        if ln.startswith("SEQRES"):
            chain = ln.ljust(12)[11]
            for name in ln[19:].split():
                counters[chain] = counters.get(chain, 0) + 1
                seq_rows.append({"pdb_strand_id": chain, "seq_id": str(counters[chain]), "mon_id": name,
                                 "entity_id": chain_to_entity.get(chain, "?")})
    _add_rows(d, "_pdbx_poly_seq_scheme", seq_rows)
    modres = []
    for ln in lines:
        if ln.startswith("MODRES"):
            ln = ln.ljust(80)
            modres.append({"auth_comp_id": ln[12:15].strip(), "auth_asym_id": ln[16],
                           "auth_seq_id": ln[18:22].strip(), "parent_comp_id": ln[24:27].strip()})
    _add_rows(d, "_pdbx_struct_mod_residue", modres)

    # ---- REMARK 350 -> assemblies ------------------------------------------ #
    d.update(_remark350_to_mmcif([ln for ln in lines if ln.startswith("REMARK 350")]))
    return d


def _remark350_to_mmcif(lines: List[str]) -> Dict[str, List[str]]:
    assemblies: List[Dict[str, str]] = []
    gens: List[Dict[str, str]] = []
    opers: Dict[str, Dict[str, str]] = {}
    oper_key: Dict[Tuple, str] = {}
    cur = None
    chains: List[str] = []
    matrices: Dict[str, Dict[int, List[float]]] = {}
    order: List[str] = []

    def flush():
        """Emit one _pdbx_struct_assembly_gen row for the current chains+operators block."""
        if cur is None or (not chains and not order):
            return
        ids = []
        for mid in order:
            rows = matrices[mid]
            if len(rows) != 3:
                continue
            flat = tuple(round(x, 6) for r in (rows[1], rows[2], rows[3]) for x in r)
            if flat not in oper_key:
                oid = str(len(oper_key) + 1)
                oper_key[flat] = oid
                row = {"id": oid, "type": "?"}
                for i in range(3):
                    for j in range(3):
                        row[f"matrix[{i + 1}][{j + 1}]"] = str(rows[i + 1][j])
                    row[f"vector[{i + 1}]"] = str(rows[i + 1][3])
                opers[oid] = row
            ids.append(oper_key[flat])
        gens.append({"assembly_id": cur["id"], "oper_expression": ",".join(ids) or "?",
                     "asym_id_list": ",".join(chains) or "?"})

    for ln in lines:
        text = ln[10:].rstrip()
        m = re.match(r"\s*BIOMOLECULE:\s*(\S+)", text)
        if m:
            flush()
            cur = {"id": m.group(1), "details": "?", "oligomeric_details": "?"}
            assemblies.append(cur)
            chains, matrices, order = [], {}, []
            continue
        if cur is None:
            continue
        m = re.match(r"\s*AUTHOR DETERMINED BIOLOGICAL UNIT:\s*(.+)", text)
        if m:
            cur["details"] = "author_defined_assembly"
            cur["oligomeric_details"] = m.group(1).strip().lower()
            continue
        m = re.match(r"\s*SOFTWARE DETERMINED QUATERNARY STRUCTURE:\s*(.+)", text)
        if m:
            if cur["details"] == "?":
                cur["details"] = "software_defined_assembly"
                cur["oligomeric_details"] = m.group(1).strip().lower()
            continue
        m = re.match(r"\s*(APPLY THE FOLLOWING TO CHAINS|AND CHAINS):\s*(.*)", text)
        if m:
            if m.group(1).startswith("APPLY") and order:
                flush()                      # a new chains/operators block begins
                chains, matrices, order = [], {}, []
            chains.extend(c.strip() for c in m.group(2).split(",") if c.strip())
            continue
        m = re.match(r"\s*BIOMT([123])\s+(\d+)\s+(\S+)\s+(\S+)\s+(\S+)\s+(\S+)", text)
        if m:
            row, mid = int(m.group(1)), m.group(2)
            if mid not in matrices:
                matrices[mid] = {}
                order.append(mid)
            matrices[mid][row] = [float(m.group(k)) for k in range(3, 7)]
    flush()
    d: Dict[str, List[str]] = {}
    _add_rows(d, "_pdbx_struct_assembly", assemblies)
    _add_rows(d, "_pdbx_struct_assembly_gen", gens)
    _add_rows(d, "_pdbx_struct_oper_list", list(opers.values()))
    if gens:
        d["_pdbkit.assembly_chain_ids"] = ["auth"]   # asym_id_list holds author chain IDs
    return d


def as_mmcif_dict(structure: Structure) -> Dict[str, List[str]]:
    """The structure's metadata as an mmCIF-style dict, whatever the input format."""
    info = source_info(structure)
    if info.mmcif is not None:
        return info.mmcif
    cache = structure.xtra.get("pdbkit_meta_cache")
    if cache is None:
        cache = pdb_header_to_mmcif(info.header_lines, structure)
        structure.xtra["pdbkit_meta_cache"] = cache
    return cache


# --------------------------------------------------------------------------- #
# Derived, user-facing metadata
# --------------------------------------------------------------------------- #

def chain_descriptions(structure: Structure) -> Dict[str, str]:
    """Author chain ID -> entity description (polymer chains)."""
    d = as_mmcif_dict(structure)
    desc = dict(zip(column(d, "_entity.id"), column(d, "_entity.pdbx_description")))
    out = {}
    for eid, strands in zip(column(d, "_entity_poly.entity_id"), column(d, "_entity_poly.pdbx_strand_id")):
        if _is_null(strands):
            continue
        for c in strands.split(","):
            if c.strip() and not _is_null(desc.get(eid)):
                out[c.strip()] = desc[eid]
    return out


def experiment(structure: Structure) -> Dict[str, str]:
    """Title, classification, organism, method, resolution, R-work and R-free."""
    d = as_mmcif_dict(structure)
    organism = first(d, "_entity_src_gen.pdbx_gene_src_scientific_name") or \
        first(d, "_entity_src_nat.pdbx_organism_scientific") or \
        first(d, "_pdbx_entity_src_syn.organism_scientific") or "NA"
    resolution = first(d, "_refine.ls_d_res_high") or first(d, "_em_3d_reconstruction.resolution") or \
        first(d, "_reflns.d_resolution_high") or "NA"
    return {
        "id": first(d, "_entry.id", source_info(structure).entry_id or str(structure.id)),
        "title": first(d, "_struct.title", "NA"),
        "classification": first(d, "_struct_keywords.pdbx_keywords", "NA"),
        "organism": organism,
        "method": "; ".join(m for m in column(d, "_exptl.method") if not _is_null(m)) or "NA",
        "resolution": resolution,
        "r_work": first(d, "_refine.ls_R_factor_R_work") or first(d, "_refine.ls_R_factor_obs") or "NA",
        "r_free": first(d, "_refine.ls_R_factor_R_free", "NA"),
    }


def seqres(structure: Structure) -> Dict[str, List[str]]:
    """Author chain -> full deposited sequence (3-letter names), incl. unobserved residues."""
    d = as_mmcif_dict(structure)
    out: Dict[str, List[str]] = {}
    strands = column(d, "_pdbx_poly_seq_scheme.pdb_strand_id")
    mons = column(d, "_pdbx_poly_seq_scheme.mon_id")
    if strands:
        seen = set()
        asym = column(d, "_pdbx_poly_seq_scheme.asym_id") or strands
        seq_ids = column(d, "_pdbx_poly_seq_scheme.seq_id") or [str(i) for i in range(len(strands))]
        for a, c, s, m in zip(asym, strands, seq_ids, mons):
            key = (a, c, s)
            if key in seen:      # microheterogeneity: keep the first monomer
                continue
            seen.add(key)
            out.setdefault(c, []).append(m)
        return out
    # mmCIF without the scheme: expand _entity_poly_seq over pdbx_strand_id
    ents = column(d, "_entity_poly_seq.entity_id")
    if ents:
        by_entity: Dict[str, List[str]] = {}
        seen = set()
        for e, n, m in zip(ents, column(d, "_entity_poly_seq.num"), column(d, "_entity_poly_seq.mon_id")):
            if (e, n) in seen:
                continue
            seen.add((e, n))
            by_entity.setdefault(e, []).append(m)
        for eid, st in zip(column(d, "_entity_poly.entity_id"), column(d, "_entity_poly.pdbx_strand_id")):
            for c in (st or "").split(","):
                if c.strip() and eid in by_entity:
                    out[c.strip()] = list(by_entity[eid])
    return out


def _expand_oper_expression(expr: str) -> List[List[str]]:
    """'1', '1,2', '(1-60)', '(1-5)(6-10)', '(X0)(1-60)' -> list of operator products."""
    groups = re.findall(r"\(([^)]*)\)", expr) if "(" in expr else [expr]

    def expand(g: str) -> List[str]:
        ids = []
        for part in g.split(","):
            part = part.strip()
            m = re.fullmatch(r"(\d+)-(\d+)", part)
            if m:
                ids.extend(str(i) for i in range(int(m.group(1)), int(m.group(2)) + 1))
            elif part:
                ids.append(part)
        return ids

    products: List[List[str]] = [[]]
    for g in groups:
        products = [p + [op] for p in products for op in expand(g)]
    return products


def _label_to_auth_chains(structure: Structure) -> Dict[str, str]:
    d = source_info(structure).mmcif or {}
    return dict(zip(column(d, "_atom_site.label_asym_id"), column(d, "_atom_site.auth_asym_id")))


def assemblies(structure: Structure) -> List[Dict[str, object]]:
    """Biological assemblies: id, details, oligomeric state, author chains and operators.

    ``is_asymmetric_unit`` is True when the assembly is exactly the deposited
    coordinates (identity operator on all chains) - the question Harms-lab
    ``pdb_oligomer`` answered.
    """
    d = as_mmcif_dict(structure)
    ids = column(d, "_pdbx_struct_assembly.id")
    if not ids:
        return []
    details = column(d, "_pdbx_struct_assembly.details") or ["?"] * len(ids)
    olig = column(d, "_pdbx_struct_assembly.oligomeric_details") or ["?"] * len(ids)
    count = column(d, "_pdbx_struct_assembly.oligomeric_count") or ["?"] * len(ids)
    ops = {}
    for i, oid in enumerate(column(d, "_pdbx_struct_oper_list.id")):
        m = np.array([[float(column(d, f"_pdbx_struct_oper_list.matrix[{r}][{c}]")[i]) for c in (1, 2, 3)]
                      for r in (1, 2, 3)])
        v = np.array([float(column(d, f"_pdbx_struct_oper_list.vector[{r}]")[i]) for r in (1, 2, 3)])
        ops[oid] = (m, v)
    auth_ids = "_pdbkit.assembly_chain_ids" in d
    label_to_auth = {} if auth_ids else _label_to_auth_chains(structure)
    all_chains = sorted({c.id for c in next(iter(structure))})
    out = []
    for i, aid in enumerate(ids):
        chains: List[str] = []
        n_ops = 0
        identity_only = True
        for a, expr, asyms in zip(column(d, "_pdbx_struct_assembly_gen.assembly_id"),
                                  column(d, "_pdbx_struct_assembly_gen.oper_expression"),
                                  column(d, "_pdbx_struct_assembly_gen.asym_id_list")):
            if a != aid:
                continue
            for asym in asyms.split(","):
                c = asym.strip() if auth_ids else label_to_auth.get(asym.strip(), asym.strip())
                if c and c not in chains:
                    chains.append(c)
            products = _expand_oper_expression(expr)
            n_ops += len(products)
            for prod in products:
                for op in prod:
                    if op in ops:
                        m, v = ops[op]
                        if not (np.allclose(m, np.eye(3), atol=1e-5) and np.allclose(v, 0, atol=1e-4)):
                            identity_only = False
        n_copies = n_ops
        out.append({
            "id": aid,
            "details": details[i] if i < len(details) else "?",
            "oligomeric_details": olig[i] if i < len(olig) else "?",
            "oligomeric_count": count[i] if i < len(count) else "?",
            "chains": chains,
            "n_operators": n_copies,
            # Harms pdb_oligomer's "wysiwyg": the assembly needs no symmetry operation
            "in_asymmetric_unit": identity_only and n_copies == 1,
            # stricter: the assembly is exactly the deposited file
            "is_asymmetric_unit": identity_only and n_copies == 1 and sorted(set(chains)) == all_chains,
        })
    return out


def fractional_matrix(structure: Structure) -> Optional[Tuple[np.ndarray, np.ndarray]]:
    """Orthogonal->fractional matrix and vector (SCALEn / _atom_sites), else from the cell."""
    d = as_mmcif_dict(structure)
    try:
        m = np.array([[float(first(d, f"_atom_sites.fract_transf_matrix[{i}][{j}]")) for j in (1, 2, 3)]
                      for i in (1, 2, 3)])
        v = np.array([float(first(d, f"_atom_sites.fract_transf_vector[{i}]", 0.0)) for i in (1, 2, 3)])
        return m, v
    except (TypeError, ValueError):
        pass
    try:
        a, b, c = (float(first(d, f"_cell.length_{x}")) for x in "abc")
        al, be, ga = (np.radians(float(first(d, f"_cell.angle_{x}"))) for x in ("alpha", "beta", "gamma"))
    except (TypeError, ValueError):
        return None
    # standard PDB orthogonalisation (a along x, b in xy plane)
    vol = a * b * c * np.sqrt(1 - np.cos(al) ** 2 - np.cos(be) ** 2 - np.cos(ga) ** 2
                              + 2 * np.cos(al) * np.cos(be) * np.cos(ga))
    orth = np.array([
        [a, b * np.cos(ga), c * np.cos(be)],
        [0, b * np.sin(ga), c * (np.cos(al) - np.cos(be) * np.cos(ga)) / np.sin(ga)],
        [0, 0, vol / (a * b * np.sin(ga))],
    ])
    return np.linalg.inv(orth), np.zeros(3)
