"""Domain, family, EC, secondary-structure and heteroatom-bonding annotation
(integrated ``mmcif_domain_family.py``).

For each structure this writes ``<PDB_ID>_domain_ss_info.csv`` with one row
per chain: title, experimental method, resolution, EC number/type/KEGG
classification, enzyme type/name, UniProt accession, family (SIMILARITY
comment), UniProt domains mapped to PDB numbering, Pfam/InterPro/CATH/SCOP
cross-references, secondary-structure order and ranges, heteroatoms and their
bonding partners.

Differences from the original script:

* **PDB-format input is fully supported.** The original reported most
  columns as N/A for ``.pdb`` files, stating that the information is
  mmCIF-only. wwPDB-format files do carry it - ``HELIX``/``SHEET`` (secondary
  structure), ``SSBOND``/``LINK`` (bonding), ``DBREF`` (UniProt mapping),
  ``COMPND EC:`` - and :func:`pdbkit.metadata.as_mmcif_dict` translates
  those records, so both formats go through the same code.
* Structures are loaded through :func:`pdbkit.load` (gzip, stdin, IDs).

All parsing and lookup functions below are the original ones, unchanged, so
mmCIF output is identical to the script's (verified by the test suite).
"""

from __future__ import annotations

import csv
import json
import os
import re
import sys
import tempfile
import urllib.error
import urllib.request
from dataclasses import dataclass, field
from typing import Dict, List, Optional, Tuple

RCSB_CIF_URL = "https://files.rcsb.org/download/{pdb_id}.cif"
UNIPROT_REST_URL = "https://rest.uniprot.org/uniprotkb/{accession}.json"

FAMILY_DOMAIN_DBS = {"Pfam", "InterPro", "CATH", "SCOP", "SCOP2", "SUPFAM", "PROSITE", "Gene3D"}

STRUCTURE_FILE_EXTENSIONS = (".cif", ".mmcif", ".pdb", ".ent")

# EC first-digit -> broad enzyme category, per the standard EC classification
EC_BROAD_CATEGORIES = {
    "1": "Oxidoreductases",
    "2": "Transferases",
    "3": "Hydrolases",
    "4": "Lyases",
    "5": "Isomerases",
    "6": "Ligases",
    "7": "Translocases",
}

CSV_FIELDNAMES = [
    "PDB_ID", "Title", "Experimental_Method", "Resolution",
    "Chain_ID", "Entity_ID", "Entity_Description",
    "EC_Number", "EC_Type", "EC_Detailed_Classification",
    "Enzyme Type", "Enzyme Name",
    "UniProt_Accession", "UniProt_Family",
    "UniProt_Domains_PDB_Numbering",
    "Family_Domain_DB_Xrefs",
    "Secondary_Structure_Order",
    "Secondary_Structure_Residue_Ranges",
    "Heteroatoms_In_Chain",
    "Heteroatom_Bonding_Partners",
]


# --------------------------------------------------------------------------- #
# Generic mmCIF dict helpers
# --------------------------------------------------------------------------- #

def _as_list(value):
    """MMCIF2Dict gives a plain string for single-row categories and a list
    for multi-row ones. Normalize everything to a list so downstream code
    doesn't have to special-case it."""
    if value is None:
        return []
    if isinstance(value, list):
        return value
    return [value]


def _column(mmcif_dict: dict, category: str, tags: List[str]) -> Dict[str, List[str]]:
    """Fetch several tags of one mmCIF category as equal-length lists,
    padding missing/unknown values with '?'."""
    columns = {}
    length = None
    for tag in tags:
        key = f"{category}.{tag}"
        col = _as_list(mmcif_dict.get(key))
        columns[tag] = col
        if col:
            length = len(col) if length is None else length

    if length is None:
        return {tag: [] for tag in tags}

    for tag in tags:
        col = columns[tag]
        if len(col) != length:
            columns[tag] = (col * length)[:length] if col else ["?"] * length
    return columns


def _rows(columns: Dict[str, List[str]]) -> List[Dict[str, str]]:
    if not columns:
        return []
    length = len(next(iter(columns.values())))
    return [{k: v[i] for k, v in columns.items()} for i in range(length)]


def _first(mmcif_dict: dict, key: str, default: str = "N/A") -> str:
    values = _as_list(mmcif_dict.get(key))
    return values[0] if values and values[0] not in (None, "?", ".") else default


# --------------------------------------------------------------------------- #
# Header info: title, experimental method, resolution
# --------------------------------------------------------------------------- #

def parse_header_info(mmcif_dict: dict) -> Dict[str, str]:
    title = _first(mmcif_dict, "_struct.title")

    methods = [m for m in _as_list(mmcif_dict.get("_exptl.method")) if m not in (None, "?", ".")]
    method = "; ".join(methods) if methods else "N/A"

    # Resolution lives in different categories depending on experiment type.
    resolution = "N/A"
    for key in ("_refine.ls_d_res_high", "_em_3d_reconstruction.resolution", "_reflns.d_resolution_high"):
        val = _first(mmcif_dict, key, default=None)
        if val:
            resolution = val
            break

    return {"title": title, "method": method, "resolution": resolution}


# --------------------------------------------------------------------------- #
# EC number / classification
# --------------------------------------------------------------------------- #

def split_ec_numbers(raw: Optional[str]) -> List[str]:
    if not raw or raw in ("?", ".", "N/A"):
        return []
    return [ec.strip() for ec in raw.split(",") if ec.strip()]


def ec_broad_category(ec: str) -> str:
    first = ec.split(".")[0].strip()
    return EC_BROAD_CATEGORIES.get(first, "Unknown")


def _is_complete_ec(ec: str) -> bool:
    """KEGG only has entries for fully-specified EC numbers (four numeric
    levels); wildcarded ones like '3.4.-.-' aren't valid lookups."""
    parts = ec.split(".")
    return len(parts) == 4 and all(re.fullmatch(r"\d+", p) for p in parts)


def fetch_ec_kegg_detail(ec: str, cache: Dict[str, str]) -> str:
    """Full classification hierarchy + enzyme name for one EC number, via
    KEGG's ENZYME database, parsed with Biopython's own Bio.KEGG module."""
    if ec in cache:
        return cache[ec]

    if not _is_complete_ec(ec):
        cache[ec] = "N/A (incomplete EC number, not a valid KEGG lookup)"
        return cache[ec]

    try:
        from Bio.KEGG import REST, Enzyme
        handle = REST.kegg_get(f"ec:{ec}")
        record = next(Enzyme.parse(handle))
        name = record.name[0].rstrip(";").strip() if record.name else ""
        hierarchy = " > ".join(c.rstrip(";").strip() for c in record.classname) if record.classname else ""
        parts = [p for p in (hierarchy, name) if p]
        detail = " | ".join(parts) if parts else "N/A"
    except Exception as exc:  # noqa: BLE001 -- KEGG lookups fail for many benign reasons
        detail = f"N/A (KEGG lookup failed: {exc})"

    cache[ec] = detail
    return detail


NON_ENZYME = "Non-enzyme"


def enzyme_type(ec_numbers: List[str]) -> str:
    """Broad enzyme class, from the **first digit** of each EC number.

    An entity with no EC number is not an enzyme; saying so explicitly is more
    useful than a blank, because a blank is ambiguous between "not an enzyme"
    and "we did not look".
    """
    if not ec_numbers:
        return NON_ENZYME
    return "; ".join(sorted({ec_broad_category(ec) for ec in ec_numbers}))


def fetch_ec_name(ec: str, cache: Dict[str, str]) -> str:
    """The accepted enzyme name for a **complete four-level** EC number.

    Kept separate from :func:`fetch_ec_kegg_detail`, which returns the whole
    class hierarchy joined with the name. Here we want the name on its own so
    it can stand in its own column.
    """
    key = f"name:{ec}"
    if key in cache:
        return cache[key]

    if not _is_complete_ec(ec):
        cache[key] = f"N/A (EC {ec} is not fully specified)"
        return cache[key]

    try:
        from Bio.KEGG import REST, Enzyme
        handle = REST.kegg_get(f"ec:{ec}")
        record = next(Enzyme.parse(handle))
        name = record.name[0].rstrip(";").strip() if record.name else ""
        cache[key] = name or "N/A"
    except Exception as exc:  # noqa: BLE001 -- KEGG is often unreachable
        cache[key] = f"N/A (KEGG lookup failed: {exc})"
    return cache[key]


def enzyme_name(ec_numbers: List[str], entity: Optional["EntityInfo"],
                skip_lookup: bool, ec_cache: Dict[str, str]) -> str:
    """Name the enzyme, preferring KEGG's accepted name for the full EC.

    Falls back to the mmCIF entity description, which for an enzyme is almost
    always its name -- so the column stays populated when KEGG is unreachable
    or lookups are switched off. The fallback is labelled, so a reader can tell
    a curated EC name from a description scraped out of the file.
    """
    if not ec_numbers:
        return NON_ENZYME

    if not skip_lookup:
        names = [fetch_ec_name(ec, ec_cache) for ec in ec_numbers]
        good = [n for n in names if n and not n.startswith("N/A")]
        if good:
            return "; ".join(dict.fromkeys(good))

    description = (entity.description if entity else "") or ""
    if description and description not in ("?", ".", "N/A"):
        return f"{description} (from mmCIF entity description)"
    return "N/A"


# --------------------------------------------------------------------------- #
# Sequence / domain-family mmCIF categories (_struct_ref, _entity, ...)
# --------------------------------------------------------------------------- #

def parse_struct_ref(mmcif_dict: dict) -> List[Dict[str, str]]:
    tags = ["id", "entity_id", "db_name", "db_code", "pdbx_db_accession",
            "pdbx_seq_one_letter_code", "pdbx_align_begin"]
    return _rows(_column(mmcif_dict, "_struct_ref", tags))


def parse_struct_ref_seq(mmcif_dict: dict) -> List[Dict[str, str]]:
    tags = ["align_id", "ref_id", "pdbx_PDB_id_code", "pdbx_strand_id",
            "seq_align_beg", "seq_align_end",
            "db_align_beg", "db_align_end",
            "pdbx_auth_seq_align_beg", "pdbx_auth_seq_align_end"]
    return _rows(_column(mmcif_dict, "_struct_ref_seq", tags))


def parse_entity(mmcif_dict: dict) -> List[Dict[str, str]]:
    tags = ["id", "type", "src_method", "pdbx_description",
            "pdbx_number_of_molecules", "pdbx_ec"]
    return _rows(_column(mmcif_dict, "_entity", tags))


def parse_entity_poly(mmcif_dict: dict) -> List[Dict[str, str]]:
    tags = ["entity_id", "type", "pdbx_seq_one_letter_code_can", "pdbx_strand_id"]
    return _rows(_column(mmcif_dict, "_entity_poly", tags))


def parse_struct_descriptor(mmcif_dict: dict) -> Dict[str, str]:
    return {
        "title": _first(mmcif_dict, "_struct.title"),
        "pdbx_descriptor": _first(mmcif_dict, "_struct.pdbx_descriptor"),
    }


def parse_assembly(mmcif_dict: dict) -> List[Dict[str, str]]:
    tags = ["id", "details", "method_details", "oligomeric_details", "oligomeric_count"]
    return _rows(_column(mmcif_dict, "_pdbx_struct_assembly", tags))


@dataclass
class EntityInfo:
    entity_id: str
    description: str = "N/A"
    poly_type: str = "N/A"
    chains: List[str] = field(default_factory=list)
    ec_numbers: List[str] = field(default_factory=list)
    db_name: Optional[str] = None
    db_code: Optional[str] = None
    accession: Optional[str] = None
    alignments: List[tuple] = field(default_factory=list)


def build_entity_table(mmcif_dict: dict) -> Dict[str, EntityInfo]:
    entities: Dict[str, EntityInfo] = {}

    for row in parse_entity(mmcif_dict):
        eid = row["id"]
        entities[eid] = EntityInfo(
            entity_id=eid,
            description=row.get("pdbx_description", "N/A") or "N/A",
            ec_numbers=split_ec_numbers(row.get("pdbx_ec")),
        )

    for row in parse_entity_poly(mmcif_dict):
        eid = row["entity_id"]
        ent = entities.setdefault(eid, EntityInfo(entity_id=eid))
        ent.poly_type = row.get("type", "N/A")
        chains_raw = row.get("pdbx_strand_id", "")
        ent.chains = [c.strip() for c in chains_raw.split(",") if c.strip()] if chains_raw not in (None, "?", ".") else []

    struct_refs = {r["id"]: r for r in parse_struct_ref(mmcif_dict)}
    for ref_id, ref in struct_refs.items():
        eid = ref["entity_id"]
        ent = entities.setdefault(eid, EntityInfo(entity_id=eid))
        if ref.get("db_name") == "UNP" or (ref.get("db_name") or "").upper().startswith("UNIPROT"):
            ent.db_name = ref.get("db_name")
            ent.db_code = ref.get("db_code")
            ent.accession = ref.get("pdbx_db_accession")

    for row in parse_struct_ref_seq(mmcif_dict):
        ref = struct_refs.get(row["ref_id"])
        if not ref:
            continue
        eid = ref["entity_id"]
        ent = entities.setdefault(eid, EntityInfo(entity_id=eid))
        try:
            ent.alignments.append((
                int(row["seq_align_beg"]), int(row["seq_align_end"]),
                int(row["db_align_beg"]), int(row["db_align_end"]),
                int(row["pdbx_auth_seq_align_beg"]), int(row["pdbx_auth_seq_align_end"]),
            ))
        except (KeyError, ValueError):
            pass

    return entities


def entity_for_chain(entities: Dict[str, EntityInfo], chain_id: str) -> Optional[EntityInfo]:
    for ent in entities.values():
        if chain_id in ent.chains:
            return ent
    return None


def map_uniprot_pos_to_pdb(entity: EntityInfo, uniprot_pos: int) -> Optional[int]:
    for _, _, db_beg, db_end, auth_beg, auth_end in entity.alignments:
        if db_beg <= uniprot_pos <= db_end:
            return auth_beg + (uniprot_pos - db_beg)
    return None


# --------------------------------------------------------------------------- #
# UniProt domain / family retrieval
# --------------------------------------------------------------------------- #

def fetch_uniprot_json(accession: str, timeout: int = 20) -> dict:
    url = UNIPROT_REST_URL.format(accession=accession)
    req = urllib.request.Request(url, headers={"Accept": "application/json"})
    with urllib.request.urlopen(req, timeout=timeout) as resp:
        return json.loads(resp.read().decode("utf-8"))


def extract_domains(uniprot_json: dict) -> List[Dict[str, str]]:
    domains = []
    for feat in uniprot_json.get("features", []):
        if feat.get("type") in ("Domain", "Region", "Repeat"):
            loc = feat.get("location", {})
            start = loc.get("start", {}).get("value")
            end = loc.get("end", {}).get("value")
            domains.append({
                "feature_type": feat.get("type"),
                "description": feat.get("description", ""),
                "start": start,
                "end": end,
            })
    return domains


def extract_family_xrefs(uniprot_json: dict) -> List[Dict[str, str]]:
    xrefs = []
    for xref in uniprot_json.get("uniProtKBCrossReferences", []):
        if xref.get("database") in FAMILY_DOMAIN_DBS:
            props = {p["key"]: p["value"] for p in xref.get("properties", [])}
            xrefs.append({
                "database": xref.get("database"),
                "id": xref.get("id"),
                "description": props.get("EntryName") or props.get("MatchStatus") or "",
            })
    return xrefs


def extract_family_comment(uniprot_json: dict) -> Optional[str]:
    for comment in uniprot_json.get("comments", []):
        if comment.get("commentType") == "SIMILARITY":
            texts = comment.get("texts", [])
            if texts:
                return texts[0].get("value")
    return None


@dataclass
class UniProtSummary:
    family: str = "N/A"
    domains_pdb_numbering: str = "N/A"
    xrefs: str = "N/A"


def fetch_uniprot_summary(entity: EntityInfo, cache: Dict[str, "UniProtSummary"]) -> UniProtSummary:
    accession = entity.accession
    if not accession:
        return UniProtSummary()
    if accession in cache:
        base = cache[accession]
    else:
        try:
            data = fetch_uniprot_json(accession)
        except Exception as exc:  # noqa: BLE001
            base = UniProtSummary(family=f"N/A (UniProt lookup failed: {exc})")
            cache[accession] = base
            return base

        family = extract_family_comment(data) or "N/A"
        xrefs = extract_family_xrefs(data)
        xref_str = "; ".join(f"{x['database']}:{x['id']}" + (f"({x['description']})" if x["description"] else "")
                              for x in xrefs) or "N/A"
        base = UniProtSummary(family=family, xrefs=xref_str)
        # domains need per-entity coordinate mapping, so store raw domains on the cache entry too
        base._raw_domains = extract_domains(data)  # type: ignore[attr-defined]
        cache[accession] = base

    raw_domains = getattr(base, "_raw_domains", [])
    parts = []
    for d in raw_domains:
        if d["start"] is None or d["end"] is None:
            continue
        pdb_start = map_uniprot_pos_to_pdb(entity, d["start"])
        pdb_end = map_uniprot_pos_to_pdb(entity, d["end"])
        label = d["description"] or d["feature_type"]
        if pdb_start and pdb_end:
            parts.append(f"{label}:{pdb_start}-{pdb_end}")
        else:
            parts.append(f"{label}:{d['start']}-{d['end']}(UniProt numbering)")
    domains_str = "; ".join(parts) if parts else "N/A"

    return UniProtSummary(family=base.family, domains_pdb_numbering=domains_str, xrefs=base.xrefs)


# --------------------------------------------------------------------------- #
# Secondary structure: _struct_conf (helix/turn) + _struct_sheet_range (sheet)
# --------------------------------------------------------------------------- #

def parse_struct_conf(mmcif_dict: dict) -> List[Dict[str, str]]:
    tags = ["id", "conf_type_id", "beg_auth_asym_id", "beg_auth_seq_id",
            "end_auth_asym_id", "end_auth_seq_id"]
    return _rows(_column(mmcif_dict, "_struct_conf", tags))


def parse_struct_sheet_range(mmcif_dict: dict) -> List[Dict[str, str]]:
    tags = ["sheet_id", "id", "beg_auth_asym_id", "beg_auth_seq_id",
            "end_auth_asym_id", "end_auth_seq_id"]
    return _rows(_column(mmcif_dict, "_struct_sheet_range", tags))


def _classify_conf_type(conf_type_id: str) -> str:
    ct = (conf_type_id or "").upper()
    if ct.startswith("HELX"):
        return "Helix"
    if ct.startswith("TURN"):
        return "Turn"
    if ct.startswith("STRN") or ct.startswith("STRAND"):
        return "Sheet"
    return "Other"


@dataclass
class SSSegment:
    chain: str
    start: int
    end: int
    kind: str  # Helix / Turn / Sheet / Coil


def collect_ss_segments(mmcif_dict: dict) -> Dict[str, List[SSSegment]]:
    """Raw annotated (non-Coil) segments per chain, from _struct_conf and
    _struct_sheet_range, sorted by starting residue number."""
    per_chain: Dict[str, List[SSSegment]] = {}

    for row in parse_struct_conf(mmcif_dict):
        chain = row.get("beg_auth_asym_id")
        try:
            start, end = int(row["beg_auth_seq_id"]), int(row["end_auth_seq_id"])
        except (KeyError, ValueError, TypeError):
            continue
        kind = _classify_conf_type(row.get("conf_type_id"))
        if chain and chain not in ("?", "."):
            per_chain.setdefault(chain, []).append(SSSegment(chain, start, end, kind))

    for row in parse_struct_sheet_range(mmcif_dict):
        chain = row.get("beg_auth_asym_id")
        try:
            start, end = int(row["beg_auth_seq_id"]), int(row["end_auth_seq_id"])
        except (KeyError, ValueError, TypeError):
            continue
        if chain and chain not in ("?", "."):
            per_chain.setdefault(chain, []).append(SSSegment(chain, start, end, "Sheet"))

    for chain in per_chain:
        per_chain[chain].sort(key=lambda s: s.start)

    return per_chain


def get_chain_residue_range(structure, chain_id: str) -> Optional[Tuple[int, int]]:
    """Min/max author residue number of *standard* (non-heteroatom) residues
    in a chain, taken directly from the parsed Biopython Structure object."""
    if structure is None:
        return None
    try:
        model = next(iter(structure))
    except StopIteration:
        return None
    for chain in model:
        if chain.id != chain_id:
            continue
        seq_ids = [res.id[1] for res in chain if res.id[0] == " "]
        if not seq_ids:
            return None
        return min(seq_ids), max(seq_ids)
    return None


def fill_coil_gaps(segments: List[SSSegment], chain_range: Optional[Tuple[int, int]], chain: str) -> List[SSSegment]:
    """Insert explicit 'Coil' segments for any residue stretch not covered
    by an annotated Helix/Turn/Sheet element, bounded by the chain's actual
    observed residue range. mmCIF only annotates helices/turns/sheets
    explicitly; everything else is implicitly coil/loop."""
    if not segments:
        if chain_range:
            lo, hi = chain_range
            return [SSSegment(chain, lo, hi, "Coil")]
        return []

    filled: List[SSSegment] = []
    lo, hi = chain_range if chain_range else (segments[0].start, segments[-1].end)

    cursor = lo
    for seg in segments:
        if seg.start > cursor:
            filled.append(SSSegment(chain, cursor, seg.start - 1, "Coil"))
        filled.append(seg)
        cursor = max(cursor, seg.end + 1)
    if cursor <= hi:
        filled.append(SSSegment(chain, cursor, hi, "Coil"))
    return filled


def format_ss_order(segments: List[SSSegment]) -> str:
    return "-".join(s.kind for s in segments) if segments else "N/A"


def format_ss_ranges(segments: List[SSSegment]) -> str:
    return "; ".join(f"{s.kind}:{s.start}-{s.end}" for s in segments) if segments else "N/A"


# --------------------------------------------------------------------------- #
# Heteroatoms + bonding partners (_struct_conn), per chain
# --------------------------------------------------------------------------- #

def parse_struct_conn(mmcif_dict: dict) -> List[Dict[str, str]]:
    tags = ["id", "conn_type_id",
            "ptnr1_auth_asym_id", "ptnr1_auth_comp_id", "ptnr1_auth_seq_id", "ptnr1_label_atom_id",
            "ptnr2_auth_asym_id", "ptnr2_auth_comp_id", "ptnr2_auth_seq_id", "ptnr2_label_atom_id"]
    return _rows(_column(mmcif_dict, "_struct_conn", tags))


def get_heteroatoms_per_chain(structure) -> Dict[str, List[Tuple[str, str]]]:
    """(resname, resseq) heteroatoms per chain, from the Biopython Structure
    object. Water (Biopython hetflag 'W') is excluded -- it isn't normally
    what's meant by a bonding-relevant heteroatom/ligand."""
    result: Dict[str, List[Tuple[str, str]]] = {}
    if structure is None:
        return result
    try:
        model = next(iter(structure))
    except StopIteration:
        return result
    for chain in model:
        hets = []
        for res in chain:
            hetflag = res.id[0]
            if hetflag not in (" ", "W"):
                hets.append((res.resname.strip(), str(res.id[1])))
        if hets:
            result[chain.id] = hets
    return result


def build_heteroatom_bonds(struct_conn_rows: List[Dict[str, str]],
                            heteroatoms_per_chain: Dict[str, List[Tuple[str, str]]]) -> Dict[str, Dict[str, List[str]]]:
    """For each chain's heteroatoms, find their bonding partner residues
    from _struct_conn (covalent/disulfide/metal-coordination/H-bond),
    regardless of which side of the connection record the heteroatom is on."""
    het_lookup = set()
    for chain, hets in heteroatoms_per_chain.items():
        for resname, resseq in hets:
            het_lookup.add((chain, resname, resseq))

    bonds: Dict[str, Dict[str, List[str]]] = {chain: {f"{n}{s}": [] for n, s in hets}
                                               for chain, hets in heteroatoms_per_chain.items()}

    for row in struct_conn_rows:
        conn_type = row.get("conn_type_id", "?")
        p1 = (row.get("ptnr1_auth_asym_id"), row.get("ptnr1_auth_comp_id"), row.get("ptnr1_auth_seq_id"))
        p2 = (row.get("ptnr2_auth_asym_id"), row.get("ptnr2_auth_comp_id"), row.get("ptnr2_auth_seq_id"))

        if p1 in het_lookup:
            chain, resname, resseq = p1
            label = f"{resname}{resseq}"
            other_chain, other_res, other_seq = p2
            if all((other_chain, other_res, other_seq)):
                bonds.setdefault(chain, {}).setdefault(label, []).append(
                    f"{other_res}{other_seq}(chain {other_chain}, {conn_type})")

        if p2 in het_lookup:
            chain, resname, resseq = p2
            label = f"{resname}{resseq}"
            other_chain, other_res, other_seq = p1
            if all((other_chain, other_res, other_seq)):
                bonds.setdefault(chain, {}).setdefault(label, []).append(
                    f"{other_res}{other_seq}(chain {other_chain}, {conn_type})")

    return bonds


def format_heteroatoms(hets: List[Tuple[str, str]]) -> str:
    return "; ".join(f"{name}{seq}" for name, seq in hets) if hets else "N/A"


def format_heteroatom_bonds(bond_map: Dict[str, List[str]]) -> str:
    if not bond_map:
        return "N/A"
    parts = []
    for het_label, partners in bond_map.items():
        if partners:
            parts.append(f"{het_label} -> {', '.join(partners)}")
        else:
            parts.append(f"{het_label} -> (no _struct_conn record found)")
    return "; ".join(parts)


# --------------------------------------------------------------------------- #
# Structure loading and per-structure processing (pdbKit)
# --------------------------------------------------------------------------- #

def discover_structure_files(folder: str) -> List[str]:
    files = []
    for name in sorted(os.listdir(folder)):
        low = name.lower()
        if low.endswith(STRUCTURE_FILE_EXTENSIONS) or low.endswith(tuple(e + ".gz" for e in STRUCTURE_FILE_EXTENSIONS)):
            files.append(os.path.join(folder, name))
    return files


def build_rows_for_structure(pdb_id: str, structure, mmcif_dict: Optional[dict],
                             skip_uniprot: bool, ec_cache: dict, uniprot_cache: dict) -> List[Dict[str, str]]:
    """One CSV row per chain. ``mmcif_dict`` is the file's own dict (mmCIF) or
    the translated PDB header (see :func:`pdbkit.metadata.as_mmcif_dict`)."""
    if mmcif_dict is None:
        from .metadata import as_mmcif_dict
        mmcif_dict = as_mmcif_dict(structure)
    header = parse_header_info(mmcif_dict)
    title, method, resolution = header["title"], header["method"], header["resolution"]
    entities = build_entity_table(mmcif_dict)
    raw_ss = collect_ss_segments(mmcif_dict)
    ss_by_chain = {
        chain: fill_coil_gaps(segs, get_chain_residue_range(structure, chain), chain)
        for chain, segs in raw_ss.items()
    }
    hets_per_chain = get_heteroatoms_per_chain(structure)
    het_bonds = build_heteroatom_bonds(parse_struct_conn(mmcif_dict), hets_per_chain)

    try:
        model = next(iter(structure))
        chain_ids = [c.id for c in model]
    except StopIteration:
        chain_ids = []

    rows = []
    for chain_id in chain_ids:
        entity = entity_for_chain(entities, chain_id)

        ec_numbers = entity.ec_numbers if entity else []
        ec_number_str = "; ".join(ec_numbers) if ec_numbers else "N/A"
        ec_type_str = "; ".join(sorted({ec_broad_category(ec) for ec in ec_numbers})) if ec_numbers else "N/A"
        if ec_numbers and not skip_uniprot:
            ec_detail_str = "; ".join(fetch_ec_kegg_detail(ec, ec_cache) for ec in ec_numbers)
        elif ec_numbers:
            ec_detail_str = "N/A (skipped, --no-uniprot also skips KEGG lookups)"
        else:
            ec_detail_str = "N/A"

        enzyme_type_str = enzyme_type(ec_numbers)
        enzyme_name_str = enzyme_name(ec_numbers, entity, skip_uniprot, ec_cache)

        if entity and entity.accession and not skip_uniprot:
            summary = fetch_uniprot_summary(entity, uniprot_cache)
        else:
            summary = UniProtSummary()

        segments = ss_by_chain.get(chain_id, [])

        rows.append({
            "PDB_ID": pdb_id,
            "Title": title,
            "Experimental_Method": method,
            "Resolution": resolution,
            "Chain_ID": chain_id,
            "Entity_ID": entity.entity_id if entity else "N/A",
            "Entity_Description": entity.description if entity else "N/A",
            "EC_Number": ec_number_str,
            "EC_Type": ec_type_str,
            "EC_Detailed_Classification": ec_detail_str,
            "Enzyme Type": enzyme_type_str,
            "Enzyme Name": enzyme_name_str,
            "UniProt_Accession": entity.accession if entity and entity.accession else "N/A",
            "UniProt_Family": summary.family,
            "UniProt_Domains_PDB_Numbering": summary.domains_pdb_numbering,
            "Family_Domain_DB_Xrefs": summary.xrefs,
            "Secondary_Structure_Order": format_ss_order(segments),
            "Secondary_Structure_Residue_Ranges": format_ss_ranges(segments),
            "Heteroatoms_In_Chain": format_heteroatoms(hets_per_chain.get(chain_id, [])),
            "Heteroatom_Bonding_Partners": format_heteroatom_bonds(het_bonds.get(chain_id, {})),
        })
    return rows


def write_csv(path: str, rows: List[Dict[str, str]]) -> None:
    with open(path, "w", newline="") as fh:
        writer = csv.DictWriter(fh, fieldnames=CSV_FIELDNAMES)
        writer.writeheader()
        for row in rows:
            writer.writerow(row)


def annotate(source, *, skip_lookups: bool = False, ec_cache: Optional[dict] = None,
             uniprot_cache: Optional[dict] = None) -> List[Dict[str, str]]:
    """Rows for one structure: a Structure, a file (.cif/.pdb, optionally .gz) or a PDB ID."""
    from .structure import load, source_info
    structure = load(source)
    info = source_info(structure)
    pdb_id = (info.entry_id or str(structure.id)).upper()
    return build_rows_for_structure(pdb_id, structure, info.mmcif, skip_lookups,
                                    {} if ec_cache is None else ec_cache,
                                    {} if uniprot_cache is None else uniprot_cache)


def run(inputs, outdir: str = ".", skip_lookups: bool = False, log=print) -> List[str]:
    """Process IDs, files or directories; write one CSV per structure. Returns paths."""
    os.makedirs(outdir, exist_ok=True)
    ec_cache: Dict[str, str] = {}
    uniprot_cache: Dict[str, UniProtSummary] = {}
    items: List[str] = []
    for inp in ([inputs] if isinstance(inputs, str) else inputs):
        if os.path.isdir(inp):
            found = discover_structure_files(inp)
            if not found:
                raise ValueError(f"No structure files ({', '.join(STRUCTURE_FILE_EXTENSIONS)}) in {inp!r}")
            items.extend(found)
        else:
            items.append(inp)
    written = []
    for item in items:
        try:
            from .structure import load, source_info
            structure = load(item)
            info = source_info(structure)
            pdb_id = (info.entry_id or str(structure.id)).upper()
            rows = build_rows_for_structure(pdb_id, structure, info.mmcif, skip_lookups, ec_cache, uniprot_cache)
            path = os.path.join(outdir, f"{pdb_id}_domain_ss_info.csv")
            write_csv(path, rows)
            written.append(path)
            log(f"  {item} -> {path}")
        except Exception as exc:  # noqa: BLE001 - keep a batch going
            if len(items) == 1:
                raise
            print(f"  [error] {item}: {exc}", file=sys.stderr)
    return written
