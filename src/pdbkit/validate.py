"""File-level validation for PDB and mmCIF (extends pdb-tools ``pdb_validate``).

* **PDB** - every ATOM/HETATM line is checked field by field against the
  wwPDB column layout (the pdb_validate checks). Lines that are not exactly
  80 characters are reported as *warnings* rather than errors: pdb_validate
  failed on them, but short lines are common in valid files and every viewer
  accepts them.
* **mmCIF** - the file must parse; ``_atom_site`` must have the columns
  viewers rely on, equal-length columns, numeric coordinates and unique atom
  ids.
* **both** - after parsing: models with different atom sets (ensemble check).
"""

from __future__ import annotations

import re
from typing import Dict, List

from .structure import PdbKitError, _read_text, detect_format, format_from_path, loads

__all__ = ["validate_file", "validate_text"]

_PDB_FIELDS = (
    ("atom serial", slice(6, 11), re.compile(r"^[\d\s]+$|^[0-9A-Za-z]{5}$")),
    ("spacer", slice(11, 12), re.compile(r"^\s$")),
    ("atom name", slice(12, 16), re.compile(r"^\s*[A-Za-z0-9'\"*+-]+\s*$")),
    ("altloc", slice(16, 17), re.compile(r"^[A-Za-z0-9 ]$")),
    ("residue name", slice(17, 20), re.compile(r"^\s*[A-Za-z0-9+-]+\s*$")),
    ("chain", slice(21, 22), re.compile(r"^[A-Za-z0-9 ]$")),
    ("residue number", slice(22, 26), re.compile(r"^\s*-?\d+\s*$")),
    ("insertion code", slice(26, 27), re.compile(r"^[A-Za-z0-9 ]$")),
    ("x", slice(30, 38), re.compile(r"^\s*-?\d+\.\d+\s*$")),
    ("y", slice(38, 46), re.compile(r"^\s*-?\d+\.\d+\s*$")),
    ("z", slice(46, 54), re.compile(r"^\s*-?\d+\.\d+\s*$")),
    ("occupancy", slice(54, 60), re.compile(r"^\s*-?\d+\.\d+\s*$|^\s*$")),
    ("B-factor", slice(60, 66), re.compile(r"^\s*-?\d+\.\d+\s*$|^\s*$")),
    ("element", slice(76, 78), re.compile(r"^[\sA-Za-z]*$")),
    ("charge", slice(78, 80), re.compile(r"^(\s*|\d[+-]|[+-]\d)$")),
)


def _check_pdb(text: str) -> Dict[str, List[str]]:
    errors, warnings = [], []
    short = long = 0
    for n, raw in enumerate(text.splitlines(), start=1):
        line = raw.rstrip("\r")
        if not line.strip():
            continue
        if len(line) < 80:
            short += 1
        elif len(line) > 80:
            long += 1
        if line[:6] in ("ATOM  ", "HETATM"):
            padded = line.ljust(80)
            if len(line) < 54:
                errors.append(f"line {n}: coordinate record too short ({len(line)} < 54 columns)")
                continue
            for name, sl, rx in _PDB_FIELDS:
                if not rx.match(padded[sl]):
                    errors.append(f"line {n}: bad {name} field {padded[sl]!r} (columns {sl.start + 1}-{sl.stop})")
                    break
    if short:
        warnings.append(f"{short} line(s) shorter than 80 columns (accepted by viewers; strict spec wants 80)")
    if long:
        errors.append(f"{long} line(s) longer than 80 columns")
    return {"errors": errors, "warnings": warnings}


def _check_cif(text: str) -> Dict[str, List[str]]:
    from Bio.PDB.MMCIF2Dict import MMCIF2Dict
    import io
    errors, warnings = [], []
    try:
        d = MMCIF2Dict(io.StringIO(text))
    except Exception as exc:  # noqa: BLE001
        return {"errors": [f"mmCIF does not parse: {exc}"], "warnings": []}
    required = ["group_PDB", "id", "type_symbol", "label_atom_id", "label_comp_id", "label_asym_id",
                "Cartn_x", "Cartn_y", "Cartn_z"]
    recommended = ["auth_asym_id", "auth_seq_id", "label_seq_id", "label_entity_id", "occupancy",
                   "B_iso_or_equiv", "pdbx_PDB_model_num"]
    cols = {k.split(".", 1)[1]: (v if isinstance(v, list) else [v])
            for k, v in d.items() if k.startswith("_atom_site.")}
    if not cols:
        return {"errors": ["no _atom_site category"], "warnings": []}
    for c in required:
        if c not in cols:
            errors.append(f"_atom_site.{c} missing")
    for c in recommended:
        if c not in cols:
            warnings.append(f"_atom_site.{c} missing (recommended)")
    lengths = {len(v) for v in cols.values()}
    if len(lengths) > 1:
        errors.append(f"_atom_site columns have different lengths: {sorted(lengths)}")
    for c in ("Cartn_x", "Cartn_y", "Cartn_z"):
        bad = sum(1 for v in cols.get(c, []) if not re.fullmatch(r"-?\d+(\.\d*)?([eE][-+]?\d+)?", v))
        if bad:
            errors.append(f"{bad} non-numeric value(s) in _atom_site.{c}")
    ids = cols.get("id", [])
    if len(set(ids)) != len(ids):
        errors.append("duplicate _atom_site.id values")
    if "_entity.id" not in d:
        warnings.append("no _entity category (viewers will guess polymer/ligand/water)")
    return {"errors": errors, "warnings": warnings}


def validate_text(text: str, fmt: str = None) -> Dict[str, object]:
    fmt = fmt or detect_format(text)
    result = _check_cif(text) if fmt == "cif" else _check_pdb(text)
    try:
        s = loads(text, fmt)
        from .edit import check_ensemble
        result["warnings"].extend(check_ensemble(s))
        result["atoms"] = sum(1 for _ in s.get_atoms())
        result["models"] = len(s)
    except PdbKitError as exc:
        result["errors"].append(str(exc))
    except Exception as exc:  # noqa: BLE001
        result["errors"].append(f"parser failed: {exc}")
    result["format"] = fmt
    result["ok"] = not result["errors"]
    return result


def validate_file(path: str) -> Dict[str, object]:
    text = _read_text(path)
    return validate_text(text, format_from_path(path) or detect_format(text))
