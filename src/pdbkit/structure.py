"""Reading, writing and classifying macromolecular structures.

Every pdbKit operation works on a Biopython :class:`~Bio.PDB.Structure.Structure`.
This module is the single place where files become structures and back, so
PDB and mmCIF input are handled by the same code path everywhere else:

* :func:`load` accepts a path (``.pdb``, ``.ent``, ``.pdbN``, ``.cif``,
  ``.mmcif``, optionally gzipped), ``"-"`` for stdin, a PDB ID or an
  AlphaFold/UniProt ID (fetched and cached), or an existing Structure.
* :func:`save` / :func:`dumps` write PDB (via Biopython's ``PDBIO``) or mmCIF
  (via pdbKit's own writer in :mod:`pdbkit.cifio`, which also writes the
  ``_entity``/``_struct_asym``/``_chem_comp`` categories viewers use).
* Each residue is classified once, at load time, as ``polymer``, ``ligand``
  or ``water``. The classification is what "all atoms", "C-alpha" and
  "heteroatoms" mean throughout the package (see :func:`residue_kind`).
"""

from __future__ import annotations

import gzip
import io
import os
import re
import sys
import warnings
from dataclasses import dataclass, field
from pathlib import Path
from typing import Dict, Iterable, List, Optional, Tuple, Union

import numpy as np
from Bio.PDB.MMCIFParser import MMCIFParser
from Bio.PDB.PDBExceptions import PDBConstructionWarning
from Bio.PDB.PDBIO import PDBIO
from Bio.PDB.PDBParser import PDBParser
from Bio.PDB.Structure import Structure

from .constants import AMINO_ACIDS, NUCLEOTIDES, WATER_NAMES

__all__ = [
    "PdbKitError", "SourceInfo", "load", "loads", "save", "dumps",
    "detect_format", "source_info", "residue_kind", "polymer_type",
    "annotate", "structure_id",
]

PathLike = Union[str, os.PathLike]
StructureLike = Union[Structure, PathLike]

PDB_EXTENSIONS = (".pdb", ".ent", ".brk")
CIF_EXTENSIONS = (".cif", ".mmcif")
_PDB_ASSEMBLY_EXT = re.compile(r"\.pdb\d+$", re.IGNORECASE)
COORD_RECORDS = ("ATOM  ", "HETATM", "ANISOU", "TER   ", "MODEL ", "ENDMDL", "END   ", "SIGATM", "SIGUIJ")

PEPTIDE_BOND_MAX = 2.0   # Angstrom; C(i)-N(i+1)
PHOSPHO_BOND_MAX = 2.0   # Angstrom; O3'(i)-P(i+1)


class PdbKitError(Exception):
    """Raised for user-facing errors (bad input, impossible output, ...)."""


@dataclass
class SourceInfo:
    """Where a structure came from, plus the raw header it was read with.

    ``header_lines`` holds every non-coordinate record of a PDB file;
    ``mmcif`` holds the full ``MMCIF2Dict`` of an mmCIF file. Metadata
    functions (:mod:`pdbkit.metadata`) read whichever one is present.
    """

    fmt: str
    path: Optional[str] = None
    entry_id: str = ""
    header_lines: List[str] = field(default_factory=list)
    mmcif: Optional[dict] = None


# --------------------------------------------------------------------------- #
# Format detection
# --------------------------------------------------------------------------- #

def _strip_gz(name: str) -> str:
    return name[:-3] if name.lower().endswith(".gz") else name


def format_from_path(path: PathLike) -> Optional[str]:
    """``'pdb'`` / ``'cif'`` from a file name, or ``None`` if the extension is unknown."""
    name = _strip_gz(os.fspath(path)).lower()
    if name.endswith(CIF_EXTENSIONS):
        return "cif"
    if name.endswith(PDB_EXTENSIONS) or _PDB_ASSEMBLY_EXT.search(name):
        return "pdb"
    return None


def detect_format(text: str) -> str:
    """Sniff PDB vs mmCIF from file content."""
    for line in text.splitlines():
        s = line.strip()
        if not s or s.startswith("#"):
            continue
        if s.startswith("data_") or s.startswith("loop_") or s.startswith("_"):
            return "cif"
        return "pdb"
    raise PdbKitError("Empty input: no structure data found.")


def normalize_format(fmt: Optional[str]) -> Optional[str]:
    if fmt is None:
        return None
    f = fmt.lower().lstrip(".")
    if f in ("pdb", "ent", "brk"):
        return "pdb"
    if f in ("cif", "mmcif", "pdbx"):
        return "cif"
    raise PdbKitError(f"Unknown structure format {fmt!r}; use 'pdb' or 'cif'.")


# --------------------------------------------------------------------------- #
# Loading
# --------------------------------------------------------------------------- #

def _read_text(path: PathLike) -> str:
    p = os.fspath(path)
    if p.lower().endswith(".gz"):
        with gzip.open(p, "rt", errors="replace") as fh:
            return fh.read()
    with open(p, "r", errors="replace") as fh:
        return fh.read()


def loads(text: str, fmt: Optional[str] = None, structure_id: str = "",
          path: Optional[str] = None) -> Structure:
    """Parse a structure from a string (PDB or mmCIF; sniffed when ``fmt`` is None)."""
    fmt = normalize_format(fmt) or detect_format(text)
    with warnings.catch_warnings():
        warnings.simplefilter("ignore", PDBConstructionWarning)
        if fmt == "cif":
            parser = MMCIFParser(QUIET=True)
            structure = parser.get_structure(structure_id or "structure", io.StringIO(text))
            mmcif = dict(parser._mmcif_dict)
            info = SourceInfo(fmt="cif", path=path, mmcif=mmcif)
            eid = mmcif.get("_entry.id")
            eid = eid[0] if isinstance(eid, list) else eid
            info.entry_id = (eid or "").upper() if eid not in (None, "?", ".") else ""
        else:
            parser = PDBParser(PERMISSIVE=True, QUIET=True)
            structure = parser.get_structure(structure_id or "structure", io.StringIO(text))
            header = [ln.rstrip("\n") for ln in text.splitlines()
                      if ln[:6].ljust(6) not in COORD_RECORDS]
            info = SourceInfo(fmt="pdb", path=path, header_lines=header)
            for ln in header:
                if ln.startswith("HEADER") and len(ln) >= 66:
                    info.entry_id = ln[62:66].strip().upper()
                    break
    if not structure_id:
        stem = Path(_strip_gz(path)).name.split(".")[0] if path else ""
        structure.id = info.entry_id or stem or "structure"
    structure.xtra["pdbkit_source"] = info
    for model in structure:        # PDBParser numbers a MODEL-less file 0, MMCIFParser 1: use 1-based
        if not getattr(model, "serial_num", None):
            model.serial_num = model.id + 1
    if len(structure) == 0 or not any(True for _ in structure.get_atoms()):
        raise PdbKitError(f"No atoms found in {path or 'input'}.")
    annotate(structure)
    return structure


def load(source: StructureLike, fmt: Optional[str] = None, *, structure_id: str = "",
         fetch_missing: bool = True, cache_dir: Optional[PathLike] = None) -> Structure:
    """Load a structure from almost anything.

    Parameters
    ----------
    source:
        A Biopython Structure (returned unchanged), a file path, ``"-"`` for
        stdin, a 4-character PDB ID (downloaded from RCSB as mmCIF), or an
        AlphaFold DB model (``AF-P69905-F1`` or ``AF-P69905``).
    fmt:
        Force ``'pdb'`` or ``'cif'`` parsing; normally inferred.
    fetch_missing:
        When ``source`` is not an existing file but looks like an ID, download it.
    cache_dir:
        Where fetched files are kept (default ``~/.cache/pdbkit``).
    """
    if isinstance(source, Structure):
        if "pdbkit_source" not in source.xtra:
            source.xtra["pdbkit_source"] = SourceInfo(fmt="pdb")
            annotate(source)
        return source
    src = os.fspath(source)
    if src == "-":
        text = sys.stdin.read()
        return loads(text, fmt, structure_id=structure_id)
    if os.path.isfile(src):
        text = _read_text(src)
        fmt = normalize_format(fmt) or format_from_path(src) or detect_format(text)
        return loads(text, fmt, structure_id=structure_id, path=src)
    if fetch_missing:
        import importlib
        _fetch = importlib.import_module("pdbkit.fetch")
        kind = _fetch.classify_identifier(src)
        if kind is not None:
            path = _fetch.fetch(src, fmt=normalize_format(fmt) or "cif",
                                outdir=cache_dir or _fetch.default_cache_dir())
            return load(path, fmt, structure_id=structure_id, fetch_missing=False)
    raise PdbKitError(f"{src!r} is not an existing file or a recognisable PDB/AlphaFold ID.")


def source_info(structure: Structure) -> SourceInfo:
    info = structure.xtra.get("pdbkit_source")
    if info is None:
        info = SourceInfo(fmt="pdb")
        structure.xtra["pdbkit_source"] = info
    return info


def structure_id(structure: Structure) -> str:
    info = source_info(structure)
    return info.entry_id or str(structure.id)


# --------------------------------------------------------------------------- #
# Residue classification
# --------------------------------------------------------------------------- #

def _near(a, b, cutoff: float) -> bool:
    return float(np.linalg.norm(a.coord - b.coord)) <= cutoff


def _linked(prev, res, nxt) -> bool:
    """True when ``res`` is covalently chained to a neighbouring residue.

    Backbone-directed on purpose: the previous residue's carbonyl C bonded to
    any atom of ``res`` (peptide bond, or an ester/amide C-terminal cap), or
    any atom of ``res`` bonded to the next residue's N (N-terminal caps such
    as ACE). Glycans or ligands attached to side chains are not "chained".
    """
    if prev is not None and "C" in prev and \
            any(_near(prev["C"], a, PEPTIDE_BOND_MAX) for a in res if a.element != "H"):
        return True
    if nxt is not None and "N" in nxt and \
            any(_near(a, nxt["N"], PEPTIDE_BOND_MAX) for a in res if a.element != "H"):
        return True
    if "P" in res and prev is not None and "O3'" in prev and _near(prev["O3'"], res["P"], PHOSPHO_BOND_MAX):
        return True
    if "O3'" in res and nxt is not None and "P" in nxt and _near(res["O3'"], nxt["P"], PHOSPHO_BOND_MAX):
        return True
    return False


def _label_seq_map(mmcif: dict) -> Dict[Tuple[str, str, str, str, str], Optional[int]]:
    """(model, auth_asym, auth_seq, ins_code, comp) -> label_seq_id (None for non-polymer)."""
    def col(tag):
        v = mmcif.get(f"_atom_site.{tag}")
        return v if isinstance(v, list) else ([v] if v is not None else [])
    label_seq = col("label_seq_id")
    if not label_seq:
        return {}
    n = len(label_seq)
    asym = col("auth_asym_id") or col("label_asym_id")
    seq = col("auth_seq_id") or col("label_seq_id")
    ins = col("pdbx_PDB_ins_code") or ["?"] * n
    comp = col("label_comp_id") or col("auth_comp_id")   # Biopython names residues by label_comp_id
    model = col("pdbx_PDB_model_num") or ["1"] * n
    out: Dict[Tuple[str, str, str, str, str], Optional[int]] = {}
    for i in range(n):
        ic = ins[i] if ins[i] not in ("?", ".") else " "
        key = (model[i], asym[i], seq[i], ic, comp[i])
        if key in out:
            continue
        ls = label_seq[i]
        out[key] = int(ls) if ls not in ("?", ".") and ls.lstrip("-").isdigit() else None
    return out


def polymer_type(residue) -> Optional[str]:
    """``'protein'``, ``'nucleic'`` or ``None`` for a polymer residue."""
    name = residue.get_resname().strip()
    if name in NUCLEOTIDES or ("P" in residue and "C1'" in residue) or ("C1'" in residue and "O4'" in residue):
        return "nucleic"
    if name in AMINO_ACIDS or "CA" in residue:
        return "protein"
    return None


def annotate(structure: Structure) -> None:
    """Classify every residue as polymer / ligand / water (stored in ``residue.xtra``).

    * **water** - residue name in :data:`pdbkit.constants.WATER_NAMES`.
    * **polymer** - for mmCIF input, any residue with a ``label_seq_id``
      (the file's own statement that it belongs to a polymer). Otherwise:
      ``ATOM`` records of standard residues, plus ``HETATM`` residues that are
      covalently chained to a neighbour (peptide C-N or phosphodiester O3'-P
      within 2.0 A) - this is how modified residues such as MSE are kept in
      the chain.
    * **ligand** - everything else (the "heteroatoms" of the CLI options).
    """
    info = structure.xtra.get("pdbkit_source")
    seqmap = _label_seq_map(info.mmcif) if info is not None and info.mmcif else {}
    for model in structure:
        mnum = str(getattr(model, "serial_num", model.id + 1))
        for chain in model:
            residues = chain.get_unpacked_list()      # all microheterogeneity variants
            for i, res in enumerate(residues):
                name = res.get_resname().strip()
                hetflag, resseq, icode = res.id
                label_seq = None
                known = False
                if seqmap:
                    key = (mnum, chain.id, str(resseq), icode, name)
                    if key in seqmap:
                        known = True
                        label_seq = seqmap[key]
                if name in WATER_NAMES or hetflag == "W":
                    kind = "water"
                elif known:
                    kind = "polymer" if label_seq is not None else "ligand"
                elif hetflag == " " and (name in AMINO_ACIDS or name in NUCLEOTIDES):
                    kind = "polymer"
                else:
                    prev = residues[i - 1] if i > 0 else None
                    nxt = residues[i + 1] if i + 1 < len(residues) else None
                    kind = "polymer" if _linked(prev, res, nxt) else "ligand"
                res.xtra["pdbkit_kind"] = kind
                res.xtra["pdbkit_poly"] = polymer_type(res) if kind == "polymer" else None
                if label_seq is not None:
                    res.xtra["label_seq_id"] = label_seq


def residue_kind(residue) -> str:
    """``'polymer'``, ``'ligand'`` or ``'water'`` (classifies on demand if needed)."""
    kind = residue.xtra.get("pdbkit_kind")
    if kind is None:
        name = residue.get_resname().strip()
        if name in WATER_NAMES or residue.id[0] == "W":
            kind = "water"
        elif name in AMINO_ACIDS or name in NUCLEOTIDES:
            kind = "polymer"
        else:
            kind = "ligand"
        residue.xtra["pdbkit_kind"] = kind
    return kind


# --------------------------------------------------------------------------- #
# Writing
# --------------------------------------------------------------------------- #

_SINGLE_CHAIN_IDS = [chr(c) for c in range(ord("A"), ord("Z") + 1)] + \
                    [chr(c) for c in range(ord("a"), ord("z") + 1)] + \
                    [str(d) for d in range(10)]


def _pdb_compatibility_problems(structure: Structure) -> List[str]:
    problems = []
    long_chains = sorted({c.id for c in structure.get_chains() if len(c.id) > 1})
    if long_chains:
        problems.append(f"chain IDs longer than 1 character: {', '.join(long_chains[:10])}")
    long_res = sorted({r.get_resname().strip() for r in iter_residues_all(structure) if len(r.get_resname().strip()) > 3})
    if long_res:
        problems.append(f"residue names longer than 3 characters: {', '.join(long_res[:10])}")
    bad_num = [r for r in iter_residues_all(structure) if not (-999 <= r.id[1] <= 9999)]
    if bad_num:
        problems.append(f"{len(bad_num)} residue number(s) outside -999..9999")
    long_atoms = sorted({a.get_id() for a in structure.get_atoms() if len(a.get_id()) > 4})
    if long_atoms:
        problems.append(f"atom names longer than 4 characters: {', '.join(long_atoms[:10])}")
    n_atoms = sum(1 for _ in iter_atoms_all_altlocs(structure))
    per_model = n_atoms / max(1, len(structure))
    if per_model > 99999:
        problems.append(f"{int(per_model)} atoms per model (PDB format limit is 99999)")
    return problems


def remap_chain_ids(structure: Structure) -> Dict[str, str]:
    """Give every chain a unique single-character ID (in place). Returns old->new."""
    mapping: Dict[str, str] = {}
    used = {c.id for c in structure.get_chains() if len(c.id) == 1}
    pool = [c for c in _SINGLE_CHAIN_IDS if c not in used]
    for chain in structure.get_chains():
        if len(chain.id) <= 1:
            continue
        if chain.id not in mapping:
            if not pool:
                raise PdbKitError("More than 62 chains: cannot map chain IDs to single characters; write mmCIF instead.")
            mapping[chain.id] = pool.pop(0)
    for model in structure:
        for chain in list(model):
            if chain.id in mapping:
                chain.id = mapping[chain.id]   # Entity.id setter re-keys the parent
    return mapping


def _dump_pdb(structure: Structure, remap_chains: bool = False) -> str:
    problems = _pdb_compatibility_problems(structure)
    if problems and remap_chains and all(p.startswith("chain IDs") for p in problems):
        structure = deep_copy(structure)
        mapping = remap_chain_ids(structure)
        print("pdbKit: remapped chain IDs for PDB output: " +
              ", ".join(f"{k}->{v}" for k, v in mapping.items()), file=sys.stderr)
        problems = _pdb_compatibility_problems(structure)
    if problems:
        raise PdbKitError("This structure cannot be written in PDB format (" + "; ".join(problems) +
                          "). Write mmCIF instead (--format cif / a .cif file name)" +
                          (", or pass --remap-chains." if any(p.startswith("chain IDs") for p in problems) else "."))
    pdbio = PDBIO()
    pdbio.set_structure(structure)
    buf = io.StringIO()
    preserve = bool(structure.xtra.get("pdbkit_preserve_serials"))
    with warnings.catch_warnings():
        warnings.simplefilter("ignore")
        pdbio.save(buf, preserve_atom_numbering=preserve)
    return _fix_ter_records(buf.getvalue(), structure, preserve)


def _fix_ter_records(text: str, structure: Structure, preserve: bool) -> str:
    """Place TER after each chain's last polymer residue and pad records to 80 columns.

    Biopython's PDBIO writes TER after the chain's last residue (often a water)
    and emits an 81-column TER line. The wwPDB format puts TER at the end of
    the polymer and has 80-column records; atom serials count the TER records.
    """
    last_poly: Dict[Tuple[int, str], Tuple[int, str]] = {}
    for model in structure:
        mnum = getattr(model, "serial_num", None) or model.id + 1
        for chain in model:
            for r in chain.get_unpacked_list():
                if residue_kind(r) == "polymer":
                    last_poly[(mnum, chain.id)] = (r.id[1], r.id[2])
    out: List[str] = []
    serial = 0
    first = next(iter(structure), None)
    mnum = (getattr(first, "serial_num", None) or 1) if first is not None else 1
    lines = [ln for ln in text.splitlines() if not ln.startswith("TER")]
    for i, ln in enumerate(lines):
        rec = ln[:6]
        if rec == "MODEL ":
            mnum = int(ln[10:14])
            serial = 0
        if rec in ("ATOM  ", "HETATM"):
            serial = int(ln[6:11]) if preserve else serial + 1
            if not preserve:
                ln = f"{ln[:6]}{serial:5d}{ln[11:]}"
            out.append(ln.ljust(80)[:80])
            chain, resseq, icode = ln[21], int(ln[22:26]), ln[26]
            if last_poly.get((mnum, chain)) == (resseq, icode):
                nxt = lines[i + 1] if i + 1 < len(lines) else ""
                same_res = nxt[:6] in ("ATOM  ", "HETATM") and nxt[21] == chain and \
                    int(nxt[22:26]) == resseq and nxt[26] == icode
                if not same_res:
                    serial += 1
                    out.append(f"TER   {serial:5d}      {ln[17:20]} {chain}{resseq:4d}{icode}".ljust(80))
            continue
        out.append(ln.ljust(80)[:80] if ln.strip() else ln)
    return "\n".join(out) + "\n"


def dumps(structure: Structure, fmt: str = "pdb", *, remap_chains: bool = False) -> str:
    """Serialise a structure to a PDB or mmCIF string."""
    fmt = normalize_format(fmt)
    if fmt == "cif":
        from .cifio import dump_cif
        return dump_cif(structure)
    return _dump_pdb(structure, remap_chains=remap_chains)


def output_format(dest: Optional[PathLike], fmt: Optional[str], structure: Optional[Structure] = None) -> str:
    """Pick the output format: explicit > file extension > input format > PDB."""
    if fmt:
        return normalize_format(fmt)
    if dest not in (None, "-"):
        f = format_from_path(dest)
        if f:
            return f
    if structure is not None:
        return source_info(structure).fmt
    return "pdb"


def save(structure: Structure, dest: Union[PathLike, io.TextIOBase, None] = "-",
         fmt: Optional[str] = None, *, remap_chains: bool = False) -> Optional[str]:
    """Write a structure to a file path, an open text handle, or stdout (``"-"``/None).

    The format is taken from ``fmt``, else the file extension, else the format
    the structure was read from. ``.gz`` paths are compressed. Returns the
    path written (or None for handles/stdout).
    """
    if hasattr(dest, "write"):
        dest.write(dumps(structure, output_format(None, fmt, structure), remap_chains=remap_chains))
        return None
    fmt = output_format(dest, fmt, structure)
    text = dumps(structure, fmt, remap_chains=remap_chains)
    if dest in (None, "-"):
        sys.stdout.write(text)
        sys.stdout.flush()
        return None
    path = os.fspath(dest)
    parent = os.path.dirname(path)
    if parent:
        os.makedirs(parent, exist_ok=True)
    if path.lower().endswith(".gz"):
        with gzip.open(path, "wt") as fh:
            fh.write(text)
    else:
        with open(path, "w") as fh:
            fh.write(text)
    return path


def deep_copy(entity):
    """``entity.copy()`` with alternate-location selections repaired.

    Biopython's ``DisorderedEntityWrapper.copy()`` copies every alternate
    location but leaves the wrapper's *selected* child pointing at the
    ORIGINAL object. Reading through the wrapper (``get_atoms()``, ``.coord``)
    then returns the original's coordinates after the copy is moved or
    edited. Every copy pdbKit makes goes through this function, which points
    each selection at the copied child with the same altloc / residue name.
    """
    new = entity.copy()
    if hasattr(entity, "xtra"):
        new.xtra = dict(entity.xtra)
    level = new.get_level()
    residues = [new] if level == "R" else list(iter_residues_all(new))
    for chain in ([new] if level == "C" else (new.get_chains() if level in ("S", "M") else [])):
        for res in chain:
            if res.is_disordered() == 2:
                sel = res.selected_child
                name = sel.get_resname() if sel is not None else None
                if name is not None and res.disordered_has_id(name):
                    res.disordered_select(name)
    for res in residues:
        for atom in res.child_list:
            if atom.is_disordered():
                sel = atom.selected_child
                label = sel.altloc if sel is not None else None
                if label is not None and atom.disordered_has_id(label):
                    atom.disordered_select(label)
    return new


def iter_residues_all(entity) -> Iterable:
    """Every residue, expanding microheterogeneity (``DisorderedResidue``) into all variants.

    Iterating a chain yields one wrapper per disordered residue, whose atoms
    are those of the *selected* variant only; this yields every variant.
    """
    level = entity.get_level()
    if level == "A":
        return
    if level == "R":
        yield from (entity.disordered_get_list() if entity.is_disordered() == 2 else [entity])
        return
    chains = [entity] if level == "C" else entity.get_chains()
    for chain in chains:
        yield from chain.get_unpacked_list()


def iter_atoms_all_altlocs(entity) -> Iterable:
    """Every atom, including all alternate locations and residue variants."""
    for residue in iter_residues_all(entity):
        yield from residue.get_unpacked_list()
