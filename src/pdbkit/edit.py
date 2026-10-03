"""Structure editing: renumbering, renaming, field setting, sorting and combining.

Every function returns a new structure (the input is not modified) so the
functions chain naturally in scripts and notebooks. Equivalent tools:

==================================  =============================================
pdbKit                              pdb-tools / Harms pdbtools
==================================  =============================================
``renumber_residues(start=...)``    pdb_reres, residue_renumber (Harms)
``renumber_residues(shift=...)``    pdb_shiftres, offset (Harms)
``renumber_atoms``                  pdb_reatom, atom_renumber (Harms)
``rename_chains`` / ``set_chain``   pdb_rplchain / pdb_chain
``rename_residues``                 pdb_rplresname
``set_bfactor`` / ``set_occupancy`` pdb_b / pdb_occ
``bfactor_from_data``               bfactor (Harms)
``set_segid`` / chain<->segid       pdb_seg, pdb_chainxseg, pdb_segxchain
``sort_structure``                  pdb_sort
``merge``                           pdb_merge
``make_ensemble`` / ``check_ens.``  pdb_mkensemble / pdb_chkensemble
``intersect``                       pdb_intersect
==================================  =============================================
"""

from __future__ import annotations

import os
from typing import Dict, Iterable, List, Mapping, Optional, Sequence, Tuple, Union

from Bio.PDB.Chain import Chain
from Bio.PDB.Model import Model
from Bio.PDB.Structure import Structure

from .select import _as_list
from .structure import (PdbKitError, SourceInfo, deep_copy, iter_atoms_all_altlocs, iter_residues_all, load,
                        residue_kind, source_info)

__all__ = ["renumber_residues", "renumber_atoms", "rename_chains", "set_chain", "rename_residues",
           "set_bfactor", "set_occupancy", "bfactor_from_data", "load_data_file", "set_segid",
           "chain_to_segid", "segid_to_chain", "sort_structure", "merge", "make_ensemble",
           "check_ensemble", "intersect"]


def _copy(structure: Structure) -> Structure:
    return deep_copy(structure)


def _set_residue_ids(chain: Chain, new_ids: List[Tuple[str, int, str]]) -> None:
    """Re-key all residues of a chain at once (avoids transient ID collisions).

    Microheterogeneity wrappers (``DisorderedResidue``) are re-keyed together
    with every variant they hold; each variant keeps its own hetero flag.
    """
    residues = list(chain)
    if len(set(new_ids)) != len(new_ids):
        raise PdbKitError(f"Renumbering would create duplicate residue IDs in chain {chain.id!r}.")
    chain.child_dict = {}
    for res, rid in zip(residues, new_ids):
        if res.is_disordered() == 2:
            for v in res.disordered_get_list():
                v._id = (v.id[0], rid[1], rid[2])
                v._reset_full_id()
            res.id = rid
        else:
            res._id = rid
            res._reset_full_id()
        chain.child_dict[rid] = res


def renumber_residues(structure: Structure, *, start: Optional[int] = None, shift: Optional[int] = None,
                      chains: Optional[Sequence[str]] = None, include_het: bool = True,
                      keep_gaps: bool = False, continuous: bool = False) -> Structure:
    """Renumber residues.

    * ``start=N`` - number residues N, N+1, ... in file order (per chain;
      ``continuous=True`` carries on across chains like pdb_reres).
      ``keep_gaps=True`` instead offsets each chain so its first residue is N
      and gaps are preserved (Harms ``pdb_residue_renumber`` default).
    * ``shift=K`` - add K to every residue number (pdb_shiftres / pdb_offset).

    ``include_het=False`` leaves ligands and waters untouched. Numbering
    restarts for every model. Insertion codes are cleared when residues are
    renumbered sequentially and kept when shifting.
    """
    if (start is None) == (shift is None):
        raise PdbKitError("Give exactly one of start= or shift=.")
    s = _copy(structure)
    wanted = _as_list(chains)
    for model in s:
        counter = start
        for chain in model:
            if wanted is not None and chain.id not in wanted:
                continue
            residues = list(chain)
            new_ids = []
            first_num = None
            if start is not None and not continuous:
                counter = start
            for res in residues:
                het, num, icode = res.id
                if not include_het and residue_kind(res) != "polymer":
                    new_ids.append(res.id)
                    continue
                if shift is not None:
                    new_ids.append((het, num + shift, icode))
                elif keep_gaps:
                    if first_num is None:
                        first_num = num
                    new_ids.append((het, start + (num - first_num), icode))
                else:
                    new_ids.append((het, counter, " "))
                    counter += 1
            _set_residue_ids(chain, new_ids)
    return s


def renumber_atoms(structure: Structure, start: int = 1) -> Structure:
    """Set atom serial numbers to start, start+1, ... (kept by the writers)."""
    s = _copy(structure)
    n = start
    for model in s:
        for atom in iter_atoms_all_altlocs(model):
            atom.serial_number = n
            n += 1
    s.xtra["pdbkit_preserve_serials"] = True
    return s


def _rechain(model: Model, chain: Chain, new_id: str, old_id: Optional[str] = None) -> None:
    """Rename a chain; merge into an existing chain of that name if there is one."""
    old_id = old_id or chain.id
    if new_id == chain.id:
        return
    if new_id in model:
        target = model[new_id]
        for res in list(chain):
            if res.id in target:
                raise PdbKitError(f"Cannot rename chain {old_id!r} to {new_id!r}: chain {new_id!r} already exists "
                                  f"and residue {res.id[1]}{res.id[2].strip()} is in both. Rename or renumber one of them first.")
            chain.detach_child(res.id)
            target.add(res)
        model.detach_child(chain.id)
    else:
        chain.id = new_id


def rename_chains(structure: Structure, mapping: Union[Mapping[str, str], str]) -> Structure:
    """Rename chains, e.g. ``{"A": "H", "B": "L"}`` or ``"A:H,B:L"``.

    All renames happen simultaneously, so swaps (``A:B,B:A``) work.
    """
    if isinstance(mapping, str):
        mapping = dict(item.split(":", 1) for item in mapping.split(",") if item.strip())
    s = _copy(structure)
    for model in s:
        chains = [(c, c.id, mapping[c.id]) for c in list(model) if c.id in mapping]
        for i, (c, _, _) in enumerate(chains):
            c.id = f"__pdbkit_tmp_{i}"
        for c, old, new in chains:
            _rechain(model, c, new, old)
    return s


def set_chain(structure: Structure, chain_id: str) -> Structure:
    """Give every residue the same chain ID (pdb_chain)."""
    s = _copy(structure)
    for model in s:
        chains = [(c, c.id) for c in model]
        for i, (c, _) in enumerate(chains):
            c.id = f"__pdbkit_tmp_{i}"
        for c, old in chains:
            _rechain(model, c, chain_id, old)
    return s


def rename_residues(structure: Structure, mapping: Union[Mapping[str, str], str]) -> Structure:
    """Rename residues by name, e.g. ``{"HSD": "HIS"}`` or ``"HSD:HIS,HIP:HIS"``."""
    if isinstance(mapping, str):
        mapping = dict(item.split(":", 1) for item in mapping.split(",") if item.strip())
    s = _copy(structure)
    for chain in s.get_chains():
        residues = list(chain)
        new_ids = []
        for res in residues:
            name = res.get_resname().strip()
            het, num, icode = res.id
            if name in mapping:
                res.resname = mapping[name]
                if het.startswith("H_"):
                    het = "H_" + mapping[name]
            new_ids.append((het, num, icode))
        _set_residue_ids(chain, new_ids)
    return s


def _atoms(structure: Structure, chains=None):
    wanted = _as_list(chains)
    for chain in structure.get_chains():
        if wanted is None or chain.id in wanted:
            for res in chain.get_unpacked_list():
                yield from res.get_unpacked_list()


def set_bfactor(structure: Structure, value: float = 10.0, chains=None) -> Structure:
    s = _copy(structure)
    for a in _atoms(s, chains):
        a.bfactor = float(value)
    return s


def set_occupancy(structure: Structure, value: float = 1.0, chains=None) -> Structure:
    s = _copy(structure)
    for a in _atoms(s, chains):
        a.occupancy = float(value)
    return s


def load_data_file(path: str, key_col: int = 0, value_col: int = 1, abs_value: bool = False) -> Dict[str, float]:
    """Two-column (or more) whitespace data file -> {key: value}. '#' lines are comments.

    Keys may be residue numbers (``42``) or chain-qualified (``A:42``).
    """
    data: Dict[str, float] = {}
    with open(path) as fh:
        for line in fh:
            if not line.strip() or line.lstrip().startswith("#"):
                continue
            fields = line.split()
            try:
                v = float(fields[value_col])
            except (IndexError, ValueError):
                continue
            data[fields[key_col]] = abs(v) if abs_value else v
    return data


def bfactor_from_data(structure: Structure, data: Union[Mapping[str, float], str], *,
                      default: float = 0.0, field: str = "bfactor") -> Structure:
    """Write per-residue values into the B-factor (or occupancy) column.

    ``data`` maps residue numbers (``"42"``) or chain-qualified numbers
    (``"A:42"``) to values; a file path is read with :func:`load_data_file`.
    Residues not in ``data`` - and all ligands/waters - get ``default``.
    (Harms ``pdb_bfactor`` wrote the text ``NA`` there, which is not a valid
    number in the B-factor column.)
    """
    if isinstance(data, (str, os.PathLike)):
        data = load_data_file(os.fspath(data))
    s = _copy(structure)
    for chain in s.get_chains():
        for res in chain.get_unpacked_list():
            num = str(res.id[1]) + (res.id[2].strip() or "")
            value = default
            if residue_kind(res) == "polymer":
                value = data.get(f"{chain.id}:{num}", data.get(num, default))
            for a in res.get_unpacked_list():
                setattr(a, "bfactor" if field == "bfactor" else "occupancy", float(value))
    return s


def set_segid(structure: Structure, segid: str = "") -> Structure:
    s = _copy(structure)
    for res in iter_residues_all(s):
        res.segid = segid[:4]
    return s


def chain_to_segid(structure: Structure) -> Structure:
    s = _copy(structure)
    for chain in s.get_chains():
        for res in chain.get_unpacked_list():
            res.segid = chain.id[:4]
    return s


def segid_to_chain(structure: Structure) -> Structure:
    """Chain ID := first character of the segment ID (pdb_segxchain)."""
    s = _copy(structure)
    for model in s:
        moves = []
        for chain in list(model):
            for res in list(chain):
                seg = (res.segid.strip() or chain.id)[:1]
                if seg != chain.id:
                    moves.append((chain, res, seg))
        for chain, res, seg in moves:
            chain.detach_child(res.id)
            if seg not in model:
                model.add(Chain(seg))
            if res.id in model[seg]:
                raise PdbKitError(f"Residue {res.id} already exists in chain {seg!r}.")
            model[seg].add(res)
        for chain in list(model):
            if len(chain) == 0:
                model.detach_child(chain.id)
    return s


def sort_structure(structure: Structure) -> Structure:
    """Sort chains by ID; residues by number/insertion code, polymer first, then ligands, then water."""
    s = _copy(structure)
    order = {"polymer": 0, "ligand": 1, "water": 2}
    for model in s:
        model.child_list.sort(key=lambda c: c.id)
        for chain in model:
            chain.child_list.sort(key=lambda r: (order[residue_kind(r)], r.id[1], r.id[2]))
    return s


def _load_all(items) -> List[Structure]:
    return [load(x) for x in items]


def merge(structures: Iterable, *, rename_duplicates: bool = False) -> Structure:
    """Combine the chains of several structures into one (first model of each).

    Duplicate chain IDs are an error unless ``rename_duplicates=True``, in
    which case the later chain gets the next free letter.
    """
    structs = _load_all(structures)
    if not structs:
        raise PdbKitError("Nothing to merge.")
    out = Structure("merged")
    model = Model(0, serial_num=1)
    out.add(model)
    free = [c for c in "ABCDEFGHIJKLMNOPQRSTUVWXYZabcdefghijklmnopqrstuvwxyz0123456789"]
    for st in structs:
        for chain in next(iter(st)):
            c = deep_copy(chain)
            if c.id in model:
                if not rename_duplicates:
                    raise PdbKitError(f"Chain {c.id!r} occurs in more than one input; "
                                      "use --rename-duplicates or rename chains first.")
                used = {x.id for x in model}
                c.id = next(x for x in free if x not in used)
            model.add(c)
    out.xtra["pdbkit_source"] = SourceInfo(fmt=source_info(structs[0]).fmt)
    return out


def make_ensemble(structures: Iterable) -> Structure:
    """Stack structures as MODEL 1..N of one ensemble (pdb_mkensemble)."""
    structs = _load_all(structures)
    out = Structure("ensemble")
    for i, st in enumerate(structs):
        m = deep_copy(next(iter(st)))
        m.id = i
        m.serial_num = i + 1
        out.add(m)
    out.xtra["pdbkit_source"] = SourceInfo(fmt=source_info(structs[0]).fmt if structs else "pdb")
    return out


def _atom_keys(model) -> List[Tuple]:
    keys = []
    for chain in model:
        for res in chain:
            for a in res.get_unpacked_list():
                keys.append((chain.id, res.id[1], res.id[2], res.get_resname(), a.get_id(), a.altloc))
    return keys


def check_ensemble(structure: Structure) -> List[str]:
    """Problems that make a multi-model file an invalid ensemble (empty list = OK)."""
    models = list(structure)
    if len(models) < 2:
        return []
    ref = _atom_keys(models[0])
    ref_set = set(ref)
    problems = []
    for m in models[1:]:
        keys = _atom_keys(m)
        ks = set(keys)
        num = getattr(m, "serial_num", m.id + 1)
        missing, extra = ref_set - ks, ks - ref_set
        if missing:
            problems.append(f"model {num}: {len(missing)} atom(s) of model 1 missing, e.g. {sorted(missing)[0]}")
        if extra:
            problems.append(f"model {num}: {len(extra)} atom(s) not in model 1, e.g. {sorted(extra)[0]}")
        if not missing and not extra and keys != ref:
            problems.append(f"model {num}: same atoms as model 1 but in a different order")
    return problems


def intersect(structures: Iterable) -> Structure:
    """Atoms common to all inputs (name, altloc, residue name/number/icode, chain); coordinates of the first."""
    from .select import extract
    structs = _load_all(structures)
    common = set(_atom_keys(next(iter(structs[0]))))
    for st in structs[1:]:
        common &= set(_atom_keys(next(iter(st))))
    s = _copy(structs[0])
    for model in list(s)[1:]:
        s.detach_child(model.id)
    model = next(iter(s))
    for chain in list(model):
        for res in list(chain):
            for a in list(res.get_unpacked_list()):
                if (chain.id, res.id[1], res.id[2], res.get_resname(), a.get_id(), a.altloc) not in common:
                    if a.get_id() in res:
                        res.detach_child(a.get_id())
            if len(res) == 0:
                chain.detach_child(res.id)
        if len(chain) == 0:
            model.detach_child(chain.id)
    if len(model) == 0:
        raise PdbKitError("The inputs have no atoms in common.")
    return s
