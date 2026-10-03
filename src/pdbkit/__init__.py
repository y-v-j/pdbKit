"""pdbKit - one toolkit for PDB and mmCIF macromolecular structure files.

Merges HADDOCK ``pdb-tools`` (2.0.0-rc1), Harms-lab ``pdbtools``, the
Het-PDB Navi parser and the mmCIF domain/family annotator into one package
that treats PDB and mmCIF input identically.

Quick start (Python / IPython / Jupyter)::

    import pdbkit as pk
    s = pk.load("4HHB")                         # fetch + parse (mmCIF)
    parts = pk.split_chains(s, atoms="ca", het="all")
    pk.save(parts["A"], "4HHB_A_ca.pdb")
    centred, com = pk.center(s, chains="A", apply_to="all")
    print(pk.fasta(s))

Command line: ``pdbkit --help``.
"""

from __future__ import annotations

__version__ = "1.0.0"

from .structure import (PdbKitError, SourceInfo, load, loads, save, dumps, detect_format,
                        residue_kind, source_info, structure_id)
from .fetch import fetch, fetch_pdb, fetch_alphafold
from .select import Selection, extract, split_chains, split_models, split_segments, chain_ids, write_parts
from .geometry import center_of_mass, center, translate, radius_of_gyration, center_asu
from .edit import (renumber_residues, renumber_atoms, rename_chains, set_chain, rename_residues,
                   set_bfactor, set_occupancy, bfactor_from_data, set_segid, chain_to_segid,
                   segid_to_chain, sort_structure, merge, make_ensemble, check_ensemble, intersect)
from .analysis import (summary, sequences, to_fasta, fasta, gaps, ligands, disulfides, contacts,
                       close_contacts, water_contacts, neighbors, torsions, param, coulomb,
                       ion_distance_histogram, distances, sasa, experiment, assemblies, to_dataframe)
from .validate import validate_file, validate_text
from . import hetnavi, domains, metadata

__all__ = [n for n in dir() if not n.startswith("_")]


def write_fasta(source, path=None, **kwargs) -> str:
    """Write FASTA for all polymer chains; default path ``<stem>.fasta`` next to the input."""
    import os
    s = load(source)
    text = fasta(s, **kwargs)
    if path is None:
        info = source_info(s)
        base = info.path or structure_id(s)
        stem = os.path.basename(base)
        for ext in (".gz", ".cif", ".mmcif", ".pdb", ".ent"):
            if stem.lower().endswith(ext):
                stem = stem[: -len(ext)]
        path = os.path.join(os.path.dirname(info.path) if info.path else ".", stem + ".fasta")
    with open(path, "w") as fh:
        fh.write(text)
    return path
