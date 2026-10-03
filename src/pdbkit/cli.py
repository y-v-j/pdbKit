"""``pdbkit`` command-line interface.

Every structure command reads PDB or mmCIF (``.gz`` too), a PDB/AlphaFold ID,
or ``-`` for stdin, and writes to ``-o FILE`` or stdout, so commands pipe::

    pdbkit fetch 1brs -o - | pdbkit extract - --chains A --atoms ca | pdbkit center - -o 1brs_A_ca_c.pdb

Run ``pdbkit <command> -h`` for the options of one command.
"""

from __future__ import annotations

import argparse
import csv
import json
import os
import shutil
import sys
from typing import List, Optional

from . import __version__
from . import analysis as an
from . import edit as ed
from .geometry import center, center_asu, center_of_mass, radius_of_gyration
from .select import Selection, describe, extract, split_chains, split_models, split_segments, write_parts
from .structure import PdbKitError, load, output_format, save, source_info, structure_id


# --------------------------------------------------------------------------- #
# helpers
# --------------------------------------------------------------------------- #

def _stem(structure, given: Optional[str] = None) -> str:
    if given:
        return given
    info = source_info(structure)
    if info.path:
        name = os.path.basename(info.path)
        for ext in (".gz", ".cif", ".mmcif", ".pdb", ".ent"):
            if name.lower().endswith(ext):
                name = name[: -len(ext)]
        return name
    return structure_id(structure)


def _emit_rows(rows, args, columns=None) -> None:
    rows = list(rows)
    out = open(args.output, "w", newline="") if getattr(args, "output", None) not in (None, "-") else sys.stdout
    fmt = getattr(args, "table", "tsv")
    try:
        if fmt == "json":
            json.dump(rows, out, indent=2, default=str)
            out.write("\n")
            return
        if not rows and not columns:
            return
        cols = columns or list(rows[0].keys())
        w = csv.DictWriter(out, fieldnames=cols, delimiter="," if fmt == "csv" else "\t",
                           extrasaction="ignore", lineterminator="\n")
        w.writeheader()
        for r in rows:
            w.writerow({k: ("" if v is None else v) for k, v in r.items()})
    finally:
        if out is not sys.stdout:
            out.close()


def _emit_obj(obj, args) -> None:
    out = open(args.output, "w") if getattr(args, "output", None) not in (None, "-") else sys.stdout
    try:
        if getattr(args, "table", "text") == "json":
            json.dump(obj, out, indent=2, default=str)
            out.write("\n")
        else:
            for k, v in obj.items():
                if isinstance(v, (list, dict)):
                    v = json.dumps(v, default=str)
                out.write(f"{k}\t{v}\n")
    finally:
        if out is not sys.stdout:
            out.close()


def _write_structure(s, args) -> None:
    dest = getattr(args, "output", None) or "-"
    path = save(s, dest, getattr(args, "format", None), remap_chains=getattr(args, "remap_chains", False))
    if path:
        print(f"wrote {path}", file=sys.stderr)


def _add_io(p, output=True, table=False) -> None:
    p.add_argument("input", help="structure file (.pdb/.cif[.gz]), PDB/AlphaFold ID, or - for stdin")
    if output:
        p.add_argument("-o", "--output", default="-", help="output file (default: stdout)")
    if not table:
        p.add_argument("-f", "--format", choices=["pdb", "cif"],
                       help="output format (default: from -o extension, else the input format)")
        p.add_argument("--remap-chains", action="store_true",
                       help="PDB output only: map multi-character chain IDs to single letters")


def _add_table(p) -> None:
    p.add_argument("--table", choices=["tsv", "csv", "json"], default="tsv", help="output table format")


def _add_selection(p, het_default="all", water_default=True) -> None:
    g = p.add_argument_group("selection")
    g.add_argument("-c", "--chains", help="keep only these chains, e.g. A or A,C")
    g.add_argument("--exclude-chains", help="drop these chains")
    g.add_argument("-a", "--atoms", default="all",
                   help="polymer atoms: all | ca | backbone | trace | comma list of names (default all)")
    g.add_argument("--exclude-atoms", help="drop atoms with these names")
    g.add_argument("--het", default=het_default,
                   help=f"ligands/heteroatoms: all | none | comma list of residue names (default {het_default})")
    wg = g.add_mutually_exclusive_group()
    wg.add_argument("--water", dest="water", action="store_true", default=water_default,
                    help="keep water" + (" (default)" if water_default else ""))
    wg.add_argument("--no-water", dest="water", action="store_false",
                    help="drop water" + ("" if water_default else " (default)"))
    g.add_argument("--no-hydrogens", action="store_true", help="drop H/D atoms")
    g.add_argument("--models", help="model numbers to keep, e.g. 1 or 1,3")
    g.add_argument("--residues", help="residue ranges: 1:10, 5:, :20, 1:100:5, 42 (comma separated)")
    g.add_argument("--exclude-residues", help="residue ranges to drop")
    g.add_argument("--resnames", help="keep only residues with these names")
    g.add_argument("--exclude-resnames", help="drop residues with these names")
    g.add_argument("--elements", help="keep only these elements")
    g.add_argument("--exclude-elements", help="drop these elements")
    g.add_argument("--segids", help="keep only these segment IDs (PDB)")
    g.add_argument("--altloc", help="resolve alternate locations: best (highest occupancy) or a label, e.g. A")
    hg = g.add_mutually_exclusive_group()
    hg.add_argument("--hetatm-only", action="store_true", help="keep only HETATM records")
    hg.add_argument("--no-hetatm", action="store_true", help="drop all HETATM records (pdb_delhetatm)")


def _selection(args) -> Selection:
    return Selection(
        chains=args.chains, exclude_chains=args.exclude_chains, atoms=args.atoms,
        exclude_atoms=args.exclude_atoms, het=args.het, water=args.water,
        hydrogens=not args.no_hydrogens, models=args.models, residues=args.residues,
        exclude_residues=args.exclude_residues, resnames=args.resnames,
        exclude_resnames=args.exclude_resnames, elements=args.elements,
        exclude_elements=args.exclude_elements, segids=args.segids, altloc=args.altloc,
        only_het_records=True if args.hetatm_only else (False if args.no_hetatm else None),
    )


def _selection_given(args) -> bool:
    return any(getattr(args, k, None) for k in (
        "chains", "exclude_chains", "exclude_atoms", "models", "residues", "exclude_residues",
        "resnames", "exclude_resnames", "elements", "exclude_elements", "segids", "altloc",
        "hetatm_only", "no_hetatm", "no_hydrogens")) or args.atoms != "all" or args.het != "all" or not args.water


# --------------------------------------------------------------------------- #
# commands
# --------------------------------------------------------------------------- #

def cmd_fetch(args):
    import tempfile
    from .fetch import fetch
    if args.output is not None and len(args.ids) > 1:
        raise PdbKitError("-o/--output takes a single ID; use -d/--outdir for several.")

    def get(ident, outdir):
        return fetch(ident, fmt=args.format, source=args.source, assembly=args.assembly,
                     outdir=outdir, overwrite=args.overwrite, pae=args.pae)

    if args.output == "-":                         # stream to stdout (for pipes)
        with tempfile.TemporaryDirectory() as tmp:
            with open(get(args.ids[0], tmp)) as fh:
                sys.stdout.write(fh.read())
        return
    if args.output is not None:                    # explicit file name
        with tempfile.TemporaryDirectory() as tmp:
            path = get(args.ids[0], tmp)
            os.makedirs(os.path.dirname(os.path.abspath(args.output)), exist_ok=True)
            shutil.move(path, args.output)          # works across filesystems (/tmp is often tmpfs)
        print(args.output)
        return
    for ident in args.ids:
        print(get(ident, args.outdir))


def cmd_convert(args):
    s = load(args.input)
    if args.format is None and args.output in (None, "-"):
        args.format = "cif" if source_info(s).fmt == "pdb" else "pdb"
    _write_structure(s, args)


def cmd_info(args):
    s = load(args.input)
    info = an.summary(s)
    info.update({f"exp_{k}": v for k, v in an.experiment(s).items() if k != "id"})
    info["input_format"] = source_info(s).fmt
    _emit_obj(info, args)


def cmd_validate(args):
    from .validate import validate_file, validate_text
    res = validate_text(sys.stdin.read()) if args.input == "-" else validate_file(args.input)
    for e in res["errors"]:
        print(f"ERROR   {e}")
    for w in res["warnings"]:
        print(f"WARNING {w}")
    print(f"{'OK' if res['ok'] else 'FAILED'}: {res['format']} file, "
          f"{res.get('models', '?')} model(s), {res.get('atoms', '?')} atoms")
    return 0 if res["ok"] else 1


def cmd_fasta(args):
    s = load(args.input)
    text = an.fasta(s, source=args.source, chains=args.chains, width=args.width)
    if args.output == "-":
        sys.stdout.write(text)
        return
    info = source_info(s)
    local = info.path and os.path.isfile(args.input)
    # default: next to a local input file (4MMK.cif -> 4MMK.fasta); cwd for a fetched ID
    folder = args.outdir or (os.path.dirname(info.path) if local else ".") or "."
    dest = args.output or os.path.join(folder, _stem(s) + ".fasta")
    os.makedirs(os.path.dirname(os.path.abspath(dest)), exist_ok=True)
    with open(dest, "w") as fh:
        fh.write(text)
    print(dest)


def cmd_extract(args):
    s = load(args.input)
    _write_structure(extract(s, _selection(args)), args)


def cmd_split(args):
    s = load(args.input)
    fmt = args.format or source_info(s).fmt
    stem = _stem(s, args.prefix)
    if args.by == "chain":
        sel = _selection(args)
        chains = sel.chains
        parts = split_chains(s, chains=chains, selection=sel)
        tag = describe(sel) if not args.no_tag else ""
    elif args.by == "model":
        parts, tag = {f"model{k}": v for k, v in split_models(s).items()}, ""
    else:
        parts, tag = {f"seg{k or 'blank'}": v for k, v in split_segments(s).items()}, ""
    for p in write_parts(parts, args.outdir, stem, fmt, tag, remap_chains=args.remap_chains):
        print(p)


def cmd_com(args):
    s = load(args.input)
    sel = _selection(args) if _selection_given(args) else None
    com = center_of_mass(s, sel, weighted=not args.geometric)
    rg = radius_of_gyration(s, sel, weighted=not args.geometric)
    obj = {"x": round(float(com[0]), 4), "y": round(float(com[1]), 4), "z": round(float(com[2]), 4),
           "weighting": "geometric" if args.geometric else "mass", "radius_of_gyration": round(rg, 4)}
    _emit_obj(obj, args)


def cmd_center(args):
    s = load(args.input)
    sel = _selection(args) if _selection_given(args) else None
    new, com = center(s, sel, weighted=not args.geometric, apply_to="all" if args.keep_all else "selection",
                      per_model=args.per_model)
    print(f"centre of {'mass' if not args.geometric else 'geometry'} before: "
          f"{com[0]:.4f} {com[1]:.4f} {com[2]:.4f}", file=sys.stderr)
    if args.output is None:
        args.output = f"{_stem(s)}_centered.{output_format(None, args.format, s)}"
    _write_structure(new, args)


def cmd_centerasu(args):
    s = load(args.input)
    new, shift = center_asu(s)
    print(f"fractional shift: {shift[0]:.0f} {shift[1]:.0f} {shift[2]:.0f}", file=sys.stderr)
    _write_structure(new, args)


def cmd_renumber(args):
    s = load(args.input)
    if args.atoms_from is not None:
        s = ed.renumber_atoms(s, args.atoms_from)
    if args.start is not None or args.shift is not None:
        s = ed.renumber_residues(s, start=args.start, shift=args.shift, chains=args.chains,
                                 include_het=not args.polymer_only, keep_gaps=args.keep_gaps,
                                 continuous=args.continuous)
    _write_structure(s, args)


def cmd_rename(args):
    s = load(args.input)
    if args.set_chain is not None:
        s = ed.set_chain(s, args.set_chain)
    if args.chain_map:
        s = ed.rename_chains(s, args.chain_map)
    if args.residue_map:
        s = ed.rename_residues(s, args.residue_map)
    if args.set_segid is not None:
        s = ed.set_segid(s, args.set_segid)
    if args.chain_to_segid:
        s = ed.chain_to_segid(s)
    if args.segid_to_chain:
        s = ed.segid_to_chain(s)
    _write_structure(s, args)


def cmd_setfield(args):
    s = load(args.input)
    if args.bfactor is not None:
        s = ed.set_bfactor(s, args.bfactor, args.chains)
    if args.occupancy is not None:
        s = ed.set_occupancy(s, args.occupancy, args.chains)
    if args.bfactor_file:
        s = ed.bfactor_from_data(s, args.bfactor_file, default=args.default,
                                 field="occupancy" if args.into_occupancy else "bfactor")
    _write_structure(s, args)


def cmd_sort(args):
    _write_structure(ed.sort_structure(load(args.input)), args)


def cmd_merge(args):
    s = ed.merge(args.inputs, rename_duplicates=args.rename_duplicates)
    _write_structure(s, args)


def cmd_ensemble(args):
    if args.check:
        problems = ed.check_ensemble(load(args.inputs[0]))
        for p in problems:
            print(p)
        print("ensemble OK" if not problems else f"{len(problems)} problem(s)")
        return 0 if not problems else 1
    _write_structure(ed.make_ensemble(args.inputs), args)


def cmd_intersect(args):
    _write_structure(ed.intersect(args.inputs), args)


def _rows_cmd(fn, columns=None):
    def run(args):
        s = load(args.input)
        _emit_rows(fn(s, args), args, columns)
    return run


def cmd_disulfide(args):
    s = load(args.input)
    res = an.disulfides(s, args.cutoff)
    rows = [dict(r, state="bonded") for r in res["bonds"]] + \
           [{"cys1": c, "cys2": "", "distance": "", "state": "free"} for c in res["free"]]
    _emit_rows(rows, args, ["cys1", "cys2", "distance", "state"])


def cmd_watercontacts(args):
    s = load(args.input)
    nb = an.water_contacts(s, args.probe)
    _emit_rows([{"residue": k, "n": len(v), "neighbors": ";".join(v)} for k, v in nb.items()], args)


def cmd_param(args):
    _emit_obj(an.param(load(args.input), args.chains, args.source), args)


def cmd_coulomb(args):
    e = an.coulomb(load(args.input), args.dielectric, args.ionic_strength, args.ph, args.temperature)
    _emit_obj({"coulomb_energy_kcal_mol": round(e, 4)}, args)


def cmd_exper(args):
    _emit_obj(an.experiment(load(args.input)), args)


def cmd_assembly(args):
    rows = an.assemblies(load(args.input))
    _emit_rows([dict(r, chains=",".join(r["chains"])) for r in rows], args)


def cmd_hetnavi(args):
    from . import hetnavi
    records, written = hetnavi.interactions(args.input, assembly=args.assembly, dat_file=args.file,
                                            outdir=args.outdir, save_raw=args.save_raw)
    print(f"Parsed {len(records)} chain record(s).")
    for p in written:
        print(p)


def cmd_domains(args):
    from . import domains
    written = domains.run(args.inputs, args.outdir, skip_lookups=args.no_lookups,
                          log=lambda m: print(m, file=sys.stderr))
    for p in written:
        print(p)


# --------------------------------------------------------------------------- #
# parser
# --------------------------------------------------------------------------- #

def build_parser() -> argparse.ArgumentParser:
    p = argparse.ArgumentParser(prog="pdbkit", description="pdbKit: tools for PDB and mmCIF structure files.",
                                epilog="Run 'pdbkit <command> -h' for help on a command.")
    p.add_argument("--version", action="version", version=f"pdbKit {__version__}")
    sub = p.add_subparsers(dest="command", metavar="<command>")
    sub.required = True

    def add(name, fn, help_, **kw):
        sp = sub.add_parser(name, help=help_, description=help_, **kw)
        sp.set_defaults(func=fn)
        return sp

    # ---- input / output ---------------------------------------------------- #
    sp = add("fetch", cmd_fetch, "download structures from RCSB PDB or AlphaFold DB")
    sp.add_argument("ids", nargs="+", help="PDB IDs (1ABC), UniProt accessions (P69905) or AF-<acc>-F1")
    sp.add_argument("-f", "--format", choices=["pdb", "cif"], default="cif", help="file format (default cif)")
    sp.add_argument("--source", choices=["auto", "pdb", "alphafold"], default="auto")
    sp.add_argument("--assembly", type=int, help="biological assembly number (RCSB)")
    sp.add_argument("-d", "--outdir", default=".", help="download directory (default .)")
    sp.add_argument("-o", "--output", default=None, help="file name for a single download; - for stdout")
    sp.add_argument("--pae", action="store_true", help="AlphaFold: also save the PAE JSON")
    sp.add_argument("--overwrite", action="store_true")

    sp = add("convert", cmd_convert, "convert PDB <-> mmCIF (also tidies: writes a clean coordinate file)")
    _add_io(sp)

    sp = add("info", cmd_info, "summary: counts, chains, ligands, method, resolution (pdb_wc + exper)")
    _add_io(sp, table=True)
    sp.add_argument("--table", choices=["text", "json"], default="text")

    sp = add("validate", cmd_validate, "check file format (pdb_validate for PDB; mmCIF consistency)")
    sp.add_argument("input")

    sp = add("fasta", cmd_fasta, "write FASTA for all protein/nucleic-acid chains (4MMK.cif -> 4MMK.fasta)")
    sp.add_argument("input")
    sp.add_argument("-o", "--output", default=None, help="output file (default <stem>.fasta; - for stdout)")
    sp.add_argument("-d", "--outdir", default=None, help="directory for the default-named file")
    sp.add_argument("-c", "--chains", help="only these chains")
    sp.add_argument("--source", choices=["auto", "seqres", "atom"], default="auto",
                    help="seqres = deposited sequence, atom = modelled residues, auto = seqres if present")
    sp.add_argument("--width", type=int, default=60, help="line width (default 60)")

    # ---- selection / splitting --------------------------------------------- #
    sp = add("extract", cmd_extract, "select/delete atoms, residues, chains, ligands, water (pdb_sel*/del*)")
    _add_io(sp)
    _add_selection(sp)

    sp = add("split", cmd_split, "split into one file per chain (with representation options), model or segment")
    sp.add_argument("input")
    sp.add_argument("--by", choices=["chain", "model", "segment"], default="chain")
    sp.add_argument("-d", "--outdir", default=".", help="output directory (default .)")
    sp.add_argument("-p", "--prefix", help="file-name prefix (default: input file stem)")
    sp.add_argument("-f", "--format", choices=["pdb", "cif"], help="output format (default: input format)")
    sp.add_argument("--no-tag", action="store_true", help="do not append the representation tag to file names")
    sp.add_argument("--remap-chains", action="store_true")
    _add_selection(sp, het_default="none", water_default=False)

    # ---- geometry ------------------------------------------------------------ #
    sp = add("com", cmd_com, "report centre of mass (whole structure or selection) and radius of gyration")
    _add_io(sp, table=True)
    sp.add_argument("--table", choices=["text", "json"], default="text")
    sp.add_argument("--geometric", action="store_true", help="unweighted geometric centre")
    _add_selection(sp)

    sp = add("center", cmd_center, "move the centre of mass to the origin and write a new file")
    _add_io(sp)
    sp.set_defaults(output=None)   # default: <stem>_centered.<ext>; '-o -' for stdout
    sp.add_argument("--geometric", action="store_true", help="use the geometric centre")
    sp.add_argument("--keep-all", action="store_true",
                    help="write the whole structure, shifted by the selection's COM (default: selection only)")
    sp.add_argument("--per-model", action="store_true", help="centre each model separately")
    _add_selection(sp)

    sp = add("centerasu", cmd_centerasu, "move the asymmetric unit into the unit cell (Harms centerasu)")
    _add_io(sp)

    # ---- editing ------------------------------------------------------------ #
    sp = add("renumber", cmd_renumber, "renumber residues (start/shift) and/or atoms")
    _add_io(sp)
    sp.add_argument("--start", type=int, help="number residues from N")
    sp.add_argument("--shift", type=int, help="add K to residue numbers")
    sp.add_argument("--keep-gaps", action="store_true", help="with --start: offset, preserving gaps")
    sp.add_argument("--continuous", action="store_true", help="with --start: do not restart for each chain (pdb_reres)")
    sp.add_argument("-c", "--chains", help="only these chains")
    sp.add_argument("--polymer-only", action="store_true", help="leave ligands and water numbers unchanged")
    sp.add_argument("--atoms-from", type=int, help="renumber atom serials from N")

    sp = add("rename", cmd_rename, "rename chains / residues / segment IDs")
    _add_io(sp)
    sp.add_argument("--chain-map", help="A:H,B:L (swaps allowed)")
    sp.add_argument("--set-chain", help="set every chain ID to this value")
    sp.add_argument("--residue-map", help="HSD:HIS,HIP:HIS")
    sp.add_argument("--set-segid", help="set segment ID")
    sp.add_argument("--chain-to-segid", action="store_true")
    sp.add_argument("--segid-to-chain", action="store_true")

    sp = add("setfield", cmd_setfield, "set B-factor / occupancy, or load per-residue values into the B column")
    _add_io(sp)
    sp.add_argument("--bfactor", type=float)
    sp.add_argument("--occupancy", type=float)
    sp.add_argument("--bfactor-file", help="two-column file: residue (42 or A:42) and value")
    sp.add_argument("--into-occupancy", action="store_true", help="put --bfactor-file values in the occupancy column")
    sp.add_argument("--default", type=float, default=0.0, help="value for residues not in the file")
    sp.add_argument("-c", "--chains", help="only these chains (constant values)")

    sp = add("sort", cmd_sort, "sort chains and residues (pdb_sort)")
    _add_io(sp)

    for name, fn, help_ in (("merge", cmd_merge, "combine chains of several files into one (pdb_merge)"),
                            ("ensemble", cmd_ensemble, "stack files as models of an ensemble, or --check one"),
                            ("intersect", cmd_intersect, "keep atoms common to all files (pdb_intersect)")):
        sp = add(name, fn, help_)
        sp.add_argument("inputs", nargs="+")
        sp.add_argument("-o", "--output", default="-")
        sp.add_argument("-f", "--format", choices=["pdb", "cif"])
        sp.add_argument("--remap-chains", action="store_true")
        if name == "merge":
            sp.add_argument("--rename-duplicates", action="store_true")
        if name == "ensemble":
            sp.add_argument("--check", action="store_true", help="check a multi-model file (pdb_chkensemble)")

    # ---- analysis ----------------------------------------------------------- #
    sp = add("gaps", _rows_cmd(lambda s, a: an.gaps(s, a.cutoff),
                               ["model", "chain", "from", "to", "type", "distance"]), "chain breaks (pdb_gap)")
    _add_io(sp, table=True); _add_table(sp)
    sp.add_argument("--cutoff", type=float, default=4.0, help="CA-CA distance (A) counted as a break")

    sp = add("ligands", _rows_cmd(lambda s, a: an.ligands(s, a.skip_boring),
                                  ["chain", "resname", "resseq", "icode", "atoms"]), "list ligands (Harms ligand)")
    _add_io(sp, table=True); _add_table(sp)
    sp.add_argument("--skip-boring", action="store_true", help="omit HOH, CA, SO4, IOD, NA, CL, GOL, PO4")

    sp = add("disulfide", cmd_disulfide, "disulfide bonds and free cysteines")
    _add_io(sp, table=True); _add_table(sp)
    sp.add_argument("--cutoff", type=float, default=3.0, help="SG-SG cutoff (A, default 3.0); Harms used 4.0")

    sp = add("contacts", _rows_cmd(lambda s, a: an.contacts(s, a.target, a.cutoff, a.water)),
             "residues within a cutoff of a target residue (Harms contact)")
    _add_io(sp, table=True); _add_table(sp)
    sp.add_argument("-t", "--target", default="TRP", help="residue name (HEM), number (42) or chain:number (A:42)")
    sp.add_argument("--cutoff", type=float, default=3.5)
    sp.add_argument("--water", action="store_true", help="include water as contacts")

    sp = add("closecontacts", _rows_cmd(lambda s, a: an.close_contacts(s, a.cutoff, a.ca_only)),
             "atom pairs closer than a cutoff (Harms closecontacts)")
    _add_io(sp, table=True); _add_table(sp)
    sp.add_argument("--cutoff", type=float, default=3.5)
    sp.add_argument("--ca-only", action="store_true")

    sp = add("watercontacts", cmd_watercontacts, "residue neighbours by vdW + water probe (Harms watercontact)")
    _add_io(sp, table=True); _add_table(sp)
    sp.add_argument("--probe", type=float, default=1.4)

    sp = add("neighbors", _rows_cmd(lambda s, a: an.neighbors(s, a.cutoff, a.seq_cutoff)),
             "CB neighbour counts per residue (Harms neighbors)")
    _add_io(sp, table=True); _add_table(sp)
    sp.add_argument("--cutoff", type=float, default=10.0)
    sp.add_argument("--seq-cutoff", type=int, default=4)

    sp = add("torsion", _rows_cmd(lambda s, a: an.torsions(s)), "backbone phi/psi (Harms torsion)")
    _add_io(sp, table=True); _add_table(sp)

    sp = add("param", cmd_param, "MW, pI, composition, charge (Harms param)")
    _add_io(sp, table=True)
    sp.add_argument("--table", choices=["text", "json"], default="text")
    sp.add_argument("-c", "--chains")
    sp.add_argument("--source", choices=["atom", "seqres"], default="atom")

    sp = add("coulomb", cmd_coulomb, "screened Coulomb energy of titratable groups (Harms coulomb)")
    _add_io(sp, table=True)
    sp.add_argument("--table", choices=["text", "json"], default="text")
    sp.add_argument("--dielectric", type=float, default=40.0)
    sp.add_argument("--ionic-strength", type=float, default=100.0, help="mM")
    sp.add_argument("--ph", type=float, default=7.0)
    sp.add_argument("--temperature", type=float, default=298.0)

    sp = add("iondist", _rows_cmd(lambda s, a: an.ion_distance_histogram(s, a.step)),
             "distance histogram of titratable groups (Harms iondist)")
    _add_io(sp, table=True); _add_table(sp)
    sp.add_argument("--step", type=float, default=0.1)

    sp = add("distances", _rows_cmd(lambda s, a: an.distances(s, a.reference, a.atom, a.chains)),
             "distances from one atom to all atoms of a name (Harms dist_filter)")
    _add_io(sp, table=True); _add_table(sp)
    sp.add_argument("-r", "--reference", required=True, help="CHAIN:RESNUM:ATOM, e.g. A:42:NZ")
    sp.add_argument("--atom", help="target atom name (default: the reference atom's name)")
    sp.add_argument("-c", "--chains")

    sp = add("sasa", _rows_cmd(lambda s, a: an.sasa(s, a.level, a.probe, a.points, a.het)),
             "solvent accessible surface area, Shrake-Rupley (replaces NACCESS-based Harms sasa)")
    _add_io(sp, table=True); _add_table(sp)
    sp.add_argument("--level", choices=["residue", "atom"], default="residue")
    sp.add_argument("--probe", type=float, default=1.4)
    sp.add_argument("--points", type=int, default=100)
    sp.add_argument("--het", action="store_true", help="include ligands in the surface")

    sp = add("exper", cmd_exper, "experiment summary: method, resolution, R, organism (Harms exper)")
    _add_io(sp, table=True)
    sp.add_argument("--table", choices=["text", "json"], default="text")

    sp = add("assembly", cmd_assembly, "biological assemblies; is the file the assembly? (Harms oligomer)")
    _add_io(sp, table=True); _add_table(sp)

    # ---- integrations ------------------------------------------------------- #
    sp = add("hetnavi", cmd_hetnavi, "Het-PDB Navi protein-ligand / subunit interactions -> CSV")
    sp.add_argument("input", help="PDB ID or a structure file (the ID is read from its header)")
    sp.add_argument("--file", help="already-downloaded Het-PDB Navi .dat file")
    sp.add_argument("--assembly", default="assembly1")
    sp.add_argument("-d", "--outdir", default=".")
    sp.add_argument("--save-raw", action="store_true", help="also keep the downloaded .dat file")

    sp = add("domains", cmd_domains, "domain / family / EC / secondary structure / bonding table per chain")
    sp.add_argument("inputs", nargs="+", help="PDB IDs, structure files (.cif/.pdb) or folders")
    sp.add_argument("-d", "--outdir", default=".")
    sp.add_argument("-q", "--no-lookups", action="store_true", help="skip live UniProt and KEGG lookups")
    return p


def main(argv: Optional[List[str]] = None) -> int:
    parser = build_parser()
    args = parser.parse_args(argv)
    try:
        rc = args.func(args)
        return int(rc or 0)
    except PdbKitError as exc:
        print(f"pdbkit: error: {exc}", file=sys.stderr)
        return 2
    except BrokenPipeError:
        return 0
    except KeyboardInterrupt:
        return 130


if __name__ == "__main__":
    sys.exit(main())
