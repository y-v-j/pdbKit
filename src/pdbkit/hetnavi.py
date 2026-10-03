"""Het-PDB Navi interaction data -> CSV tables (integrated ``hetpdbnavi_parser.py``).

Het-PDB Navi (https://hetpdbnavi.nagahama-i-bio.ac.jp/) publishes, for every
PDB entry, a flat file of protein-heteroatom and protein-protein contacts.
This module downloads (or reads) that file and writes:

1. ``<PDBID>_<CODE>_<RESNUM>_<CHAIN>_interaction.csv`` - one per heteroatom
   instance (the site's "Select hetatm" dropdown);
2. ``<PDBID>_all_het_interactions.csv`` - all protein-heteroatom contacts;
3. ``<PDBID>_subunit_interactions.csv`` - chain-chain interface contacts;
4. ``<PDBID>_disorder.csv`` - unresolved residue ranges.

The CSV layout is unchanged from the original script, and for data files
that include hydrogen-bond counts the outputs are byte-identical to it.
Fixed here: Het-PDB Navi also writes ``#HET`` lines *without* the trailing
hydrogen-bond count (most entries); the original parser required that field
and silently dropped every such contact. They are now kept, with an empty
``num_hbonds``. Column names are the original ones; note that ``res_start``
is the residue's position in the sequence and ``res_end`` its author
(PDB) residue number. Additions: :func:`interactions`, which
accepts a PDB ID *or* a structure file/Structure (the ID is read from the
file's header), and the option to keep the raw ``.dat`` file.
"""

from __future__ import annotations

import csv
import os
import re
import sys
import urllib.request
from collections import defaultdict

DOWNLOAD_URL_TEMPLATE = "https://hetpdbnavi.nagahama-i-bio.ac.jp/gethetpdbdata/{pdbid}/{assembly}"

HET_FIELDS = [
    "res_start", "res_end", "res_name", "res_chain",
    "het_index", "het_resnum", "het_code", "het_chain",
    "contact_type", "num_contacts", "num_hbonds",
]

SUB_FIELDS = [
    "res1_start", "res1_end", "res1_name", "res1_chain",
    "res2_start", "res2_end", "res2_name", "res2_chain",
    "interaction_type", "contact_type", "num_contacts",
]


def fetch_data(pdbid: str, assembly: str = "assembly1", timeout: int = 60) -> str:
    """Download the raw .dat interaction file for a PDB ID from Het-PDB Navi.

    Uses pdbKit's downloader: retries transient network errors (timeouts,
    resets) and raises a clear PdbKitError instead of a traceback.
    """
    import urllib.error
    from .fetch import _http_get
    from .structure import PdbKitError
    url = DOWNLOAD_URL_TEMPLATE.format(pdbid=pdbid.upper(), assembly=assembly)
    try:
        return _http_get(url, timeout=timeout, retries=3).decode("utf-8", errors="replace")
    except urllib.error.HTTPError as exc:
        raise PdbKitError(f"Het-PDB Navi returned HTTP {exc.code} for {pdbid.upper()} {assembly}.") from exc


def split_fields(line: str):
    """Split a tab-delimited FT line and strip padding whitespace from each field."""
    return [f.strip() for f in line.rstrip("\n").split("\t")]


def parse_records(text: str):
    """
    Parse the flat-file text into a list of records.
    Each record corresponds to one protein chain and is a dict with:
        id, pdbid, chain, het_contacts (list of dict), sub_contacts (list of dict),
        disorder (list of (start, end))
    """
    records = []
    current = None

    for raw_line in text.splitlines():
        if not raw_line.strip():
            continue
        fields = split_fields(raw_line)
        tag = fields[0]

        if tag == "ID":
            # e.g. ID  9VMHA  STANDARD;  PRT;  1127  AA.
            entry_id = fields[1]
            pdbid, chain = entry_id[:4], entry_id[4:]
            current = {
                "id": entry_id,
                "pdbid": pdbid,
                "chain": chain,
                "het_contacts": [],
                "sub_contacts": [],
                "disorder": [],
            }
            records.append(current)

        elif tag == "FT" and current is not None:
            subtag = fields[1]
            data = fields[2:]

            # Het-PDB Navi writes #HET lines with or without the trailing hydrogen-bond
            # count; the original parser required it and silently dropped the others
            if subtag == "#HET" and len(data) >= len(HET_FIELDS) - 1:
                row = dict(zip(HET_FIELDS, data[: len(HET_FIELDS)]))
                row.setdefault("num_hbonds", "")
                row["protein_chain"] = current["chain"]
                current["het_contacts"].append(row)

            elif subtag == "#SUB" and len(data) >= len(SUB_FIELDS):
                row = dict(zip(SUB_FIELDS, data[: len(SUB_FIELDS)]))
                row["protein_chain"] = current["chain"]
                current["sub_contacts"].append(row)

            elif subtag == "DISORDER" and len(data) >= 2:
                current["disorder"].append((data[0], data[1]))

        elif tag == "//":
            current = None

    return records


def heteroatom_key(row):
    """Unique key for a heteroatom instance, matching the site's dropdown label."""
    return (row["het_code"], row["het_resnum"], row["het_chain"])


def heteroatom_label(key):
    code, resnum, chain = key
    return f"{code}:{resnum}:{chain}"


def write_csv(path, fieldnames, rows):
    with open(path, "w", newline="") as fh:
        writer = csv.DictWriter(fh, fieldnames=fieldnames)
        writer.writeheader()
        for row in rows:
            writer.writerow(row)
    return path


def export(records, outdir, pdbid):
    os.makedirs(outdir, exist_ok=True)
    written = []

    # ---- 1. Per-heteroatom interaction tables (mirrors the dropdown on the site) ----
    grouped = defaultdict(list)
    for rec in records:
        for row in rec["het_contacts"]:
            grouped[heteroatom_key(row)].append(row)

    het_out_fields = [
        "protein_chain", "res_start", "res_end", "res_name",
        "contact_type", "num_contacts", "num_hbonds",
        "het_code", "het_resnum", "het_chain", "het_index",
    ]

    all_het_rows = []
    for key, rows in sorted(grouped.items()):
        label = heteroatom_label(key)
        safe_label = re.sub(r"[^A-Za-z0-9]+", "_", label).strip("_")
        fname = f"{pdbid}_{safe_label}_interaction.csv"
        fpath = os.path.join(outdir, fname)
        ordered_rows = [{k: r[k] for k in het_out_fields} for r in rows]
        write_csv(fpath, het_out_fields, ordered_rows)
        written.append(fpath)

        for r in rows:
            r_out = dict(r)
            r_out["heteroatom"] = label
            all_het_rows.append(r_out)

    # ---- 2. Combined table of all protein-heteroatom contacts ----
    combined_fields = ["heteroatom"] + het_out_fields
    combined_path = os.path.join(outdir, f"{pdbid}_all_het_interactions.csv")
    write_csv(
        combined_path,
        combined_fields,
        [{k: r[k] for k in combined_fields} for r in all_het_rows],
    )
    written.append(combined_path)

    # ---- 3. Protein-protein (subunit interface) contacts ----
    sub_rows = []
    for rec in records:
        sub_rows.extend(rec["sub_contacts"])
    if sub_rows:
        sub_fields = ["protein_chain"] + SUB_FIELDS
        sub_path = os.path.join(outdir, f"{pdbid}_subunit_interactions.csv")
        write_csv(sub_path, sub_fields, [{k: r[k] for k in sub_fields} for r in sub_rows])
        written.append(sub_path)

    # ---- 4. Disordered (unmodeled) residue ranges ----
    disorder_rows = []
    for rec in records:
        for start, end in rec["disorder"]:
            disorder_rows.append({"chain": rec["chain"], "start": start, "end": end})
    if disorder_rows:
        dis_path = os.path.join(outdir, f"{pdbid}_disorder.csv")
        write_csv(dis_path, ["chain", "start", "end"], disorder_rows)
        written.append(dis_path)

    return written


def resolve_pdb_id(source) -> str:
    """PDB ID from an ID string, a structure file, or a loaded Structure."""
    from Bio.PDB.Structure import Structure
    from .structure import load, source_info
    if isinstance(source, str) and re.fullmatch(r"[0-9][A-Za-z0-9]{3}", source.strip()) \
            and not os.path.exists(source):
        return source.strip().upper()
    s = source if isinstance(source, Structure) else load(source, fetch_missing=False)
    from .structure import PdbKitError
    pid = source_info(s).entry_id
    if not re.fullmatch(r"[0-9][A-Za-z0-9]{3}", pid or ""):
        raise PdbKitError(f"{source!s} has no 4-character PDB ID (entry ID {pid or 'none'!r}). Het-PDB Navi covers "
                          "PDB entries only - predicted models (AlphaFold) have no interaction data there. "
                          "Pass a PDB ID if the file is a copy of one.")
    return pid.upper()


def interactions(source, *, assembly: str = "assembly1", dat_file: str = None, outdir: str = None,
                 save_raw: bool = False):
    """Parse Het-PDB Navi data for a PDB ID / structure.

    Returns ``(records, written_paths)``; CSVs are written only when
    ``outdir`` is given. ``dat_file`` skips the download.
    """
    pdbid = resolve_pdb_id(source)
    if dat_file:
        with open(dat_file, "r", errors="replace") as fh:
            text = fh.read()
    else:
        text = fetch_data(pdbid, assembly)
    records = parse_records(text)
    if not records:
        from .structure import PdbKitError
        raise PdbKitError(f"No Het-PDB Navi records for {pdbid} - check the ID/file.")
    written = []
    if outdir is not None:
        written = export(records, outdir, pdbid)
        if save_raw and not dat_file:
            raw = os.path.join(outdir, f"{pdbid}_{assembly}.dat")
            with open(raw, "w") as fh:
                fh.write(text)
            written.append(raw)
    return records, written


def het_contacts_table(records):
    """Flat list of protein-heteroatom contact rows (with a 'heteroatom' label)."""
    rows = []
    for rec in records:
        for row in rec["het_contacts"]:
            r = dict(row)
            r["heteroatom"] = heteroatom_label(heteroatom_key(row))
            rows.append(r)
    return rows
