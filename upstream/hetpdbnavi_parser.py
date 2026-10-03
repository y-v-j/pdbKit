#!/usr/bin/env python3
"""
hetpdbnavi_parser.py
=====================

Downloads (or reads) the interaction data file for a PDB entry from
Het-PDB Navi (https://hetpdbnavi.nagahama-i-bio.ac.jp/) and converts the
per-heteroatom contact tables shown in the "Interaction" tab of the
web page into CSV files -- one per heteroatom instance, matching the
entries in the "Select hetatm" dropdown (e.g. "ADP:1201:A").

The site's raw data (obtained from the "Download interaction data" link,
i.e. https://hetpdbnavi.nagahama-i-bio.ac.jp/gethetpdbdata/<PDBID>/assembly1)
is a SwissProt/EMBL-style flat file with one record per protein chain,
separated by "//". Within each record, "FT" lines encode contacts:

  FT  #HET  <res_start> <res_end> <res_name> <res_chain> \
            <het_index> <het_resnum> <het_code> <het_chain> \
            <contact_type> <num_contacts> <num_hbonds>

  FT  #SUB  <res1_start> <res1_end> <res1_name> <res1_chain> \
            <res2_start> <res2_end> <res2_name> <res2_chain> \
            <interaction_type> <contact_type> <num_contacts>

  FT  DISORDER  <start>  <end>

This script parses those "FT" lines and writes:

  1. <PDBID>_<CODE>_<RESNUM>_<CHAIN>_interaction.csv
        one file per heteroatom instance (protein-heteroatom contacts) --
        equivalent to selecting that entry in the "Select hetatm" dropdown
        on the Interaction tab.
  2. <PDBID>_all_het_interactions.csv
        all protein-heteroatom contacts combined, with a "Heteroatom"
        column identifying which ligand instance each row belongs to.
  3. <PDBID>_subunit_interactions.csv
        protein-protein (chain-chain) interface contacts, if present.
  4. <PDBID>_disorder.csv
        residue ranges not resolved in the structure, per chain.

Usage
-----
    # Download directly from Het-PDB Navi and write CSVs to ./out
    python hetpdbnavi_parser.py 9VMH --outdir out

    # Use an already-downloaded .dat file instead of downloading
    python hetpdbnavi_parser.py 9VMH --file 9vmh.dat --outdir out

    # Use a different biological assembly
    python hetpdbnavi_parser.py 9VMH --assembly assembly2 --outdir out
"""

import argparse
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


def fetch_data(pdbid: str, assembly: str = "assembly1", timeout: int = 30) -> str:
    """Download the raw .dat interaction file for a PDB ID from Het-PDB Navi."""
    url = DOWNLOAD_URL_TEMPLATE.format(pdbid=pdbid.upper(), assembly=assembly)
    req = urllib.request.Request(url, headers={"User-Agent": "Mozilla/5.0"})
    with urllib.request.urlopen(req, timeout=timeout) as resp:
        return resp.read().decode("utf-8", errors="replace")


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

            if subtag == "#HET" and len(data) >= len(HET_FIELDS):
                row = dict(zip(HET_FIELDS, data[: len(HET_FIELDS)]))
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


def main():
    parser = argparse.ArgumentParser(description="Convert Het-PDB Navi interaction data to CSV.")
    parser.add_argument("pdbid", help="4-character PDB ID, e.g. 9VMH")
    parser.add_argument("--file", help="Path to an already-downloaded .dat file "
                                        "(skips downloading from the web).")
    parser.add_argument("--assembly", default="assembly1",
                         help="Biological assembly to download (default: assembly1).")
    parser.add_argument("--outdir", default=".", help="Directory to write CSV files into.")
    args = parser.parse_args()

    pdbid = args.pdbid.upper()

    if args.file:
        with open(args.file, "r", errors="replace") as fh:
            text = fh.read()
    else:
        try:
            text = fetch_data(pdbid, args.assembly)
        except Exception as exc:  # noqa: BLE001
            print(f"Failed to download data for {pdbid}: {exc}", file=sys.stderr)
            sys.exit(1)

    records = parse_records(text)
    if not records:
        print("No records parsed -- check the input file/PDB ID.", file=sys.stderr)
        sys.exit(1)

    written = export(records, args.outdir, pdbid)

    print(f"Parsed {len(records)} chain record(s) for {pdbid}.")
    print("Wrote:")
    for path in written:
        print(f"  {path}")


if __name__ == "__main__":
    main()
