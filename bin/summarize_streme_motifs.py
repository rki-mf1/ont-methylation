#!/usr/bin/env python3
"""
Parse a streme.xml into a clean, sorted per-motif summary TSV -- instead of
having to read the raw streme.txt by hand.

Stdlib only (no pandas) -- this runs inside the MEME Suite container, which
doesn't have Python data-science packages installed.

Usage: python summarize_streme_motifs.py --xml streme.xml --output motifs_summary.tsv
"""
import argparse
import csv
import xml.etree.ElementTree as ET


FIELDS = ["rank", "consensus", "width", "n_sites", "p_value", "distance_to_center", "is_palindromic"]


def parse_args():
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument("--xml", required=True, help="streme.xml")
    p.add_argument("--output", required=True, help="Output summary TSV")
    return p.parse_args()


def main():
    args = parse_args()
    root = ET.parse(args.xml).getroot()

    rows = []
    for motif in root.iter("motif"):
        motif_id = motif.get("id", "")
        rank, _, consensus = motif_id.partition("-")
        rows.append({
            "rank": int(rank) if rank.isdigit() else None,
            "consensus": consensus or motif_id,
            "width": motif.get("width"),
            "n_sites": motif.get("train_pos_count"),
            "p_value": motif.get("train_pvalue"),
            "distance_to_center": motif.get("train_dtc"),
            "is_palindromic": motif.get("is_palindromic"),
        })

    rows.sort(key=lambda r: (r["rank"] is None, r["rank"]))

    with open(args.output, "w", newline="") as fh:
        writer = csv.DictWriter(fh, fieldnames=FIELDS, delimiter="\t")
        writer.writeheader()
        writer.writerows(rows)

    print(f"Motif summary -> {args.output} ({len(rows)} motifs)")
    for r in rows:
        print(f"  #{r['rank']} {r['consensus']}  width={r['width']} n_sites={r['n_sites']} "
              f"p={r['p_value']} dist_to_center={r['distance_to_center']} palindromic={r['is_palindromic']}")
    if not rows:
        print("No motifs found.")


if __name__ == "__main__":
    main()
