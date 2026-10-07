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


FIELDS = ["core", "p_value", "n_sites", "consensus", "width", "distance_to_center", "is_palindromic", "streme_rank"]

# streme extends a motif's window in both directions to maximize its statistical score,
# which routinely pads a short real motif with low-information degenerate IUPAC positions
# (e.g. the real motif AGCTGC reported as AGCTGCBSSV) -- stripping edge positions that
# aren't a plain A/C/G/T recovers the actual conserved core for a human to read at a glance.
AMBIGUOUS_IUPAC = set("NRYSWKMBDHV")


def trim_to_core(consensus):
    chars = list(consensus)
    start = 0
    while start < len(chars) and chars[start] in AMBIGUOUS_IUPAC:
        start += 1
    end = len(chars)
    while end > start and chars[end - 1] in AMBIGUOUS_IUPAC:
        end -= 1
    core = "".join(chars[start:end])
    return core or consensus


def parse_args():
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument("--xml", required=True, help="streme.xml")
    p.add_argument("--output", required=True, help="Output summary TSV")
    return p.parse_args()


def main():
    args = parse_args()
    root = ET.parse(args.xml).getroot()

    MIN_CORE_LEN = 4  # a 1-3nt "core" after trimming is too short to mean anything on its
    # own (e.g. just "A") -- drop it from this clean table rather than show a non-motif; the
    # full streme report (published separately, see discover_motifs' streme_raw output) still
    # has it for anyone who wants to look.

    rows = []
    for motif in root.iter("motif"):
        motif_id = motif.get("id", "")
        rank, _, consensus = motif_id.partition("-")
        consensus = consensus or motif_id
        core = trim_to_core(consensus)
        if len(core) < MIN_CORE_LEN:
            continue
        rows.append({
            "streme_rank": int(rank) if rank.isdigit() else None,
            "core": core,
            "consensus": consensus,
            "width": motif.get("width"),
            "n_sites": motif.get("train_pos_count"),
            "p_value": motif.get("train_pvalue"),
            "distance_to_center": motif.get("train_dtc"),
            "is_palindromic": motif.get("is_palindromic"),
        })

    # streme's own "rank" is its discovery/erasure order, not sorted by significance (a
    # weak rank-2 hit can have a far worse p-value than rank-5) -- sort by p-value so the
    # strongest, most trustworthy motifs are always at the top.
    rows.sort(key=lambda r: (r["p_value"] is None, float(r["p_value"]) if r["p_value"] is not None else 0))

    with open(args.output, "w", newline="") as fh:
        writer = csv.DictWriter(fh, fieldnames=FIELDS, delimiter="\t")
        writer.writeheader()
        writer.writerows(rows)

    print(f"Motif summary -> {args.output} ({len(rows)} motifs)")
    for r in rows:
        print(f"  {r['core']}  (full={r['consensus']}, width={r['width']}, n_sites={r['n_sites']}, "
              f"p={r['p_value']}, dist_to_center={r['distance_to_center']}, palindromic={r['is_palindromic']})")
    if not rows:
        print("No motifs found.")


if __name__ == "__main__":
    main()
