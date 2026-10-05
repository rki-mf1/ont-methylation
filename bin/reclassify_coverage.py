#!/usr/bin/env python3
"""
Rewrite a modkit pileup bedMethyl file so that "coverage" means total reads
covering a site (including below-threshold/ambiguous and other/mismatched
reads), not just modkit's own confidently-classified subset -- while keeping
the bedMethyl format itself unchanged, so it can be fed straight into
`modkit dmr pair` (or anything else) exactly as before.

Ambiguous reads (below-threshold + other-bases) are folded into the
"canonical/unmodified" count. Modified-read counts are untouched. This
changes what modkit's own statistics are computed from, without touching
modkit itself.

Accepts plain or gzip-compressed input (detected from the .gz extension).
Output is always plain text (bgzip it afterward if needed).

Usage: python reclassify_coverage.py input.bed[.gz] output.bed
"""
import argparse
import gzip
import sys


def parse_args():
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument("input", help="Input modkit pileup bedMethyl file (plain or .gz)")
    p.add_argument("output", help="Output reclassified bedMethyl file (plain text)")
    return p.parse_args()


def main():
    args = parse_args()

    opener = gzip.open if args.input.endswith(".gz") else open
    n_lines = 0
    with opener(args.input, "rt") as fin, open(args.output, "w") as fout:
        for line in fin:
            f = line.rstrip("\n").split("\t")
            n_mod = int(f[11])
            n_unmod = int(f[12])
            n_other_mod = int(f[13])
            n_below_thresh = int(f[15])
            n_other_bases = int(f[16])

            new_unmod = n_unmod + n_below_thresh + n_other_bases
            new_valid_cov = n_mod + new_unmod + n_other_mod
            new_pct = round(100 * n_mod / new_valid_cov, 2) if new_valid_cov else 0.0

            f[4] = str(new_valid_cov)   # score column mirrors Nvalid_cov
            f[9] = str(new_valid_cov)
            f[10] = f"{new_pct:.2f}"
            f[12] = str(new_unmod)
            f[15] = "0"
            f[16] = "0"

            fout.write("\t".join(f) + "\n")
            n_lines += 1

    print(f"Reclassified {n_lines} lines -> {args.output}", file=sys.stderr)


if __name__ == "__main__":
    main()
