#!/usr/bin/env python3
"""
Filter and sort modkit dmr pair output.
Usage: python filter_dmr.py input.bed output.bed [--min-coverage 10] [--max-pvalue 0.05] [--min-effect 0.5] [--min-score 5.0] [--chroms 1 2 ...]
"""

import argparse
import pandas as pd

def parse_args():
    parser = argparse.ArgumentParser(description="Filter and sort modkit dmr pair output")
    parser.add_argument("input", help="Input dmr pair BED file")
    parser.add_argument("output", help="Output filtered BED file")
    parser.add_argument("--min-coverage", type=int, default=10,
                        help="Minimum read coverage for both samples (default: 10)")
    parser.add_argument("--max-pvalue", type=float, default=0.05,
                        help="Maximum p-MAP value (default: 0.05)")
    parser.add_argument("--min-effect", type=float, default=0.5,
                        help="Minimum absolute effect size / pct_modified difference (default: 0.5)")
    parser.add_argument("--min-score", type=float, default=None,
                        help="Minimum absolute modkit score (optional, off by default). "
                             "Useful as an alternative/complement to --max-pvalue, since p-value "
                             "underflows to 0 for many sites at high coverage and can't rank them; "
                             "score stays continuous. E.g. the Y. pestis modkit-dmr paper uses "
                             "|score| > 5.0 as their significance cutoff.")
    parser.add_argument("--chroms", nargs="+", default=None,
                        help="Only keep these contigs/chromosomes e.g. --chroms 1 (default: keep all)")
    return parser.parse_args()

def main():
    args = parse_args()

    # Column names based on modkit dmr pair output schema
    col_names = [
        "chrom", "start", "end", "name", "score", "strand",
        "a_counts", "a_total", "b_counts", "b_total",
        "a_pct_mod_raw", "b_pct_mod_raw",
        "a_pct_modified", "b_pct_modified",
        "p_value", "effect_size",
        "cohen_h", "cohen_h_low", "cohen_h_high"
    ]

    print(f"Reading {args.input}...")
    df = pd.read_csv(
        args.input,
        sep="\t",
        header=None,
        names=col_names,
        comment="#"
    )

    print(f"Total positions: {len(df)}")

    # Filter by chromosome/contig if specified
    if args.chroms:
        df = df[df["chrom"].astype(str).isin([str(c) for c in args.chroms])]
        print(f"Positions after contig filter ({', '.join(args.chroms)}): {len(df)}")

    # Convert numeric columns
    for col in ["a_total", "b_total", "a_pct_modified", "b_pct_modified",
                "p_value", "effect_size", "score",
                "cohen_h", "cohen_h_low", "cohen_h_high"]:
        df[col] = pd.to_numeric(df[col], errors="coerce")

    neg = df["cohen_h"] < 0
    df.loc[neg, ["cohen_h_low", "cohen_h_high"]] = -df.loc[neg, ["cohen_h_high", "cohen_h_low"]].values

    mask = (
        (df["a_total"] >= args.min_coverage) &
        (df["b_total"] >= args.min_coverage) &
        (df["p_value"] < args.max_pvalue) &
        (df["effect_size"].abs() >= args.min_effect)
    )
    if args.min_score is not None:
        mask &= df["score"].abs() >= args.min_score

    filtered = df[mask].copy()

    print(f"Positions after filtering: {len(filtered)}")

    # Sort by absolute score descending (most different first)
    filtered["abs_score"] = filtered["score"].abs()
    filtered = filtered.sort_values("abs_score", ascending=False).drop(columns="abs_score")

    filtered.to_csv(args.output, sep="\t", index=False, header=True)
    print(f"Saved to {args.output}")


if __name__ == "__main__":
    main()
