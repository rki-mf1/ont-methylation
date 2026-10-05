#!/usr/bin/env python3
"""
Volcano-style plot for a single modkit dmr pair comparison: effect size (x) vs
|score| (y) -- |score| is used instead of -log10(p_value) because modkit's
p-value routinely underflows to exactly 0.0 at high coverage, which would
otherwise collapse many of the strongest hits onto the same y position.

Usage:
    python dmr_volcano_plot.py --raw dmr_A_B.bed --filtered dmr_A_B_filtered.tsv \
        --output volcano.png --title "A vs B" [--min-coverage 10] [--min-effect 0.5] \
        [--max-background 200000]
"""
import argparse
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import pandas as pd

RAW_COLUMNS = [
    "chrom", "start", "end", "name", "score", "strand",
    "a_counts", "a_total", "b_counts", "b_total",
    "a_pct_mod_raw", "b_pct_mod_raw",
    "a_pct_modified", "b_pct_modified",
    "p_value", "effect_size",
    "cohen_h", "cohen_h_low", "cohen_h_high"
]


def parse_args():
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument("--raw", required=True, help="Raw modkit dmr pair output (.bed, no header)")
    p.add_argument("--filtered", required=True, help="filter_dmr.py output (.tsv, with header) -- defines which sites are highlighted")
    p.add_argument("--output", required=True, help="Output PNG path")
    p.add_argument("--title", default="", help="Plot title")
    p.add_argument("--min-coverage", type=int, default=10,
                    help="Coverage floor for points shown in the background (default: 10)")
    p.add_argument("--min-effect", type=float, default=0.5,
                    help="Effect-size threshold to draw as dashed guide lines (default: 0.5)")
    p.add_argument("--max-background", type=int, default=200000,
                    help="Downsample non-selected points above this count, for speed/file size (default: 200000)")
    return p.parse_args()


def main():
    args = parse_args()

    df = pd.read_csv(args.raw, sep="\t", header=None, names=RAW_COLUMNS, comment="#")
    for col in ["a_total", "b_total", "effect_size", "score"]:
        df[col] = pd.to_numeric(df[col], errors="coerce")
    df = df[(df.a_total >= args.min_coverage) & (df.b_total >= args.min_coverage)].copy()
    df["abs_score"] = df["score"].abs()

    selected = pd.read_csv(args.filtered, sep="\t")
    selected_keys = set(zip(selected["chrom"], selected["start"]))
    df["is_selected"] = [k in selected_keys for k in zip(df["chrom"], df["start"])]

    background = df[~df.is_selected]
    if len(background) > args.max_background:
        background = background.sample(args.max_background, random_state=0)
    hits = df[df.is_selected]

    fig, ax = plt.subplots(figsize=(7, 6))
    ax.scatter(background.effect_size, background.abs_score, s=3, color="lightgray", alpha=0.4,
               label=f"not selected (n={len(df) - len(hits)})")
    ax.scatter(hits.effect_size, hits.abs_score, s=6, color="crimson", alpha=0.7,
               label=f"selected by filter_dmr.py (n={len(hits)})")
    ax.axvline(args.min_effect, color="black", linestyle="--", linewidth=0.7)
    ax.axvline(-args.min_effect, color="black", linestyle="--", linewidth=0.7)
    ax.set_xlim(-1.05, 1.05)
    ax.set_xlabel("effect size (a_pct_modified - b_pct_modified)")
    ax.set_ylabel("|score|  (modkit's magnitude-of-difference statistic)")
    ax.set_title(args.title or "DMR: effect size vs |score|")
    ax.legend(loc="upper left", fontsize=8)
    plt.tight_layout()
    plt.savefig(args.output, dpi=150, bbox_inches="tight")
    print(f"Saved {args.output}  (n={len(df)} plotted, {len(hits)} selected)")


if __name__ == "__main__":
    main()
