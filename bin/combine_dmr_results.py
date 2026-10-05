#!/usr/bin/env python3
"""
Combine multiple filter_dmr.py outputs (one per pairwise comparison, all for the
same modification/base) into:
  - one long-format table with a `comparison` column
  - a per-site overlap summary (how many comparisons flag each site, and which)
  - an overlap plot: a schematic 2/3-circle Venn when there are 2 or 3
    comparisons (the only cases a circle-Venn stays readable), otherwise a bar
    chart of "sites shared by exactly k comparisons" (works for any N).

Usage:
    python combine_dmr_results.py --labels A_vs_B A_vs_C --tsvs a.tsv b.tsv \
        --base A --outdir .
"""
import argparse
import os
import pandas as pd
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.patches import Circle


def parse_args():
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument("--labels", nargs="+", required=True, help="One label per comparison, same order as --tsvs")
    p.add_argument("--tsvs", nargs="+", required=True, help="filter_dmr.py output TSVs, same order as --labels")
    p.add_argument("--base", required=True, help="Modification base these comparisons were run for (e.g. A or C)")
    p.add_argument("--outdir", default=".", help="Output directory")
    return p.parse_args()


def draw_venn2(ax, sizes, labels):
    only_a, only_b, both = sizes
    ax.add_patch(Circle((-0.5, 0), 1.2, alpha=0.4, color="tab:blue"))
    ax.add_patch(Circle((0.5, 0), 1.2, alpha=0.4, color="tab:orange"))
    ax.text(-0.9, 0, str(only_a), ha="center", va="center", fontsize=13)
    ax.text(0.9, 0, str(only_b), ha="center", va="center", fontsize=13)
    ax.text(0, 0, str(both), ha="center", va="center", fontsize=13)
    ax.text(-0.5, 1.5, labels[0], ha="center", fontsize=10)
    ax.text(0.5, 1.5, labels[1], ha="center", fontsize=10)
    ax.set_xlim(-2, 2); ax.set_ylim(-2, 2); ax.axis("off")


def draw_venn3(ax, sizes, labels):
    # sizes keyed as: only_a, only_b, only_c, ab, ac, bc, abc
    only_a, only_b, only_c, ab, ac, bc, abc = sizes
    centers = [(-0.5, 0.4), (0.5, 0.4), (0, -0.5)]
    colors = ["tab:blue", "tab:orange", "tab:green"]
    for c, col in zip(centers, colors):
        ax.add_patch(Circle(c, 1.0, alpha=0.35, color=col))
    ax.text(-0.9, 0.7, str(only_a), ha="center", va="center", fontsize=12)
    ax.text(0.9, 0.7, str(only_b), ha="center", va="center", fontsize=12)
    ax.text(0, -1.1, str(only_c), ha="center", va="center", fontsize=12)
    ax.text(0, 0.9, str(ab), ha="center", va="center", fontsize=12)
    ax.text(-0.55, -0.25, str(ac), ha="center", va="center", fontsize=12)
    ax.text(0.55, -0.25, str(bc), ha="center", va="center", fontsize=12)
    ax.text(0, 0.15, str(abc), ha="center", va="center", fontsize=12, fontweight="bold")
    ax.text(-0.9, 1.6, labels[0], ha="center", fontsize=10)
    ax.text(0.9, 1.6, labels[1], ha="center", fontsize=10)
    ax.text(0, -1.7, labels[2], ha="center", fontsize=10)
    ax.set_xlim(-2, 2); ax.set_ylim(-2.2, 2.2); ax.axis("off")


def main():
    args = parse_args()
    os.makedirs(args.outdir, exist_ok=True)

    if len(args.labels) != len(args.tsvs):
        raise SystemExit(f"--labels ({len(args.labels)}) and --tsvs ({len(args.tsvs)}) must have the same length")

    tables = {}
    for label, path in zip(args.labels, args.tsvs):
        df = pd.read_csv(path, sep="\t")
        df.insert(0, "comparison", label)
        tables[label] = df

    combined = pd.concat(tables.values(), ignore_index=True)
    combined_path = os.path.join(args.outdir, f"dmr_combined_{args.base}.tsv")
    combined.to_csv(combined_path, sep="\t", index=False)
    print(f"Combined table -> {combined_path}  ({len(combined)} rows across {len(tables)} comparisons)")

    site_sets = {label: set(zip(df["chrom"], df["start"])) for label, df in tables.items()}
    all_sites = set().union(*site_sets.values()) if site_sets else set()

    rows = []
    for site in all_sites:
        members = [label for label, s in site_sets.items() if site in s]
        rows.append({"chrom": site[0], "start": site[1],
                     "n_comparisons": len(members), "comparisons": ",".join(members)})
    summary = pd.DataFrame(rows).sort_values(["n_comparisons", "chrom", "start"], ascending=[False, True, True])
    summary_path = os.path.join(args.outdir, f"dmr_overlap_{args.base}_summary.tsv")
    summary.to_csv(summary_path, sep="\t", index=False)
    print(f"Overlap summary -> {summary_path}  ({len(summary)} unique sites)")

    plot_path = os.path.join(args.outdir, f"dmr_overlap_{args.base}.png")
    n = len(args.labels)

    if n == 2:
        a, b = site_sets[args.labels[0]], site_sets[args.labels[1]]
        sizes = (len(a - b), len(b - a), len(a & b))
        fig, ax = plt.subplots(figsize=(6, 6))
        draw_venn2(ax, sizes, args.labels)
        ax.set_title(f"6mA/5mC DMR site overlap ({args.base})")
        plt.tight_layout()
        plt.savefig(plot_path, dpi=150, bbox_inches="tight")
    elif n == 3:
        a, b, c = (site_sets[l] for l in args.labels)
        only_a = len(a - b - c); only_b = len(b - a - c); only_c = len(c - a - b)
        ab = len((a & b) - c); ac = len((a & c) - b); bc = len((b & c) - a)
        abc = len(a & b & c)
        fig, ax = plt.subplots(figsize=(6, 6))
        draw_venn3(ax, (only_a, only_b, only_c, ab, ac, bc, abc), args.labels)
        ax.set_title(f"6mA/5mC DMR site overlap ({args.base})")
        plt.tight_layout()
        plt.savefig(plot_path, dpi=150, bbox_inches="tight")
    else:
        # a circle-Venn stops being readable past 3 sets; fall back to a
        # bar chart of "sites shared by exactly k comparisons", which stays
        # informative for any number of comparisons.
        counts = summary["n_comparisons"].value_counts().sort_index()
        fig, ax = plt.subplots(figsize=(6, 5))
        ax.bar(counts.index.astype(str), counts.values, color="tab:blue", alpha=0.7)
        ax.set_xlabel("number of comparisons a site is significant in")
        ax.set_ylabel("number of sites")
        ax.set_title(f"DMR site reproducibility across {n} comparisons ({args.base})")
        plt.tight_layout()
        plt.savefig(plot_path, dpi=150, bbox_inches="tight")

    print(f"Overlap plot -> {plot_path}")


if __name__ == "__main__":
    main()
