#!/usr/bin/env python3
"""
Combine multiple filter_dmr.py outputs (one per pairwise comparison, all for the
same modification/base) into:
  - one long-format table with a `comparison` column
  - a per-site overlap summary (how many comparisons flag each site, and which)
  - an overlap plot: a schematic 2/3-circle Venn when there are 2 or 3
    comparisons (the only cases a circle-Venn stays readable), an UpSet plot
    for 4 or more, and a bar chart of "sites shared by exactly k comparisons"
    as the fallback (a single comparison, or nothing to intersect).

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
    p.add_argument("--base", required=True, help="Modification these comparisons were run for (e.g. 6mA, 5mC, 4mC)")
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


def draw_upset(site_sets, labels, base, plot_path, max_combos=30):
    """UpSet plot in plain matplotlib (no extra dependency): one bar per exact
    combination of comparisons a site is significant in, with a dot matrix underneath
    and the per-comparison totals on the left. Returns False if there is nothing to draw."""
    members = {}
    for label in labels:
        for site in site_sets[label]:
            members.setdefault(site, []).append(label)
    combo_counts = {}
    for labs in members.values():
        key = tuple(l for l in labels if l in labs)
        combo_counts[key] = combo_counts.get(key, 0) + 1
    if not combo_counts:
        return False

    # biggest bars first; ties broken by fewer sets, then label order, so the layout is stable
    ranked = sorted(combo_counts.items(), key=lambda kv: (-kv[1], len(kv[0]), [labels.index(l) for l in kv[0]]))
    shown = ranked[:max_combos]
    n_sets, n_combos = len(labels), len(shown)

    fig = plt.figure(figsize=(max(7, 0.45 * n_combos + 4), 4.5 + 0.35 * n_sets))
    # three columns: per-comparison totals | room for the comparison labels | main plot
    gs = fig.add_gridspec(2, 3, width_ratios=[1, 1.4, max(3, 0.45 * n_combos)],
                          height_ratios=[3, max(1.2, 0.35 * n_sets)], hspace=0.05, wspace=0.03)
    ax_bar = fig.add_subplot(gs[0, 2])
    ax_mat = fig.add_subplot(gs[1, 2], sharex=ax_bar)
    ax_tot = fig.add_subplot(gs[1, 0], sharey=ax_mat)

    xs = range(n_combos)
    ax_bar.bar(xs, [c for _, c in shown], color="tab:blue", alpha=0.8)
    for x, (_, c) in zip(xs, shown):
        ax_bar.text(x, c, str(c), ha="center", va="bottom", fontsize=8)
    ax_bar.set_ylabel("sites in exactly this combination")
    ax_bar.set_ylim(0, max(c for _, c in shown) * 1.12)
    ax_bar.tick_params(axis="x", bottom=False, labelbottom=False)
    for side in ("top", "right"):
        ax_bar.spines[side].set_visible(False)
    title = f"DMR site overlap across {n_sets} comparisons ({base})"
    if len(ranked) > n_combos:
        title += f"\ntop {n_combos} of {len(ranked)} combinations"
    ax_bar.set_title(title)

    # one row per comparison, first label on top
    ypos = {l: n_sets - 1 - i for i, l in enumerate(labels)}
    for y in ypos.values():
        ax_mat.axhspan(y - 0.5, y + 0.5, color="0.95" if y % 2 else "white", zorder=0)
    for x, (combo, _) in zip(xs, shown):
        for l in labels:
            ax_mat.scatter(x, ypos[l], s=40, color="black" if l in combo else "0.82", zorder=3)
        if len(combo) > 1:
            ys = [ypos[l] for l in combo]
            ax_mat.plot([x, x], [min(ys), max(ys)], color="black", lw=2, zorder=2)
    ax_mat.set_yticks([ypos[l] for l in labels])
    ax_mat.set_ylim(-0.5, n_sets - 0.5)
    ax_mat.set_xlim(-0.7, n_combos - 0.3)
    ax_mat.tick_params(axis="x", bottom=False, labelbottom=False)
    ax_mat.set_yticklabels(labels, fontsize=9)
    ax_mat.tick_params(axis="y", left=False, pad=6)
    for side in ("top", "right", "left", "bottom"):
        ax_mat.spines[side].set_visible(False)

    totals = [len(site_sets[l]) for l in labels]
    ax_tot.barh([ypos[l] for l in labels], totals, color="tab:gray", alpha=0.8)
    for l, t in zip(labels, totals):
        ax_tot.text(t * 1.04, ypos[l], str(t), va="center", ha="right", fontsize=8)
    ax_tot.set_xlim(max(totals) * 1.4, 0)
    ax_tot.tick_params(axis="y", left=False, labelleft=False)
    ax_tot.set_xlabel("sites per comparison", fontsize=8)
    for side in ("top", "left", "right"):
        ax_tot.spines[side].set_visible(False)

    fig.savefig(plot_path, dpi=150, bbox_inches="tight")
    plt.close(fig)
    return True


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

    # gene name is purely a function of genomic position (independent of which
    # comparison found the site), so any one table that has it is enough to look it up --
    # only present when annotate_dmr.py's output (not the bare filter_dmr.py one) was
    # passed in, i.e. when --gff3 was given to the DMR flow.
    gene_by_site = {}
    for df in tables.values():
        if "gene" not in df.columns:
            continue
        for chrom, start, gene in zip(df["chrom"], df["start"], df["gene"]):
            gene_by_site.setdefault((chrom, start), gene)

    rows = []
    for site in all_sites:
        members = [label for label, s in site_sets.items() if site in s]
        rows.append({"chrom": site[0], "start": site[1], "gene": gene_by_site.get(site, ""),
                     "n_comparisons": len(members), "comparisons": ",".join(members)})
    # pd.DataFrame([]) (no DMR sites in any comparison, e.g. a rare modification like 4mC
    # under the strict coverage mode) has no columns at all -- explicit columns keep
    # sort_values/value_counts below working on an empty result instead of KeyError'ing.
    summary = pd.DataFrame(rows, columns=["chrom", "start", "gene", "n_comparisons", "comparisons"]) \
        .sort_values(["n_comparisons", "chrom", "start"], ascending=[False, True, True])
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
        ax.set_title(f"DMR site overlap ({args.base})")
        plt.tight_layout()
        plt.savefig(plot_path, dpi=150, bbox_inches="tight")
    elif n == 3:
        a, b, c = (site_sets[l] for l in args.labels)
        only_a = len(a - b - c); only_b = len(b - a - c); only_c = len(c - a - b)
        ab = len((a & b) - c); ac = len((a & c) - b); bc = len((b & c) - a)
        abc = len(a & b & c)
        fig, ax = plt.subplots(figsize=(6, 6))
        draw_venn3(ax, (only_a, only_b, only_c, ab, ac, bc, abc), args.labels)
        ax.set_title(f"DMR site overlap ({args.base})")
        plt.tight_layout()
        plt.savefig(plot_path, dpi=150, bbox_inches="tight")
    elif n >= 4 and draw_upset(site_sets, args.labels, args.base, plot_path):
        # a circle-Venn stops being readable past 3 sets, an UpSet plot doesn't
        pass
    else:
        # single comparison, or no sites at all: bar chart of "sites shared by exactly
        # k comparisons", which stays informative for any number of comparisons.
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
