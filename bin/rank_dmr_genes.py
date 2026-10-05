#!/usr/bin/env python3
"""
Rank genes by how many distinct DMR sites (loci) they contain, across every
comparison and modification base fed in.

Usage: python rank_dmr_genes.py --annotated a1.tsv a2.tsv ... --output ranked_genes.tsv
"""
import argparse
import pandas as pd


def parse_args():
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument("--annotated", nargs="+", required=True,
                   help="One or more annotate_dmr.py output TSVs (any mix of comparisons/bases)")
    p.add_argument("--output", required=True, help="Output ranked TSV")
    return p.parse_args()


def main():
    args = parse_args()

    dfs = [pd.read_csv(path, sep="\t") for path in args.annotated]
    combined = pd.concat(dfs, ignore_index=True)

    genic = combined[combined["region_type"] == "genic"].copy()
    genic["site_id"] = genic["chrom"].astype(str) + ":" + genic["start"].astype(str)

    # Bakta (and most annotators) reuse generic names like "Phage tail protein" for many
    # distinct genes scattered across the genome -- group by the actual locus (name +
    # coordinates), not by name alone, or unrelated genes that share a name get merged
    # into one fake "super-gene" with a misleading single coordinate span.
    per_gene = (
        genic.groupby(["gene", "chrom", "gene_start", "gene_end"])
        .agg(
            n_dmr_loci=("site_id", "nunique"),      # distinct genomic positions
            n_hits=("site_id", "size"),             # total rows (>n_dmr_loci if a site recurs across comparisons/bases)
            mean_abs_effect=("effect_size", lambda x: x.abs().mean()),
            max_abs_score=("score", lambda x: x.abs().max()),
        )
        .reset_index()
        .sort_values(["n_dmr_loci", "max_abs_score"], ascending=[False, False])
    )

    per_gene.to_csv(args.output, sep="\t", index=False)
    print(f"Per-locus gene ranking -> {args.output} ({len(per_gene)} distinct gene loci with at least one DMR site)")
    print(per_gene.head(20).to_string(index=False))

    # Separate question: is this *name* (e.g. a generic annotation like "Phage tail
    # protein" that Bakta assigns to many distinct genes genome-wide) disproportionately
    # associated with DMR sites overall? Reported honestly as a name-level count, with
    # how many distinct loci that name spans -- not as a single fake gene span.
    per_name = (
        per_gene.groupby("gene")
        .agg(
            n_distinct_loci_with_this_name=("gene_start", "nunique"),
            n_dmr_loci_total=("n_dmr_loci", "sum"),
            n_hits_total=("n_hits", "sum"),
        )
        .reset_index()
        .sort_values("n_dmr_loci_total", ascending=False)
    )
    name_output = args.output.replace(".tsv", "_by_name.tsv")
    per_name.to_csv(name_output, sep="\t", index=False)
    print(f"\nBy gene *name* (may span multiple distinct loci) -> {name_output} ({len(per_name)} names)")
    print(per_name.head(20).to_string(index=False))


if __name__ == "__main__":
    main()
