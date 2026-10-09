#!/usr/bin/env python3
"""
Annotate the output of filter_dmr.py with gene names, sequence context,
and write a MEME-ready FASTA for motif discovery.

Usage:
    python annotate_dmr.py \
        --dmr <filtered_bed> \
        --gff3 <bakta_gff3> \
        --fasta <referece_file> \
        --output <output_file> \
        --meme   <output_meme> \
        --window 25
"""
import argparse
import pandas as pd
from Bio import SeqIO
from Bio.Seq import Seq
from utils import parse_gff3_gene_names


def get_sequence_context(chrom, pos, genome, window, strand):
    """
    Return window bp either side of pos (0-based), centered exactly on pos. Pads with N
    at edges.

    The reference FASTA only has the + strand. When the real modified base was called on
    the - strand, the + strand shows its complement at that exact position (a real A reads
    as T, etc.) -- reverse-complementing the window whenever strand is "-" puts every
    sequence in the orientation of the strand that's actually modified, so "the center
    base" always means the real modified base, not its complement.
    """
    seq = genome.get(str(chrom))
    if seq is None:
        return "N" * (2 * window + 1)
    left  = max(0, pos - window)
    right = min(len(seq), pos + window + 1)
    context = ("N" * max(0, window - pos)) + str(seq[left:right]) + ("N" * max(0, pos + window + 1 - len(seq)))
    if strand == "-":
        context = str(Seq(context).reverse_complement())
    return context


def find_genes(chrom, pos, annotation):
    """Return all gene tuples overlapping pos (0-based BED -> 1-based GFF3)."""
    p = pos + 1
    return [g for g in annotation if str(g[0]) == str(chrom) and g[1] <= p <= g[2]]


def parse_args():
    p = argparse.ArgumentParser()
    p.add_argument("--dmr",    required=True, help="Output of filter_dmr.py (TSV with header)")
    p.add_argument("--gff3",   required=True, help="GFF3 annotation file")
    p.add_argument("--fasta",  required=True, help="Reference FASTA")
    p.add_argument("--output", required=True, help="Annotated output TSV")
    p.add_argument("--meme",   required=True, help="MEME-ready FASTA output")
    p.add_argument("--window", type=int, default=25,
                   help="bp each side of methylated base for sequence context (default: 25 -- "
                        "too narrow (~10bp) gives streme too little context to tell apart "
                        "distinct short motifs and it blends them into one blurry composite)")
    return p.parse_args()


def main():
    args = parse_args()

    print("Loading FASTA...")
    genome = {r.id: r.seq for r in SeqIO.parse(args.fasta, "fasta")}

    print("Parsing GFF3...")
    annotation = parse_gff3_gene_names(args.gff3)
    print(f"  {len(annotation)} features loaded")

    # read the header from filter_dmr.py
    print(f"Reading {args.dmr}...")
    df = pd.read_csv(args.dmr, sep="\t")
    print(f"  {len(df)} positions")

    rows = []
    for _, rec in df.iterrows():
        chrom = str(rec["chrom"])
        pos   = int(rec["start"])
        seq   = get_sequence_context(chrom, pos, genome, args.window, rec["strand"])
        hits  = find_genes(chrom, pos, annotation)

        if hits:
            for g in hits:
                rows.append({**rec, "gene": g[4], "gene_start": g[1],
                             "gene_end": g[2], "gene_strand": g[3],
                             "region_type": "genic", "sequence_context": seq})
        else:
            rows.append({**rec, "gene": "intergenic", "gene_start": None,
                         "gene_end": None, "gene_strand": None,
                         "region_type": "intergenic", "sequence_context": seq})

    # pd.DataFrame([]) (zero DMR sites survived filtering, e.g. a rare modification like
    # 4mC under the strict coverage mode) has no columns at all, not even the expected
    # ones -- explicit columns keep every downstream access (region_type, etc.) working on
    # an empty result instead of KeyError'ing.
    out_columns = list(df.columns) + ["gene", "gene_start", "gene_end", "gene_strand",
                                       "region_type", "sequence_context"]
    out = pd.DataFrame(rows, columns=out_columns)
    out[["gene_start", "gene_end"]] = out[["gene_start", "gene_end"]].astype("Int64")
    out.to_csv(args.output, sep="\t", index=False)
    print(f"Annotated table -> {args.output}  ({len(out)} rows)")

    # MEME FASTA - one entry per unique position, simple numeric header
    seen = set()
    meme_i = 0
    with open(args.meme, "w") as fh:
        for _, row in out.iterrows():
            key = (row["chrom"], row["start"])
            if key in seen:
                continue
            seen.add(key)
            meme_i += 1
            fh.write(f">seq{meme_i}\n")
            fh.write(row["sequence_context"] + "\n")
    print(f"MEME FASTA       -> {args.meme}  ({meme_i} sequences)")

    # Save intergenic positions separately
    intergenic_out = args.output.replace(".tsv", "_intergenic.tsv")
    intergenic = out[out["region_type"] == "intergenic"].drop_duplicates(["chrom", "start"])
    intergenic.to_csv(intergenic_out, sep="\t", index=False)
    print(f"Intergenic table -> {intergenic_out}  ({len(intergenic)} positions)")

    # summary
    print("\nTop 10 hits:")
    cols = ["chrom", "start", "score", "a_pct_modified", "b_pct_modified",
            "effect_size", "p_value", "gene", "region_type"]
    print(out[cols].head(10).to_string(index=False))
    print("\nGenic / intergenic breakdown:")
    print(out.drop_duplicates(["chrom", "start"])["region_type"].value_counts().to_string())


if __name__ == "__main__":
    main()
