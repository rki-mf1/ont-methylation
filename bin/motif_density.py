#!/usr/bin/env python3
"""
motif_density.py
Compute genome-wide motif density and peak annotation from a motifs file + FASTA.
"""
import argparse
import os
import sys
import numpy as np
import pandas as pd
from Bio import SeqIO
from Bio.Seq import Seq
from Bio.SeqUtils import nt_search
from scipy.signal import find_peaks
from utils import parse_gff3_gene_names, select_contigs


def parse_args():
    p = argparse.ArgumentParser()
    p.add_argument("--fasta",                required=True)
    p.add_argument("--gff3",                 required=True)
    p.add_argument("--motif_file",           required=True)
    p.add_argument("--sample",               required=True)
    p.add_argument("--outdir",               required=True)
    p.add_argument("--top_n",                type=int,   default=30)
    p.add_argument("--window_size",          type=int,   default=500)
    p.add_argument("--step_size",            type=int,   default=10)
    p.add_argument("--smoothing_window",     type=int,   default=100)
    p.add_argument("--enrichment_threshold", type=float, default=3.0)
    p.add_argument("--min_absolute_count",   type=int,   default=2)
    p.add_argument("--chromosome_only",      action="store_true", help="Only analyse the largest contig")
    return p.parse_args()


def deduplicate_motifs(motif_df):
    seen, kept = set(), []
    for _, row in motif_df.iterrows():
        motif = row["motif"].upper()
        if motif in seen:
            continue
        rc = str(Seq(motif).reverse_complement()).upper()
        seen.update([motif, rc])
        palindrome = (rc == motif)
        kept.append(dict(motif=motif, mod_code=row["mod_code"],
                         is_palindrome=palindrome, rc_motif=rc))
        print(f"  Keeping: {motif}  (rc={rc}, palindrome={palindrome})")
    return pd.DataFrame(kept)


def find_motif_positions(motif, genome_seq):
    seq = str(genome_seq)
    positions = nt_search(seq, motif)[1:]
    rc = str(Seq(motif).reverse_complement())
    if rc.upper() != motif.upper():
        positions += nt_search(seq, rc)[1:]
    return np.array(sorted(set(positions)))


def compute_genome_density(mod_positions, genome_length, window_size, step_size):
    genome_positions = np.arange(0, genome_length, step_size)
    counts = np.array([
        np.sum((mod_positions >= pos) & (mod_positions < pos + window_size))
        for pos in genome_positions
    ])
    return genome_positions + window_size // 2, counts


def annotate_peaks(peak_positions, peak_heights, peak_counts, annotation):
    rows = []
    for pos, height, count in zip(peak_positions, peak_heights, peak_counts):
        overlapping = [g for g in annotation if g[1] <= pos <= g[2]]
        if overlapping:
            for g in overlapping:
                rows.append(dict(Peak_Position=int(pos), Density=round(float(height), 6),
                                 Count=int(count), Gene=g[4], Gene_Start=g[1],
                                 Gene_End=g[2], Strand=g[3]))
        else:
            rows.append(dict(Peak_Position=int(pos), Density=round(float(height), 6),
                             Count=int(count), Gene="intergenic", Gene_Start=None,
                             Gene_End=None, Strand=None))
    if not rows:
        return pd.DataFrame(columns=["Peak_Position", "Density", "Count", "Gene", "Gene_Start", "Gene_End", "Strand"])
    return pd.DataFrame(rows).sort_values("Density", ascending=False)


def process_contig(positions, contig_length, annotation, args):
    genome_pos, raw_counts = compute_genome_density(
        positions, contig_length, args.window_size, args.step_size)
    density = raw_counts / args.window_size
    density_smooth = (pd.Series(density)
                      .rolling(window=args.smoothing_window, center=True)
                      .mean().bfill().ffill().values)

    expected          = (len(positions) / contig_length) * args.window_size
    enrichment_based  = (expected * args.enrichment_threshold) / args.window_size
    absolute_floor    = args.min_absolute_count / args.window_size
    min_peak_height   = max(enrichment_based, absolute_floor)
    min_peak_distance = args.window_size // args.step_size
    print(f"  Expected/window: {expected:.2f}  |  threshold: {min_peak_height:.5f}")

    peaks_idx, _ = find_peaks(density, height=min_peak_height, distance=min_peak_distance)
    print(f"  Peaks found: {len(peaks_idx)}")
    peaks_df = annotate_peaks(genome_pos[peaks_idx], density[peaks_idx], raw_counts[peaks_idx], annotation)
    density_df = pd.DataFrame({"Position": genome_pos, "Density_raw": density, "Density_smooth": density_smooth})
    return density_df, peaks_df


def process_motif(motif, contigs, reference_genome, annotation, args):
    densities, peaks = [], []
    for contig in contigs:
        positions = find_motif_positions(motif, reference_genome[contig])
        print(f"\n--- {motif} / {contig} ---  occurrences: {len(positions)}")
        if len(positions) < 50:
            print("  Skipping — too few occurrences.")
            continue
        contig_genes = [g for g in annotation if g[0] == contig]
        density_df, peaks_df = process_contig(positions, len(reference_genome[contig]), contig_genes, args)
        density_df.insert(0, "Contig", contig)
        peaks_df.insert(0, "Contig", contig)
        densities.append(density_df)
        peaks.append(peaks_df)

    if not densities:
        return None
    return pd.concat(densities, ignore_index=True), pd.concat(peaks, ignore_index=True)


def save_results(results, outdir):
    os.makedirs(outdir, exist_ok=True)
    for label, (density_df, peaks_df) in results.items():
        safe = label.replace(":", "_")
        density_df.to_csv(os.path.join(outdir, f"density_{safe}.csv"), index=False)
        peaks_df.to_csv(os.path.join(outdir, f"peaks_{safe}.csv"), index=False)
        print(f"  Saved {label}")


def main():
    args = parse_args()

    annotation = parse_gff3_gene_names(args.gff3)
    reference_genome = {r.id: r.seq for r in SeqIO.parse(args.fasta, "fasta")}
    contigs = select_contigs(reference_genome, args.chromosome_only)
    print(f"Contigs: {', '.join(contigs)}")

    motif_df = pd.read_csv(args.motif_file, sep=r'\s+')
    motif_df.columns = motif_df.columns.str.strip()
    print("\nDeduplicating motifs:")
    motif_df_dedup = deduplicate_motifs(motif_df)
    print(f"Kept {len(motif_df_dedup)} motifs after deduplication")

    results = {}
    for _, row in motif_df_dedup.iterrows():
        result = process_motif(row["motif"], contigs, reference_genome, annotation, args)
        if result is not None:
            results[row["motif"]] = result

    if not results:
        print("No results produced — exiting.")
        sys.exit(1)

    save_results(results, args.outdir)


if __name__ == "__main__":
    main()