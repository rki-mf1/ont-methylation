#!/usr/bin/env python3
"""
methylation_density.py
Compute genome-wide methylation density and peak annotation from a modkit pileup BED.
"""
import argparse
import os
import sys
import numpy as np
import pandas as pd
from Bio import SeqIO
from scipy.signal import find_peaks
from utils import parse_gff3_gene_names, read_modkit, select_contigs


MOD_CODES = {
    "6mA": "a",
    "5mC": "m",
    "4mC": "21839",
}


def parse_args():
    p = argparse.ArgumentParser()
    p.add_argument("--modkit_bed",                 required=True)
    p.add_argument("--fasta",                      required=True)
    p.add_argument("--gff3",                       required=True)
    p.add_argument("--sample",                     required=True)
    p.add_argument("--outdir",                     required=True)
    p.add_argument("--modifications",              default="6mA,5mC,4mC")
    p.add_argument("--percent_modified_threshold", type=float, default=0.5)
    p.add_argument("--min_coverage",               type=int,   default=10)
    p.add_argument("--top_n",                      type=int,   default=30)
    p.add_argument("--window_size",                type=int,   default=500)
    p.add_argument("--step_size",                  type=int,   default=10)
    p.add_argument("--smoothing_window",           type=int,   default=100)
    p.add_argument("--enrichment_threshold",       type=float, default=3.0)
    p.add_argument("--min_sites",                  type=int,   default=50)
    p.add_argument("--chromosome_only",            action="store_true", help="Only analyse the largest contig")
    return p.parse_args()


def compute_genome_density(mod_positions, genome_length, window_size, step_size):
    genome_positions = np.arange(0, genome_length, step_size)
    counts = np.array([
        np.sum((mod_positions >= pos) & (mod_positions < pos + window_size))
        for pos in genome_positions
    ])
    return genome_positions + window_size // 2, counts


PEAKS_COLUMNS = ["Peak_Position", "Density", "Count", "Gene", "Gene_Start", "Gene_End", "Strand"]


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
        return pd.DataFrame(columns=PEAKS_COLUMNS)
    return pd.DataFrame(rows).sort_values("Density", ascending=False)


def process_contig(label, args, mod_positions, contig_length, annotation):
    genome_pos, raw_counts = compute_genome_density(
        mod_positions, contig_length, args.window_size, args.step_size)
    density = raw_counts / args.window_size
    density_smooth = (pd.Series(density)
                      .rolling(window=args.smoothing_window, center=True)
                      .mean().bfill().ffill().values)

    expected = (len(mod_positions) / contig_length) * args.window_size
    min_peak_height = (expected * args.enrichment_threshold) / args.window_size
    min_peak_distance = args.window_size // args.step_size
    print(f"  Expected/window: {expected:.1f}  |  threshold: {min_peak_height:.5f}")

    peaks_idx, _ = find_peaks(density, height=min_peak_height, distance=min_peak_distance)
    print(f"  Peaks found: {len(peaks_idx)}")
    peaks_df = annotate_peaks(genome_pos[peaks_idx], density[peaks_idx], raw_counts[peaks_idx], annotation)
    density_df = pd.DataFrame({"Position": genome_pos, "Density_smooth": density_smooth})
    return density_df, peaks_df


def process_modification(label, args, contigs, reference_genome, annotation):
    mod_code = MOD_CODES.get(label)
    if mod_code is None:
        print(f"WARNING: unknown modification '{label}', skipping.")
        return None

    mod = read_modkit(args.modkit_bed, args.percent_modified_threshold, mod_code, args.min_coverage)
    mod["Contig"] = mod["Contig"].astype(str).str.strip()

    densities, peaks = [], []
    for contig in contigs:
        positions = mod.loc[mod["Contig"] == contig, "Position"].values
        print(f"\n--- {label} / {contig} ---  sites: {len(positions)}")
        if len(positions) < args.min_sites:
            print(f"  Skipping — too few sites (< {args.min_sites}).")
            continue
        contig_genes = [g for g in annotation if g[0] == contig]
        density_df, peaks_df = process_contig(label, args, positions, len(reference_genome[contig]), contig_genes)
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
        density_df.to_csv(os.path.join(outdir, f"density_{label}.csv"), index=False)
        peaks_df.to_csv(os.path.join(outdir, f"peaks_{label}.csv"), index=False)
        print(f"  Saved {label}")


def main():
    args = parse_args()

    annotation = parse_gff3_gene_names(args.gff3)
    reference_genome = {r.id: r.seq for r in SeqIO.parse(args.fasta, "fasta")}
    contigs = select_contigs(reference_genome, args.chromosome_only)
    print(f"Contigs: {', '.join(contigs)}")

    requested = [m.strip() for m in args.modifications.split(",")]
    results = {}
    for label in requested:
        result = process_modification(label, args, contigs, reference_genome, annotation)
        if result is not None:
            results[label] = result

    if not results:
        print("No results produced — exiting.")
        sys.exit(1)

    save_results(results, args.outdir)


if __name__ == "__main__":
    main()