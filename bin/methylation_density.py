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
from utils import parse_gff3_gene_names, read_modkit


MOD_CODES = {
    "6mA": "a",
    "5mC": "m",
    "4mC": 21839,
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
    p.add_argument("--top_n",                      type=int,   default=30)
    p.add_argument("--window_size",                type=int,   default=500)
    p.add_argument("--step_size",                  type=int,   default=10)
    p.add_argument("--smoothing_window",           type=int,   default=100)
    p.add_argument("--enrichment_threshold",       type=float, default=3.0)
    return p.parse_args()


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
    return pd.DataFrame(rows).sort_values("Density", ascending=False)


def process_modification(label, args, modkit_file, largest_contig, genome_length, annotation):
    mod_code = MOD_CODES.get(label)
    if mod_code is None:
        print(f"WARNING: unknown modification '{label}', skipping.")
        return None

    mod = read_modkit(modkit_file, args.percent_modified_threshold, mod_code)
    mod["Contig"] = mod["Contig"].astype(str).str.strip()
    mod = mod[mod["Contig"] == largest_contig].copy()

    print(f"\n--- {label} ---  sites: {len(mod)}")
    if len(mod) < 50:
        print("  Skipping — too few sites.")
        return None

    mod_positions = mod["Position"].values
    genome_pos, raw_counts = compute_genome_density(
        mod_positions, genome_length, args.window_size, args.step_size)
    density = raw_counts / args.window_size
    density_smooth = (pd.Series(density)
                      .rolling(window=args.smoothing_window, center=True)
                      .mean().bfill().ffill().values)

    expected = (len(mod_positions) / genome_length) * args.window_size
    min_peak_height = (expected * args.enrichment_threshold) / args.window_size
    min_peak_distance = args.window_size // args.step_size
    print(f"  Expected/window: {expected:.1f}  |  threshold: {min_peak_height:.5f}")

    peaks_idx, _ = find_peaks(density, height=min_peak_height, distance=min_peak_distance)
    peak_positions = genome_pos[peaks_idx]
    peak_heights   = density[peaks_idx]
    peak_counts    = raw_counts[peaks_idx]
    print(f"  Peaks found: {len(peaks_idx)}")

    peaks_df = annotate_peaks(peak_positions, peak_heights, peak_counts, annotation)

    return genome_pos, density_smooth, peaks_df, peak_positions, peak_heights


def save_results(results, outdir):
    os.makedirs(outdir, exist_ok=True)
    for label, (genome_pos, density_smooth, peaks_df, _, _) in results.items():
        pd.DataFrame({"Position": genome_pos, "Density_smooth": density_smooth}) \
            .to_csv(os.path.join(outdir, f"density_{label}.csv"), index=False)
        peaks_df.to_csv(
            os.path.join(outdir, f"peaks_{label}.csv"), index=False)
        print(f"  Saved {label}")


def main():
    args = parse_args()

    annotation = parse_gff3_gene_names(args.gff3)
    reference_genome = {r.id: r.seq for r in SeqIO.parse(args.fasta, "fasta")}
    largest_contig = max(reference_genome, key=lambda k: len(reference_genome[k]))
    genome_length  = len(reference_genome[largest_contig])
    print(f"Contig: {largest_contig}  length: {genome_length}")
    annotation = [g for g in annotation if g[0] == largest_contig]

    requested = [m.strip() for m in args.modifications.split(",")]
    results = {}
    for label in requested:
        result = process_modification(
            label, args, args.modkit_bed,
            largest_contig, genome_length, annotation)
        if result is not None:
            results[label] = result

    if not results:
        print("No results produced — exiting.")
        sys.exit(1)

    save_results(results, args.outdir)


if __name__ == "__main__":
    main()