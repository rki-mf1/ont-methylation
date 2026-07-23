#!/usr/bin/env python3
"""
combine_peaks.py
Merge peak tables from methylation (and motif) results into a single wide CSV.
"""
import argparse
import glob
import os
import sys
import pandas as pd


def parse_args():
    p = argparse.ArgumentParser()
    p.add_argument("--meth_peaks",  nargs="+", required=True)
    p.add_argument("--motif_peaks", nargs="+", default=[])
    p.add_argument("--outdir",      default=".")
    return p.parse_args()


INDEX_COLS = ["Gene", "Gene_Start", "Gene_End", "Strand"]

def load_peak_files(filepaths, source_type):
    dfs = []
    for filepath in filepaths:
        label = os.path.basename(filepath).replace("peaks_", "").replace(".csv", "")
        df = pd.read_csv(filepath)
        if df.empty:
            continue
        df["Signal"]      = label
        df["Source_type"] = source_type
        dfs.append(df)
        print(f"  Loaded {len(df)} rows from {os.path.basename(filepath)}")
    return dfs


def main():
    args = parse_args()
    all_dfs = load_peak_files(args.meth_peaks, "methylation")
    if args.motif_peaks:
        all_dfs += load_peak_files(args.motif_peaks, "motif")

    if not all_dfs:
        print("No tables found — exiting.")
        sys.exit(1)

    os.makedirs(args.outdir, exist_ok=True)

    combined = pd.concat(all_dfs, ignore_index=True)
    print(f"\nTotal rows: {len(combined)}  |  Unique genes: {combined['Gene'].nunique()}")

    # Long format
    long_file = os.path.join(args.outdir, "all_peaks_long.csv")
    combined.to_csv(long_file, index=False)

    meth_combined  = combined[combined["Source_type"] == "methylation"]
    motif_combined = combined[combined["Source_type"] == "motif"]

    # Wide pivot — methylation
    pivot_meth = pd.DataFrame(columns=INDEX_COLS)
    meth_signal_cols = []
    if not meth_combined.empty:
        pivot_meth = meth_combined.pivot_table(
            index=INDEX_COLS, columns="Signal", values="Count", aggfunc="sum"
        ).reset_index()
        pivot_meth.columns.name = None
        meth_signal_cols = [c for c in pivot_meth.columns if c not in INDEX_COLS]

    # Wide pivot — motifs
    pivot_motif = pd.DataFrame(columns=INDEX_COLS)
    motif_signal_cols = []
    if not motif_combined.empty:
        pivot_motif = motif_combined.pivot_table(
            index=INDEX_COLS, columns="Signal", values="Count", aggfunc="sum"
        ).reset_index()
        pivot_motif.columns.name = None
        motif_signal_cols = [c for c in pivot_motif.columns if c not in INDEX_COLS]
        pivot_motif[motif_signal_cols] = pivot_motif[motif_signal_cols].fillna(0).astype(int)

    # Merge
    if not pivot_meth.empty and not pivot_motif.empty:
        pivot = pd.merge(pivot_meth, pivot_motif, on=INDEX_COLS, how="outer")
    elif not pivot_meth.empty:
        pivot = pivot_meth
    else:
        pivot = pivot_motif

    if meth_signal_cols:
        pivot[meth_signal_cols] = pivot[meth_signal_cols].fillna(0).astype(int)
    if motif_signal_cols:
        pivot[motif_signal_cols] = pivot[motif_signal_cols].fillna(0).astype(int)

    if meth_signal_cols:
        pivot["Total_meth"] = pivot[meth_signal_cols].sum(axis=1)
        pivot = pivot.sort_values("Total_meth", ascending=False)

    wide_file = os.path.join(args.outdir, "all_peaks_combined.csv")
    pivot.to_csv(wide_file, index=False)
    print(f"Saved wide → {wide_file}")


if __name__ == "__main__":
    main()