#!/usr/bin/env python3
"""
circular_plot.py
Generate circular methylation density + gene peak plot.
"""
import argparse
import pandas as pd
from pycirclize import Circos
import numpy as np


def parse_args():
    p = argparse.ArgumentParser()
    p.add_argument("--peaks",        required=True, help="peaks_<mod>.csv file")
    p.add_argument("--density",      required=True, help="density_<mod>.csv file")
    p.add_argument("--modification", required=True, help="Modification label e.g. 6mA")
    p.add_argument("--species",      required=True, help="Sample/species label")
    p.add_argument("--top_n",        type=int, default=30)
    p.add_argument("--output",       required=True, help="Output PNG filename")
    return p.parse_args()


def main():
    args = parse_args()

    # --- Load and clean peaks ---
    df_peaks = pd.read_csv(args.peaks)
    df_peaks = df_peaks[df_peaks["Gene"] != "intergenic"]
    df_peaks = df_peaks[pd.to_numeric(df_peaks["Gene_Start"], errors="coerce").notnull()]
    df_peaks = df_peaks[pd.to_numeric(df_peaks["Gene_End"],   errors="coerce").notnull()]
    df_peaks = df_peaks.drop_duplicates(subset="Gene").head(args.top_n)
    df_peaks["Start"] = df_peaks["Gene_Start"].astype(int)
    df_peaks["End"]   = df_peaks["Gene_End"].astype(int)

    # --- Load density ---
    df_density = pd.read_csv(args.density)
    genome_length = df_density["Position"].max() + 10000

    # --- Build Circos plot ---
    circos = Circos(sectors={"1": genome_length})
    sector = circos.sectors[0]

    # Density track
    track_density = sector.add_track((45, 90))
    track_density.axis()

    pos         = df_density["Position"].values
    density     = df_density["Density_smooth"].values
    avg_density = density.mean()

    track_density.fill_between(pos, density, y2=avg_density, color="#6A9EE8", alpha=1)
    above = density > avg_density
    track_density.fill_between(pos[above], density[above], y2=avg_density, color="#C0392B", alpha=1)

    track_density.xticks_by_interval(
        interval=500_000,
        outer=False,
        label_formatter=lambda v: "" if v == 0 else f"{v/1e6:.1f} Mb",
        label_orientation="vertical",
        line_kws=dict(ec="grey"),
        label_size=4,
    )

    # Gene track
    track_genes = sector.add_track((92, 100))
    track_genes.axis(fc="#EEEEEE", ec="none")

    for _, row in df_peaks.iterrows():
        mid = (row["Start"] + row["End"]) / 2
        track_genes.bar([mid], [1], width=5000, color="#C0392B", alpha=1, ec="none")
        track_genes.annotate(mid, row["Gene"], label_size=5)

    circos.text(f"{args.species}\n{args.modification} methylation", size=8, r=0, color="black")

    fig = circos.plotfig(figsize=(6, 6))
    fig.savefig(args.output, dpi=300, bbox_inches="tight")
    print(f"Saved: {args.output}")


if __name__ == "__main__":
    main()