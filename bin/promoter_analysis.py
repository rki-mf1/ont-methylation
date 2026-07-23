#!/usr/bin/env python3

import argparse
import pandas as pd
import numpy as np
from Bio import SeqIO
from Bio.Seq import Seq
from Bio.SeqUtils import nt_search
import os
from utils import parse_gff3_gene_names, read_modkit

def parse_args():
    parser = argparse.ArgumentParser(description="Promoter methylation analysis")
    parser.add_argument("--fasta",                       required=True)
    parser.add_argument("--modkit_bed",                  required=True)
    parser.add_argument("--gff3",                        required=True)
    parser.add_argument("--motif_file",                  required=True)
    parser.add_argument("--sample",                      required=True)
    parser.add_argument("--outdir",                      required=True)
    parser.add_argument("--promoter_window",             type=int,   default=250)
    parser.add_argument("--percent_modified_promoter",  type=float, default=0.3)
    return parser.parse_args()


def deduplicate_motifs(motif_df):
    seen = set()
    kept = []
    for _, row in motif_df.iterrows():
        motif = row["motif"].upper()
        if motif in seen:
            continue
        rc = str(Seq(motif).reverse_complement()).upper()
        seen.add(motif)
        seen.add(rc)
        kept.append({
            "motif":    motif,
            "mod_code": row["mod_code"],
            "rc_motif": rc,
        })
    return pd.DataFrame(kept)

def find_motif_positions(motif, genome_seq):
    return np.array(nt_search(str(genome_seq), motif)[1:])

def get_promoter_region(start, end, strand, window, genome_length):
    if strand == "+":
        return max(0, start - window), start
    else:
        return end, min(genome_length, end + window)

def is_motif_methylated(motif_start, motif_len, methylated_set):
    return any(p in methylated_set for p in range(motif_start, motif_start + motif_len))

def main():
    args = parse_args()
    os.makedirs(args.outdir, exist_ok=True)

    reference_genome = {record.id: record.seq for record in SeqIO.parse(args.fasta, "fasta")}
    largest_contig   = max(reference_genome, key=lambda k: len(reference_genome[k]))
    genome_length    = len(reference_genome[largest_contig])
    genome_seq       = reference_genome[largest_contig]
    print(f"Using contig: {largest_contig}")

    annotation = parse_gff3_gene_names(args.gff3)
    annotation = [g for g in annotation if g[0] == largest_contig]

    motif_df       = pd.read_csv(args.motif_file, sep=r'\s+')
    motif_df_dedup = deduplicate_motifs(motif_df)
    print(f"Kept {len(motif_df_dedup)} motifs after deduplication")

    motif_positions = {}
    for _, row in motif_df_dedup.iterrows():
        positions = find_motif_positions(row["motif"], genome_seq)
        print(f"  {row['motif']}: {len(positions)} occurrences")
        motif_positions[row["motif"]] = positions

    mod_code_map  = motif_df_dedup.set_index("motif")["mod_code"].to_dict()
    motif_len_map = {row["motif"]: len(row["motif"]) for _, row in motif_df_dedup.iterrows()}

    all_positions = {}
    for mod_code in motif_df_dedup["mod_code"].unique():
        print(f"\nLoading modkit data for mod_code: {mod_code}")
        mod = read_modkit(args.modkit_bed, 0, mod_code)
        mod["Contig"] = mod["Contig"].astype(str).str.strip()
        mod = mod[mod["Contig"] == largest_contig].copy()
        all_positions[mod_code] = (
            mod.groupby("Position")["Percent_modified"]
            .max() # here taking the highest value if both strands have methylated bases
            .to_dict()
        )
        print(f"  Total positions loaded: {len(all_positions[mod_code])}")
    methylated_sets = {
        mod_code: set(pos for pos, pct in pos_dict.items() if pct >= args.percent_modified_promoter)
        for mod_code, pos_dict in all_positions.items()
    }

    rows = []
    for contig, start, end, strand, gene in annotation:
        prom_start, prom_end = get_promoter_region(
            start, end, strand, args.promoter_window, genome_length
        )

        row = {
            "Gene":       gene,
            "Gene_Start": start,
            "Gene_End":   end,
            "Strand":     strand,
            "Prom_Start": prom_start,
            "Prom_End":   prom_end,
        }

        total_methylated   = 0
        total_unmethylated = 0
        total_hits         = 0

        for motif, positions in motif_positions.items():
            mod_code  = mod_code_map[motif]
            motif_len = motif_len_map[motif]
            meth_set  = methylated_sets[mod_code]

            window_hits = positions[(positions >= prom_start) & (positions < prom_end)]

            n_methylated   = 0
            n_unmethylated = 0
            n_not_covered  = 0

            for hit_pos in window_hits:
                motif_seq = str(genome_seq[hit_pos:hit_pos + motif_len])
                print(f"\n  Motif: {motif} | genomic seq: {motif_seq} | start: {hit_pos} | end: {hit_pos + motif_len}")
                for p in range(hit_pos, hit_pos + motif_len):
                    in_all   = p in all_positions[mod_code]
                    in_meth  = p in methylated_sets[mod_code]
                    pct      = all_positions[mod_code].get(p, None)
                    print(f"    pos {p}: covered={in_all} | methylated={in_meth} | percent_modified={pct}")


                covered = any(
                    p in all_positions[mod_code]
                    for p in range(hit_pos, hit_pos + motif_len)
                )
                if not covered:
                    n_not_covered += 1
                elif is_motif_methylated(hit_pos, motif_len, meth_set):
                    n_methylated += 1
                else:
                    n_unmethylated += 1

            row[f"{motif}_hits"]         = len(window_hits)
            row[f"{motif}_methylated"]   = n_methylated
            row[f"{motif}_unmethylated"] = n_unmethylated
            row[f"{motif}_not_covered"]  = n_not_covered

            total_methylated   += n_methylated
            total_unmethylated += n_unmethylated
            total_hits         += len(window_hits)

        row["Total_hits"]         = total_hits
        row["Total_methylated"]   = total_methylated
        row["Total_unmethylated"] = total_unmethylated

        rows.append(row)

    df_promoter = pd.DataFrame(rows)
    df_promoter = df_promoter[df_promoter["Total_hits"] > 0]

    for motif in motif_positions.keys():
        motif_cols = [
            "Gene", "Gene_Start", "Gene_End", "Strand",
            "Prom_Start", "Prom_End",
            f"{motif}_hits",
            f"{motif}_methylated",
            f"{motif}_unmethylated",
            f"{motif}_not_covered",
        ]

        df_motif = (
            df_promoter[df_promoter[f"{motif}_hits"] > 0][motif_cols]
            .copy()
            .sort_values([f"{motif}_unmethylated", f"{motif}_hits"], ascending=[False, False])
        )

        safe_motif = motif.replace("/", "_")
        out_file   = os.path.join(args.outdir, f"promoter_{safe_motif}.csv")
        df_motif.to_csv(out_file, index=False)
        print(f"Saved → {out_file}")

        print(f"\n===== {motif} — unmethylated in promoter =====")
        unmeth = df_motif[df_motif[f"{motif}_unmethylated"] > 0]
        print(f"Genes with unmethylated {motif}: {len(unmeth)}")
        print(unmeth.to_string(index=False))

        print(f"\n===== {motif} — methylated only in promoter =====")
        meth_only = df_motif[df_motif[f"{motif}_unmethylated"] == 0]
        print(f"Genes with only methylated {motif}: {len(meth_only)}")
        print(meth_only.to_string(index=False))


if __name__ == "__main__":
    main()