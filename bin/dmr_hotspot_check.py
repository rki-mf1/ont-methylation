#!/usr/bin/env python3
"""
Flag DMR sites that sit on a sequence variant or a basecalling-error hotspot, using the
aligned reads of both samples. Writes the full table with QC columns, the passing sites
only, and a one-line summary.
"""

import argparse
import csv
import re
import subprocess
from collections import defaultdict

INDEL = re.compile(r"[+-](\d+)")


def parse_args():
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument("--dmr", required=True, help="DMR sites table (filter_dmr.py or annotate_dmr.py output)")
    p.add_argument("--bam-a", required=True)
    p.add_argument("--bam-b", required=True)
    p.add_argument("--fasta", required=True)
    p.add_argument("--label", required=True, help="Comparison label for the summary, e.g. bc21_vs_bc22_5mC")
    p.add_argument("--out-checked", required=True)
    p.add_argument("--out-pass", required=True)
    p.add_argument("--out-summary", required=True)
    p.add_argument("--window", type=int, default=5, help="bp each side used for the error rate")
    p.add_argument("--variant-window", type=int, default=10, help="bp each side searched for variants")
    p.add_argument("--min-mapq", type=int, default=20)
    p.add_argument("--min-depth", type=int, default=10)
    p.add_argument("--variant-frac", type=float, default=0.3)
    p.add_argument("--max-error", type=float, default=0.15)
    p.add_argument("--max-error-diff", type=float, default=0.05)
    return p.parse_args()


def parse_bases(s):
    """Per-strand counts from one mpileup bases string."""
    c = {"+": defaultdict(int), "-": defaultdict(int)}
    i = 0
    while i < len(s):
        ch = s[i]
        if ch == "^":
            i += 2
            continue
        if ch == "$":
            i += 1
            continue
        if ch in "+-":
            m = INDEL.match(s, i)
            n = int(m.group(1))
            seq = s[m.end():m.end() + n]
            if ch == "+":
                c["+" if seq[:1].isupper() else "-"]["ins"] += 1
            i = m.end() + n
            continue
        if ch == ".":
            c["+"]["match"] += 1
        elif ch == ",":
            c["-"]["match"] += 1
        elif ch == "*":
            c["+"]["del"] += 1
        elif ch == "#":
            c["-"]["del"] += 1
        elif ch in "ACGTN":
            c["+"][ch] += 1
        elif ch in "acgtn":
            c["-"][ch.upper()] += 1
        i += 1
    return c


def pileup(bam, fasta, regions_bed, min_mapq):
    cmd = ["samtools", "mpileup", "-B", "-Q", "0", "-d", "0", "--reverse-del",
           "-q", str(min_mapq), "-f", fasta, "-l", regions_bed, bam]
    out = {}
    proc = subprocess.run(cmd, capture_output=True, text=True, check=True)
    for line in proc.stdout.splitlines():
        f = line.split("\t")
        out[(f[0], int(f[1]) - 1)] = (f[2].upper(), parse_bases(f[4]))
    return out


def site_metrics(pile, chrom, pos, strand, window, variant_window, min_depth):
    errs, alt_best, alt_desc = [], 0.0, ""
    for p in range(pos - variant_window, pos + variant_window + 1):
        if (chrom, p) not in pile:
            continue
        ref, c = pile[(chrom, p)]
        s = c[strand]
        depth = s["match"] + s["del"] + sum(s[b] for b in "ACGTN")
        if abs(p - pos) <= window and depth >= min_depth:
            mism = sum(s[b] for b in "ACGTN" if b != ref)
            errs.append((mism + s["del"] + s["ins"]) / depth)
        both = {k: c["+"][k] + c["-"][k] for k in ["match", "del", "A", "C", "G", "T", "N"]}
        total = sum(both.values())
        if total >= min_depth:
            for b in "ACGT":
                if b != ref and both[b] / total > alt_best:
                    alt_best, alt_desc = both[b] / total, f"{p + 1}:{ref}>{b}"
            if both["del"] / total > alt_best:
                alt_best, alt_desc = both["del"] / total, f"{p + 1}:{ref}>del"
    err = sum(errs) / len(errs) if errs else float("nan")
    return err, alt_best, alt_desc


def main():
    args = parse_args()
    with open(args.dmr) as fh:
        reader = csv.DictReader(fh, delimiter="\t")
        fields = reader.fieldnames or []
        rows = list(reader)

    sites = sorted({(r["chrom"], int(r["start"]), r["strand"]) for r in rows})
    regions = "regions.bed"
    with open(regions, "w") as fh:
        for chrom, pos, _ in sites:
            w = max(args.window, args.variant_window)
            fh.write(f"{chrom}\t{max(0, pos - w)}\t{pos + w + 1}\n")

    piles = {}
    if sites:
        piles = {"a": pileup(args.bam_a, args.fasta, regions, args.min_mapq),
                 "b": pileup(args.bam_b, args.fasta, regions, args.min_mapq)}

    qc = {}
    for chrom, pos, strand in sites:
        m = {k: site_metrics(piles[k], chrom, pos, strand, args.window, args.variant_window, args.min_depth) for k in "ab"}
        reasons = []
        variant = [f"{k}:{m[k][2]}({m[k][1]:.0%})" for k in "ab" if m[k][1] > args.variant_frac]
        if variant:
            reasons.append("variant")
        ea, eb = m["a"][0], m["b"][0]
        if (ea == ea and eb == eb) and (max(ea, eb) > args.max_error or abs(ea - eb) > args.max_error_diff):
            reasons.append("error_hotspot")
        qc[(chrom, pos, strand)] = {
            "local_error_a": f"{ea:.3f}", "local_error_b": f"{eb:.3f}",
            "max_alt_frac_a": f"{m['a'][1]:.3f}", "max_alt_frac_b": f"{m['b'][1]:.3f}",
            "variant_detail": ";".join(variant),
            "qc_status": "+".join(reasons) if reasons else "pass",
        }

    out_fields = fields + ["local_error_a", "local_error_b", "max_alt_frac_a", "max_alt_frac_b",
                           "variant_detail", "qc_status"]
    with open(args.out_checked, "w", newline="") as fc, open(args.out_pass, "w", newline="") as fp:
        wc = csv.DictWriter(fc, out_fields, delimiter="\t")
        wp = csv.DictWriter(fp, out_fields, delimiter="\t")
        wc.writeheader()
        wp.writeheader()
        for r in rows:
            r.update(qc[(r["chrom"], int(r["start"]), r["strand"])])
            wc.writerow(r)
            if r["qc_status"] == "pass":
                wp.writerow(r)

    status = [q["qc_status"] for q in qc.values()]
    with open(args.out_summary, "w") as fh:
        fh.write("comparison\tsites\tpass\tvariant\terror_hotspot\tvariant_and_error_hotspot\n")
        fh.write(f"{args.label}\t{len(status)}\t{status.count('pass')}\t{status.count('variant')}\t"
                 f"{status.count('error_hotspot')}\t{status.count('variant+error_hotspot')}\n")
    print(f"{args.label}: {len(status)} sites, {status.count('pass')} pass")


if __name__ == "__main__":
    main()
