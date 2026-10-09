#!/usr/bin/env python3
"""
Center-aware motif detection for DMR sites -- an alternative to streme for the
centered training sequences that annotate_dmr.py writes with --meme.

Every sequence is centered on the real modified base, so a motif only has to be
looked for at a fixed position relative to that base. The search runs in rounds:

  1. For every width (--minw..--maxw) and every offset that keeps the modified
     base inside the window, each exact k-mer seen there is used as a seed.
     Each seed is refined one position at a time into IUPAC codes (e.g. A -> W,
     or -> N), keeping a change only if it makes the motif more enriched.
  2. The most enriched motif is reported, and the sequences it explains are
     removed.
  3. Repeat on the remaining sequences until the best motif explains fewer than
     the minimum number of sites, or --max-motifs motifs have been reported. The
     minimum is --min-fraction of all input sequences (default 5%), but never
     less than --min-sites (default 5).

Enrichment is a log-likelihood ratio of the observed hit count against the count
expected by chance, where chance comes from the base composition of the sequence
flanks. Ranking by enrichment instead of raw hit count matters: a vague pattern
like CNGS matches nearly every sequence, but it is also expected to match many
by chance, so it loses to a tight motif like CCWGG.

There is no significance cut-off. Every motif is reported with the number of
sequences it explains and the number it would explain by chance. The share of
all sequences (pct_total) is the most reliable junk signal: the chance estimate
is per motif and ignores how many patterns were tried, so on shuffled sequences
junk motifs still look enriched, but each one only covers ~1% of sequences.

Usage: python detect_motifs.py --fasta meme.fasta --output motifs.tsv
"""
import argparse
import csv
import math
from collections import Counter

BASES = "ACGT"
# Single bases, two-base IUPAC codes and N. Three-base codes are left out on
# purpose: they are almost as vague as N and make motifs hard to read.
IUPAC = {
    "A": "A", "C": "C", "G": "G", "T": "T",
    "R": "AG", "Y": "CT", "S": "CG", "W": "AT", "K": "GT", "M": "AC",
    "N": "ACGT",
}
FIELDS = ["rank", "motif", "offset", "modified_base_pos", "sites", "sites_total",
          "pct_total", "expected_by_chance", "total_sequences"]


def read_fasta(path):
    seqs, cur = [], []
    with open(path) as fh:
        for line in fh:
            line = line.strip()
            if line.startswith(">"):
                if cur:
                    seqs.append("".join(cur).upper())
                cur = []
            elif line:
                cur.append(line)
    if cur:
        seqs.append("".join(cur).upper())
    return seqs


def popcount(x):
    return bin(x).count("1")


class Search:
    def __init__(self, seqs, maxw):
        self.seqs = seqs
        self.n = len(seqs)
        self.all_mask = (1 << self.n) - 1
        self.centers = [len(s) // 2 for s in seqs]
        self.center_base = Counter(s[c] for s, c in zip(seqs, self.centers)).most_common(1)[0][0]

        # base_masks[(pos, base)] = bitmask of sequences with that base at
        # position pos relative to the modified base
        self.base_masks = {}
        for pos in range(-(maxw - 1), maxw):
            for b in BASES:
                self.base_masks[(pos, b)] = 0
            for i, (s, c) in enumerate(zip(seqs, self.centers)):
                j = c + pos
                if 0 <= j < len(s) and s[j] in BASES:
                    self.base_masks[(pos, s[j])] |= 1 << i
        self.code_cache = {}

        # Background base composition, from the flanks outside the motif search
        # range (or from every non-center position if the sequences are short)
        comp = Counter()
        for s, c in zip(seqs, self.centers):
            flank = s[:max(0, c - maxw)] + s[c + maxw:]
            comp.update(b for b in flank if b in BASES)
        if sum(comp.values()) == 0:
            for s, c in zip(seqs, self.centers):
                comp.update(b for j, b in enumerate(s) if j != c and b in BASES)
        total = sum(comp.values()) or 1
        self.freq = {b: (comp[b] / total if total else 0.25) for b in BASES}

    def code_mask(self, pos, code):
        key = (pos, code)
        if key not in self.code_cache:
            if code == "N":
                m = self.all_mask
            else:
                m = 0
                for b in IUPAC[code]:
                    m |= self.base_masks[(pos, b)]
            self.code_cache[key] = m
        return self.code_cache[key]

    def matches(self, pattern, offset, active):
        m = active
        for i, code in enumerate(pattern):
            m &= self.code_mask(offset + i, code)
            if not m:
                break
        return m

    def chance_prob(self, pattern, offset):
        """Probability that a random site with the modified base at the center
        matches the pattern, from the flank base composition."""
        p = 1.0
        for i, code in enumerate(pattern):
            if offset + i != 0:
                p *= sum(self.freq[b] for b in IUPAC[code])
        return p

    def score(self, pattern, offset, active, n_active):
        k = popcount(self.matches(pattern, offset, active))
        p0 = self.chance_prob(pattern, offset)
        if k == 0 or n_active == 0 or p0 <= 0:
            return 0.0, k
        p = k / n_active
        if p <= p0:
            return 0.0, k
        llr = k * math.log(p / p0)
        if k < n_active and p0 < 1:
            llr += (n_active - k) * math.log((1 - p) / (1 - p0))
        return llr, k

    def best_motif(self, active, minw, maxw, min_seed):
        n_active = popcount(active)
        active_idx = [i for i in range(self.n) if active >> i & 1]
        best = None  # (llr, -width, pattern, offset)
        seen = set()
        for w in range(minw, maxw + 1):
            for offset in range(-(w - 1), 1):
                seeds = Counter()
                for i in active_idx:
                    s, c = self.seqs[i], self.centers[i]
                    lo, hi = c + offset, c + offset + w
                    if lo >= 0 and hi <= len(s):
                        kmer = s[lo:hi]
                        if all(b in BASES for b in kmer):
                            seeds[kmer] += 1
                for seed, cnt in seeds.items():
                    if cnt < min_seed or seed[-offset] != self.center_base:
                        continue
                    pattern = list(seed)
                    cur, _ = self.score(pattern, offset, active, n_active)
                    # hill-climb: change one position's code at a time, keep the
                    # single best improvement, repeat until nothing improves
                    while True:
                        step = None
                        for i in range(w):
                            if offset + i == 0:
                                continue  # the modified base itself stays fixed
                            for code in IUPAC:
                                if code == pattern[i]:
                                    continue
                                trial = pattern[:i] + [code] + pattern[i + 1:]
                                s_, _ = self.score(trial, offset, active, n_active)
                                if s_ > cur + 1e-9 and (step is None or s_ > step[0]):
                                    step = (s_, trial)
                        if step is None:
                            break
                        cur, pattern = step
                    # trim N padding from both ends
                    off = offset
                    while pattern and pattern[0] == "N" and off != 0:
                        pattern = pattern[1:]
                        off += 1
                    while pattern and pattern[-1] == "N" and off + len(pattern) - 1 != 0:
                        pattern = pattern[:-1]
                    key = ("".join(pattern), off)
                    if key in seen:
                        continue
                    seen.add(key)
                    cand = (cur, -len(pattern), "".join(pattern), off)
                    if best is None or cand[:2] > best[:2]:
                        best = cand
        return best


def main():
    p = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    p.add_argument("--fasta", required=True, help="Centered training sequences (annotate_dmr.py's --meme output)")
    p.add_argument("--output", required=True, help="Output TSV, one row per motif")
    p.add_argument("--minw", type=int, default=4, help="Minimum motif width (default 4)")
    p.add_argument("--maxw", type=int, default=10, help="Maximum motif width (default 10)")
    p.add_argument("--min-sites", type=int, default=5,
                   help="Absolute floor for the minimum number of sites per motif (default 5)")
    p.add_argument("--min-fraction", type=float, default=0.05,
                   help="Stop when the best motif explains fewer than this fraction of all input "
                        "sequences (default 0.05); --min-sites still applies as a floor")
    p.add_argument("--max-motifs", type=int, default=5, help="Maximum number of motifs to report (default 5)")
    args = p.parse_args()

    seqs = read_fasta(args.fasta)
    rows = []
    if seqs:
        search = Search(seqs, args.maxw)
        n = search.n
        active = search.all_mask
        min_sites = max(args.min_sites, math.ceil(args.min_fraction * n))
        # seeds need at least 2 copies so that a single sequence can't seed a motif
        min_seed = min(2, min_sites)
        while active and len(rows) < args.max_motifs:
            best = search.best_motif(active, args.minw, args.maxw, min_seed)
            if best is None:
                break
            _, _, motif, offset = best
            pattern = list(motif)
            hit = search.matches(pattern, offset, active)
            sites = popcount(hit)
            if sites < min_sites:
                break
            sites_total = popcount(search.matches(pattern, offset, search.all_mask))
            rows.append({
                "rank": len(rows) + 1,
                "motif": motif,
                "offset": offset,
                "modified_base_pos": -offset + 1,
                "sites": sites,
                "sites_total": sites_total,
                "pct_total": round(100 * sites_total / n, 1),
                "expected_by_chance": round(n * search.chance_prob(pattern, offset), 2),
                "total_sequences": n,
            })
            active &= ~hit

    with open(args.output, "w", newline="") as fh:
        writer = csv.DictWriter(fh, fieldnames=FIELDS, delimiter="\t")
        writer.writeheader()
        writer.writerows(rows)

    if not seqs:
        print("No sequences -- nothing to do.")
        return
    print(f"Motif detection -> {args.output} ({len(seqs)} sequences, modified base {search.center_base}, "
          f"minimum {min_sites} sites per motif)")
    for r in rows:
        print(f"  {r['rank']}. {r['motif']}  modified base at position {r['modified_base_pos']}: "
              f"{r['sites']} new sites, {r['sites_total']}/{r['total_sequences']} total "
              f"({r['pct_total']}%), {r['expected_by_chance']} expected by chance")
    print(f"  {popcount(active)} sequences not explained by any reported motif")


if __name__ == "__main__":
    main()
