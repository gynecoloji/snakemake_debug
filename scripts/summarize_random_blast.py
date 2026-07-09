#!/usr/bin/env python3
"""Summarize remote-BLAST hits of the random unaligned reads into a species table:
best hit per read -> 'Genus species' -> read count. Surfaces organisms that Kraken2
(DB-limited) and the motif screen miss. Reads a blastn -outfmt 6 file with columns
qseqid pident length evalue bitscore stitle."""
import argparse, collections, os, re


def species_of(stitle):
    """Crude 'Genus species' from a subject title; drop PREDICTED:/UNVERIFIED: prefixes."""
    t = re.sub(r"^(PREDICTED|UNVERIFIED|TPA[_:]?\w*):\s*", "", stitle.strip())
    toks = t.split()
    if len(toks) >= 2 and toks[0][:1].isalpha():
        return f"{toks[0]} {toks[1]}"
    return toks[0] if toks else "unknown"


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--raw", required=True)
    ap.add_argument("--out-species", required=True)
    ap.add_argument("--out-hits", required=True)
    a = ap.parse_args()

    best = {}  # qseqid -> (bitscore, pident, length, evalue, stitle)
    if os.path.exists(a.raw):
        with open(a.raw) as f:
            for line in f:
                p = line.rstrip("\n").split("\t")
                if len(p) < 6:
                    continue
                q, pid, ln, ev, bit, stitle = p[0], p[1], p[2], p[3], p[4], p[5]
                try:
                    b = float(bit)
                except ValueError:
                    continue
                if q not in best or b > best[q][0]:
                    best[q] = (b, pid, ln, ev, stitle)

    spc = collections.Counter()
    rows = []
    for q, (b, pid, ln, ev, stitle) in best.items():
        sp = species_of(stitle)
        spc[sp] += 1
        rows.append((q, pid, ln, ev, f"{b:.0f}", sp, stitle))

    for path in (a.out_species, a.out_hits):
        d = os.path.dirname(path)
        if d:
            os.makedirs(d, exist_ok=True)

    with open(a.out_hits, "w") as o:
        o.write("query\tpident\tlength\tevalue\tbitscore\tspecies\tsubject\n")
        for r in sorted(rows, key=lambda x: x[5]):
            o.write("\t".join(r) + "\n")

    total = max(sum(spc.values()), 1)
    with open(a.out_species, "w") as o:
        o.write("species\treads\tpct_of_hits\n")
        for sp, c in spc.most_common():
            o.write(f"{sp}\t{c}\t{c/total*100:.1f}\n")

    print(f"[04] {len(best)} random reads with hits -> {len(spc)} species -> {a.out_species}")


if __name__ == "__main__":
    main()
