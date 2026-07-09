#!/usr/bin/env python3
"""Pool each sample's top-N unaligned sequences, sum counts across samples,
emit the M most abundant distinct sequences as a FASTA for BLAST."""
import argparse, os, collections

def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--in", dest="inp", nargs="+", required=True, help="per-sample *.topseq.tsv")
    ap.add_argument("--n-per-sample", type=int, default=6)
    ap.add_argument("--n-total", type=int, default=20)
    ap.add_argument("--out", required=True)
    a = ap.parse_args()

    g = collections.Counter()
    for f in a.inp:
        rows = [ln.rstrip("\n").split("\t") for ln in open(f)][1:]  # skip header
        picked = 0
        for r in rows:                       # rows already sorted by count desc
            if len(r) != 3: continue
            count, _frac, seq = r
            g[seq] += int(count); picked += 1
            if picked >= a.n_per_sample: break

    os.makedirs(os.path.dirname(a.out), exist_ok=True)
    with open(a.out, "w") as o:
        for i, (seq, c) in enumerate(g.most_common(a.n_total), 1):
            o.write(f">top{i:02d}_totalcount={c}_len={len(seq)}\n{seq}\n")
    print(f"[pool] wrote {a.out} ({min(a.n_total, len(g))} sequences)")

if __name__ == "__main__":
    main()
