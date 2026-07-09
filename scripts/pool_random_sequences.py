#!/usr/bin/env python3
"""Reservoir-sample N random reads per sample from the unaligned FASTQs and pool
them into one FASTA. Unlike pool_top_sequences (most-duplicated), a random sample
is unbiased, so a comprehensive remote BLAST of it surfaces organisms that Kraken2
(DB-limited) and the fixed motif screen both miss. Seeded -> reproducible."""
import argparse, gzip, os, random


def sample_fastq(path, n, rng):
    """Reservoir-sample up to n sequences from a fastq(.gz)."""
    op = gzip.open if path.endswith(".gz") else open
    res, i = [], 0
    with op(path, "rt") as f:
        while True:
            h = f.readline()
            if not h:
                break
            seq = f.readline().strip(); f.readline(); f.readline()
            if not seq or seq == "*":
                continue
            if len(res) < n:
                res.append(seq)
            else:
                j = rng.randint(0, i)
                if j < n:
                    res[j] = seq
            i += 1
    return res


def sample_name(path, suffix=".unaligned.fq.gz"):
    b = os.path.basename(path)
    return b[:-len(suffix)] if b.endswith(suffix) else os.path.splitext(b)[0]


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--n-per-sample", type=int, default=25)
    ap.add_argument("--seed", type=int, default=100)
    ap.add_argument("--out", required=True)
    ap.add_argument("--in", dest="inp", nargs="+", required=True)
    a = ap.parse_args()
    os.makedirs(os.path.dirname(a.out), exist_ok=True)
    rng = random.Random(a.seed)
    n = 0
    with open(a.out, "w") as o:
        for path in sorted(a.inp):
            s = sample_name(path)
            for i, seq in enumerate(sample_fastq(path, a.n_per_sample, rng)):
                o.write(f">random_{s}_{i}\n{seq}\n")
                n += 1
    print(f"[04] wrote {n} random unaligned reads -> {a.out}")


if __name__ == "__main__":
    main()
