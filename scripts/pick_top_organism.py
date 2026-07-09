#!/usr/bin/env python3
"""Pick the single most-supported NON-HOST organism from the analysis, combining
Kraken2 (top_species_overall.tsv) and random-BLAST (random_blast_species.tsv).
A credibility threshold skips low-count noise (e.g. Kraken2's stray Borrelia in a
sample whose real contaminant, Drosophila, is only visible via BLAST)."""
import argparse, csv, os

# substrings that are not fetchable single organisms / are common false positives
NOISE = ["synthetic", "vector", "uncultured", "unclassified", "construct",
         "orfeome", "plasmid", " clone", "mag:", "environmental", "unknown"]


def _rd(p):
    if not os.path.exists(p):
        return []
    with open(p) as f:
        return list(csv.DictReader(f, delimiter="\t"))


def _num(x):
    try:
        return float(x)
    except (TypeError, ValueError):
        return 0.0


def pick_organism(kraken_tsv, blast_tsv, host_label, min_frac=0.05):
    """Return (name, score, source) for the top non-host organism, or (None, 0, '')."""
    host = (host_label or "host").lower()

    def ok(name):
        n = name.lower()
        if host in n or "human" in n or n.startswith("homo "):
            return False
        return not any(tok in n for tok in NOISE) and len(n.split()) >= 2

    cand, src = {}, {}
    # Kraken2: fraction of total classified reads across samples
    kr = _rd(kraken_tsv)
    ktot = sum(_num(r.get("total_reads_across_samples")) for r in kr) or 1.0
    for r in kr:
        sp = r.get("species", "")
        if ok(sp):
            frac = _num(r.get("total_reads_across_samples")) / ktot
            if frac > cand.get(sp, 0):
                cand[sp], src[sp] = frac, "kraken2"
    # random-BLAST: fraction of hits
    for r in _rd(blast_tsv):
        sp = r.get("species", "")
        if ok(sp):
            frac = _num(r.get("pct_of_hits")) / 100.0
            if frac > cand.get(sp, 0):
                cand[sp], src[sp] = frac, ("blast" if src.get(sp) is None else "kraken2+blast")
    for name, score in sorted(cand.items(), key=lambda kv: -kv[1]):
        if score >= min_frac:
            return name, score, src[name]
    return None, 0.0, ""


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--kraken", required=True)
    ap.add_argument("--blast", default="")
    ap.add_argument("--host", default="host")
    ap.add_argument("--min-fraction", type=float, default=0.05)
    a = ap.parse_args()
    name, score, source = pick_organism(a.kraken, a.blast, a.host, a.min_fraction)
    if name:
        print(f"{name}\t{score:.3f}\t{source}")
    else:
        print("NONE\t0\t")


if __name__ == "__main__":
    main()
