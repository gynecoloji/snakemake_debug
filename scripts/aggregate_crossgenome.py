#!/usr/bin/env python3
"""Step 09 cross-genome aggregation: combine per-sample cross-genome checks into
crossgenome_summary.tsv — of the sampled MULTIMAPPED reads, the % mapping to BOTH the host
and spike-in genome (and the codominant subset = equally-good hit in each genome)."""
import argparse, csv, os


def rd(p):
    if not os.path.exists(p):
        return []
    with open(p) as f:
        return list(csv.DictReader(f, delimiter="\t"))


def num(x):
    try:
        return float(x)
    except (TypeError, ValueError):
        return 0.0


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--out-summary", required=True)
    ap.add_argument("--in", dest="inp", nargs="+", required=True)
    a = ap.parse_args()
    d = os.path.dirname(a.out_summary)
    if d:
        os.makedirs(d, exist_ok=True)

    rows = []
    for f in sorted(a.inp):
        for x in rd(f):
            tot = num(x.get("n_reads")) or 1.0
            both = num(x.get("both"))
            codom = num(x.get("both_codominant"))
            rows.append({
                "sample": x.get("sample", ""),
                "n_reads": int(tot),
                "pct_both": round(both / tot * 100, 3),            # maps to BOTH genomes
                "pct_both_codominant": round(codom / tot * 100, 3),  # equally-good in each = ambiguous
                "pct_host_only": round(num(x.get("host_only")) / tot * 100, 3),
                "pct_spike_only": round(num(x.get("spike_only")) / tot * 100, 3),
                "pct_neither": round(num(x.get("neither")) / tot * 100, 3),
            })
    rows.sort(key=lambda r: -r["pct_both"])

    cols = ["sample", "n_reads", "pct_both", "pct_both_codominant",
            "pct_host_only", "pct_spike_only", "pct_neither"]
    with open(a.out_summary, "w") as o:
        o.write("\t".join(cols) + "\n")
        for r in rows:
            o.write("\t".join(str(r[c]) for c in cols) + "\n")

    print(f"[09x] wrote {a.out_summary} ({len(rows)} samples)")


if __name__ == "__main__":
    main()
