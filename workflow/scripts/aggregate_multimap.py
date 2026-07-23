#!/usr/bin/env python3
"""Step 09 aggregation: combine per-sample multimapped-read summaries into
multimap_summary.tsv (per-sample % host / spike-in / multimapped + spike-in
fraction), and pool the per-sample multimap loci into top_multimap_loci.tsv
(which contigs the multimappers concentrate at — the 'why' for non-spike-in runs)."""
import argparse, csv, collections, os


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
    ap.add_argument("--out-loci", required=True)
    ap.add_argument("--in", dest="inp", nargs="+", required=True)          # per-sample .multimap.tsv
    ap.add_argument("--in-loci", dest="loci", nargs="*", default=[])       # per-sample .loci.tsv
    a = ap.parse_args()
    for p in (a.out_summary, a.out_loci):
        d = os.path.dirname(p)
        if d:
            os.makedirs(d, exist_ok=True)

    rows = []
    for f in sorted(a.inp):
        for x in rd(f):
            total = num(x.get("total_mapped")) or 1.0
            mm = num(x.get("multimapped"))
            host = num(x.get("host_reads"))
            spk = num(x.get("spikein_reads"))
            host_mm = num(x.get("host_multimapped"))
            spk_mm = num(x.get("spikein_multimapped"))
            geno = host + spk
            rows.append({
                "sample": x.get("sample", ""),
                "method": x.get("method", ""),
                "total_mapped": int(total),
                "pct_multimapped": round(mm / total * 100, 3),   # of all mapped reads
                "pct_host": round(host / geno * 100, 3) if geno else 0.0,      # host vs spike-in
                "pct_spikein": round(spk / geno * 100, 3) if geno else 0.0,    # genome split (=spike-in fraction)
                # multimapped DISTRIBUTION: fraction of EACH genome's reads that multimap
                "pct_host_multimapped": round(host_mm / host * 100, 3) if host else 0.0,
                "pct_spikein_multimapped": round(spk_mm / spk * 100, 3) if spk else 0.0,
            })
    rows.sort(key=lambda r: -r["pct_multimapped"])

    with open(a.out_summary, "w") as o:
        cols = ["sample", "method", "total_mapped", "pct_multimapped", "pct_host", "pct_spikein",
                "pct_host_multimapped", "pct_spikein_multimapped"]
        o.write("\t".join(cols) + "\n")
        for r in rows:
            o.write("\t".join(str(r[c]) for c in cols) + "\n")

    loci = collections.Counter()
    for f in sorted(a.loci):
        for x in rd(f):
            loci[x.get("contig", "")] += int(num(x.get("multimapped_reads")))
    tot = max(sum(loci.values()), 1)
    with open(a.out_loci, "w") as o:
        o.write("contig\tmultimapped_reads\tpct_of_multimapped\n")
        for contig, n in loci.most_common(25):
            o.write(f"{contig}\t{n}\t{n/tot*100:.2f}\n")

    print(f"[09] wrote {a.out_summary} ({len(rows)} samples) + {a.out_loci} ({len(loci)} contigs)")


if __name__ == "__main__":
    main()
