#!/usr/bin/env python3
"""Step 08: aggregate per-sample `samtools idxstats` of unaligned reads mapped to a
custom FASTA into (a) per-sample overall % mapping and (b) a sample x sequence
read-count matrix (which custom sequence each read hits)."""
import argparse, os


def sample_name(path):
    b = os.path.basename(path)
    return b[:-len(".idxstats")] if b.endswith(".idxstats") else os.path.splitext(b)[0]


def parse_idxstats(path):
    """Return (ordered [(seqname, mapped)], unmapped_count)."""
    seqs, unmapped = [], 0
    with open(path) as f:
        for line in f:
            p = line.rstrip("\n").split("\t")
            if len(p) < 4:
                continue
            name, _length, mapped, unmap = p[0], p[1], int(p[2]), int(p[3])
            if name == "*":
                unmapped += unmap
            else:
                seqs.append((name, mapped))
                unmapped += unmap  # mate-unmapped segments assigned to a ref (SE: ~0)
    return seqs, unmapped


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--out-mapping", required=True)
    ap.add_argument("--out-perseq", required=True)
    ap.add_argument("--in", dest="inp", nargs="+", required=True)
    a = ap.parse_args()
    for p in (a.out_mapping, a.out_perseq):
        d = os.path.dirname(p)
        if d:
            os.makedirs(d, exist_ok=True)

    rows, seq_order = [], []
    for path in sorted(a.inp):
        s = sample_name(path)
        seqs, unmapped = parse_idxstats(path)
        if not seq_order:
            seq_order = [name for name, _ in seqs]
        mapped_by = dict(seqs)
        mapped = sum(mapped_by.values())
        total = mapped + unmapped
        rows.append({"sample": s, "total": total, "mapped": mapped,
                     "pct": (mapped / total * 100) if total else 0.0,
                     "by_seq": mapped_by})

    rows.sort(key=lambda r: -r["pct"])

    with open(a.out_mapping, "w") as o:
        o.write("sample\ttotal_unaligned_reads\treads_mapping_custom\tpct_mapping_custom\n")
        for r in rows:
            o.write(f"{r['sample']}\t{r['total']}\t{r['mapped']}\t{r['pct']:.4f}\n")

    with open(a.out_perseq, "w") as o:
        o.write("sample\t" + "\t".join(seq_order) + "\n")
        for r in rows:
            o.write(r["sample"] + "\t" + "\t".join(str(r["by_seq"].get(name, 0)) for name in seq_order) + "\n")

    print(f"[07] wrote {a.out_mapping} + {a.out_perseq} "
          f"({len(rows)} samples, {len(seq_order)} custom sequences)")


if __name__ == "__main__":
    main()
