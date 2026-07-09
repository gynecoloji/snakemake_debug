#!/usr/bin/env python3
"""Step 02+03 for ONE sample: sample unaligned reads -> FASTQ, top-sequence
tally, and motif signature counts. samtools + stdlib only."""
import argparse, subprocess, os, re, gzip, collections

def rc(s):
    return s.translate(str.maketrans("ACGTN", "TGCAN"))[::-1]

# --- motif signatures (edit/extend here to add rRNA, vector, mito, ...) -----
BACT_16S = ["GGGCGTAAAGC", "TTACGCCC"]          # conserved bacterial 16S rRNA
ADAPTER  = ["AGATCGGAAGAGC"]                     # Illumina TruSeq adapter start
HOST_REP = ["GGCTGAGGCAGG", "TGTAATCCCAGC",      # human Alu consensus segments
            "GGCCGGGCGCGGTGG", "GAGGCTGAGGCAGGAGAA"]
def has(seq, motifs): return any(m in seq or rc(m) in seq for m in motifs)
def is_poly(seq):     return bool(re.search(r"A{20,}|T{20,}", seq))

def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--align", required=True)
    ap.add_argument("--sample", required=True)
    ap.add_argument("--flag", type=int, default=4)
    ap.add_argument("--unmapped-cap", type=int, default=300000)
    ap.add_argument("--records-cap", type=int, default=3000000)
    ap.add_argument("--out-fastq", required=True)
    ap.add_argument("--out-topseq", required=True)
    ap.add_argument("--out-sig", required=True)
    a = ap.parse_args()
    for p in (a.out_fastq, a.out_topseq, a.out_sig):
        os.makedirs(os.path.dirname(p), exist_ok=True)

    p = subprocess.Popen(["samtools", "view", a.align], stdout=subprocess.PIPE, text=True, bufsize=1 << 20)
    seqcount = collections.Counter()
    c = dict(bact_16S=0, illumina_adapter=0, polyA=0, host_repeat=0, other=0)
    scanned = kept = 0
    fq = gzip.open(a.out_fastq, "wt")
    try:
        for line in p.stdout:
            scanned += 1
            if scanned > a.records_cap: break
            f = line.split("\t", 11)
            if len(f) < 11: continue
            if int(f[1]) & a.flag:
                name, seq, qual = f[0], f[9], f[10]
                if seq == "*" or not seq: continue
                fq.write(f"@{name}\n{seq}\n+\n{qual}\n")
                seqcount[seq] += 1
                if   has(seq, BACT_16S): c["bact_16S"] += 1
                elif has(seq, ADAPTER):  c["illumina_adapter"] += 1
                elif is_poly(seq):       c["polyA"] += 1
                elif has(seq, HOST_REP): c["host_repeat"] += 1
                else:                    c["other"] += 1
                kept += 1
                if kept >= a.unmapped_cap: break
    finally:
        fq.close(); p.stdout.close()
        try: p.terminate()
        except Exception: pass

    with open(a.out_topseq, "w") as o:
        o.write("count\tfrac_of_sampled\tsequence\n")
        for seq, n in seqcount.most_common(30):
            o.write(f"{n}\t{n/max(kept,1):.4f}\t{seq}\n")

    tot = max(kept, 1)
    with open(a.out_sig, "w") as o:   # single data row (no header; aggregator adds it)
        o.write(f"{a.sample}\t{kept}\t{scanned}\t{c['bact_16S']/tot:.4f}\t"
                f"{c['illumina_adapter']/tot:.4f}\t{c['polyA']/tot:.4f}\t"
                f"{c['host_repeat']/tot:.4f}\t{c['other']/tot:.4f}\n")
    print(f"[extract] {a.sample}: kept {kept} unaligned reads (scanned {scanned})")

if __name__ == "__main__":
    main()
