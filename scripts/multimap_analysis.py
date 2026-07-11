#!/usr/bin/env python3
"""Step 09 for ONE sample: investigate MULTIMAPPED reads in the raw aligner output.

Multimapper detection is aligner-aware (auto):
  - NH tag present  -> NH:i > 1              (HISAT2 / STAR, RNA-seq)
  - else XS tag     -> XS:i present           (raw bowtie2; --xs-strict: require XS == AS, an
                                               equally-best 2nd hit). XS-present reproduces
                                               bowtie2's own ">1 time" alignment rate.
  - else            -> MAPQ <= --mapq-max     (processed BAM with tags stripped: MAPQ survives)

Two modes:
  * non-spike-in: scan primary alignments (capped by --records-cap) -> % multimapped +
    the top contigs where multimappers concentrate (explains WHY they multimap).
  * spike-in: reads are aligned to a combined host+spike-in genome (spike-in contigs carry
    --spikein-prefix). Reports, over ALL reads:
      - host vs spike-in genome split (the normalization signal),
      - overall % multimapped (aligner-aware, so bowtie2 XS matches its own >1-time rate),
      - the multimapped DISTRIBUTION: what % of host reads and of spike-in reads multimap,
        and the top contigs the multimappers land on.
    host/spike-in totals come from `samtools idxstats` (exact over all reads, index-based, so
    spike-in contigs at the end of a coordinate-sorted BAM are never missed); the multimapped
    subset is tallied by streaming only the flagged reads (fast). A SAM / unindexed input
    falls back to a single full scan.
samtools + stdlib only."""
import argparse, subprocess, os, re, collections

NH = re.compile(r"\tNH:i:(\d+)")
AS = re.compile(r"\tAS:i:(-?\d+)")
XS = re.compile(r"\tXS:i:(-?\d+)")


def popen(cmd):
    return subprocess.Popen(cmd, stdout=subprocess.PIPE, stderr=subprocess.DEVNULL,
                            text=True, bufsize=1 << 20)


def run(cmd):
    return subprocess.run(cmd, stdout=subprocess.PIPE, stderr=subprocess.DEVNULL, text=True).stdout


def detect_method(align, forced):
    """Peek the head of the primary alignments to choose the detector (nh/xs/mapq)."""
    if forced in ("nh", "xs", "mapq"):
        return forced
    p = popen(["samtools", "view", "-F", "0x900", align])
    has_nh = has_xs = False
    for _ in range(100000):
        ln = p.stdout.readline()
        if not ln:
            break
        if not has_nh and "\tNH:i:" in ln:
            has_nh = True
        if not has_xs and "\tXS:i:" in ln:
            has_xs = True
        if has_nh:
            break
    p.stdout.close()
    try:
        p.terminate()
    except Exception:
        pass
    return "nh" if has_nh else ("xs" if has_xs else "mapq")


def mm_expr(method, mapq_max, xs_strict):
    """samtools filter expression selecting the MULTIMAPPED reads for this method."""
    if method == "nh":
        return "[NH]>1"
    if method == "xs":
        return "[XS]==[AS]" if xs_strict else "[XS]"
    return f"mapq<={mapq_max}"


def idxstats_split(align, prefix):
    """(host, spikein) mapped-read totals over ALL reads; (0,0) if idxstats fails / unindexed."""
    host = spk = 0
    for ln in run(["samtools", "idxstats", align]).splitlines():
        f = ln.split("\t")
        if len(f) < 3 or f[0] == "*":
            continue
        m = int(f[2]) if f[2].isdigit() else 0
        if f[0].startswith(prefix):
            spk += m
        else:
            host += m
    return host, spk


def stream_mm(align, expr, prefix):
    """Stream only the multimapped primary reads; tally by genome + per-contig."""
    p = popen(["samtools", "view", "-@", "4", "-F", "0x904", "-e", expr, align])
    host_mm = spk_mm = 0
    loci = collections.Counter()
    for ln in p.stdout:
        rname = ln.split("\t", 3)[2]
        if rname.startswith(prefix):
            spk_mm += 1
        else:
            host_mm += 1
        loci[rname] += 1
    p.stdout.close()
    p.wait()
    return host_mm, spk_mm, loci


def scan(align, method, mapq_max, xs_strict, prefix, cap):
    """Full per-read scan (non-spike-in / unindexed fallback). cap=0 -> no cap.
    Counts genome totals AND the multimapped subset, split by genome."""
    p = popen(["samtools", "view", "-@", "4", "-F", "0x904", align])
    total = mm = host = spk = host_mm = spk_mm = 0
    loci = collections.Counter()
    for ln in p.stdout:
        if cap and total >= cap:
            break
        f = ln.split("\t", 6)
        rname, mapq = f[2], int(f[4])
        if method == "nh":
            m = NH.search(ln); is_mm = (int(m.group(1)) > 1) if m else False
        elif method == "xs":
            if xs_strict:
                mx, ma = XS.search(ln), AS.search(ln)
                is_mm = bool(mx and ma and int(mx.group(1)) == int(ma.group(1)))
            else:
                is_mm = "\tXS:i:" in ln
        else:
            is_mm = mapq <= mapq_max
        total += 1
        spike_read = rname.startswith(prefix)
        if spike_read:
            spk += 1
        else:
            host += 1
        if is_mm:
            mm += 1; loci[rname] += 1
            if spike_read:
                spk_mm += 1
            else:
                host_mm += 1
    p.stdout.close()
    try:
        p.terminate()
    except Exception:
        pass
    return host, spk, mm, total, host_mm, spk_mm, loci


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--align", required=True)
    ap.add_argument("--sample", required=True)
    ap.add_argument("--spikein", action="store_true")
    ap.add_argument("--spikein-prefix", default="spikein_")
    ap.add_argument("--method", default="auto")
    ap.add_argument("--xs-strict", action="store_true")
    ap.add_argument("--mapq-max", type=int, default=10)
    ap.add_argument("--records-cap", type=int, default=3000000)
    ap.add_argument("--out", required=True)
    ap.add_argument("--out-loci", required=True)
    a = ap.parse_args()
    for p in (a.out, a.out_loci):
        os.makedirs(os.path.dirname(p), exist_ok=True)

    if a.spikein:
        host, spk = idxstats_split(a.align, a.spikein_prefix)
        method = detect_method(a.align, a.method)
        if host + spk > 0:                       # indexed BAM: exact totals + streamed mm subset
            expr = mm_expr(method, a.mapq_max, a.xs_strict)
            host_mm, spk_mm, loci = stream_mm(a.align, expr, a.spikein_prefix)
            total, mm = host + spk, host_mm + spk_mm
            method = f"{method}(idxstats)"
        else:                                    # SAM / unindexed: one full scan (no cap)
            host, spk, mm, total, host_mm, spk_mm, loci = scan(
                a.align, method, a.mapq_max, a.xs_strict, a.spikein_prefix, 0)
    else:
        method = detect_method(a.align, a.method)
        host, spk, mm, total, host_mm, spk_mm, loci = scan(
            a.align, method, a.mapq_max, a.xs_strict, a.spikein_prefix, a.records_cap)

    with open(a.out, "w") as o:
        o.write("sample\tmethod\ttotal_mapped\tmultimapped\thost_reads\tspikein_reads\t"
                "host_multimapped\tspikein_multimapped\n")
        o.write(f"{a.sample}\t{method}\t{total}\t{mm}\t{host}\t{spk}\t{host_mm}\t{spk_mm}\n")
    with open(a.out_loci, "w") as o:
        o.write("sample\tcontig\tmultimapped_reads\tpct_of_multimapped\n")
        for contig, n in loci.most_common(20):
            o.write(f"{a.sample}\t{contig}\t{n}\t{n/max(mm,1)*100:.2f}\n")

    msg = (f"[09] {a.sample}: method={method} total={total} multimapped={mm} "
           f"({mm/max(total,1)*100:.2f}%)")
    if a.spikein:
        msg += (f" | host={host} (mm {host_mm}={host_mm/max(host,1)*100:.2f}%)"
                f" spikein={spk} (mm {spk_mm}={spk_mm/max(spk,1)*100:.2f}%)")
    print(msg)


if __name__ == "__main__":
    main()
