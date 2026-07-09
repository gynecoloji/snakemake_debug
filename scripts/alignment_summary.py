#!/usr/bin/env python3
"""Step 01: per-sample alignment rate via `samtools flagstat` (aligner-agnostic)."""
import argparse, glob, os, re, subprocess
from concurrent.futures import ThreadPoolExecutor

def rate(path):
    out = subprocess.run(["samtools", "flagstat", path], capture_output=True, text=True).stdout
    total = pct = None
    for ln in out.splitlines():
        m = re.match(r"(\d+) \+ \d+ in total", ln)
        if m: total = int(m.group(1))
        m = re.match(r"\d+ \+ \d+ primary mapped \(([\d.]+)%", ln)
        if m: pct = float(m.group(1))
    if pct is None:                                  # older samtools fallback
        for ln in out.splitlines():
            m = re.search(r"mapped \(([\d.]+)%", ln)
            if m: pct = float(m.group(1)); break
    return total, pct

def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--align-dir", required=True)
    ap.add_argument("--suffix", required=True)
    ap.add_argument("--threads", type=int, default=8)
    ap.add_argument("--out", required=True)
    a = ap.parse_args()

    files = sorted(glob.glob(f"{a.align_dir}/*{a.suffix}"))
    samples = [os.path.basename(f)[:-len(a.suffix)] for f in files]
    with ThreadPoolExecutor(max_workers=a.threads) as ex:
        res = list(ex.map(rate, files))

    rows = sorted(zip(samples, res), key=lambda r: (r[1][1] if r[1][1] is not None else 999))
    os.makedirs(os.path.dirname(a.out), exist_ok=True)
    with open(a.out, "w") as o:
        o.write("sample\talignment_rate_pct\ttotal_reads\n")
        for s, (tot, pct) in rows:
            o.write(f"{s}\t{pct if pct is not None else 'NA'}\t{tot if tot is not None else 'NA'}\n")
    print(f"[01] wrote {a.out} ({len(samples)} samples)")

if __name__ == "__main__":
    main()
