#!/usr/bin/env python3
"""Parse bowtie2 logs -> % of unaligned reads mapping to the contaminant genome."""
import argparse, os, re

def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--in", dest="inp", nargs="+", required=True, help="*.bt2.log files")
    ap.add_argument("--out", required=True)
    a = ap.parse_args()
    os.makedirs(os.path.dirname(a.out), exist_ok=True)

    rows = []
    for log in a.inp:
        smp = os.path.basename(log).replace(".bt2.log", "")
        t = open(log).read()
        mt = re.search(r"(\d+) reads; of these", t)
        mr = re.search(r"([\d.]+)% overall alignment rate", t)
        if mt and mr:
            rows.append((smp, int(mt.group(1)), float(mr.group(1))))
    with open(a.out, "w") as o:
        o.write("sample\tsampled_unaligned_reads\tpct_reads_mapping_contaminant\n")
        for smp, tot, pct in sorted(rows, key=lambda r: -r[2]):
            o.write(f"{smp}\t{tot}\t{pct:.2f}\n")
    print(f"[05] wrote {a.out} ({len(rows)} samples)")

if __name__ == "__main__":
    main()
