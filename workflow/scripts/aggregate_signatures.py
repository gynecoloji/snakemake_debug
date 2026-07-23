#!/usr/bin/env python3
"""Concatenate per-sample signature rows into one table (with header)."""
import argparse, os

def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--in", dest="inp", nargs="+", required=True)
    ap.add_argument("--out", required=True)
    a = ap.parse_args()
    os.makedirs(os.path.dirname(a.out), exist_ok=True)
    rows = []
    for f in a.inp:
        rows.append(open(f).read().rstrip("\n"))
    rows = [r for r in rows if r]
    # sort by bact_16S fraction desc (col index 3)
    rows.sort(key=lambda r: -float(r.split("\t")[3]))
    with open(a.out, "w") as o:
        o.write("sample\tsampled_unaligned\trecords_scanned\tbact_16S\t"
                "illumina_adapter\tpolyA\thost_repeat\tother\n")
        o.write("\n".join(rows) + "\n")
    print(f"[02] wrote {a.out} ({len(rows)} samples)")

if __name__ == "__main__":
    main()
