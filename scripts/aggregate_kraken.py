#!/usr/bin/env python3
"""Aggregate Kraken2 per-sample reports: %classified, %bacteria, top species."""
import argparse, os, collections

def parse(rp):
    classified = unclass = bacteria = 0
    species = []
    for line in open(rp):
        c = line.rstrip("\n").split("\t")
        if len(c) < 6: continue
        clade = int(c[1]); rank = c[3]; name = c[5].strip()
        if   rank == "U":                       unclass = clade
        elif rank == "R" and name == "root":    classified = clade
        elif rank == "D" and name == "Bacteria": bacteria = clade
        elif rank == "S":                        species.append((clade, name))
    tot = classified + unclass
    species.sort(reverse=True)
    return dict(total=max(tot, 1), classified=classified, bacteria=bacteria, top=species[:5])

def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--in", dest="inp", nargs="+", required=True, help="*.kreport files")
    ap.add_argument("--out-summary", required=True)
    ap.add_argument("--out-species", required=True)
    a = ap.parse_args()
    os.makedirs(os.path.dirname(a.out_summary), exist_ok=True)

    per = {}; overall = collections.Counter()
    for rp in a.inp:
        smp = os.path.basename(rp).replace(".kreport", "")
        d = per[smp] = parse(rp)
        for reads, name in d["top"]: overall[name] += reads

    with open(a.out_summary, "w") as o:
        o.write("sample\ttotal_reads\tpct_classified\tpct_bacteria\ttop_species\n")
        for smp in sorted(per, key=lambda s: -per[s]["bacteria"]/per[s]["total"]):
            d = per[smp]; t = d["total"]
            top = "; ".join(f"{n}({r})" for r, n in d["top"][:3])
            o.write(f"{smp}\t{t}\t{100*d['classified']/t:.1f}\t{100*d['bacteria']/t:.1f}\t{top}\n")
    with open(a.out_species, "w") as o:
        o.write("species\ttotal_reads_across_samples\n")
        for name, c in overall.most_common(25):
            o.write(f"{name}\t{c}\n")
    print(f"[06] wrote {a.out_summary} + {a.out_species}")

if __name__ == "__main__":
    main()
