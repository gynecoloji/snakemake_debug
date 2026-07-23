#!/usr/bin/env python3
"""Generate a generic 00_SUMMARY.md from whatever step outputs exist. Prints to stdout."""
import argparse, csv, os

def rd(p):
    with open(p) as f: return list(csv.DictReader(f, delimiter="\t"))

def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--out-dir", required=True)
    ap.add_argument("--host", default="host")
    a = ap.parse_args(); O = a.out_dir
    L = []
    P = L.append
    P("# Low-alignment / contamination diagnostics — summary\n")

    # alignment
    pa = f"{O}/01_alignment_summary/alignment_summary.tsv"
    low = []
    if os.path.exists(pa):
        r = [x for x in rd(pa) if x["alignment_rate_pct"] not in ("", "NA")]
        rates = [(x["sample"], float(x["alignment_rate_pct"])) for x in r]
        if rates:
            hi = max(v for _, v in rates); lo = min(v for _, v in rates)
            low = [s for s, v in rates if v < hi - 10]
            P(f"- **Samples:** {len(rates)}  |  **alignment rate:** {lo:.1f}%–{hi:.1f}%")
            if low:
                P(f"- **Low-alignment samples (>10% below best):** {', '.join(low)}")
            else:
                P("- No sample is >10% below the best — alignment looks uniform.")
    # kraken top species
    ps = f"{O}/06_kraken2/top_species_overall.tsv"
    if os.path.exists(ps):
        sp = rd(ps)
        if sp:
            top = sp[0]
            P(f"- **Top organism in unaligned reads (Kraken2):** {top['species']} "
              f"({top['total_reads_across_samples']} reads across samples)")
            P("\n### Most abundant organisms in the unaligned pool")
            P("| species | reads |\n|---|---|")
            for x in sp[:8]:
                P(f"| {x['species']} | {x['total_reads_across_samples']} |")
    # contaminant genome
    pc = f"{O}/05_contaminant_genome_screen/contaminant_alignment.tsv"
    if os.path.exists(pc):
        c = rd(pc)
        vals = [float(x["pct_reads_mapping_contaminant"]) for x in c]
        if vals:
            P(f"\n- **Contaminant-genome mapping of unaligned reads:** "
              f"{min(vals):.1f}%–{max(vals):.1f}% across samples")

    # custom sequences (construct/transgene)
    pcu = f"{O}/07_custom_sequences/custom_mapping.tsv"
    if os.path.exists(pcu):
        cu = rd(pcu)
        vals = [float(x["pct_mapping_custom"]) for x in cu]
        if vals:
            P(f"- **Custom-sequence (construct) mapping of unaligned reads:** "
              f"{min(vals):.1f}%–{max(vals):.1f}% across samples")

    # multimapped reads (step 09)
    pmm = f"{O}/09_multimapped_reads/multimap_summary.tsv"
    if os.path.exists(pmm):
        m = rd(pmm)
        vals = [float(x["pct_multimapped"]) for x in m]
        if vals:
            P(f"- **Multimapped reads:** {min(vals):.1f}%–{max(vals):.1f}% of mapped reads across samples")
        spk = [float(x.get("pct_spikein", 0) or 0) for x in m]
        if any(v > 0 for v in spk):
            P(f"- **Spike-in genome fraction:** {min(spk):.1f}%–{max(spk):.1f}% of mapped reads "
              f"(the normalization signal)")
    pxg = f"{O}/09_multimapped_reads/crossgenome_summary.tsv"
    if os.path.exists(pxg):
        x = rd(pxg)
        b = [float(r["pct_both"]) for r in x]
        c = [float(r["pct_both_codominant"]) for r in x]
        if b:
            P(f"- **Multimapped reads mapping to BOTH genomes:** {min(b):.2f}%–{max(b):.2f}% "
              f"({min(c):.2f}%–{max(c):.2f}% codominant/ambiguous) — cross-genome mapping that can "
              f"distort spike-in normalization")

    # interpretation hint
    P("\n### Reading this")
    P("- Low alignment confined to a subset of samples + a dominant organism in the "
      "unaligned reads = **contamination** (see that organism above).")
    P(f"- Unaligned reads that are mostly {a.host} repeats/uniquely-unplaceable are **benign**.")
    P("- All samples low + reads match the host = wrong reference / strandedness / index.")

    # folder guide
    P("\n### Outputs")
    P("| folder | contents |\n|---|---|")
    P("| `01_alignment_summary/` | per-sample alignment rate |")
    P("| `unaligned_reads/` | sampled unaligned reads + top-sequence tallies |")
    P("| `02_sequence_signatures/` | motif composition of unaligned reads |")
    P("| `03_blast_top_sequences/` | top sequences (+ BLAST identity if run) |")
    P("| `05_contaminant_genome_screen/` | reads mapping to the suspect genome |")
    P("| `06_kraken2/` | taxonomic classification |")
    P("| `07_custom_sequences/` | unaligned reads mapped to a custom fasta (per-sequence) |")
    P("| `09_multimapped_reads/` | multimapped-read investigation (spike-in split, per-genome rate, cross-genome both-genomes %) |")
    P("| `plots/` | figures |")

    print("\n".join(L))

if __name__ == "__main__":
    main()
