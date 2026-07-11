#!/usr/bin/env python3
"""Generic diagnostic plots. Reads whatever step TSVs exist under --out-dir,
skips plots whose inputs are missing. matplotlib Agg -> PNG.
Palette is colorblind-safe (validated): human=blue, contaminant/fly=orange,
poly/satellite=purple, bacterial/other-organism=green, neutral=grey."""
import argparse, csv, os, re
import matplotlib; matplotlib.use("Agg")
import matplotlib.pyplot as plt

# --- validated CVD-safe palette (see dataviz palette validator) --------------
HUMAN, FLY, SAT, NOISE, GREY = "#0072B2", "#D55E00", "#CC79A7", "#009E73", "#c3c5c1"
SURF, INK, MUTED = "#fcfcfb", "#0b0b0b", "#898781"
plt.rcParams.update({"figure.facecolor": SURF, "axes.facecolor": SURF, "savefig.facecolor": SURF,
    "text.color": INK, "axes.labelcolor": INK, "xtick.color": MUTED, "ytick.color": MUTED,
    "axes.edgecolor": MUTED, "font.size": 9, "axes.titlesize": 11, "axes.titleweight": "bold"})


def rd(p):
    with open(p) as f: return list(csv.DictReader(f, delimiter="\t"))
def style(ax):
    ax.spines["top"].set_visible(False); ax.spines["right"].set_visible(False); ax.set_axisbelow(True)
def short(s): return re.sub(r"^GSF\d+-", "", s)


def parse_fasta_counts(path):
    out = {}
    if not os.path.exists(path): return out
    cur = None
    with open(path) as f:
        for line in f:
            line = line.rstrip("\n")
            if line.startswith(">"):
                cur = line[1:].split("_", 1)[0]
                m = re.search(r"totalcount=(\d+)", line)
                out[cur] = {"n": int(m.group(1)) if m else 0, "seq": ""}
            elif cur is not None:
                out[cur]["seq"] += line.strip()
    return out
def is_lowcomplex(seq):
    if len(seq) < 12: return True
    k = {seq[i:i+3] for i in range(len(seq)-2)}
    return len(k)/(len(seq)-2) < 0.5


def main():
    ap = argparse.ArgumentParser(); ap.add_argument("--out-dir", required=True)
    O = ap.parse_args().out_dir
    os.makedirs(f"{O}/plots", exist_ok=True)
    made = []

    # 1) alignment rate per sample (outliers highlighted)
    p = f"{O}/01_alignment_summary/alignment_summary.tsv"
    if os.path.exists(p):
        r = [x for x in rd(p) if x["alignment_rate_pct"] not in ("", "NA")]
        r.sort(key=lambda x: float(x["alignment_rate_pct"]))
        vals = [float(x["alignment_rate_pct"]) for x in r]
        hi = max(vals) if vals else 100
        cols = [FLY if v < hi - 10 else HUMAN for v in vals]
        fig, ax = plt.subplots(figsize=(9, max(3, 0.32*len(r))))
        ax.barh([short(x["sample"]) for x in r], vals, color=cols)
        ax.set_xlabel("alignment rate (%)"); ax.set_xlim(0, 100)
        ax.set_title("Alignment rate per sample (orange = >10% below best)"); style(ax)
        fig.tight_layout(); fig.savefig(f"{O}/plots/alignment_rate.png", dpi=140); plt.close(fig)
        made.append("alignment_rate.png")

    # 2) signature composition (stacked)
    p = f"{O}/02_sequence_signatures/signature_fractions.tsv"
    if os.path.exists(p):
        r = rd(p)
        cats = [("bact_16S", NOISE), ("host_repeat", HUMAN), ("illumina_adapter", FLY),
                ("polyA", SAT), ("other", GREY)]
        fig, ax = plt.subplots(figsize=(9, max(3, 0.32*len(r)))); left = [0.0]*len(r)
        for key, col in cats:
            v = [float(x[key])*100 for x in r]
            ax.barh([short(x["sample"]) for x in r], v, left=left, color=col, label=key, edgecolor=SURF, linewidth=0.5)
            left = [a+b for a, b in zip(left, v)]
        ax.set_xlabel("composition of unaligned reads (%)"); ax.set_xlim(0, 100)
        ax.set_title("Unaligned-read motif signatures"); style(ax)
        ax.legend(loc="lower right", frameon=False, fontsize=8)
        fig.tight_layout(); fig.savefig(f"{O}/plots/signature_composition.png", dpi=140); plt.close(fig)
        made.append("signature_composition.png")

    # 3) BLAST identity of the top unaligned sequences (step 03)
    fa = f"{O}/03_blast_top_sequences/top_unaligned.fasta"
    counts = parse_fasta_counts(fa)
    if counts:
        hits = {}
        for x in ([] if not os.path.exists(f"{O}/03_blast_top_sequences/blast_best_hits.tsv")
                  else rd(f"{O}/03_blast_top_sequences/blast_best_hits.tsv")):
            hits.setdefault(x.get("query", "").split("_", 1)[0], x.get("subject", ""))
        def cat(tid, info):
            subj = hits.get(tid, "").lower()
            if subj:
                if "mitochond" in subj or "d-loop" in subj: return HUMAN, "human mtDNA"
                if "drosophila" in subj: return FLY, "Drosophila"
                if "homo sapiens" in subj or "human" in subj: return HUMAN, "Homo sapiens"
                return NOISE, "other hit"
            return (SAT, "repeat/no hit") if is_lowcomplex(info["seq"]) else (GREY, "no hit")
        items = sorted(counts.items(), key=lambda kv: kv[1]["n"])
        labs = [k for k, _ in items]; vals = [v["n"] for _, v in items]
        cols = [cat(k, v)[0] for k, v in items]
        fig, ax = plt.subplots(figsize=(9, max(3, 0.34*len(items))))
        ax.barh(labs, vals, color=cols)
        ax.set_xlabel("duplicate count across samples")
        ax.set_title("Top unaligned sequences by BLAST identity"); style(ax)
        seen = {}
        for k, v in items:
            c, name = cat(k, v)
            if name not in seen: seen[name] = c
        handles = [plt.Rectangle((0, 0), 1, 1, color=c) for c in seen.values()]
        ax.legend(handles, list(seen.keys()), loc="lower right", frameon=False, fontsize=8)
        fig.tight_layout(); fig.savefig(f"{O}/plots/blast_top_sequences.png", dpi=140); plt.close(fig)
        made.append("blast_top_sequences.png")

    # 3b) random-read BLAST species (step 04)
    p = f"{O}/04_blast_random/random_blast_species.tsv"
    if os.path.exists(p):
        r = rd(p)[:12][::-1]
        if r:
            labs = [x["species"] for x in r]; vals = [int(x["reads"]) for x in r]
            cols = [HUMAN if re.search(r"homo sapiens", x["species"], re.I)
                    else (FLY if re.search(r"drosophila", x["species"], re.I) else NOISE) for x in r]
            fig, ax = plt.subplots(figsize=(9, max(3, 0.36*len(r))))
            ax.barh(labs, vals, color=cols)
            ax.set_xlabel("random unaligned reads with this best nt hit")
            ax.set_title("Random-read BLAST: species discovered"); style(ax)
            fig.tight_layout(); fig.savefig(f"{O}/plots/random_blast_species.png", dpi=140); plt.close(fig)
            made.append("random_blast_species.png")

    # 4) Kraken bacterial fraction
    p = f"{O}/06_kraken2/kraken_summary.tsv"
    if os.path.exists(p):
        r = rd(p); r.sort(key=lambda x: float(x["pct_bacteria"]))
        fig, ax = plt.subplots(figsize=(9, max(3, 0.32*len(r))))
        ax.barh([short(x["sample"]) for x in r], [float(x["pct_bacteria"]) for x in r], color=NOISE)
        ax.set_xlabel("% of unaligned reads = Bacteria (Kraken2)")
        ax.set_title("Bacterial contamination of unaligned reads"); style(ax)
        fig.tight_layout(); fig.savefig(f"{O}/plots/bacterial_fraction.png", dpi=140); plt.close(fig)
        made.append("bacterial_fraction.png")

    # 5) Kraken top classified species overall (step 06)
    p = f"{O}/06_kraken2/top_species_overall.tsv"
    if os.path.exists(p):
        r = rd(p)[:10][::-1]
        if r:
            labs = [x["species"] for x in r]
            vals = [max(int(float(x["total_reads_across_samples"])), 1) for x in r]
            cols = [HUMAN if re.search(r"homo sapiens", x["species"], re.I) else NOISE for x in r]
            fig, ax = plt.subplots(figsize=(9, max(3, 0.36*len(r))))
            ax.barh(labs, vals, color=cols)
            ax.set_xscale("log"); ax.set_xlabel("reads classified across samples (log scale)")
            ax.set_title("Top classified species in unaligned reads"); style(ax)
            fig.tight_layout(); fig.savefig(f"{O}/plots/kraken_top_species.png", dpi=140); plt.close(fig)
            made.append("kraken_top_species.png")

    # 6) contaminant genome mapping
    p = f"{O}/05_contaminant_genome_screen/contaminant_alignment.tsv"
    if os.path.exists(p):
        r = rd(p); r.sort(key=lambda x: float(x["pct_reads_mapping_contaminant"]))
        fig, ax = plt.subplots(figsize=(9, max(3, 0.32*len(r))))
        ax.barh([short(x["sample"]) for x in r], [float(x["pct_reads_mapping_contaminant"]) for x in r], color=FLY)
        ax.set_xlabel("% of unaligned reads mapping to contaminant genome")
        ax.set_title("Contaminant genome mapping (bowtie2)"); style(ax)
        fig.tight_layout(); fig.savefig(f"{O}/plots/contaminant_mapping.png", dpi=140); plt.close(fig)
        made.append("contaminant_mapping.png")

    # 6b) custom-sequence (construct) mapping (step 07)
    pm = f"{O}/07_custom_sequences/custom_mapping.tsv"
    if os.path.exists(pm):
        r = rd(pm); r.sort(key=lambda x: float(x["pct_mapping_custom"]))
        fig, ax = plt.subplots(figsize=(9, max(3, 0.32*len(r))))
        ax.barh([short(x["sample"]) for x in r], [float(x["pct_mapping_custom"]) for x in r], color=SAT)
        ax.set_xlabel("% of unaligned reads mapping to custom sequences")
        ax.set_title("Custom-sequence (construct) mapping"); style(ax)
        fig.tight_layout(); fig.savefig(f"{O}/plots/custom_mapping.png", dpi=140); plt.close(fig)
        made.append("custom_mapping.png")
        pp = f"{O}/07_custom_sequences/custom_per_sequence.tsv"
        rows = rd(pp) if os.path.exists(pp) else []
        seqcols = [c for c in rows[0].keys() if c != "sample"] if rows else []
        if seqcols:
            names = [re.sub(r"_5'.*$", "", c) for c in seqcols]
            pal = [SAT, FLY, NOISE, HUMAN]
            labs = [short(x["sample"]) for x in rows]
            fig, ax = plt.subplots(figsize=(9, max(3, 0.32*len(rows)))); left = [0.0]*len(rows)
            tots = [max(sum(int(x[cc]) for cc in seqcols), 1) for x in rows]
            for i, c in enumerate(seqcols):
                v = [int(x[c])/t*100 for x, t in zip(rows, tots)]
                ax.barh(labs, v, left=left, color=pal[min(i, 3)], label=names[i], edgecolor=SURF, linewidth=0.5)
                left = [a+b for a, b in zip(left, v)]
            ax.set_xlabel("composition of construct-mapped reads (%)"); ax.set_xlim(0, 100)
            ax.set_title("Which custom sequence the reads hit"); style(ax)
            ax.legend(loc="lower right", frameon=False, fontsize=8)
            fig.tight_layout(); fig.savefig(f"{O}/plots/custom_per_sequence.png", dpi=140); plt.close(fig)
            made.append("custom_per_sequence.png")

    # 6c) top-organism (auto-fetched RefSeq) genome screen (step 08)
    pto = f"{O}/08_top_organism/top_organism_alignment.tsv"
    if os.path.exists(pto):
        r = [x for x in rd(pto) if x.get("organism") and x["organism"] != "NONE"]
        if r:
            org = r[0]["organism"]
            r.sort(key=lambda x: float(x["pct_reads_mapping_top_organism"]))
            fig, ax = plt.subplots(figsize=(9, max(3, 0.32*len(r))))
            ax.barh([short(x["sample"]) for x in r],
                    [float(x["pct_reads_mapping_top_organism"]) for x in r], color=NOISE)
            ax.set_xlabel(f"% of unaligned reads mapping to {org}")
            ax.set_title(f"Top-organism genome screen: {org}"); style(ax)
            fig.tight_layout(); fig.savefig(f"{O}/plots/top_organism_mapping.png", dpi=140); plt.close(fig)
            made.append("top_organism_mapping.png")

    # 7) whole-library composition: host / top contamination (bowtie2) / others
    #    Use the precise bowtie2 genome alignment — step 08 auto top-organism if
    #    present, else the step-05 contaminant genome — and lump the rest as "others".
    pa = f"{O}/01_alignment_summary/alignment_summary.tsv"
    if os.path.exists(pa):
        aln = {x["sample"]: float(x["alignment_rate_pct"]) for x in rd(pa)
               if x["alignment_rate_pct"] not in ("", "NA")}
        top, label = {}, "contaminant"
        p09 = f"{O}/08_top_organism/top_organism_alignment.tsv"
        p05 = f"{O}/05_contaminant_genome_screen/contaminant_alignment.tsv"
        if os.path.exists(p09):
            r9 = [x for x in rd(p09) if x.get("organism") and x["organism"] != "NONE"]
            if r9:
                top = {x["sample"]: float(x["pct_reads_mapping_top_organism"]) for x in r9}
                label = r9[0]["organism"]
        if not top and os.path.exists(p05):
            top = {x["sample"]: float(x["pct_reads_mapping_contaminant"]) for x in rd(p05)}
            label = "contaminant genome"
        smps = sorted(aln, key=lambda s: aln[s])
        if smps:
            host = [aln[s] for s in smps]
            contam = [max(100-aln[s], 0.0) * min(top.get(s, 0)/100, 1.0) for s in smps]
            oth = [max(100-aln[s], 0.0) - c for s, c in zip(smps, contam)]
            labs = [short(s) for s in smps]
            fig, ax = plt.subplots(figsize=(9, max(3, 0.32*len(smps))))
            ax.barh(labs, host, color=HUMAN, label="host (aligned)")
            ax.barh(labs, contam, left=host, color=FLY, label=f"top contamination: {label}")
            ax.barh(labs, oth, left=[h+c for h, c in zip(host, contam)], color=GREY, label="others")
            ax.set_xlabel("share of library (%)"); ax.set_xlim(0, 100)
            ax.set_title("Whole-library composition"); style(ax)
            ax.legend(loc="lower right", frameon=False, fontsize=8)
            fig.tight_layout(); fig.savefig(f"{O}/plots/library_composition.png", dpi=140); plt.close(fig)
            made.append("library_composition.png")

    # 8) multimapped reads (step 09)
    p = f"{O}/09_multimapped_reads/multimap_summary.tsv"
    spike = False
    if os.path.exists(p):
        r = rd(p)
        if r:
            spike = any(float(x.get("pct_spikein", 0) or 0) > 0 for x in r)
            r.sort(key=lambda x: float(x["pct_multimapped"]))
            labs = [short(x["sample"]) for x in r]
            fig, ax = plt.subplots(figsize=(9, max(3, 0.32*len(r))))
            if spike:
                host = [float(x["pct_host"]) for x in r]
                spk = [float(x["pct_spikein"]) for x in r]
                ax.barh(labs, host, color=HUMAN, label="host genome")
                ax.barh(labs, spk, left=host, color=FLY, label="spike-in genome")
                ax.set_xlabel("host vs spike-in genome split (%)"); ax.set_xlim(0, 100)
                ax.set_title("Spike-in genome split (multimapped % in report)")
                ax.legend(loc="lower right", frameon=False, fontsize=8)
            else:
                ax.barh(labs, [float(x["pct_multimapped"]) for x in r], color=SAT)
                ax.set_xlabel("% of mapped reads that are multimapped")
                ax.set_title("Multimapped reads per sample")
            style(ax); fig.tight_layout()
            fig.savefig(f"{O}/plots/multimap_composition.png", dpi=140); plt.close(fig)
            made.append("multimap_composition.png")
            if spike:   # multimapper distribution BY genome: host mm% vs spike-in mm%
                y = list(range(len(r)))
                hmm = [float(x.get("pct_host_multimapped", 0) or 0) for x in r]
                smm = [float(x.get("pct_spikein_multimapped", 0) or 0) for x in r]
                fig, ax = plt.subplots(figsize=(9, max(3, 0.42*len(r))))
                ax.barh([yi-0.2 for yi in y], hmm, height=0.4, color=HUMAN, label="host reads")
                ax.barh([yi+0.2 for yi in y], smm, height=0.4, color=FLY, label="spike-in reads")
                ax.set_yticks(y); ax.set_yticklabels(labs)
                ax.set_xlabel("% of that genome's reads that are multimapped")
                ax.set_title("Multimapping rate within each genome")
                ax.legend(loc="lower right", frameon=False, fontsize=8)
                style(ax); fig.tight_layout()
                fig.savefig(f"{O}/plots/multimap_by_genome.png", dpi=140); plt.close(fig)
                made.append("multimap_by_genome.png")
        pl = f"{O}/09_multimapped_reads/top_multimap_loci.tsv"
        rr = rd(pl) if os.path.exists(pl) else []
        if rr:
            rr = rr[:15][::-1]
            fig, ax = plt.subplots(figsize=(9, max(3, 0.34*len(rr))))
            ax.barh([x["contig"] for x in rr], [int(x["multimapped_reads"]) for x in rr], color=NOISE)
            ax.set_xlabel("multimapped reads (pooled across samples)")
            ax.set_title("Where multimappers concentrate (top contigs)"); style(ax)
            fig.tight_layout(); fig.savefig(f"{O}/plots/multimap_loci.png", dpi=140); plt.close(fig)
            made.append("multimap_loci.png")
        xgp = f"{O}/09_multimapped_reads/crossgenome_summary.tsv"
        xg = rd(xgp) if os.path.exists(xgp) else []
        if xg:
            xg.sort(key=lambda x: float(x.get("pct_both", 0) or 0))
            labs = [short(x["sample"]) for x in xg]
            both = [float(x.get("pct_both", 0) or 0) for x in xg]
            codom = [float(x.get("pct_both_codominant", 0) or 0) for x in xg]
            fig, ax = plt.subplots(figsize=(9, max(3, 0.34*len(xg))))
            ax.barh(labs, both, color=SAT, label="maps to both genomes")
            ax.barh(labs, codom, color=FLY, label="codominant (ambiguous)")
            ax.set_xlabel("% of multimapped reads mapping to BOTH genomes")
            ax.set_title("Cross-genome multimapped reads")
            ax.legend(loc="lower right", frameon=False, fontsize=8); style(ax)
            fig.tight_layout(); fig.savefig(f"{O}/plots/multimap_crossgenome.png", dpi=140); plt.close(fig)
            made.append("multimap_crossgenome.png")
            # re-alignment outcome composition: host-only / spike-only / both / neither (stacked to 100%)
            xg.sort(key=lambda x: float(x.get("pct_host_only", 0) or 0))
            labs = [short(x["sample"]) for x in xg]
            ho = [float(x.get("pct_host_only", 0) or 0) for x in xg]
            so = [float(x.get("pct_spike_only", 0) or 0) for x in xg]
            bo = [float(x.get("pct_both", 0) or 0) for x in xg]
            ne = [float(x.get("pct_neither", 0) or 0) for x in xg]
            l1 = ho
            l2 = [a + b for a, b in zip(ho, so)]
            l3 = [a + b for a, b in zip(l2, bo)]
            fig, ax = plt.subplots(figsize=(9, max(3, 0.34*len(xg))))
            ax.barh(labs, ho, color=HUMAN, label="host-only")
            ax.barh(labs, so, left=l1, color=FLY, label="spike-only")
            ax.barh(labs, bo, left=l2, color=SAT, label="both genomes")
            ax.barh(labs, ne, left=l3, color=GREY, label="neither")
            ax.set_xlim(0, 100)
            ax.set_xlabel("% of multimapped reads (re-aligned per genome)")
            ax.set_title("Re-alignment outcome of multimapped reads")
            ax.legend(loc="upper center", bbox_to_anchor=(0.5, -0.28), ncol=4,
                      frameon=False, fontsize=8); style(ax)
            fig.tight_layout()
            fig.savefig(f"{O}/plots/multimap_crossgenome_composition.png", dpi=140, bbox_inches="tight")
            plt.close(fig)
            made.append("multimap_crossgenome_composition.png")

    print(f"[plots] wrote: {', '.join(made) if made else '(none)'}")

if __name__ == "__main__":
    main()
