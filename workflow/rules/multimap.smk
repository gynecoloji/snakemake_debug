# Step 09 — MULTIMAPPED-read investigation (SPIKE-IN SEQUENCING ONLY): the host
# vs spike-in split and per-genome multimap rate (multimap_analysis), plus the
# cross-genome check of how many multimapped reads map to BOTH genomes
# (crossgenome_analysis). Gated by MM_ON / MM_XGENOME. Shared config, constants
# and target lists live in common.smk.


# Aggregate target for the multimap stage. Run it alone with:
#   snakemake --cores N multimap_all
rule multimap_all:
    input:
        MULTIMAP_TARGETS,


# ---- 09: investigate MULTIMAPPED reads in the raw aligner output (mode-aware) ----
# Detection auto: NH:i>1 (HISAT2) -> XS present (bowtie2) -> MAPQ<=max (fallback).
# Spike-in mode: exact host vs spike-in split (idxstats) + the multimapped DISTRIBUTION
# (what % of each genome's reads multimap, and the top contigs they land on).
# Reads MM_ALIGN_DIR (may differ from the unaligned-analysis align_dir).
rule multimap_analysis:
    input: MM_ALIGN_DIR + "/{sample}" + MM_SUF
    output:
        summ = f"{OUT}/09_multimapped_reads/persample/{{sample}}.multimap.tsv",
        loci = f"{OUT}/09_multimapped_reads/persample/{{sample}}.loci.tsv",
    params: spike=("--spikein" if MM_SPIKE else ""), prefix=MM_PREFIX, method=MM_METHOD,
            xs=("--xs-strict" if MM_XSSTR else ""), mapq=MM_MAPQ, rcap=RCAP
    shell:
        "python " + SCRIPTS + "/multimap_analysis.py --align {input} --sample {wildcards.sample} "
        "{params.spike} --spikein-prefix {params.prefix:q} --method {params.method} {params.xs} "
        "--mapq-max {params.mapq} --records-cap {params.rcap} "
        "--out {output.summ} --out-loci {output.loci}"

rule aggregate_multimap:
    input:
        summ = expand(f"{OUT}/09_multimapped_reads/persample/{{s}}.multimap.tsv", s=MM_SAMPLES),
        loci = expand(f"{OUT}/09_multimapped_reads/persample/{{s}}.loci.tsv", s=MM_SAMPLES),
    output:
        summ = f"{OUT}/09_multimapped_reads/multimap_summary.tsv",
        loci = f"{OUT}/09_multimapped_reads/top_multimap_loci.tsv",
    shell:
        "python " + SCRIPTS + "/aggregate_multimap.py --out-summary {output.summ} "
        "--out-loci {output.loci} --in {input.summ} --in-loci {input.loci}"

# ---- 09 cross-genome: of the multimapped reads, how many map to BOTH genomes ----
# Re-aligns a sample of multimapped reads to the host-only and spike-in-only indexes.
rule crossgenome_analysis:
    input: MM_ALIGN_DIR + "/{sample}" + MM_SUF
    output: f"{OUT}/09_multimapped_reads/crossgenome/{{sample}}.crossgenome.tsv"
    params: host=MM_HOST_IDX, spk=MM_SPK_IDX, prefix=MM_PREFIX, frac=MM_XFRAC,
            seed=MM_XSEED, cap=MM_XCAP, margin=MM_XMARGIN
    threads: THREADS
    shell:
        "python " + SCRIPTS + "/crossgenome_analysis.py --align {input} --sample {wildcards.sample} "
        "--host-index {params.host:q} --spikein-index {params.spk:q} "
        "--spikein-prefix {params.prefix:q} --fraction {params.frac} --seed {params.seed} "
        "--cap {params.cap} --codominant-margin {params.margin} --threads {threads} --out {output}"

rule aggregate_crossgenome:
    input: expand(f"{OUT}/09_multimapped_reads/crossgenome/{{s}}.crossgenome.tsv", s=MM_SAMPLES)
    output: f"{OUT}/09_multimapped_reads/crossgenome_summary.tsv"
    shell:
        "python " + SCRIPTS + "/aggregate_crossgenome.py --out-summary {output} --in {input}"
