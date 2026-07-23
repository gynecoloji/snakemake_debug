# Unaligned-read foundation: per-sample alignment rate (01), extraction of the
# unaligned reads (a bounded sample for the motif/BLAST/Kraken2 steps, plus the
# full set for the bowtie2 genome screens), and the sequence-signature screen
# (02). Shared config, constants and target lists live in common.smk.


# Aggregate target for the unaligned/signature stage. Run it alone with:
#   snakemake --cores N unaligned_all
rule unaligned_all:
    input:
        UNALIGNED_TARGETS,


# ---- 01: alignment rate per sample (aligner-agnostic, via samtools flagstat)
rule alignment_summary:
    input: expand(ALIGN_DIR + "/{s}" + SUF, s=SAMPLES)
    output: f"{OUT}/01_alignment_summary/alignment_summary.tsv"
    params: d=ALIGN_DIR, suf=SUF
    threads: THREADS
    shell:
        "python " + SCRIPTS + "/alignment_summary.py "
        "--align-dir {params.d} --suffix {params.suf} --threads {threads} --out {output}"

# ---- extract unaligned reads (bounded sample); tally top seqs; score signatures (02)
rule extract_unaligned:
    input: ALIGN_DIR + "/{sample}" + SUF
    output:
        fq  = f"{OUT}/unaligned_reads/{{sample}}.unaligned.fq.gz",
        top = f"{OUT}/unaligned_reads/{{sample}}.topseq.tsv",
        sig = f"{OUT}/02_sequence_signatures/persample/{{sample}}.sig.tsv",
    params: flag=FLAG, ucap=UCAP, rcap=RCAP
    shell:
        "python " + SCRIPTS + "/extract_unaligned.py --align {input} --sample {wildcards.sample} "
        "--flag {params.flag} --unmapped-cap {params.ucap} --records-cap {params.rcap} "
        "--out-fastq {output.fq} --out-topseq {output.top} --out-sig {output.sig}"

# ---- full (uncapped) unaligned-read extraction for the bowtie2 genome screens ------
# The bowtie2 steps (05 contaminant, 07 custom, 08 top-organism) align ALL unaligned
# reads, not the bounded sample used by the motif/BLAST/Kraken2 steps. Like the sampled
# extraction, this feeds the screens but is not itself a report section.
rule extract_all_unaligned:
    input: ALIGN_DIR + "/{sample}" + SUF
    output: f"{OUT}/unaligned_reads/{{sample}}.all_unaligned.fq.gz"
    params: flag=FLAG
    threads: THREADS
    shell:
        "samtools view -@ {threads} -f {params.flag} -F 0x900 {input} | "
        "awk -F'\\t' '$10!=\"*\"{{print \"@\"$1\"\\n\"$10\"\\n+\\n\"$11}}' | gzip > {output}"

rule aggregate_signatures:
    input: expand(f"{OUT}/02_sequence_signatures/persample/{{s}}.sig.tsv", s=SAMPLES)
    output: f"{OUT}/02_sequence_signatures/signature_fractions.tsv"
    shell: "python " + SCRIPTS + "/aggregate_signatures.py --out {output} --in {input}"
