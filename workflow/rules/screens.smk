# Genome / taxonomic screens of the unaligned reads: a suspect contaminant
# genome (05, bowtie2), Kraken2 classification (06), a custom construct/vector
# fasta (07, bowtie2 + idxstats), and the auto-picked top non-host organism
# (08, Kraken2 + random-BLAST -> NCBI RefSeq -> bowtie2). Steps 06/05/07 are
# gated by their config paths; step 08 needs internet (login node). Shared
# config, constants and target lists live in common.smk.


# Aggregate target for the screens stage. Run it alone with:
#   snakemake --cores N screens_all      (step 08 needs a login node)
rule screens_all:
    input:
        SCREENS_TARGETS,


# ---- 08: auto-pick the top non-host organism (Kraken2 + random-BLAST), fetch its
#          RefSeq genome from NCBI, and align unaligned reads to it. Needs internet
#          -> run on a login node:  snakemake ... --cores 4 auto_ref_screen ----------
rule auto_ref_screen:
    input:
        kraken = f"{OUT}/06_kraken2/top_species_overall.tsv",
        fqs = expand(f"{OUT}/unaligned_reads/{{s}}.all_unaligned.fq.gz", s=SAMPLES),
    output:
        tsv  = f"{OUT}/08_top_organism/top_organism_alignment.tsv",
        pick = f"{OUT}/08_top_organism/picked_organism.txt",
    params: outdir=OUT, host=HOST, minf=AUTO_MINF
    threads: THREADS
    shell:
        "python " + SCRIPTS + "/auto_ref_screen.py --out-dir {params.outdir} "
        "--host {params.host:q} --min-fraction {params.minf} --threads {threads} "
        "--out-tsv {output.tsv} --out-pick {output.pick} --fastqs {input.fqs}"

# ---- 06: Kraken2 taxonomic classification (optional) ------------------------
rule kraken:
    input: f"{OUT}/unaligned_reads/{{sample}}.unaligned.fq.gz"
    output: f"{OUT}/06_kraken2/reports/{{sample}}.kreport"
    params: db=KDB
    threads: THREADS
    shell:
        "kraken2 --db {params.db} --threads {threads} --gzip-compressed --use-names "
        "--report {output} --output /dev/null {input}"

rule aggregate_kraken:
    input: expand(f"{OUT}/06_kraken2/reports/{{s}}.kreport", s=SAMPLES)
    output:
        summ = f"{OUT}/06_kraken2/kraken_summary.tsv",
        spp  = f"{OUT}/06_kraken2/top_species_overall.tsv",
    shell:
        "python " + SCRIPTS + "/aggregate_kraken.py "
        "--out-summary {output.summ} --out-species {output.spp} --in {input}"

# ---- 05: align unaligned reads to the suspect genome (optional) -------------
rule contaminant_index:
    input: REF if REF else []
    output: f"{OUT}/05_contaminant_genome_screen/ref/contaminant.1.bt2"
    params: pfx=f"{OUT}/05_contaminant_genome_screen/ref/contaminant"
    threads: THREADS
    shell: "bowtie2-build --threads {threads} {input} {params.pfx} > /dev/null 2>&1"

rule contaminant_align:
    input:
        fq  = f"{OUT}/unaligned_reads/{{sample}}.all_unaligned.fq.gz",
        idx = f"{OUT}/05_contaminant_genome_screen/ref/contaminant.1.bt2",
    output: f"{OUT}/05_contaminant_genome_screen/bt2_logs/{{sample}}.bt2.log"
    params: pfx=f"{OUT}/05_contaminant_genome_screen/ref/contaminant"
    threads: THREADS
    shell:
        "bowtie2 -x {params.pfx} -U {input.fq} -p {threads} --very-sensitive-local "
        "2> {output} > /dev/null"

rule aggregate_contaminant:
    input: expand(f"{OUT}/05_contaminant_genome_screen/bt2_logs/{{s}}.bt2.log", s=SAMPLES)
    output: f"{OUT}/05_contaminant_genome_screen/contaminant_alignment.tsv"
    shell: "python " + SCRIPTS + "/aggregate_contaminant.py --out {output} --in {input}"

# ---- 07: map unaligned reads to a CUSTOM fasta (e.g. a transgene/vector
#          construct), with a per-sequence breakdown via samtools idxstats -------
rule custom_index:
    input: CUSTOM if CUSTOM else []
    output: f"{OUT}/07_custom_sequences/ref/custom.1.bt2"
    params: pfx=f"{OUT}/07_custom_sequences/ref/custom"
    threads: THREADS
    shell: "bowtie2-build --threads {threads} {input} {params.pfx} > /dev/null 2>&1"

rule custom_align:
    input:
        fq  = f"{OUT}/unaligned_reads/{{sample}}.all_unaligned.fq.gz",
        idx = f"{OUT}/07_custom_sequences/ref/custom.1.bt2",
    output: f"{OUT}/07_custom_sequences/idxstats/{{sample}}.idxstats"
    params: pfx=f"{OUT}/07_custom_sequences/ref/custom",
            bam=f"{OUT}/07_custom_sequences/idxstats/{{sample}}.bam"
    threads: THREADS
    shell:
        "bowtie2 -x {params.pfx} -U {input.fq} -p {threads} --very-sensitive-local "
        "2> {output}.log | samtools sort -@ {threads} -o {params.bam} - ; "
        "samtools index {params.bam} ; samtools idxstats {params.bam} > {output} ; "
        "rm -f {params.bam} {params.bam}.bai"

rule aggregate_custom:
    input: expand(f"{OUT}/07_custom_sequences/idxstats/{{s}}.idxstats", s=SAMPLES)
    output:
        mp = f"{OUT}/07_custom_sequences/custom_mapping.tsv",
        ps = f"{OUT}/07_custom_sequences/custom_per_sequence.tsv",
    shell:
        "python " + SCRIPTS + "/aggregate_custom.py "
        "--out-mapping {output.mp} --out-perseq {output.ps} --in {input}"
