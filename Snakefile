# ============================================================================
# Low-alignment / contamination diagnostics workflow
# Runs entirely inside debug_tools.sif:
#     apptainer exec --cleanenv -B <project> debug_tools.sif \
#         snakemake -s Snakefile --configfile config.yaml --cores 8
# See README.md for the two-phase (setup on login node / run via srun) usage.
# ============================================================================
import os

configfile: "config.yaml"

PROJ    = config["project_dir"]
def A(p):  # resolve relative to project_dir
    return p if os.path.isabs(p) else os.path.normpath(os.path.join(PROJ, p))

ALIGN_DIR = A(config["input"]["align_dir"])
SUF       = config["input"]["align_suffix"]
FLAG      = config["input"]["unmapped_flag"]
OUT       = A(config["output_dir"])
HOST      = config["host_label"]
UCAP      = config["sampling"]["unmapped_cap"]
RCAP      = config["sampling"]["records_cap"]
KDB       = config["kraken"]["db"]
REF       = config["contaminant"]["reference_fasta"]
CUSTOM    = config.get("custom", {}).get("fasta", "")
BLAST_ON  = bool(config["blast"]["enabled"])
RANDOM_ON = bool(config["blast"].get("random", False))
NRAND     = int(config["blast"].get("n_random", 25))
AUTO_ON   = bool(config.get("auto_ref", {}).get("enabled", False))
AUTO_MINF = float(config.get("auto_ref", {}).get("min_fraction", 0.05))
THREADS   = config["threads"]
SCRIPTS   = os.path.join(os.path.dirname(os.path.abspath(workflow.snakefile)), "scripts")

SAMPLES, = glob_wildcards(os.path.join(ALIGN_DIR, "{sample}" + SUF))
SAMPLES  = sorted(s for s in SAMPLES if "/" not in s)
if not SAMPLES:
    raise WorkflowError(f"No '{SUF}' files found in {ALIGN_DIR}")

wildcard_constraints:
    sample = r"[^/]+"

# ---- aggregate outputs (some conditional on config) ------------------------
AGG = [
    f"{OUT}/01_alignment_summary/alignment_summary.tsv",
    f"{OUT}/02_sequence_signatures/signature_fractions.tsv",
    f"{OUT}/03_blast_top_sequences/top_unaligned.fasta",
]
if KDB: AGG.append(f"{OUT}/06_kraken2/kraken_summary.tsv")
if REF: AGG.append(f"{OUT}/05_contaminant_genome_screen/contaminant_alignment.tsv")
if CUSTOM: AGG.append(f"{OUT}/07_custom_sequences/custom_mapping.tsv")

FINAL = AGG + [f"{OUT}/plots/.done", f"{OUT}/00_SUMMARY.md", f"{OUT}/report.html"]
if BLAST_ON:  FINAL.append(f"{OUT}/03_blast_top_sequences/blast_best_hits.tsv")
if RANDOM_ON: FINAL.append(f"{OUT}/04_blast_random/random_blast_species.tsv")
if AUTO_ON:   FINAL.append(f"{OUT}/08_top_organism/top_organism_alignment.tsv")

rule all:
    input: FINAL

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

rule pool_top_sequences:
    input: expand(f"{OUT}/unaligned_reads/{{s}}.topseq.tsv", s=SAMPLES)
    output: f"{OUT}/03_blast_top_sequences/top_unaligned.fasta"
    params: nper=6, ntot=20
    shell:
        "python " + SCRIPTS + "/pool_top_sequences.py --n-per-sample {params.nper} "
        "--n-total {params.ntot} --out {output} --in {input}"

# ---- 03: remote BLAST (optional; needs internet -> run on a login node) -----
rule blast:
    input: f"{OUT}/03_blast_top_sequences/top_unaligned.fasta"
    output: f"{OUT}/03_blast_top_sequences/blast_best_hits.tsv"
    params: db=config["blast"]["db"], mts=config["blast"]["max_target_seqs"]
    shell:
        "blastn -query {input} -db {params.db} -remote "
        "-outfmt '6 qseqid pident length evalue bitscore stitle' -max_target_seqs {params.mts} "
        "> {output}.raw 2> {output}.log ; "
        "{{ printf 'query\\tpident\\tlength\\tevalue\\tbitscore\\tsubject\\n'; "
        "awk -F'\\t' '!seen[$1]++' {output}.raw ; }} > {output}"

# ---- 04: random unaligned reads -> remote BLAST (find species Kraken2/motifs miss)
# Random (unbiased) sampling, unlike pool_top_sequences (most-duplicated). Needs
# internet -> run on a login node, same as the `blast` rule above.
rule pool_random_sequences:
    input: expand(f"{OUT}/unaligned_reads/{{s}}.unaligned.fq.gz", s=SAMPLES)
    output: f"{OUT}/04_blast_random/random_unaligned.fasta"
    params: n=NRAND, seed=100
    shell:
        "python " + SCRIPTS + "/pool_random_sequences.py --n-per-sample {params.n} "
        "--seed {params.seed} --out {output} --in {input}"

rule blast_random:
    input: f"{OUT}/04_blast_random/random_unaligned.fasta"
    output:
        spp  = f"{OUT}/04_blast_random/random_blast_species.tsv",
        hits = f"{OUT}/04_blast_random/random_blast_best_hits.tsv",
    params: db=config["blast"]["db"], mts=config["blast"]["max_target_seqs"]
    shell:
        # remote BLAST is flaky (NCBI throttles/drops); retry with backoff, then
        # always summarize whatever came back so a bad connection degrades to an
        # empty/partial table instead of failing the whole run.
        "for i in 1 2 3 4 5; do "
        "  blastn -query {input} -db {params.db} -remote "
        "    -outfmt '6 qseqid pident length evalue bitscore stitle' -max_target_seqs {params.mts} "
        "    > {output.hits}.raw 2> {output.hits}.log && break "
        "  || {{ echo \"[blast_random] attempt $i failed; retry in 45s\" >> {output.hits}.log; sleep 45; }}; "
        "done; "
        "python " + SCRIPTS + "/summarize_random_blast.py --raw {output.hits}.raw "
        "--out-species {output.spp} --out-hits {output.hits}"

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

# ---- plots + summary --------------------------------------------------------
rule plots:
    input: AGG
    output: touch(f"{OUT}/plots/.done")
    params: outdir=OUT
    shell: "MPLCONFIGDIR=/tmp/mpl python " + SCRIPTS + "/make_plots.py --out-dir {params.outdir}"

rule summary:
    input: AGG + [f"{OUT}/plots/.done"]
    output: f"{OUT}/00_SUMMARY.md"
    params: outdir=OUT, host=HOST
    shell:
        "python " + SCRIPTS + "/write_summary.py --out-dir {params.outdir} "
        "--host {params.host:q} > {output}"

# ---- interactive self-contained HTML report (one section per step) ----------
rule report:
    input: AGG + ([f"{OUT}/03_blast_top_sequences/blast_best_hits.tsv"] if BLAST_ON else [])
    output: f"{OUT}/report.html"
    params: outdir=OUT, host=HOST,
            clabel=(os.path.basename(REF) if REF else "contaminant genome")
    shell:
        "python " + SCRIPTS + "/build_report.py --out-dir {params.outdir} "
        "--out {output} --host {params.host:q} --contaminant-label {params.clabel:q}"
