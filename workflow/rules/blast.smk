# Remote BLAST identification: pool the most-duplicated unaligned sequences (03)
# and a random unbiased sample (04), then BLAST both against NCBI. These rules
# need internet, so run them on a login node. Shared config, constants and
# target lists live in common.smk.


# Aggregate target for the BLAST stage. Run it alone (login node) with:
#   snakemake --cores N blast_all
rule blast_all:
    input:
        BLAST_TARGETS,


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
