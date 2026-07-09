# SOP: Diagnosing a Low / Variable Alignment Rate

A reusable, tool-agnostic workflow for finding **why some sequencing samples align poorly** to their
reference. Written for RNA-seq (HISAT2/STAR) but applies to any short-read DNA/RNA alignment.

**Core idea:** the reads that *don't* align are the evidence. Localize the problem → extract the
unaligned reads → identify what they are with three independent methods → quantify.

---

## 0. Parameters — set these per project

```bash
PROJECT=/path/to/project
SAMDIR=$PROJECT/results/align        # per-sample aligned BAM/SAM
OUT=$PROJECT/debug                   # where results go (subfolders below)
HOST_REF=GRCh38                      # the reference you aligned to (human/mouse/...)
READLEN=150                          # read length (for context)
SAMPLE_CAP=300000                    # max unaligned reads to sample per file
```

Create the layout:
```bash
mkdir -p $OUT/{01_alignment_summary,unaligned_reads,02_sequence_signatures,\
03_blast_top_sequences,05_contaminant_genome_screen,06_kraken2,_containers,_logs}
```

---

## 1. Build a diagnostics environment (once)

The aligner image rarely contains the classification tools. Package them once. On HPC without
root/`subuid`, build unprivileged with a `proot`-backed `--fakeroot`:

```bash
# fetch a static proot if fakeroot is unavailable
curl -sL -o proot https://proot.gitlab.io/proot/bin/proot && chmod +x proot
export PATH=$PWD:$PATH
apptainer build --fakeroot debug_tools.sif debug_tools.def
```

`debug_tools.def`:
```
Bootstrap: docker
From: condaforge/miniforge3:latest
%post
    mamba install -y -c conda-forge -c bioconda \
        seqtk bowtie2 bwa samtools blast kraken2 krakentools entrez-direct \
        python=3.10 pandas numpy matplotlib
    mamba clean -a -y
%environment
    export PATH=/opt/conda/bin:$PATH
```
> Fallback if you can't build: `apptainer pull` ready biocontainers (e.g. `docker://staphb/kraken2`,
> `docker://ncbi/blast`) — pulling also produces runnable `.sif` files.

---

## 2. The workflow

```
              aligner summaries (aln-rate) + QC reports (fastp/FastQC)
                                 │
  [01] alignment_summary ────────┤  rank aln-rate per sample; correlate with QC
                                 │  Q: does low aln track a CONDITION/batch, or read quality?
                                 ▼
  [--] unaligned_reads   ── extract unaligned reads (SAM flag 4), sample N/file -> FASTQ
                            + tally the most frequent sequences
                                 │
  [02] sequence_signatures ──────┤  motif-screen the unaligned reads (see cheat-sheet)
                                 │  Q: are they foreign, or normal host repeats?
                                 ▼
        ┌────────────────────────┼────────────────────────┐
        ▼                        ▼                         ▼
 [03] BLAST              [06] Kraken2            [05] contaminant_genome
 top seqs vs nt          taxonomic profile       align unaligned reads to the
 = IDENTITY              = FULL BREAKDOWN         suspect genome = QUANTIFY
        └────────────────────────┼─────────────────────────┘
                                 ▼
                    [00] SUMMARY — synthesis + recommendation
```

### Step 01 — Localize
Tabulate per sample: alignment rate, duplication, adapter %, GC, Q30, read length. **Interpretation:**
- Low aln correlates with **condition/batch** → suspect contamination or a transgene/vector.
- Low aln correlates with **low Q30 / short reads** → quality/trimming problem.
- **All** samples low → wrong reference, wrong strandedness, or a broken index.

### Extract the unaligned reads (unnumbered — feeds the screens)
```bash
for bam in $SAMDIR/*.bam; do
  s=$(basename $bam .bam)
  samtools view -f 4 $bam | head -n $SAMPLE_CAP \
    | awk '{print "@"$1"\n"$10"\n+\n"$11}' | gzip > $OUT/unaligned_reads/$s.unaligned.fq.gz
  samtools view -f 4 $bam | head -n $SAMPLE_CAP | cut -f10 \
    | sort | uniq -c | sort -rn | head -30 > $OUT/unaligned_reads/$s.topseq.txt
done
```
A **clean** sample has no dominant sequence. A **problem** sample has a few sequences at thousands of copies.

### Step 02 — Signature screen (cheat-sheet)
Score each unaligned read against known motifs:

| Signature | Motif / pattern | Meaning |
|---|---|---|
| Bacterial 16S rRNA | `GGGCGTAAAGCG` | microbial contamination (e.g. Mycoplasma) |
| Illumina adapter | `AGATCGGAAGAGC` | under-trimmed adapter / short inserts |
| Poly-A / poly-T | `A{20,}` / `T{20,}` | degradation, internal priming |
| Host repeat (Alu/LINE) | `TGTAATCCCAGC`, `GGCTGAGGCAGG` | normal unplaced host reads (benign) |
| Vector / transgene | (BLAST it) | plasmid/lentiviral construct, GFP, WPRE, selection marker |

### Steps 03–06 — Identify & quantify (three independent methods)
```bash
# 03  IDENTITY: BLAST the most abundant sequences (needs internet)
blastn -remote -db nt -query $OUT/03_blast_top_sequences/top_unaligned.fasta \
  -outfmt "6 qseqid pident length evalue bitscore stitle" -max_target_seqs 5 \
  > $OUT/03_blast_top_sequences/blast_nt_results.tsv

# 06  BREAKDOWN: Kraken2 taxonomic classification (download a standard-8 DB once)
kraken2 --db $KRAKEN_DB --threads 16 --gzip-compressed --use-names \
  --report $OUT/06_kraken2/$s.kreport --output /dev/null $OUT/unaligned_reads/$s.unaligned.fq.gz

# 05  QUANTIFY: once identified, align unaligned reads to the suspect genome
#     (download via NCBI datasets API by accession), then:
bowtie2-build $REF.fa idx
bowtie2 -x idx -U $OUT/unaligned_reads/$s.unaligned.fq.gz --very-sensitive-local \
  2> $s.bt2.log > /dev/null   # parse "% overall alignment rate" from the log
```
**Why three:** BLAST answers *what is the top read*, Kraken2 answers *what's the whole mix*, genome
alignment answers *exactly how much*. Convergence across all three makes the conclusion robust — no
single tool's blind spot can mislead. (Note: bowtie2-to-genome usually quantifies higher than Kraken2,
which leaves fast-evolving coding reads "unclassified.")

---

## 3. Orchestration principles (HPC)

- **Never run heavy compute on the login node** — send extraction / Kraken2 / alignment to the scheduler
  (`srun`/`sbatch`); the login node is only for internet-dependent steps (BLAST, downloads).
- **Parallelize independent work** — build the container, download the DB, and extract reads at once;
  chain only the dependent steps.
- **Bound the input** — sampling ≤`SAMPLE_CAP` unaligned reads per file (and capping records scanned)
  keeps every step fast while remaining statistically ample for percentages.
- **Save every script + log** so the workflow re-runs reproducibly against the container.

---

## 4. Common causes this workflow distinguishes

| Finding in unaligned reads | Likely cause | Fix |
|---|---|---|
| One organism's genome (bacteria/Mycoplasma) | cell-culture / sample contamination | confirm (qPCR), re-culture clean, note confound |
| Plasmid / vector / reporter sequence | transgene construct not in reference | add construct to reference; count separately |
| Host rRNA | insufficient rRNA depletion | rRNA-deplete; is benign for DE if modest |
| Adapter / very short reads | under-trimming / short inserts | re-trim; check library prep |
| Poly-A, strong 3′ bias | RNA degradation (low RIN) | flag/exclude; check RIN |
| Reads *do* match host but unplaced | normal repeats | benign — no action |
| **All** samples low, reads match host | wrong reference / strandedness / bad index | fix config; rebuild index |

---

## 5. Adaptation checklist

- [ ] Set the **§0 parameters** (paths, host reference, read length).
- [ ] Point the extraction step at your **BAM/SAM** location and naming.
- [ ] Adjust the **§0 sample→condition** grouping (regex) for your naming scheme.
- [ ] Pick the right **Kraken2 DB** (standard-8 for a quick screen; PlusPF/standard-16 for more taxa).
- [ ] In step 05, set `$REF` to the **genome your BLAST/Kraken step identified** (fetch by accession via
      the NCBI datasets API).
- [ ] Match scheduler flags (`-A account`, `-p partition`) and the container runtime (`apptainer`/`singularity`)
      to your cluster.
