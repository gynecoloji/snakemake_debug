# Technical documentation

Step-by-step documentation of the low-alignment / contamination diagnostics
workflow. For installation, container usage, and the full narrative, see the
top-level [`README.md`](../README.md); for every configuration parameter, see
[`config/README.md`](../config/README.md) and the schema
[`workflow/schemas/config.schema.yaml`](schemas/config.schema.yaml).

## Layout

This is a standard-layout Snakemake workflow. Everything — Snakemake and every
tool — runs inside the single Apptainer image `containers/debug_tools.sif`, so
there are no per-rule conda environments.

```
workflow/Snakefile           entry point: min_version, configfile, include: rules, rule all
workflow/rules/common.smk    shared setup: config validation, constants, sample discovery, target lists
workflow/rules/unaligned.smk 01 alignment summary + unaligned-read extraction + 02 signatures
workflow/rules/blast.smk     03 top-sequence BLAST + 04 random-read BLAST (login node)
workflow/rules/screens.smk   05 contaminant + 06 Kraken2 + 07 custom + 08 top-organism screens
workflow/rules/multimap.smk  09 multimapped-read investigation (spike-in only)
workflow/rules/report.smk    plots + 00_SUMMARY.md + report.html
workflow/scripts/            argparse-based analysis scripts invoked by the rules
```

`common.smk` is included first, so every constant and helper it defines
(`OUT`, `SAMPLES`, `AGG`, `FINAL`, the `MM_*` flags, …) is visible to the rules
in the other modules. It validates `config` against the schema before deriving
any values.

## Overview

A single `snakemake -s workflow/Snakefile --configfile config/config.yaml` run
builds the unaligned/multimapped diagnostics and the synthesis report in
dependency order (one unified DAG). The default target (`rule all`) resolves the
conditional `FINAL` list; each stage also has an aggregator target for subset
runs.

## Steps

| step | rule(s) | output | question | tool |
|---|---|---|---|---|
| 01 | `alignment_summary` | `01_alignment_summary/alignment_summary.tsv` | which samples are low, and by how much? | samtools flagstat |
| — | `extract_unaligned`, `extract_all_unaligned` | `unaligned_reads/` | pull the unmapped reads (bounded sample + full set) | samtools |
| 02 | `aggregate_signatures` | `02_sequence_signatures/signature_fractions.tsv` | foreign or host? (16S / adapter / polyA / Alu) | motif screen |
| 03 | `pool_top_sequences`, `blast` | `03_blast_top_sequences/` | identity of the most-duplicated unaligned seqs | blastn -remote |
| 04 | `pool_random_sequences`, `blast_random` | `04_blast_random/` | species Kraken2 / motifs miss (random reads) | blastn -remote |
| 05 | `contaminant_index`, `contaminant_align`, `aggregate_contaminant` | `05_contaminant_genome_screen/` | quantify vs a suspect genome | bowtie2 |
| 06 | `kraken`, `aggregate_kraken` | `06_kraken2/` | full taxonomic breakdown | Kraken2 |
| 07 | `custom_index`, `custom_align`, `aggregate_custom` | `07_custom_sequences/` | map to a custom construct (per sequence) | bowtie2 + idxstats |
| 08 | `auto_ref_screen` | `08_top_organism/` | auto-pick top non-host organism → align to its RefSeq genome | pick + NCBI + bowtie2 |
| 09 | `multimap_analysis`, `aggregate_multimap`, `crossgenome_analysis`, `aggregate_crossgenome` | `09_multimapped_reads/` | **(spike-in only)** host vs spike-in split, per-genome multimap rate, top loci, % mapping to both genomes | samtools + NH/XS/MAPQ; bowtie2 |
| — | `plots`, `summary`, `report` | `plots/`, `00_SUMMARY.md`, `report.html` | figures + synthesis + interactive report | matplotlib / Python |

Steps 03/04/05/06/07/08/09 are optional and enabled by their config fields.
Steps 03, 04 and 08 need internet, so run them on a login node; step 09 runs
offline and is spike-in only (requires `multimap.enabled` **and**
`multimap.spikein`). See the top-level README for the two-phase workflow.

## Stage aggregator targets

Each rule module defines a convenience `*_all` target that builds just that
stage's outputs (respecting the config gates), mirroring the standard-layout
`atacseq_all` / `qc_all` idiom:

- `unaligned_all` — 01 + 02
- `blast_all` — 03 + 04 (login node)
- `screens_all` — 05 + 06 + 07 + 08 (08 needs a login node)
- `multimap_all` — 09
- `report_all` — plots + summary + report (pulls the whole pipeline via `AGG`)

The canonical target remains `rule all` (`input: FINAL`, `default_target: True`).
