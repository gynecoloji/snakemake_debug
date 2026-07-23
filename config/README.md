# Configuration

This workflow is configured through a single file in this directory:

- `config.yaml` — all workflow parameters (validated at the start of every run
  against [`workflow/schemas/config.schema.yaml`](../workflow/schemas/config.schema.yaml)).

Unlike a sample-sheet workflow, samples are **discovered automatically** by
globbing `input.align_dir` for files ending in `input.align_suffix` — there is no
`samples.csv`. The part of each filename before the suffix becomes the sample
label used throughout the reports.

All relative paths resolve under `project_dir` unless they are absolute.

## Input

Point `input.align_dir` at a directory of per-sample **raw** aligner output
(SAM/BAM that still contains the unmapped reads — not a post-filtered
unique-mapped BAM, which has the reads to diagnose removed):

```
<project_dir>/results/hisat2/     # <- input.align_dir
├── sampleA.sam                   # sample label = "sampleA"
├── sampleB.sam
└── sampleC.sam
```

```yaml
input:
  align_dir:    "results/hisat2"  # relative to project_dir (or absolute)
  align_suffix: ".sam"            # ".sam" or ".bam" — must match your files
  unmapped_flag: 4                # SAM flag for unmapped reads (almost always 4)
```

## Parameters

| key | meaning |
|---|---|
| `project_dir` | absolute path; all relative paths resolve against it |
| `input.align_dir` / `align_suffix` / `unmapped_flag` | raw aligner output dir, suffix, and the SAM unmapped flag |
| `output_dir` | where results go (under `project_dir`, or absolute) |
| `host_label` | host organism name (reporting only) |
| `sampling.unmapped_cap` / `records_cap` | bound reads sampled / records scanned per file (sampled steps 02/03/04/06) |
| `kraken.db` | Kraken2 DB dir — empty ⇒ skip step 06 |
| `contaminant.reference_fasta` | suspect-genome FASTA — empty ⇒ skip step 05 |
| `custom.fasta` | custom construct/vector FASTA — empty ⇒ skip step 07 |
| `blast.enabled` / `random` / `n_random` / `db` / `max_target_seqs` | remote BLAST of the top (03) and random (04) unaligned reads (login node) |
| `auto_ref.enabled` / `min_fraction` | step 08 — auto-pick the top non-host organism and align to its RefSeq genome (login node) |
| `multimap.*` | step 09 — multimapped-read investigation (**spike-in only**; requires `enabled` **and** `spikein`) |
| `threads` | cores per multi-threaded rule |

Optional steps switch **off** when their path/flag is left empty. See the
top-level [`README.md`](../README.md) for the full narrative and the two-phase
(login-node setup / compute-node run) usage, and
[`workflow/documentation.md`](../workflow/documentation.md) for the per-step
technical reference.

## Running

```bash
# dry run (see the DAG without executing)
apptainer exec containers/debug_tools.sif \
    snakemake -s workflow/Snakefile --configfile config/config.yaml -n

# full run
apptainer exec --cleanenv -B "$PWD" -B <project_dir> containers/debug_tools.sif \
    snakemake -s workflow/Snakefile --configfile config/config.yaml --cores 8
```

Run a subset with a stage aggregator target: `unaligned_all`, `blast_all`
(login node), `screens_all`, `multimap_all`, or `report_all`.
