# snakemake_debug — low-alignment / contamination diagnostics

[![License: MIT](https://img.shields.io/badge/License-MIT-yellow.svg)](LICENSE)
[![Snakemake](https://img.shields.io/badge/Snakemake-workflow-039475.svg)](https://snakemake.github.io)
[![CI](https://github.com/gynecoloji/snakemake_debug/actions/workflows/ci.yml/badge.svg)](https://github.com/gynecoloji/snakemake_debug/actions/workflows/ci.yml)

A **reusable, containerized Snakemake workflow** that finds *why some samples align poorly*
to their reference — most often contamination (bacteria/Mycoplasma), a transgene/vector,
rRNA, adapter, or degradation. Drop it into any short-read project by editing `config/config.yaml`.

Everything — Snakemake **and** every tool (samtools, seqtk, bowtie2, BLAST+, Kraken2, …) —
lives in **one Apptainer image**, so the whole workflow runs as:

```
apptainer exec debug_tools.sif snakemake -s workflow/Snakefile ...
```

No host Python, no host Snakemake, no per-tool modules. It uses the standard
Snakemake layout: the workflow lives in `workflow/` and its config in `config/`.

---

## What it does

See **[`docs/WORKFLOW.md`](docs/WORKFLOW.md)** for the workflow diagram. Steps 02–08 take the
**unaligned** reads (from the extraction steps) and identify what they are a different way; step 09
instead investigates the **multimapped** reads. The two extraction steps (marked `—`) feed the
unaligned-read screens but are not report sections.

> **Spike-in vs normal sequencing:** steps **01–08 apply to any short-read library** (spike-in or
> normal). **Step 09 is spike-in-only** — it splits reads into host vs spike-in genome and checks
> cross-genome mapping, which only makes sense for a combined host+spike-in alignment. It is optional
> and off by default.

| step | folder | question | tool |
|---|---|---|---|
| 01 | `01_alignment_summary/` | which samples are low, and by how much? | `samtools flagstat` |
| — | `unaligned_reads/` | extract unaligned reads (bounded sample) + most-duplicated seqs | samtools |
| 02 | `02_sequence_signatures/` | are they foreign or host? (16S/adapter/polyA/Alu) | motif screen |
| 03 | `03_blast_top_sequences/` | *identity* of the top sequences | `blastn -remote` |
| 04 | `04_blast_random/random_*` | species Kraken2/motifs miss (random reads) | `blastn -remote` |
| 05 | `05_contaminant_genome_screen/` | *quantify* vs a suspect genome | bowtie2 |
| 06 | `06_kraken2/` | *full taxonomic breakdown* | Kraken2 |
| — | `unaligned_reads/*.all_unaligned.fq.gz` | extract **all** unaligned reads for the bowtie2 screens (05/07/08) | `samtools view -f 4` |
| 07 | `07_custom_sequences/` | map to a custom construct (per sequence) | bowtie2 + idxstats |
| 08 | `08_top_organism/` | auto-pick top non-host organism → align to its RefSeq genome | pick + NCBI + bowtie2 |
| 09 | `09_multimapped_reads/` | **(spike-in only)** multimapped reads: host vs spike-in split & per-genome rate, top loci, and the % mapping to **both** genomes | samtools + NH/XS/MAPQ; bowtie2 |
| — | `plots/`, `report.html`, `00_SUMMARY.md` | figures + interactive report + synthesis | matplotlib / html |

Steps 03/04/05/06/07/08/09 are **optional** — enabled by their config fields. Steps 03, 04 and 08
need internet, so run them on a **login node**; step 09 runs **offline** and is **spike-in only**
(requires `multimap.enabled` **and** `multimap.spikein`; enabling it for a normal library is skipped
with a warning — see [Step 09](#step-09-multimapped-reads-spike-in-only)).

---

## Requirements

You do **not** install Snakemake, Python, or any bioinformatics tool yourself — they all ship
inside the container. You only need:

| requirement | where | notes |
|---|---|---|
| **Apptainer** (≥ 1.1, formerly Singularity) | login **and** compute nodes | on HPC usually `module load apptainer` |
| **Internet access** | login node only | for `setup.sh` (builds container, downloads DBs) and optional remote BLAST |
| **SLURM** | optional | `run.slurm` is provided; you can also run interactively with `srun`/`apptainer` |
| **Disk** | ~1 GB image + 8–16 GB Kraken2 DB | the container + DB are built once and reused across projects |
| **RAM** | ≥ 10 GB | Kraken2 standard-8 DB needs ~10 GB; `run.slurm` requests 24 GB |

Your **input** is one directory of per-sample aligner output (`.sam`/`.bam`) that **still
contains the unmapped reads** (see [Input layout](#input-layout)).

> **No LLM or Claude Code needed to run this.** Everything — including `report.html` — is plain
> Python plus the tools in the container. `build_report.py` (Python stdlib only) and
> `make_plots.py` (matplotlib) read the step TSVs and emit a **deterministic, self-contained**
> report offline: the same inputs always produce the same HTML, on any machine. Claude Code was
> used only to *develop* the scripts, never at runtime.

---

## Quickstart

```bash
git clone <this-repo-url> snakemake_debug
cd snakemake_debug

# 1. one-time setup on a LOGIN NODE (internet): build container + Kraken2 DB
./setup.sh                                  # add --ref-accession GCF_xxx to also fetch a genome

# 2. edit config/config.yaml  (project_dir, input.align_dir, and paste the kraken.db /
#    contaminant.reference_fasta paths that setup.sh printed)

# 3. run on a compute node
sbatch run.slurm                            # edit -A/-p first
#    ...or interactively:
#    srun -c 8 --mem 24G apptainer exec --cleanenv -B "$PWD" -B <project_dir> \
#         containers/debug_tools.sif snakemake -s workflow/Snakefile --configfile config/config.yaml --cores 8
```

Results land in the `output_dir` from `config/config.yaml`; start at `00_SUMMARY.md`.

> **Tip:** do a dry run first to see the DAG without executing anything:
> ```bash
> apptainer exec containers/debug_tools.sif \
>     snakemake -s workflow/Snakefile --configfile config/config.yaml -n
> ```

---

## Input layout

Point `input.align_dir` at a directory of per-sample aligner output. One file per sample,
each named `<sample><suffix>`; the sample name becomes the label in every report.

```
<project_dir>/results/hisat2/          # <- input.align_dir
├── sampleA.sam                        # sample name = "sampleA"
├── sampleB.sam
└── sampleC.sam
```

with, in `config/config.yaml`:

```yaml
input:
  align_dir:    "results/hisat2"   # relative to project_dir (or an absolute path)
  align_suffix: ".sam"             # ".sam" or ".bam" — must match your files
  unmapped_flag: 4                 # SAM flag for unmapped reads (almost always 4)
```

> **Important:** point `align_dir` at the aligner's **raw** output. A post-filtered
> unique-mapped BAM has the unaligned reads *removed*, so there'd be nothing to diagnose.

---

## Understanding your results

Everything is written under `output_dir`. Read them in this order:

0. **`report.html`** — open this in a browser. A self-contained, interactive report with one
   section per step (alignment, signatures, BLAST, contaminant/fly screen, Kraken2, whole-library
   composition), colorblind-safe and theme-aware. Static PNG versions are in `plots/`.
1. **`00_SUMMARY.md`** — start here. A human-readable synthesis: which samples are low,
   what the signatures/BLAST/Kraken2 point at, and the likely cause.
2. **`plots/`** — figures for a quick visual scan: `alignment_rate.png`, `signature_composition.png`,
   and (when the optional steps run) `bacterial_fraction.png`, `contaminant_mapping.png`, `library_composition.png`,
   and the step-09 set (`multimap_composition.png`, `multimap_by_genome.png`, `multimap_loci.png`,
   `multimap_crossgenome.png`, `multimap_crossgenome_composition.png`).
3. **`01_alignment_summary/alignment_summary.tsv`** — mapped % per sample.
4. **`02_sequence_signatures/signature_fractions.tsv`** — per-sample fraction of unaligned reads
   matching bacterial 16S / Illumina adapter / polyA / host-repeat motifs.
5. **`03_blast_top_sequences/blast_best_hits.tsv`** *(if BLAST enabled)* — the actual identity of
   the most-duplicated unaligned sequences.
6. **`05_contaminant_genome_screen/contaminant_alignment.tsv`** *(if a reference given)* — % of
   unaligned reads that map to your suspect genome.
7. **`06_kraken2/kraken_summary.tsv`** + `top_species_overall.tsv` *(if a DB given)* — full
   taxonomic breakdown of the unaligned reads.

8. **`09_multimapped_reads/`** *(if `multimap.enabled`)* — `multimap_summary.tsv` (% multimapped,
   host vs spike-in split, per-genome multimapping rate), `top_multimap_loci.tsv`, and
   `crossgenome_summary.tsv` *(if `multimap.cross_genome`)* — % of multimapped reads mapping to
   **both** genomes. See [Step 09](#step-09-multimapped-reads-spike-in-only).

Raw per-sample unaligned reads are kept in `unaligned_reads/<sample>.unaligned.fq.gz` if you
want to investigate further by hand.

---

## The two phases (why setup is separate)

Compute nodes usually have **no internet**, so anything that downloads or queries the web
runs in phase 1 on the login node:

- **Phase 1 (login node):** `setup.sh` builds the container + downloads the Kraken2 DB / reference genome.
  Remote **BLAST** (step 03) also needs internet — run it on the login node:
  ```bash
  apptainer exec --cleanenv -B "$PWD" -B <project_dir> containers/debug_tools.sif \
      snakemake -s workflow/Snakefile --configfile config/config.yaml --cores 4 blast
  ```
  (set `blast.enabled: true` in config first)

  To also **BLAST random (unbiased) unaligned reads** — the best way to find species that
  Kraken2's DB and the motif screen miss (how the Drosophila contamination was first found) —
  set `blast.random: true` (and `blast.n_random`) and run the `blast_random` target on the login node:
  ```bash
  apptainer exec --cleanenv -B "$PWD" -B <project_dir> containers/debug_tools.sif \
      snakemake -s workflow/Snakefile --configfile config/config.yaml --cores 4 blast_random
  ```
  Results: `04_blast_random/random_blast_species.tsv` (species → read count) — also shown in `report.html`.
- **Phase 2 (compute node):** everything else runs offline in the container via `run.slurm`.

### `setup.sh` options

```bash
./setup.sh                                   # build container + Kraken2 standard-8 DB
./setup.sh --kraken-db-size 16               # use the 16 GB DB instead
./setup.sh --no-kraken                       # container only (skip step 06)
./setup.sh --ref-accession GCF_900476065.1   # also fetch a genome for step 05
./setup.sh --drosophila                      # also fetch the Drosophila melanogaster dm6 genome (step 05)
```

`setup.sh` is **idempotent** — it skips anything already present, so it's safe to re-run.
When it finishes it prints the exact `kraken.db` / `contaminant.reference_fasta` paths to paste
into `config/config.yaml`.

---

## Step 09: multimapped reads (spike-in only)

Independent of the unaligned-read screens (01–08), step 09 investigates the reads that aligned to
**multiple** locations of a **combined host + spike-in genome**. It is **spike-in only** and
**optional** — enable with both `multimap.enabled: true` **and** `multimap.spikein: true`. Enabling
it for a normal (non-spike-in) library is skipped with a warning, since the host/spike-in split and
cross-genome check have no meaning without a spike-in genome. (Steps 01–08 run for any library.)

- **Detection is aligner-aware (auto):** `NH:i > 1` (HISAT2/STAR) → `XS` present (bowtie2 — this
  reproduces bowtie2's own ">1 time" alignment rate) → `MAPQ ≤ mapq_max` (fallback when tags are
  stripped). Override with `multimap.method`.
- **Point it at the raw alignment** (multimappers not yet removed). Step 09 can read a **different**
  file than steps 01–08 via `multimap.align_dir` / `multimap.align_suffix` (empty ⇒ inherit `input.*`).
- **Outputs** (`09_multimapped_reads/`): `multimap_summary.tsv` (% multimapped, host vs spike-in
  split, and what fraction of each genome's reads multimap) and `top_multimap_loci.tsv` (the contigs
  the multimappers land on — usually chrM and repeats).

### Host vs spike-in split

Spike-in contigs are identified by a prefix (`multimap.spikein_prefix`, e.g. `spikein_`).
`samtools idxstats` gives the exact host vs spike-in split over all reads; the report shows the
split, the per-genome multimapping rate, and the top loci the multimappers land on.

### Cross-genome check — `multimap.cross_genome: true`

Of the multimapped reads, how many map to **both** genomes? These are the ambiguous reads that can
distort spike-in normalization. The combined BAM only keeps each read's best location (and exhaustive
`bowtie2 -a` is intractable on repeats), so this **re-aligns** the multimapped reads separately to
**host-only** and **spike-in-only** indexes and classifies each read/fragment as
**host-only / spike-only / both / neither**, plus **codominant** (an equally-good hit in each genome
= genuinely ambiguous). Paired-end BAMs are re-aligned paired-end (fragment-level, matching a typical
spike-in pipeline), which is auto-detected.

This is the one step-09 feature that needs a bit of setup — **single-genome bowtie2 indexes**:

```yaml
multimap:
  enabled:        true
  spikein:        true
  spikein_prefix: "spikein_"
  cross_genome:   true
  host_index:     "ref/HG38/genome"   # bowtie2 index basename, host only
  spikein_index:  "ref/DM6/genome"    # bowtie2 index basename, spike-in only
```

Build the two indexes once (inside the container), from the same FASTAs used for the combined genome:

```bash
apptainer exec containers/debug_tools.sif bowtie2-build --threads 16 host.fa    ref/HG38/genome
apptainer exec containers/debug_tools.sif bowtie2-build --threads 16 spikein.fa ref/DM6/genome
```

All multimapped reads are used by default (`cross_fraction: 1.0`, `cross_cap: 0`); set
`cross_fraction < 1` (or `cross_cap > 0`) to sample uniformly on very large inputs. Output:
`09_multimapped_reads/crossgenome_summary.tsv` (+ the `multimap_crossgenome*.png` plots and a report
card). This step is compute-heavy (it re-aligns millions of repeat reads) — give the job plenty of
cores; e.g. this project's 6 samples took ~11 min single-end / ~78 min paired-end on 32 cores.

---

## config/config.yaml reference

Every parameter is also documented, with its type, in the schema
[`workflow/schemas/config.schema.yaml`](workflow/schemas/config.schema.yaml),
which validates `config/config.yaml` at the start of every run.

| key | meaning |
|---|---|
| `project_dir` | absolute path; all relative paths resolve against it |
| `input.align_dir` | dir of **raw** aligner output (SAM/BAM that still contains unmapped reads) |
| `input.align_suffix` | `.sam` or `.bam` |
| `input.unmapped_flag` | SAM flag for unmapped reads (usually `4`) |
| `output_dir` | where results go |
| `host_label` | host organism (reporting only) |
| `sampling.unmapped_cap` / `records_cap` | bound reads sampled / records scanned per file |
| `kraken.db` | Kraken2 DB dir — empty ⇒ skip step 06 |
| `contaminant.reference_fasta` | suspect genome FASTA — empty ⇒ skip step 05 |
| `custom.fasta` | custom FASTA (transgene/vector construct) — unaligned reads mapped to it with a per-sequence breakdown (step 07); empty ⇒ skip |
| `auto_ref.enabled` | step 08 — auto-pick the top non-host organism (Kraken2 + random-BLAST), fetch its RefSeq genome, align unaligned reads to it (login node). Run: `snakemake … auto_ref_screen` |
| `multimap.enabled` + `multimap.spikein` | step 09 (**spike-in only**) — investigate multimapped reads; **both** required, enabling without `spikein` is skipped with a warning |
| `multimap.spikein_prefix` | contigs starting with this prefix are the spike-in genome |
| `multimap.method` / `xs_strict` / `mapq_max` | override the multimapper detector (`auto`/`nh`/`xs`/`mapq`) |
| `multimap.align_dir` / `align_suffix` | step-09 input override (the raw combined BAM); empty ⇒ inherit `input.*` |
| `multimap.cross_genome` | of the multimapped reads, % mapping to **both** genomes (re-aligns to single-genome indexes) |
| `multimap.host_index` / `spikein_index` | bowtie2 index basenames (host-only / spike-in-only) for the cross-genome check |
| `multimap.cross_fraction` / `cross_cap` | `1.0` / `0` ⇒ all reads; `<1` / `>0` ⇒ sample for speed |
| `blast.enabled` | run remote BLAST (login node) |
| `blast.db` / `blast.max_target_seqs` | remote BLAST database (e.g. `nt`) and hit cap |
| `threads` | cores per rule |

---

## Adapting to another project — checklist

- [ ] `./setup.sh` (once per machine; the container + DB are reusable across projects)
- [ ] set `project_dir`, `input.align_dir`, `input.align_suffix`
- [ ] set `host_label` (e.g. `Mus musculus`)
- [ ] leave `kraken.db`/`contaminant.reference_fasta` empty for a first pass; the run tells you
      the suspect organism, then `./setup.sh --ref-accession <GCF_...>` and re-run for step 05
- [ ] edit `run.slurm` `-A`/`-p` (and `debug_tools.def` tool list if you want more tools)
- [ ] extend the motif set in `workflow/scripts/extract_unaligned.py` for non-bacterial suspects (vector, mito, rRNA)
- [ ] *(optional, spike-in only)* for step 09's cross-genome check, build single-genome bowtie2
      indexes and set `multimap.host_index` / `multimap.spikein_index` (see [Step 09](#step-09-multimapped-reads-spike-in-only))

---

## Troubleshooting

| symptom | likely cause / fix |
|---|---|
| `No '.sam' files found in <dir>` | `align_dir`/`align_suffix` don't match your files; check the path resolves under `project_dir` |
| Nothing to diagnose / 0 unaligned reads | you pointed at a filtered unique-mapped BAM — use the **raw** aligner output instead |
| `apptainer: command not found` | `module load apptainer` (or your cluster's equivalent) on **both** login and compute nodes |
| BLAST step hangs or errors on a compute node | remote BLAST needs internet — run the `blast` target on the **login** node (see [two phases](#the-two-phases-why-setup-is-separate)) |
| step 09 `cross_genome` errors / no index found | build the single-genome bowtie2 indexes and set `multimap.host_index` / `spikein_index` (see [Step 09](#step-09-multimapped-reads-spike-in-only)) |
| step 09 finds 0 multimapped reads | you pointed at a filtered BAM with multimappers already removed — use the **raw** alignment (`multimap.align_dir`) |
| Kraken2 killed / OOM | the standard-8 DB needs ~10 GB RAM; raise `--mem` in `run.slurm` or use `--no-kraken` |
| container build fails with permission errors | `setup.sh` builds unprivileged via proot-backed `--fakeroot`; ensure the fetched `containers/proot` is executable |
| a run was interrupted | re-run — `run.slurm` passes `--rerun-incomplete`; Snakemake resumes from where it stopped |

---

## Files

```
snakemake_debug/                # standard Snakemake layout (workflow/ + config/)
├── config/
│   ├── config.yaml            # <- edit this
│   └── README.md              # per-parameter config reference
├── workflow/
│   ├── Snakefile              # entry point: min_version, configfile, include: rules
│   ├── rules/                 # common.smk + topical modules
│   │                          #   (unaligned / blast / screens / multimap / report)
│   ├── schemas/
│   │   └── config.schema.yaml # config validation + parameter docs
│   ├── scripts/               # generic, argparse-based analysis scripts
│   └── documentation.md       # per-step technical reference
├── setup.sh                   # phase-1: build container + download DB/ref (login node)
├── run.slurm                  # phase-2: run everything in the container (compute node)
├── containers/
│   ├── debug_tools.def        # container recipe (Snakemake + all tools)
│   └── debug_tools.sif        # built by setup.sh (git-ignored)
├── refs/                      # Kraken2 DB + reference genomes (created by setup.sh, git-ignored)
├── docs/METHOD.md             # the conceptual SOP (tool-agnostic)
├── .test/                     # catalog test case (stubbed inputs; CI dry-runs it)
├── images/rulegraph.svg       # rule graph (catalog "tube map")
├── .github/workflows/         # CI (dry-run + lint) and release-please
├── CITATION.cff               # citation metadata (release-please-managed)
├── CHANGELOG.md               # release history (release-please-managed)
├── CONTRIBUTING.md            # contribution guide
├── CODE_OF_CONDUCT.md         # Contributor Covenant
├── version.txt                # current version (release-please-managed)
├── LICENSE                    # MIT
└── results/                   # outputs (git-ignored)
```

See `docs/METHOD.md` for the reasoning behind each step and the common-causes table,
and `workflow/documentation.md` for the per-step technical reference.

---

## Rule graph

The full rule dependency graph (the Snakemake Workflow Catalog "tube map") is in
[`images/rulegraph.svg`](images/rulegraph.svg), rendered from the catalog test
case in [`.test/`](.test/):

```bash
snakemake -s workflow/Snakefile -d .test --rulegraph -c 1 | dot -Tsvg > images/rulegraph.svg
```

---

## Contributing, releases & citation

- **Contributing:** see [`CONTRIBUTING.md`](CONTRIBUTING.md) and the
  [`CODE_OF_CONDUCT.md`](CODE_OF_CONDUCT.md). CI (`.github/workflows/ci.yml`) dry-runs
  the workflow over [`.test/`](.test/) and lints it on every push/PR.
- **Releases:** commits follow [Conventional Commits](https://www.conventionalcommits.org);
  [release-please](https://github.com/googleapis/release-please) automates the version,
  [`CHANGELOG.md`](CHANGELOG.md), and GitHub Releases.
- **Citation:** citation metadata is in [`CITATION.cff`](CITATION.cff) — GitHub renders a
  "Cite this repository" button from it.

---

## License

Released under the [MIT License](LICENSE). © 2026 gynecoloji.
