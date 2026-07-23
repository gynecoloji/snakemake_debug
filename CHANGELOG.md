# Changelog

All notable changes to this project are documented in this file.

The format is based on [Keep a Changelog](https://keepachangelog.com/en/1.1.0/),
and this project adheres to [Semantic Versioning](https://semver.org/spec/v2.0.0.html).
Going forward this file is maintained automatically by
[release-please](https://github.com/googleapis/release-please) from Conventional
Commit messages — do not edit it by hand.

## [0.1.0] - 2026-07-22

### Added

- Standard-layout Snakemake workflow: `workflow/Snakefile` with topical rule
  modules (`common`, `unaligned`, `blast`, `screens`, `multimap`, `report`),
  `workflow/schemas/config.schema.yaml` config validation, and `config/config.yaml`.
- Nine diagnostic steps — per-sample alignment rate, unaligned-read motif
  signatures, remote BLAST of the top and of random unaligned reads,
  contaminant- and custom-genome screens, Kraken2 classification, an
  auto-selected top-non-host-organism screen, and a spike-in multimapped-read
  investigation (host vs spike-in split, per-genome rate, top loci, and the
  cross-genome both-genomes check) — plus figures, a Markdown summary, and a
  self-contained interactive HTML report.
- Single-image Apptainer container (`containers/debug_tools.def`) bundling
  Snakemake and every tool, with `setup.sh` (build + optional Kraken2 DB /
  reference download) and `run.slurm`.
- Snakemake Workflow Catalog metadata, a `.test/` catalog case with a rendered
  `images/rulegraph.svg`, a CI dry-run + lint workflow, and release-please
  version/changelog automation.
