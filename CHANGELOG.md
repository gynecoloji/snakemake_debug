# Changelog

All notable changes to this project are documented in this file.

The format is based on [Keep a Changelog](https://keepachangelog.com/en/1.1.0/),
and this project adheres to [Semantic Versioning](https://semver.org/spec/v2.0.0.html).
Going forward this file is maintained automatically by
[release-please](https://github.com/googleapis/release-please) from Conventional
Commit messages — do not edit it by hand.

## [0.2.3](https://github.com/gynecoloji/snakemake_debug/compare/v0.2.2...v0.2.3) (2026-07-30)


### Documentation

* embed rule graph image in README.md ([2eb7086](https://github.com/gynecoloji/snakemake_debug/commit/2eb7086564fdb98202e450897f6ecc6ac719d7a6))

## [0.2.2](https://github.com/gynecoloji/snakemake_debug/compare/v0.2.1...v0.2.2) (2026-07-30)


### Documentation

* snakevision tube map ([f9a985d](https://github.com/gynecoloji/snakemake_debug/commit/f9a985da8344825cb3396a02dfce6e384103a2ff))

## [0.2.1](https://github.com/gynecoloji/snakemake_debug/compare/v0.2.0...v0.2.1) (2026-07-30)


### Documentation

* add Zenodo DOI ([8d86df3](https://github.com/gynecoloji/snakemake_debug/commit/8d86df34db73e59d097639025cf9b179a2452ab5))

## [0.2.0](https://github.com/gynecoloji/snakemake_debug/compare/v0.1.0...v0.2.0) (2026-07-29)


### Added

* standard Snakemake layout, Apptainer container, and catalog template ([81dee9f](https://github.com/gynecoloji/snakemake_debug/commit/81dee9f21241fc072145c2484f0f1049f5648648))

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
