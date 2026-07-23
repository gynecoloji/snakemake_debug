# Contributing

Thanks for your interest in improving this workflow! This guide covers how to
report problems, set up a development environment, run the checks, and propose
changes.

By participating, you agree to abide by the [Code of Conduct](CODE_OF_CONDUCT.md).

## Reporting issues

Open a [GitHub issue](https://github.com/gynecoloji/snakemake_debug/issues) with:

- what you ran (the exact `snakemake` command / target),
- what you expected vs. what happened,
- the relevant log(s) from `logs/`, and your OS + Apptainer/Snakemake version.

## Development setup

Everything — Snakemake and every tool — ships inside one Apptainer image, so you
do not install a per-rule software stack yourself:

```bash
git clone https://github.com/gynecoloji/snakemake_debug.git
cd snakemake_debug

# build the container (login node, needs internet); add --no-kraken to skip the DB
module load apptainer      # or your cluster's equivalent
./setup.sh --no-kraken
```

Your input is a directory of raw per-sample aligner output (`.sam`/`.bam` that
still contains unmapped reads). Configuration is validated against
[`workflow/schemas/config.schema.yaml`](workflow/schemas/config.schema.yaml) on
every run — that schema is the single source of truth for parameters.

## Running the workflow

```bash
# dry run (validates config + builds the DAG) — no container needed for -n if
# you have Snakemake on PATH; otherwise run it inside the image:
apptainer exec containers/debug_tools.sif \
    snakemake -s workflow/Snakefile --configfile config/config.yaml -n

# full run
apptainer exec --cleanenv -B "$PWD" -B <project_dir> containers/debug_tools.sif \
    snakemake -s workflow/Snakefile --configfile config/config.yaml --cores 8
```

Run a subset with a stage aggregator target: `unaligned_all`, `blast_all`
(login node), `screens_all`, `multimap_all`, or `report_all`.

## Checks

CI (`.github/workflows/ci.yml`) builds the full DAG via a dry run over the
stubbed [`.test/`](.test/) inputs (which also validates the config against the
schema) and lints the workflow. Before opening a PR, please make sure a dry run
passes:

```bash
snakemake -s workflow/Snakefile -d .test -n --cores 1
```

## Commit messages & releases

This repo uses **[Conventional Commits](https://www.conventionalcommits.org)**
and [release-please](https://github.com/googleapis/release-please) to automate
versioning, the changelog, and releases. Prefix your commits:

| Prefix | Effect |
|---|---|
| `feat: …` | new feature → minor version bump, listed under *Added* |
| `fix: …` | bug fix → patch bump, under *Fixed* |
| `feat!: …` or a `BREAKING CHANGE:` footer | major bump |
| `docs:` / `refactor:` / `perf:` / `chore:` | grouped in the changelog; no release on their own |

On merge to `main`, release-please opens/updates a "release PR"; merging that PR
tags the version, publishes a GitHub Release, and updates `CHANGELOG.md`,
`CITATION.cff`, and `version.txt`. You do **not** edit the changelog or version
numbers by hand.

## Pull requests

1. Branch from `main` and make your change.
2. Run `snakemake -s workflow/Snakefile -d .test -n --cores 1`.
3. Open a PR with a clear, Conventional-Commit-style title and description.
