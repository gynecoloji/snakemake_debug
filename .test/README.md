# Catalog test case (`.test/`)

This directory lets the [Snakemake Workflow Catalog](https://snakemake.github.io/snakemake-workflow-catalog/)
render the workflow's **tube map** (rule graph) and lets CI check that the DAG
builds and the config validates. Both run the workflow with `-d .test`:

```bash
# rule graph (what the catalog renders; see ../images/rulegraph.svg)
snakemake -s workflow/Snakefile -d .test --rulegraph -c 1

# dry run (what CI runs: DAG builds + config validates against the schema)
snakemake -s workflow/Snakefile -d .test -n -c 1
```

Neither command **executes** any rule, so the inputs here are tiny/empty
placeholders — nothing reads their contents:

- `config/config.yaml` — a copy of the config with **every optional step turned
  on**, so all 25 rules appear in the graph. Paths resolve under `.test/`.
- `align_stub/sample1.sam` — one placeholder alignment (a bare SAM header) so the
  sample glob is non-empty.
- `ref/contaminant_stub.fa`, `ref/custom_stub.fa` — empty placeholder FASTAs (the
  inputs of the two `bowtie2-build` rules).

`kraken.db`, `multimap.host_index` and `multimap.spikein_index` are passed as
**parameters** (not tracked inputs), so they only need to be non-empty strings.

**This is not an end-to-end integration test.** Turning it into one would require
a real miniature combined host+spike-in genome and reads that actually run
through samtools/bowtie2/BLAST/Kraken2 inside the container — a larger, separate
effort.
