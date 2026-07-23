# Shared setup for every rule module (unaligned / blast / screens / multimap /
# report): imports, config validation, the project-dir path helper, all
# config-derived constants, sample discovery, wildcard constraints, and the
# conditional aggregate/final target lists. Included first by workflow/Snakefile,
# so every name defined here is visible to the rules in the other included files.
import os
from snakemake.utils import validate

# ── Config validation ───────────────────────────────────────────────────
# Validate against workflow/schemas/config.schema.yaml (structure + types). The
# catalog renders the parameter table from the same schema. Path is relative to
# this file (workflow/rules/).
validate(config, "../schemas/config.schema.yaml")

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
# Step 09 (multimapped-read investigation) is SPIKE-IN ONLY: it reports the host vs spike-in
# split + the cross-genome "both genomes" check, which only make sense for a combined host+spike-in
# alignment. It requires BOTH multimap.enabled AND multimap.spikein; enabled without spikein is
# skipped with a warning. Steps 01-08 are library-agnostic (spike-in or normal sequencing).
MM_ENABLED = bool(config.get("multimap", {}).get("enabled", False))
MM_SPIKE   = bool(config.get("multimap", {}).get("spikein", False))
MM_ON      = MM_ENABLED and MM_SPIKE
if MM_ENABLED and not MM_SPIKE:
    logger.warning("[config] step 09 is spike-in only; multimap.enabled is set but multimap.spikein "
                   "is not -> skipping step 09 (steps 01-08 still run).")
MM_PREFIX = config.get("multimap", {}).get("spikein_prefix", "spikein_")
MM_METHOD = config.get("multimap", {}).get("method", "auto")
MM_XSSTR  = bool(config.get("multimap", {}).get("xs_strict", False))
MM_MAPQ   = int(config.get("multimap", {}).get("mapq_max", 10))
# Step 09 can read a DIFFERENT alignment than the unaligned-analysis steps (01-08): point it at
# the raw combined host+spike-in BAM (multimappers intact). Empty -> inherit input.align_dir/suffix.
MM_ALIGN_DIR = A(config.get("multimap", {}).get("align_dir", "") or config["input"]["align_dir"])
MM_SUF       = config.get("multimap", {}).get("align_suffix", "") or SUF
# Step 09 cross-genome check (spike-in): of the multimapped reads, how many map to BOTH genomes.
# Needs SEPARATE host-only and spike-in-only bowtie2 indexes (re-aligns the reads to each).
MM_HOST_IDX = A(config.get("multimap", {}).get("host_index", "")) if config.get("multimap", {}).get("host_index", "") else ""
MM_SPK_IDX  = A(config.get("multimap", {}).get("spikein_index", "")) if config.get("multimap", {}).get("spikein_index", "") else ""
# Enable only when BOTH indexes are set — otherwise skip gracefully (a from-scratch run still
# produces the report) instead of hard-failing bowtie2 on an empty index path.
_XG_REQ    = bool(config.get("multimap", {}).get("cross_genome", False))
MM_XGENOME = _XG_REQ and bool(MM_HOST_IDX) and bool(MM_SPK_IDX)
if _XG_REQ and not MM_XGENOME:
    logger.warning("[config] multimap.cross_genome is on but host_index/spikein_index are not "
                   "both set -> skipping the cross-genome check (rest of step 09 still runs).")
MM_XFRAC   = float(config.get("multimap", {}).get("cross_fraction", 1.0))   # 1.0 = ALL reads
MM_XSEED   = int(config.get("multimap", {}).get("cross_seed", 13))
MM_XCAP    = int(config.get("multimap", {}).get("cross_cap", 0))             # 0 = no cap
MM_XMARGIN = int(config.get("multimap", {}).get("cross_codominant_margin", 5))
THREADS   = config["threads"]
# workflow.basedir is the directory of the MAIN Snakefile (workflow/), constant
# regardless of which included .smk reads it — unlike workflow.snakefile, which
# from inside an include resolves to the included file (workflow/rules/).
SCRIPTS   = os.path.join(workflow.basedir, "scripts")

SAMPLES, = glob_wildcards(os.path.join(ALIGN_DIR, "{sample}" + SUF))
SAMPLES  = sorted(s for s in SAMPLES if "/" not in s)
if not SAMPLES:
    raise WorkflowError(f"No '{SUF}' files found in {ALIGN_DIR}")

# Step 09 sample list: from its own align dir if decoupled, else the main SAMPLES.
if MM_ON and (MM_ALIGN_DIR != ALIGN_DIR or MM_SUF != SUF):
    MM_SAMPLES, = glob_wildcards(os.path.join(MM_ALIGN_DIR, "{sample}" + MM_SUF))
    MM_SAMPLES  = sorted(s for s in MM_SAMPLES if "/" not in s)
    if not MM_SAMPLES:
        raise WorkflowError(f"multimap: no '{MM_SUF}' files found in {MM_ALIGN_DIR}")
else:
    MM_SAMPLES = SAMPLES

wildcard_constraints:
    sample = r"[^/]+"

# ---- aggregate outputs (some conditional on config) ------------------------
# AGG is the set of per-step summary tables; it is also the direct input to the
# plots/summary/report rules, so it is kept as a single list here.
AGG = [
    f"{OUT}/01_alignment_summary/alignment_summary.tsv",
    f"{OUT}/02_sequence_signatures/signature_fractions.tsv",
    f"{OUT}/03_blast_top_sequences/top_unaligned.fasta",
]
if KDB: AGG.append(f"{OUT}/06_kraken2/kraken_summary.tsv")
if REF: AGG.append(f"{OUT}/05_contaminant_genome_screen/contaminant_alignment.tsv")
if CUSTOM: AGG.append(f"{OUT}/07_custom_sequences/custom_mapping.tsv")
if MM_ON: AGG.append(f"{OUT}/09_multimapped_reads/multimap_summary.tsv")   # MM_ON => spike-in
if MM_ON and MM_XGENOME: AGG.append(f"{OUT}/09_multimapped_reads/crossgenome_summary.tsv")

FINAL = AGG + [f"{OUT}/plots/.done", f"{OUT}/00_SUMMARY.md", f"{OUT}/report.html"]
if BLAST_ON:  FINAL.append(f"{OUT}/03_blast_top_sequences/blast_best_hits.tsv")
if RANDOM_ON: FINAL.append(f"{OUT}/04_blast_random/random_blast_species.tsv")
if AUTO_ON:   FINAL.append(f"{OUT}/08_top_organism/top_organism_alignment.tsv")

# ---- per-stage convenience target lists (for the *_all aggregator rules) ----
# These mirror the entries in AGG/FINAL, sliced by rule module, so a subset of
# the pipeline can be requested by target name (e.g. `snakemake screens_all`).
# The canonical target remains `rule all` (input: FINAL).
UNALIGNED_TARGETS = [
    f"{OUT}/01_alignment_summary/alignment_summary.tsv",
    f"{OUT}/02_sequence_signatures/signature_fractions.tsv",
]
BLAST_TARGETS = [f"{OUT}/03_blast_top_sequences/top_unaligned.fasta"]
if BLAST_ON:  BLAST_TARGETS.append(f"{OUT}/03_blast_top_sequences/blast_best_hits.tsv")
if RANDOM_ON: BLAST_TARGETS.append(f"{OUT}/04_blast_random/random_blast_species.tsv")
SCREENS_TARGETS = []
if REF:     SCREENS_TARGETS.append(f"{OUT}/05_contaminant_genome_screen/contaminant_alignment.tsv")
if KDB:     SCREENS_TARGETS.append(f"{OUT}/06_kraken2/kraken_summary.tsv")
if CUSTOM:  SCREENS_TARGETS.append(f"{OUT}/07_custom_sequences/custom_mapping.tsv")
if AUTO_ON: SCREENS_TARGETS.append(f"{OUT}/08_top_organism/top_organism_alignment.tsv")
MULTIMAP_TARGETS = []
if MM_ON:                MULTIMAP_TARGETS.append(f"{OUT}/09_multimapped_reads/multimap_summary.tsv")
if MM_ON and MM_XGENOME: MULTIMAP_TARGETS.append(f"{OUT}/09_multimapped_reads/crossgenome_summary.tsv")
REPORT_TARGETS = [f"{OUT}/plots/.done", f"{OUT}/00_SUMMARY.md", f"{OUT}/report.html"]
