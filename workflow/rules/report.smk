# Synthesis outputs: figures (plots), the Markdown synthesis (00_SUMMARY.md),
# and the self-contained interactive HTML report (one section per step). All
# three consume the AGG per-step summary tables, so requesting report_all pulls
# the whole diagnostic pipeline. Shared config, constants and target lists live
# in common.smk.


# Aggregate target for the report stage (plots + summary + HTML report). Because
# these depend on the AGG tables, this effectively builds the whole pipeline
# except the internet-only extras (blast / random / auto). Run it with:
#   snakemake --cores N report_all
rule report_all:
    input:
        REPORT_TARGETS,


# ---- plots + summary --------------------------------------------------------
rule plots:
    input: AGG
    output: touch(f"{OUT}/plots/.done")
    params: outdir=OUT
    shell: "MPLCONFIGDIR=/tmp/mpl python " + SCRIPTS + "/make_plots.py --out-dir {params.outdir}"

rule summary:
    input: AGG + [f"{OUT}/plots/.done"]
    output: f"{OUT}/00_SUMMARY.md"
    params: outdir=OUT, host=HOST
    shell:
        "python " + SCRIPTS + "/write_summary.py --out-dir {params.outdir} "
        "--host {params.host:q} > {output}"

# ---- interactive self-contained HTML report (one section per step) ----------
rule report:
    input: AGG + ([f"{OUT}/03_blast_top_sequences/blast_best_hits.tsv"] if BLAST_ON else [])
    output: f"{OUT}/report.html"
    params: outdir=OUT, host=HOST,
            clabel=(os.path.basename(REF) if REF else "contaminant genome")
    shell:
        "python " + SCRIPTS + "/build_report.py --out-dir {params.outdir} "
        "--out {output} --host {params.host:q} --contaminant-label {params.clabel:q}"
