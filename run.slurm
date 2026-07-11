#!/bin/bash
# ============================================================================
# Run the diagnostics workflow via SLURM. Everything executes INSIDE the
# container (Snakemake + all tools). Edit -A / -p for your cluster.
# Kraken2 standard-8 needs ~10 GB RAM.
#     sbatch run.slurm
# ============================================================================
#SBATCH --job-name=aln_debug
#SBATCH --cpus-per-task=16
#SBATCH --mem=48G
#SBATCH --time=08:00:00
#SBATCH --output=logs/debug_%j.out
##SBATCH --account=YOUR_ACCOUNT
##SBATCH --partition=YOUR_PARTITION
set -euo pipefail
HERE="${SLURM_SUBMIT_DIR:-$(cd "$(dirname "$0")" && pwd)}"   # under SLURM, $0 is the spooled copy; use the submit dir
module load apptainer/1.4.1-u7wa 2>/dev/null || true

# project_dir from config (simple parse, no yaml dependency)
PROJECT=$(grep -E '^project_dir:' "$HERE/config.yaml" | head -1 \
  | sed -E 's/^project_dir:[[:space:]]*//; s/[[:space:]]*#.*$//; s/^"//; s/"[[:space:]]*$//; s/[[:space:]]*$//')

apptainer exec --cleanenv \
  -B "$HERE" -B "${PROJECT:-$HERE}" \
  --env MPLCONFIGDIR=/tmp/mpl \
  "$HERE/containers/debug_tools.sif" \
  snakemake -s "$HERE/Snakefile" --configfile "$HERE/config.yaml" \
    --cores "${SLURM_CPUS_PER_TASK:-8}" --rerun-incomplete -p

echo "done -> results in the output_dir set in config.yaml"
