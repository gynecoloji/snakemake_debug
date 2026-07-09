#!/bin/bash
# ============================================================================
# One-time setup — run on a LOGIN NODE (needs internet).
# Builds the container and downloads the Kraken2 DB (and optionally a
# contaminant reference genome). Idempotent: skips anything already present.
#
#   ./setup.sh                          # build container + Kraken2 standard-8 DB
#   ./setup.sh --kraken-db-size 16      # use the 16 GB DB instead
#   ./setup.sh --no-kraken              # container only
#   ./setup.sh --ref-accession GCF_900476065.1   # also fetch a genome for step 05
#   ./setup.sh --drosophila             # also fetch the Drosophila melanogaster dm6 genome
# ============================================================================
set -euo pipefail
HERE=$(cd "$(dirname "$0")" && pwd)
MODULE=${APPTAINER_MODULE:-apptainer/1.4.1-u7wa}
KSIZE=8; ACC=""; DMEL=0
while [[ $# -gt 0 ]]; do case "$1" in
  --kraken-db-size) KSIZE="$2"; shift 2;;
  --no-kraken)      KSIZE=0; shift;;
  --ref-accession)  ACC="$2"; shift 2;;
  --drosophila)     DMEL=1; shift;;
  *) echo "unknown option: $1"; exit 1;;
esac; done

module load "$MODULE" 2>/dev/null || true
command -v apptainer >/dev/null || { echo "ERROR: apptainer not on PATH (module load it)"; exit 1; }
export APPTAINER_CACHEDIR=${APPTAINER_CACHEDIR:-$HERE/.cache}
export APPTAINER_TMPDIR=${APPTAINER_TMPDIR:-$HERE/.tmp}
mkdir -p "$APPTAINER_CACHEDIR" "$APPTAINER_TMPDIR" "$HERE/refs"

# 1) container (unprivileged, proot-backed fakeroot) --------------------------
SIF=$HERE/containers/debug_tools.sif
if [[ ! -f $SIF ]]; then
  if [[ ! -x $HERE/containers/proot ]]; then
    echo "[setup] fetching static proot"
    curl -sL -o "$HERE/containers/proot" https://proot.gitlab.io/proot/bin/proot
    chmod +x "$HERE/containers/proot"
  fi
  export PATH="$HERE/containers:$PATH"
  echo "[setup] building $SIF ..."
  apptainer build --fakeroot "$SIF" "$HERE/containers/debug_tools.def"
else
  echo "[setup] container present: $SIF"
fi

# 2) Kraken2 DB ---------------------------------------------------------------
if [[ $KSIZE != 0 ]]; then
  DB=$HERE/refs/kraken_db; mkdir -p "$DB"
  if [[ ! -f $DB/hash.k2d ]]; then
    KEY=$(curl -s "https://genome-idx.s3.amazonaws.com/?list-type=2&prefix=kraken/k2_standard_$(printf %02d "$KSIZE")gb" \
          | grep -oE "kraken/k2_standard_$(printf %02d "$KSIZE")gb_[0-9]+\.tar\.gz" | sort | tail -1)
    [[ -n $KEY ]] || { echo "ERROR: could not find a k2_standard_${KSIZE}gb DB"; exit 1; }
    echo "[setup] downloading Kraken2 DB: $KEY"
    curl -sL "https://genome-idx.s3.amazonaws.com/$KEY" -o "$DB/db.tar.gz"
    tar -xzf "$DB/db.tar.gz" -C "$DB" && rm "$DB/db.tar.gz"
  else
    echo "[setup] Kraken2 DB present: $DB"
  fi
  echo "[setup]   -> set  kraken.db: \"$DB\"  in config.yaml"
fi

# 3) optional contaminant reference genome ------------------------------------
if [[ -n $ACC ]]; then
  RF=$HERE/refs/$ACC.fa
  if [[ ! -f $RF ]]; then
    echo "[setup] fetching genome $ACC from NCBI datasets"
    curl -sL "https://api.ncbi.nlm.nih.gov/datasets/v2alpha/genome/accession/$ACC/download?include_annotation_type=GENOME_FASTA" -o "$HERE/refs/$ACC.zip"
    python3 -c "import zipfile; zipfile.ZipFile('$HERE/refs/$ACC.zip').extractall('$HERE/refs/dl_$ACC')"
    cat "$HERE"/refs/dl_$ACC/ncbi_dataset/data/"$ACC"/*.fna > "$RF"
  fi
  echo "[setup]   -> set  contaminant.reference_fasta: \"$RF\"  in config.yaml"
fi

# 4) convenience: Drosophila melanogaster dm6 genome (UCSC) ------------------
if [[ $DMEL == 1 ]]; then
  DM=$HERE/refs/dm6.fa
  if [[ ! -f $DM ]]; then
    echo "[setup] fetching Drosophila melanogaster dm6 genome (UCSC)"
    curl -sL "https://hgdownload.soe.ucsc.edu/goldenPath/dm6/bigZips/dm6.fa.gz" -o "$DM.gz"
    gunzip -f "$DM.gz"
  else
    echo "[setup] dm6 genome present: $DM"
  fi
  echo "[setup]   -> set  contaminant.reference_fasta: \"$DM\"  in config.yaml"
fi
echo "[setup] done."
