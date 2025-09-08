#!/bin/bash
set -euo pipefail

PROJECT="${PROJECT:-$(cd "$(dirname "$0")/.." && pwd)}"
RAW="$PROJECT/data/raw"
PROC="$PROJECT/data/processed"

mkdir -p "$RAW" "$PROC"

# If subsampled FASTQs already exist, skip download/subsampling
if compgen -G "$PROC/*_1.fastq.subsampled" > /dev/null; then
  echo "Detected existing subsampled FASTQs in $PROC; skipping fetch/subsample."
  exit 0
fi

# --- Fetch SRX26680000 (example; left commented) ---
# echo "Downloading SRX26680000 FASTQs..."
# fastq-dl SRX26680000 -o "$RAW" --split-files --gzip

# Optional: create small subsamples (1M pairs) to data/processed/
# Requires seqtk
if command -v seqtk >/dev/null 2>&1; then
  echo "Creating 1M-pair subsamples in $PROC ..."
  for R1 in "$RAW"/*/*_1.fastq.gz "$RAW"/*/*R1*.fastq.gz 2>/dev/null; do
    [[ -f "$R1" ]] || continue
    R2="${R1/_1.fastq.gz/_2.fastq.gz}"
    [[ -f "$R2" ]] || R2="${R1/R1/R2}"
    base=$(basename "$R1" | sed 's/_1\.fastq\.gz//; s/_R1\.fastq\.gz//')
    seqtk sample -s100 "$R1" 1000000 > "$PROC/${base}_1.fastq.subsampled"
    seqtk sample -s100 "$R2" 1000000 > "$PROC/${base}_2.fastq.subsampled"
  done
  echo "Processed subsamples in: $PROC"
else
  echo "NOTE: seqtk not found; skipping subsampling. Provide subsamples in $PROC manually." 
fi
