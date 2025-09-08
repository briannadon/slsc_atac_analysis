#!/bin/bash
set -euo pipefail

PROJECT="${PROJECT:-$(cd "$(dirname "$0")/.." && pwd)}"
REFDIR="$PROJECT/refs"
PROC="$PROJECT/data/processed"
mkdir -p "$REFDIR"

# --- Build hg38  ---
cd "$REFDIR"

# Prefer existing downloaded reference in data/processed if available
if [ ! -f hg38.fa ] && [ -f "$PROC/hg38.fa.gz" ]; then
  echo "Using existing hg38.fa.gz from $PROC"
  ln -sf "$PROC/hg38.fa.gz" hg38.fa.gz
fi

if [ ! -f hg38.fa.bwt ]; then
  echo "Building BWA index for hg38..."
  bwa index hg38.fa
fi

if [ ! -f hg38.chrom.sizes ]; then
  echo "Creating chrom.sizes..."
  faidx hg38.fa
  cut -f1,2 hg38.fa.fai > hg38.chrom.sizes
fi

# Optional: blacklist
if [ ! -f hg38.blacklist.bed ]; then
  echo "Downloading ENCODE blacklist (hg38)..."
  curl -L -o hg38.blacklist.bed https://raw.githubusercontent.com/Boyle-Lab/Blacklist/master/lists/hg38-blacklist.v2.bed.gz
  gunzip hg38.blacklist.bed.gz
fi

echo " hg38 reference ready."
