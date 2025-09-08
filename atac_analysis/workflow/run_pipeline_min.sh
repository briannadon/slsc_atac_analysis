#!/usr/bin/env bash
# Minimal ATAC-seq pipeline (post-refs) with robust resume behavior.
# Steps: fastp → bwa mem → name-sort → fixmate → coord-sort → markdup → MACS → FRiP
set -euo pipefail

# -------------------------- config --------------------------------
PROJECT="${PROJECT:-$(cd "$(dirname "$0")/.." && pwd)}"
THREADS="${THREADS:-8}"
FORCE="${FORCE:-0}"   # 1 = redo even if outputs exist

REFDIR="$PROJECT/refs"
SUB="$PROJECT/data/processed"
WORK="$PROJECT/work"
ALIGN="$PROJECT/align"
PEAKS="$PROJECT/peaks"
QC="$PROJECT/qc"

# Default samples and autodiscovery from processed folder if not provided
SAMPLES=(SRX26680000)
if [[ -n "${SAMPLES:-}" && "$SAMPLES" != *"("*")"* ]]; then
  read -r -a SAMPLES <<< "$SAMPLES"
fi

# If user didn't override samples, attempt autodiscovery of SRR/SRX bases
if [[ "${SAMPLES[*]}" == "SRX26680000" ]]; then
  mapfile -t FOUND <<< "$(
    {
      for f in "$SUB"/*_1.fastq.subsampled; do [[ -f "$f" ]] && basename "$f" | sed 's/_1\.fastq\.subsampled$//' || true; done
      for f in "$SUB"/*_R1.sub.fastq.gz; do [[ -f "$f" ]] && basename "$f" | sed 's/_R1\.sub\.fastq\.gz$//' || true; done
    } | sort -u
  )"
  if [[ ${#FOUND[@]} -gt 0 ]]; then
    SAMPLES=("${FOUND[@]}")
    echo "Auto-detected samples: ${SAMPLES[*]}"
  fi
fi

mkdir -p "$WORK" "$ALIGN" "$PEAKS" "$QC"

# -------------------------- tool checks ---------------------------------------
need(){ command -v "$1" >/dev/null 2>&1 || { echo "ERROR: '$1' not found"; exit 1; }; }
need fastp; need bwa; need samtools; need bedtools

# macs3 or macs2
MACS_BIN=""
if command -v macs3 >/dev/null 2>&1; then MACS_BIN="macs3";
elif command -v macs2 >/dev/null 2>&1; then MACS_BIN="macs2";
else echo "ERROR: need macs3 or macs2 installed"; exit 1; fi

# Map sample → genome 
genome_for() { echo "hg38"; }

# -------------------------- preflight refs -----------------------------------
# If refs are missing, but data/processed/hg38.fa.gz exists, prepare refs here
prepare_refs_if_needed() {
  local ref_fa="$REFDIR/hg38.fa"
  local ref_fa_gz_src="$SUB/hg38.fa.gz"

  if [[ -f "${ref_fa}.bwt" ]]; then
    return 0
  fi

  mkdir -p "$REFDIR"

  if [[ ! -f "$ref_fa" ]]; then
    if [[ -f "$ref_fa_gz_src" ]]; then
      echo "Preparing refs from $ref_fa_gz_src"
      ln -sf "$ref_fa_gz_src" "$REFDIR/hg38.fa.gz"
      gzip -cd "$REFDIR/hg38.fa.gz" > "$ref_fa"
    elif [[ -f "$REFDIR/hg38.fa.gz" ]]; then
      echo "Decompressing existing $REFDIR/hg38.fa.gz"
      gzip -cd "$REFDIR/hg38.fa.gz" > "$ref_fa"
    else
      echo "ERROR: Missing refs. Provide $SUB/hg38.fa.gz or run 02_build_refs.sh"; exit 1
    fi
  fi

  if [[ ! -f "${ref_fa}.bwt" ]]; then
    echo "Building BWA index for hg38..."
    bwa index "$ref_fa"
  fi

  if [[ ! -f "$REFDIR/hg38.fa.fai" ]]; then
    echo "Indexing FASTA and creating chrom.sizes..."
    faidx "$ref_fa"
    cut -f1,2 "$REFDIR/hg38.fa.fai" > "$REFDIR/hg38.chrom.sizes"
  fi

  # Blacklist optional; if present, it will be used for filtering
}

# ------------------------------- helpers --------------------------------------
call_and_filter() {
  local ACC="$1"
  local genome_size_flag="$2"    # hs (hg38) / mm (mm10)
  local blacklist="$3"

  local BAM="$ALIGN/${ACC}.dedup.bam"
  local OUTPFX="$PEAKS/${ACC}"
  local RAWNP="${OUTPFX}_peaks.narrowPeak"
  local FILT="${PEAKS}/${ACC}.peaks.filt.bed"

  if [[ ! -f "$BAM" ]]; then
    echo "$ACC: skip (no BAM)"
    return 0
  fi

  if [[ "$FORCE" -eq 1 || ! -f "$RAWNP" ]]; then
    echo "$ACC: calling peaks"
    $MACS_BIN callpeak -t "$BAM" -f BAMPE -g "$genome_size_flag" -n "$ACC" --outdir "$PEAKS" \
      --pvalue 1e-3 --nomodel --shift -100 --extsize 200 --keep-dup all
  else
    echo "$ACC: peaks exist -> skip callpeak"
  fi

  if [[ -f "$RAWNP" && -f "$blacklist" ]]; then
    if [[ "$FORCE" -eq 1 || ! -f "$FILT" ]]; then
      bedtools intersect -v -a "$RAWNP" -b "$blacklist" > "$FILT"
      echo "$ACC: wrote blacklist-filtered peaks -> $FILT"
    else
      echo "$ACC: filtered peaks exist -> skip"
    fi
  fi
}

# === [1/4] fastp (trim/QC) ===
echo "=== [1/4] fastp (trim/QC) ==="
for ACC in "${SAMPLES[@]}"; do
  # Prefer new subsampled naming; fallback to older R1/R2 .sub.fastq.gz
  IN1="$SUB/${ACC}_1.fastq.subsampled"
  IN2="$SUB/${ACC}_2.fastq.subsampled"
  if [[ ! -f "$IN1" || ! -f "$IN2" ]]; then
    IN1="$SUB/${ACC}_R1.sub.fastq.gz"
    IN2="$SUB/${ACC}_R2.sub.fastq.gz"
  fi
  OUT1="$WORK/${ACC}_R1.trim.fastq.gz"
  OUT2="$WORK/${ACC}_R2.trim.fastq.gz"
  HTML="$QC/${ACC}.fastp.html"
  JSON="$QC/${ACC}.fastp.json"

  [[ -f "$IN1" && -f "$IN2" ]] || { echo "Missing $IN1 or $IN2"; exit 1; }

  if [[ "$FORCE" -eq 1 || ! -f "$OUT1" || ! -f "$OUT2" ]]; then
    echo "fastp -> $ACC"
    fastp -w "$THREADS" -i "$IN1" -I "$IN2" -o "$OUT1" -O "$OUT2" -h "$HTML" -j "$JSON" >/dev/null
  else
    echo "fastp -> $ACC (skip: outputs exist)"
  fi
done

# === [2/4] bwa mem → fixmate chain ===
echo "=== [2/4] bwa mem → fixmate chain ==="
# Ensure references exist (auto-prepare from data/processed if needed)
prepare_refs_if_needed
for ACC in "${SAMPLES[@]}"; do
  REFBASE="$(genome_for "$ACC")"
  REFPREFIX="$REFDIR/${REFBASE}.fa"  # match bwa index prefix from 02_build_refs.sh
  [[ -f "${REFPREFIX}.bwt" ]] || { echo "Missing BWA index for $REFBASE at ${REFPREFIX}.*"; exit 1; }

  TRIM1="$WORK/${ACC}_R1.trim.fastq.gz"
  TRIM2="$WORK/${ACC}_R2.trim.fastq.gz"
  NSORT="$ALIGN/${ACC}.nsort.bam"
  FIXM="$ALIGN/${ACC}.fixmate.bam"
  SORTB="$ALIGN/${ACC}.sorted.bam"
  DEDUP="$ALIGN/${ACC}.dedup.bam"

  if [[ "$FORCE" -eq 1 ]]; then
    rm -f "$NSORT" "$FIXM" "$SORTB" "$SORTB.bai" "$DEDUP" "$DEDUP.bai" 2>/dev/null || true
  fi

  if [[ -f "$DEDUP" && -f "$DEDUP.bai" ]]; then
    echo "$ACC: dedup BAM exists -> skip alignment"
    continue
  fi

  echo "$ACC: bwa mem"
  bwa mem -t "$THREADS" "$REFPREFIX" "$TRIM1" "$TRIM2" \
    | samtools view -bh -F 0x4 - \
    | samtools sort -n -@ "$THREADS" -o "$NSORT" -

  echo "$ACC: fixmate"
  samtools fixmate -m "$NSORT" "$FIXM"

  echo "$ACC: position sort"
  samtools sort -@ "$THREADS" -o "$SORTB" "$FIXM"
  samtools index "$SORTB"

  echo "$ACC: markdup (remove PCR dups)"
  samtools markdup -r -@ "$THREADS" "$SORTB" "$DEDUP"
  samtools index "$DEDUP"

  rm -f "$NSORT" "$FIXM" "$SORTB" "$SORTB.bai" 2>/dev/null || true
done

# === [3/4] macs3 peaks + blacklist ===
echo "=== [3/4] macs3 peaks + blacklist ==="
for ACC in "${SAMPLES[@]}"; do
  # hg38 only for now (PDX/mm10 can be added later if needed)
  call_and_filter "$ACC" hs "$REFDIR/hg38.blacklist.bed"
done

# === [4/4] FRiP ===
echo "=== [4/4] FRiP ==="
FRIP="$QC/frip.tsv"
echo -e "sample\treads_in_peaks\tproper_pairs\tfrip" > "$FRIP"

for ACC in "${SAMPLES[@]}"; do
  BAM="$ALIGN/${ACC}.dedup.bam"
  NP="${PEAKS}/${ACC}_peaks.narrowPeak"

  if [[ ! -f "$BAM" || ! -f "$NP" ]]; then
    echo "$ACC: skip FRiP (missing BAM or peaks)"
    continue
  fi

  PP=$(samtools view -c -f 0x2 -F 0x904 "$BAM")
  if [[ "$PP" -eq 0 ]]; then
    echo -e "${ACC}\t0\t0\t0.0" >> "$FRIP"
    continue
  fi

  QNAME="$ALIGN/${ACC}.qname.bam"
  FRAG="$ALIGN/${ACC}.frags.bed"
  samtools view -u -f 0x2 -F 0x904 "$BAM" \
    | samtools sort -n -@ "$THREADS" -o "$QNAME" -
  bedtools bamtobed -bedpe -i "$QNAME" \
    | awk '($1!="." && $2>=0 && $3>$2){print $1"\t"$2"\t"$3}' > "$FRAG"

  INPK=$(bedtools intersect -u -a "$FRAG" -b "$NP" | wc -l | tr -d ' ')
  FRIP_VAL=$(python3 - <<PY
pp=int("$PP"); ip=int("$INPK")
print(0.0 if pp==0 else round(ip/pp,4))
PY
)
  echo -e "${ACC}\t${INPK}\t${PP}\t${FRIP_VAL}" >> "$FRIP"
done

echo "FRiP -> $FRIP"
echo "Done ✅"
