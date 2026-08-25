#!/bin/bash
# Download + align village-dog / wild-canid pilot samples against CanFam3.1,
# WITHOUT ever writing full raw fastq to disk (the pilot's raw .sra data,
# 114.6GB total, does not fit alongside everything else on this machine's
# disk if downloaded simultaneously -- streaming straight into alignment
# and keeping only the much smaller BAM output is the only way to fit all
# 12 samples). Same pattern used successfully in the belyaev-fox pipeline:
# prefetch -> fasterq-dump --stdout -> bwa mem -> samtools view -F 4 ->
# samtools sort, with the .sra deleted only after a verified (flagstat,
# not just quickcheck) successful alignment.
#
# LibraryLayout varies per sample here (unlike the fox pilot): some wild
# canid runs are SINGLE-end, not PAIRED. fasterq-dump/bwa mem invocation
# branches on the LibraryLayout column in the sample sheet accordingly.
#
# Lessons carried over from the fox pipeline (see belyaev_foxes/scripts
# in the fennec_data repo for the original writeups):
#   - prefetch: retry up to 3x, clearing only a stale .sra.lock between
#     attempts -- never wipe the whole working dir, since prefetch
#     resumes/dedups partial or complete downloads on its own
#   - fasterq-dump: --size-check off (its default heuristic throws
#     false-positive "disk-limit exeeded!" errors)
#   - verification: samtools flagstat (mapped reads > 0), not just
#     quickcheck, before deleting the source .sra
#   - the .sra is only deleted once the sample is fully DONE, so a
#     downstream failure (e.g. fasterq-dump/bwa error) doesn't force a
#     wasted redownload on retry
#
# Region-restriction (added after 4 full-genome BAMs, ~5GB each, ate
# enough disk headroom to make every subsequent sample's fasterq-dump
# fail with storage-exhausted -- the accumulating BAM floor, not
# per-sample size, was the real constraint): immediately after
# alignment, each BAM is subset to the 246 candidate domestication
# regions from Pendleton et al. 2018 (BMC Biology, Table S5/Additional
# File 1; 10.81Mb total, see metadata/pendleton2018_246CDR.bed) and the
# full-genome BAM is discarded. This cut per-sample footprint from
# ~4-5GB to ~20-25MB (measured, not estimated) and removes the
# accumulating-floor problem entirely.
set -uo pipefail

source ~/miniforge3/etc/profile.d/conda.sh
conda activate village-dogs-genomics

cd "$(dirname "$0")/.."

SAMPLE_SHEET="${1:-metadata/pilot6_samples.tsv}"
REF=reference/reference.fna
BED=metadata/pendleton2018_246CDR.bed
SRA_TMP=sra_tmp
FQD_TMP=fqd_tmp
mkdir -p aligned "$SRA_TMP" "$FQD_TMP"
MASTER_LOG=aligned/align.master.log

log() { echo "[$(date '+%Y-%m-%d %H:%M:%S')] $*" | tee -a "$MASTER_LOG"; }

align_one_sample() {
  local RUN=$1 NAME=$2 CAT=$3 LAYOUT=$4
  local BAM="aligned/${RUN}.sorted.bam"

  local ATTEMPT
  for ATTEMPT in 1 2 3; do
    log "[$RUN/$NAME/$CAT] prefetch attempt $ATTEMPT"
    rm -f "${SRA_TMP}/${RUN}/${RUN}.sra.lock"
    prefetch "$RUN" -O "$SRA_TMP" >> "$MASTER_LOG" 2>&1 && break
    log "[$RUN/$NAME/$CAT] prefetch attempt $ATTEMPT failed"
    if [ "$ATTEMPT" = 3 ]; then
      log "[$RUN/$NAME/$CAT] GIVING UP on prefetch after 3 attempts"
      return 1
    fi
  done

  if [ ! -f "${SRA_TMP}/${RUN}/${RUN}.sra" ]; then
    log "[$RUN/$NAME/$CAT] no .sra file after prefetch, skipping"
    return 1
  fi

  log "[$RUN/$NAME/$CAT] stream fasterq-dump | bwa mem ($LAYOUT) | samtools sort"
  if [ "$LAYOUT" = "PAIRED" ]; then
    fasterq-dump --stdout --split-spot --skip-technical --size-check off \
        -t "$FQD_TMP" -e 4 "${SRA_TMP}/${RUN}/${RUN}.sra" \
        2> "aligned/${RUN}.fqd.log" \
      | bwa mem -p -t 4 -R "@RG\tID:${RUN}\tSM:${RUN}\tPL:ILLUMINA" "$REF" - \
        2> "aligned/${RUN}.bwa.log" \
      | samtools view -b -F 4 - \
      | samtools sort -o "$BAM" -
  else
    # --concatenate-reads: some runs reported as SINGLE in runinfo
    # actually have mixed/partial-pair spot structure that makes
    # fasterq-dump want a 3-way split (forward/reverse/unpaired),
    # which conflicts with plain --stdout (single stream only) --
    # e.g. SRR2149863/Dingo failed instantly with "requested mode
    # (FASTQ split 3) would produce multiple files". Concatenating
    # whole spots into one stream sidesteps that regardless of the
    # actual underlying layout.
    fasterq-dump --stdout --concatenate-reads --skip-technical --size-check off \
        -t "$FQD_TMP" -e 4 "${SRA_TMP}/${RUN}/${RUN}.sra" \
        2> "aligned/${RUN}.fqd.log" \
      | bwa mem -t 4 -R "@RG\tID:${RUN}\tSM:${RUN}\tPL:ILLUMINA" "$REF" - \
        2> "aligned/${RUN}.bwa.log" \
      | samtools view -b -F 4 - \
      | samtools sort -o "$BAM" -
  fi

  rm -rf "${FQD_TMP:?}"/*

  if ! samtools quickcheck -v "$BAM" 2>> "$MASTER_LOG"; then
    log "[$RUN/$NAME/$CAT] FAILED quickcheck"
    rm -f "$BAM"
    return 1
  fi

  local MAPPED
  MAPPED=$(samtools flagstat "$BAM" | head -1 | awk '{print $1}')
  if [ -z "$MAPPED" ] || [ "$MAPPED" -eq 0 ]; then
    log "[$RUN/$NAME/$CAT] FAILED: 0 mapped reads in BAM"
    rm -f "$BAM"
    return 1
  fi

  # Subset to the 246 candidate domestication regions and discard the
  # full-genome BAM -- see header comment. Verified full-genome BAM
  # (MAPPED>0 above) is what we subset from, so a bad region-extraction
  # here is a real failure, not a false negative.
  local REGION_BAM="aligned/${RUN}.region.bam"
  if ! samtools view -b -L "$BED" "$BAM" -o "$REGION_BAM" 2>> "$MASTER_LOG"; then
    log "[$RUN/$NAME/$CAT] FAILED region extraction"
    rm -f "$BAM" "$REGION_BAM"
    return 1
  fi
  mv -f "$REGION_BAM" "$BAM"

  samtools index "$BAM"
  rm -rf "${SRA_TMP:?}/${RUN}"
  local REGION_MAPPED
  REGION_MAPPED=$(samtools flagstat "$BAM" | head -1 | awk '{print $1}')
  log "[$RUN/$NAME/$CAT] DONE: $MAPPED mapped reads genome-wide, $REGION_MAPPED in the 246 CDR regions"
  return 0
}

tail -n +2 "$SAMPLE_SHEET" | while IFS=$'\t' read -r RUN BIOSAMPLE NAME BREED CAT SIZE_MB BASES LAYOUT; do
  BAM="aligned/${RUN}.sorted.bam"

  if [ -f "$BAM" ] && samtools quickcheck -v "$BAM" 2>/dev/null; then
    MAPPED=$(samtools flagstat "$BAM" 2>/dev/null | head -1 | awk '{print $1}')
    if [ -n "$MAPPED" ] && [ "$MAPPED" -gt 0 ]; then
      log "[$RUN/$NAME/$CAT] already present and verified ($MAPPED mapped reads), skipping"
      continue
    fi
  fi

  SAMPLE_ATTEMPT_OK=0
  for SAMPLE_ATTEMPT in 1 2; do
    if align_one_sample "$RUN" "$NAME" "$CAT" "$LAYOUT"; then
      SAMPLE_ATTEMPT_OK=1
      break
    fi
    log "[$RUN/$NAME/$CAT] sample-level attempt $SAMPLE_ATTEMPT failed, retrying whole sample"
  done

  if [ "$SAMPLE_ATTEMPT_OK" -eq 0 ]; then
    log "[$RUN/$NAME/$CAT] *** GIVING UP after repeated failures -- needs manual attention ***"
    # Once a sample has exhausted both sample-level attempts, no further
    # retry will reuse its .sra automatically -- keeping it around is
    # just dead weight (this orphaned-sra accumulation was directly
    # responsible for a disk-exhaustion cascade during the run: three
    # abandoned .sra files, ~10GB, sat unused while later samples starved
    # for temp space). Safe to delete; a manual retry re-downloads it.
    rm -rf "${SRA_TMP:?}/${RUN}"
  fi
done

log "=== Village-dog/wild-canid alignment run finished (sample sheet: $SAMPLE_SHEET) ==="
