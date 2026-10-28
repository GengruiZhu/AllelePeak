#!/bin/bash
# ==============================================================================
# ATAC-Quant v4 | Step 00b: alignment + mark PCR duplicates + call peaks (sensitive parameters variant)
#
# This file is a copy of bin/00_map_call.sh; the only difference is that bowtie2 uses a
# "more sensitive + no truncation of ties" parameter combination, to fix the problem of
# "missing tied-best alignments" (background and the upper bound of the gain are in
# calibration/calibration_report_missed_ties.md, diagnostics are in diagnostics/):
#   1. the default --bt2-k is raised from 20 to 100 (less -k truncation);
#   2. --very-sensitive is added (higher search sensitivity, fewer misses);
#   3. --no-1mm-upfront is added (bowtie2 >= 2.3.0; appended automatically when supported,
#      silently omitted otherwise), so that seeds carrying 1 mismatch are also examined.
#   All other lines are identical to 00_map_call.sh. Note: this variant has not been
#   validated on the server; before routine use, it is advisable to first run
#   diagnostics/02_remap_subset.sh under the same parameter combination for a small-scale
#   control.
#
# Run it with a new --outdir or a new --prefix (otherwise a BAM already present under the
# same prefix is reused).
#
# Usage (same as 00_map_call.sh):
#   00b_map_call_sensitive.sh --r1 R1.fq.gz --r2 R2.fq.gz --index /path/bt2/genome --outdir DIR \
#                  --prefix SAMPLE --threads 48 [--bt2-k 100] [--maxins 1000] \
#                  [--gsize 1.0e10] [--qvalue 0.05] [--macs2-env ATAC] [--macs2 macs2]
#   Peak calling only (alignment BAM already available):
#   00b_map_call_sensitive.sh --bam existing.bam --outdir DIR --prefix SAMPLE --threads 48 ...
#   Alignment only, no peak calling: add --skip-peaks
# ==============================================================================
set -euo pipefail

R1=""; R2=""; IDX=""; OUTDIR=""; PREFIX=""; THREADS=8
BT2_K=100; MAXINS=1000; GSIZE="1.0e10"; QVAL="0.05"
IN_BAM=""; SKIP_PEAKS=0; MACS2_ENV=""; MACS2_BIN="macs2"; SORT_MEM="1G"
ASSAY="atac"; MACS2_EXTRA=""        # v4.2: kept in sync with 00_map_call.sh

while [ $# -gt 0 ]; do
    case "$1" in
        --r1) R1="$2"; shift 2;;
        --r2) R2="$2"; shift 2;;
        --index) IDX="$2"; shift 2;;
        --outdir) OUTDIR="$2"; shift 2;;
        --prefix) PREFIX="$2"; shift 2;;
        --threads) THREADS="$2"; shift 2;;
        --bt2-k) BT2_K="$2"; shift 2;;
        --maxins) MAXINS="$2"; shift 2;;
        --gsize) GSIZE="$2"; shift 2;;
        --qvalue) QVAL="$2"; shift 2;;
        --bam) IN_BAM="$2"; shift 2;;
        --skip-peaks) SKIP_PEAKS=1; shift 1;;
        --macs2-env) MACS2_ENV="$2"; shift 2;;
        --macs2) MACS2_BIN="$2"; shift 2;;
        --sort-mem) SORT_MEM="$2"; shift 2;;
        --assay) ASSAY="$2"; shift 2;;
        --macs2-extra) MACS2_EXTRA="$2"; shift 2;;
        *) echo "[FATAL] unknown argument: $1"; exit 1;;
    esac
done

[ -z "$OUTDIR" ] && { echo "[FATAL] --outdir is required"; exit 1; }
[ -z "$PREFIX" ] && { echo "[FATAL] --prefix is required"; exit 1; }
case "$ASSAY" in atac|chip) ;; *) echo "[FATAL] --assay must be atac or chip, got: $ASSAY"; exit 1;; esac
mkdir -p "$OUTDIR"
TMP="$OUTDIR/tmp_${PREFIX}"
mkdir -p "$TMP"
now() { date "+%Y-%m-%d %H:%M:%S"; }
need() { command -v "$1" >/dev/null 2>&1 || { echo "[FATAL] cannot find $1"; exit 1; }; }

# Decide paired-end/single-end (sample 20000 records, look at bit 0x1 of the flag). For single-end, -f 2 (proper pair) selects an empty set, and MACS2 needs -f BAM.
detect_pe() {
    local n
    n=$(samtools view -@ 2 "$1" 2>/dev/null | head -n 20000 \
        | awk '{ if ($2 % 2 == 1) c++ } END { print c+0 }' || true)
    if [ "${n:-0}" -gt 100 ]; then echo 1; else echo 0; fi
}

# ---- sensitive-parameter compatibility probe: quietly omit --no-1mm-upfront on older bowtie2 ----
BOWTIE2_EXTRA="--very-sensitive"
if command -v bowtie2 >/dev/null 2>&1 && bowtie2 --help 2>&1 | grep -q -- "--no-1mm-upfront"; then
    BOWTIE2_EXTRA="$BOWTIE2_EXTRA --no-1mm-upfront"
fi

need samtools
MARKED="$OUTDIR/${PREFIX}.marked.bam"

# ------------------------------------------------------------------ 1. alignment
if [ -z "$IN_BAM" ]; then
    need bowtie2
    [ -z "$R1" ] || [ -z "$R2" ] || [ -z "$IDX" ] && { echo "[FATAL] alignment requires --r1 --r2 --index"; exit 1; }
    if [ ! -f "${IDX}.1.bt2" ] && [ ! -f "${IDX}.1.bt2l" ]; then
        echo "[FATAL] bowtie2 index not found: ${IDX}.1.bt2(l); run bowtie2-build first"; exit 1
    fi
    RAW="$OUTDIR/${PREFIX}.raw.bam"
    if [ ! -f "${RAW}.bai" ]; then
        echo "[$(now)] [00b] bowtie2 -k ${BT2_K} ${BOWTIE2_EXTRA} (keep multi-mappings for the EM; sensitive-parameter fix)"
        bowtie2 -p "$THREADS" -x "$IDX" -1 "$R1" -2 "$R2" \
            -X "$MAXINS" --no-mixed --no-discordant $BOWTIE2_EXTRA -k "$BT2_K" \
            2> "$OUTDIR/${PREFIX}.bowtie2.log" | \
        samtools view -@ 4 -b - | \
        samtools sort -@ "$THREADS" -m "$SORT_MEM" -T "$TMP/sort" -o "$RAW" -
        samtools index -@ "$THREADS" "$RAW"
    fi
    # -------------------------------------------------------------- 2. mark duplicates
    if [ ! -f "${MARKED}.bai" ]; then
        if command -v sambamba >/dev/null 2>&1; then
            echo "[$(now)] [00b] sambamba markdup (mark only, do not remove)"
            sambamba markdup -t "$THREADS" --tmpdir="$TMP" "$RAW" "$MARKED"
        else
            echo "[$(now)] [00b] samtools markdup (mark only, do not remove)"
            samtools collate -@ "$THREADS" -O -u "$RAW" "$TMP/collate" | \
                samtools fixmate -@ "$THREADS" -m -u - - | \
                samtools sort -@ "$THREADS" -m "$SORT_MEM" -T "$TMP/sort2" -u - | \
                samtools markdup -@ "$THREADS" -T "$TMP/md" - "$MARKED"
        fi
        samtools index -@ "$THREADS" "$MARKED"
    fi
else
    [ -f "$IN_BAM" ] || { echo "[FATAL] cannot find --bam $IN_BAM"; exit 1; }
    MARKED="$IN_BAM"
    [ -f "${IN_BAM}.bai" ] || samtools index -@ "$THREADS" "$IN_BAM"
    echo "[$(now)] [00b] using the existing BAM: $IN_BAM"
fi

# Path of the quantification BAM, read by the driver
echo "$MARKED" > "$OUTDIR/${PREFIX}.quant_bam.txt"

# ------------------------------------------------------------------ 3. call peaks
if [ "$SKIP_PEAKS" -eq 0 ]; then
    # Already existing peak file (narrowPeak / broadPeak / bed, any of them)
    PEAKS=""
    for cand in "$OUTDIR/${PREFIX}_peaks.narrowPeak" "$OUTDIR/${PREFIX}_peaks.broadPeak" \
                "$OUTDIR/${PREFIX}_peaks.bed"; do
        if [ -f "$cand" ]; then PEAKS="$cand"; break; fi
    done
    if [ -z "$PEAKS" ]; then
        PRIM="$OUTDIR/${PREFIX}.primary.bam"
        # ---- single-end/paired-end: for single-end, -f 2 selects an empty set, and MACS2 needs -f BAM ----
        PE=$(detect_pe "$MARKED")
        if [ "$PE" = "1" ]; then
            PFLAG="-f 2"; PFMT="BAMPE"
        else
            PFLAG=""; PFMT="BAM"
            echo "[$(now)] [00b] [WARNING] single-end data detected"
            echo "[$(now)] [00b]        - call peaks with -f BAM and without -f 2 (adding it selects an empty set)"
            echo "[$(now)] [00b]        - [IMPORTANT] step 03 quantification infers the fragment span from read1 TLEN; single-end degenerates to"
            echo "[$(now)] [00b]          the read's own span (~read length), which makes the EffLen effective length wrong and the fragment-vs-feature"
            echo "[$(now)] [00b]          overlap assignment biased. The conclusion [may not be trustworthy]; paired-end data is recommended."
        fi
        if [ ! -f "${PRIM}.bai" ]; then
            echo "[$(now)] [00b] extracting the primary+non-duplicate subset for peak calling (PE=${PE}, assay=${ASSAY})"
            # -F 0xD04: drop unmapped(0x4) secondary(0x100) duplicate(0x400) supplementary(0x800)
            # shellcheck disable=SC2086
            samtools view -@ "$THREADS" -b $PFLAG -F 0xD04 -o "$PRIM" "$MARKED"
            samtools index -@ "$THREADS" "$PRIM"
        fi
        # ---- ATAC (Tn5) needs the +4/-5 Tn5 offset correction; ChIP (antibody enrichment) has no such physical process, so it must be skipped ----
        if [ "$ASSAY" = "atac" ]; then
            need alignmentSieve
            SHIFT="$OUTDIR/${PREFIX}.primary.shifted.bam"
            if [ ! -f "${SHIFT}.bai" ]; then
                echo "[$(now)] [00b] alignmentSieve --ATACshift (+4/-5)   [assay=atac]"
                alignmentSieve -b "$PRIM" -o "$TMP/shift.unsorted.bam" --ATACshift -p "$THREADS"
                samtools sort -@ "$THREADS" -m "$SORT_MEM" -T "$TMP/sort3" -o "$SHIFT" "$TMP/shift.unsorted.bam"
                samtools index -@ "$THREADS" "$SHIFT"
                rm -f "$TMP/shift.unsorted.bam"
            fi
            PEAKBAM="$SHIFT"
        else
            echo "[$(now)] [00b] assay=chip: skipping the Tn5 (+4/-5) offset correction (antibody enrichment has no such physical process)"
            PEAKBAM="$PRIM"
        fi
        echo "[$(now)] [00b] MACS2 callpeak (-f ${PFMT}, q=${QVAL}, g=${GSIZE}) ${MACS2_EXTRA}"
        if [ -n "$MACS2_ENV" ]; then
            set +u
            # shellcheck disable=SC1091
            source activate "$MACS2_ENV" || conda activate "$MACS2_ENV"
            set -u
        fi
        # shellcheck disable=SC2086
        "$MACS2_BIN" callpeak -t "$PEAKBAM" -f "$PFMT" -g "$GSIZE" -n "$PREFIX" \
            --outdir "$OUTDIR" -q "$QVAL" --keep-dup all ${MACS2_EXTRA}
    else
        echo "[$(now)] [00b] peak file already exists, skipping: $PEAKS"
    fi
fi

rm -rf "$TMP"
echo "[$(now)] [00b] done: $PREFIX"
