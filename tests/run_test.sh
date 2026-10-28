#!/bin/bash
# ==============================================================================
# ATAC-Quant v4 self-test: simulate a small 7-copy polyploid ATAC dataset (with ground
# truth), run the full pipeline and compare the result against the ground truth.
#   When samtools + pysam are available, the SAM is converted to a real BAM and the
#   pysam code path is tested (recommended); otherwise a mock BAM (pure Python) is used
#   and only the algorithm is verified.
# Usage: bash tests/run_test.sh [workdir]      (roughly 1-3 minutes)
# Pass criterion: every check prints a [PASS] line, the last line contains
#   "ALL CHECKS PASSED", and the exit status is 0. Any failure prints a [FAIL] line and
#   exits non-zero.
# ==============================================================================
set -euo pipefail
HERE="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
ROOT="$(dirname "$HERE")"
PY="${PYTHON:-python3}"
W="${1:-$HERE/_work}"
rm -rf "$W"
mkdir -p "$W"

NCHK=0
pass() { NCHK=$((NCHK + 1)); echo "[PASS] $1"; }
fail() { echo "[FAIL] $1"; echo "RESULT: FAILED after $NCHK passing check(s)"; exit 1; }

# count lines of a file (number of records, independent of a trailing newline)
nlines() { awk 'END { print NR }' "$1" 2>/dev/null || echo 0; }

# check_file <path> <label> <min lines>
check_file() {
    local p="$1"
    local label="$2"
    local min="${3:-1}"
    local n
    if [ ! -f "$p" ]; then
        fail "$label: expected file does not exist: $p"
    fi
    n=$(nlines "$p")
    if [ "${n:-0}" -lt "$min" ]; then
        fail "$label: $p has ${n:-0} line(s), expected at least $min"
    fi
    pass "$label: $p (${n} line(s))"
}

# ---------------------------------------------------------------- 1. simulation
if "$PY" "$HERE/simulate_mini.py" -o "$W/sim" --depth 2 --sample SIM; then
    pass "simulation exited 0"
else
    fail "simulation did not exit 0"
fi

for f in genome.fa genes.gff3 clusters.tsv ancestry.tsv peaks.bed \
         truth_fragments.tsv truth_elements.tsv samples.tsv; do
    check_file "$W/sim/$f" "simulation artifact" 1
done
# the simulation must also have produced a usable alignment file
check_file "$W/sim/sim.sam" "simulation alignment file (SAM)" 2

# ---------------------------------------------------------------- 2. BAM or mock BAM
if command -v samtools >/dev/null 2>&1 && "$PY" -c "import pysam" 2>/dev/null; then
    echo "[test] samtools + pysam detected: testing the pysam code path with a real BAM"
    samtools sort -o "$W/sim/sim.bam" "$W/sim/sim.sam"
    samtools index "$W/sim/sim.bam"
    sed "s#sim\.sam#sim.bam#" "$W/sim/samples.tsv" > "$W/samples.tsv"
    export ATACQ_MOCK_BAM=0
    check_file "$W/sim/sim.bam" "real BAM built by samtools" 1
    check_file "$W/sim/sim.bam.bai" "real BAM index" 1
    if grep -q "sim\.bam" "$W/samples.tsv"; then
        pass "samples.tsv points the BAM column at sim.bam"
    else
        fail "samples.tsv was not rewritten to reference sim.bam"
    fi
else
    echo "[test] samtools/pysam not detected: testing with a mock BAM (algorithm only, the pysam path is not exercised)"
    cp "$W/sim/samples.tsv" "$W/samples.tsv"
    export ATACQ_MOCK_BAM=1
    check_file "$W/samples.tsv" "samples.tsv copied for the mock BAM path" 2
fi

# ---------------------------------------------------------------- 3. pipeline
if "$PY" "$ROOT/atacquant_core.py" --samples "$W/samples.tsv" \
    --genome "$W/sim/genome.fa" --gff "$W/sim/genes.gff3" \
    --clusters "$W/sim/clusters.tsv" --ancestry "$W/sim/ancestry.tsv" \
    -o "$W/out" -t 4 --steps peaks,features,quant,matrix,qc \
    --bootstrap 5 --quant-extra "--buckets 8" \
    --lineage-map "Officinarum=SOF,Robustum=SRO,Spontaneum=SSP"; then
    pass "atacquant_core.py exited 0"
else
    fail "atacquant_core.py did not exit 0"
fi

# ---------------------------------------------------------------- 4. output artifacts
FEAT="$W/out/features/features.tsv"
QUANT="$W/out/quant/SIM/SIM.quant.tsv"
STATS="$W/out/quant/SIM/SIM.quant_stats.json"
TRUTH_FRAG="$W/sim/truth_fragments.tsv"

check_file "$FEAT" "feature table" 2
check_file "$QUANT" "quantification table" 2
check_file "$W/out/features/features_stats.json" "feature statistics" 1
check_file "$STATS" "quantification statistics (JSON)" 1
check_file "$W/out/matrix/feature_counts.tsv" "count matrix (per feature)" 2
check_file "$W/out/matrix/group_counts.tsv" "count matrix (per group)" 2
check_file "$W/out/qc/qc_report.md" "QC report" 1

# every quantified feature must carry a finite, non-negative EM count
if "$PY" - "$QUANT" <<'PYCHECK'
import sys

path = sys.argv[1]
with open(path, encoding="utf-8") as f:
    header = f.readline().rstrip("\n").split("\t")
    try:
        idx = header.index("EM_Count")
    except ValueError:
        sys.stderr.write("features: the EM_Count column is missing from %s\n" % path)
        sys.exit(1)
    n = 0
    for line in f:
        cols = line.rstrip("\n").split("\t")
        if len(cols) <= idx:
            continue
        try:
            v = float(cols[idx])
        except ValueError:
            sys.stderr.write("features: non-numeric EM_Count %r\n" % cols[idx])
            sys.exit(1)
        if v != v or v < 0:
            sys.stderr.write("features: EM_Count out of range: %r\n" % cols[idx])
            sys.exit(1)
        n += 1
print("[info] EM_Count is finite and non-negative for %d feature(s)" % n)
sys.exit(0 if n > 0 else 1)
PYCHECK
then
    pass "EM_Count is numeric, finite and non-negative for every feature"
else
    fail "the quantification table failed the numeric EM_Count check"
fi

# ---------------------------------------------------------------- 5. mass balance
# EM allocation over features plus the background component must add up to the number of
# non-duplicate, aligned fragments in the ground truth (tolerance 1 fragment), unless the
# run was configured with outside=ignore (no background component), in which case only the
# internal EM conservation flag is required to hold.
if "$PY" - "$STATS" "$TRUTH_FRAG" <<'PYCHECK'
import json
import sys

stats_path, truth_path = sys.argv[1], sys.argv[2]

with open(stats_path, encoding="utf-8") as f:
    st = json.load(f)


def num(key):
    if key not in st:
        sys.stderr.write("mass balance: key %r is missing from the statistics file\n" % key)
        sys.exit(1)
    try:
        return float(st[key])
    except (TypeError, ValueError):
        sys.stderr.write("mass balance: %r is not numeric: %r\n" % (key, st[key]))
        sys.exit(1)


feat_sum = num("em_sum_features")
bg_sum = num("em_background")
total = feat_sum + bg_sum
print("[info] EM total: features %.3f + background %.3f = %.3f" % (feat_sum, bg_sum, total))

n_frag = 0
with open(truth_path, encoding="utf-8") as f:
    f.readline()
    for line in f:
        cols = line.rstrip("\n").split("\t")
        if len(cols) < 7:
            continue
        dup, aligned = cols[5], cols[6]
        if dup == "1" or aligned == "0":
            continue
        n_frag += 1
print("[info] non-duplicate aligned fragments in the ground truth = %d" % n_frag)

params = st.get("params") or {}
if params.get("outside") == "ignore":
    ok = bool((st.get("conservation_check") or {}).get("ok"))
    print("[info] outside=ignore: no background component is built, so only the internal "
          "EM conservation flag applies: %s" % ok)
    sys.exit(0 if ok else 1)

diff = abs(total - n_frag)
print("[info] |EM total - non-duplicate fragments| = %.3f (tolerance 1.0) -> %s"
      % (diff, "ok" if diff <= 1.0 else "mismatch"))
sys.exit(0 if diff <= 1.0 else 1)
PYCHECK
then
    pass "mass balance holds (EM feature total + background equals the fragment count)"
else
    fail "mass balance does not hold"
fi

# ---------------------------------------------------------------- 6. ground-truth evaluation
if "$PY" "$HERE/eval_truth.py" --sim "$W/sim" --out "$W/out" --sample SIM; then
    pass "eval_truth.py exited 0 (agreement with the ground truth accepted)"
else
    fail "eval_truth.py did not exit 0"
fi

echo "RESULT: ${NCHK} checks passed, 0 failed"
echo "ALL CHECKS PASSED"
