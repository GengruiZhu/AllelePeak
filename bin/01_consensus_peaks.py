#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
ATAC-Quant v4 | Step 01: merge peaks from multiple samples -> one unified peak set

When doing differential analysis across samples, all samples must be counted on the
same feature set (v1 called peaks per sample, so peak IDs did not correspond to each
other). This step takes the union of each sample's narrowPeak and merges overlapping
intervals.

  - Only 1 sample: pass through as-is (keeps the original peak IDs, so the result can
    be compared against v1)
  - Multiple samples: merge overlaps (gap <= --merge-gap); --min-samples can filter
    low-support peaks

Output:
  consensus_peaks.bed          chrom start end name nSamples
  consensus_peaks.members.tsv  name -> source sample:original peak ID
"""
import os
import sys
import argparse
from collections import defaultdict

sys.path.insert(0, os.path.join(os.path.dirname(os.path.abspath(__file__)), "..", "lib"))
from atacq_common import log, read_peaks, die  # noqa: E402


def main():
    ap = argparse.ArgumentParser(description="ATAC-Quant v4 Step 01: consensus peaks")
    ap.add_argument("--peaks", nargs="+", required=True, metavar="SAMPLE:PATH")
    ap.add_argument("-o", "--outdir", required=True)
    ap.add_argument("--merge-gap", type=int, default=0)
    ap.add_argument("--min-samples", type=int, default=1)
    ap.add_argument("--prefix", default="cpeak")
    args = ap.parse_args()

    os.makedirs(args.outdir, exist_ok=True)
    out_bed = os.path.join(args.outdir, "consensus_peaks.bed")
    out_mem = os.path.join(args.outdir, "consensus_peaks.members.tsv")

    items = []
    for spec in args.peaks:
        if ":" not in spec:
            die("--peaks format must be SAMPLE:PATH, got {}".format(spec))
        s, path = spec.split(":", 1)
        pk = read_peaks(path)
        log("[01] {}: {:,} peaks ({})".format(s, len(pk), path))
        items.append((s, pk))

    if len(items) == 1:
        s, pk = items[0]
        with open(out_bed, "w") as fb, open(out_mem, "w") as fm:
            fm.write("PeakID\tMembers\n")
            for c, st, en, name in sorted(pk, key=lambda x: (x[0], x[1], x[2])):
                fb.write("{}\t{}\t{}\t{}\t1\n".format(c, st, en, name))
                fm.write("{}\t{}:{}\n".format(name, s, name))
        log("[01] single sample, passing through {:,} peaks unchanged -> {}".format(len(pk), out_bed))
        return

    by_chrom = defaultdict(list)
    for s, pk in items:
        for c, st, en, name in pk:
            by_chrom[c].append((st, en, s, name))

    n_out = n_drop = 0
    with open(out_bed, "w") as fb, open(out_mem, "w") as fm:
        fm.write("PeakID\tMembers\n")
        for c in sorted(by_chrom):
            iv = sorted(by_chrom[c])
            cur = None
            def flush(cur):
                nonlocal n_out, n_drop
                st, en, mem = cur
                samples = {m[0] for m in mem}
                if len(samples) < args.min_samples:
                    n_drop += 1
                    return
                n_out += 1
                name = "{}_{:07d}".format(args.prefix, n_out)
                fb.write("{}\t{}\t{}\t{}\t{}\n".format(c, st, en, name, len(samples)))
                fm.write("{}\t{}\n".format(name, ",".join("{}:{}".format(a, b) for a, b in mem)))
            for st, en, s, name in iv:
                if cur is not None and st <= cur[1] + args.merge_gap:
                    cur[1] = max(cur[1], en)
                    cur[2].append((s, name))
                else:
                    if cur is not None:
                        flush(cur)
                    cur = [st, en, [(s, name)]]
            if cur is not None:
                flush(cur)
    log("[01] merged into {:,} peaks (support samples < {}, dropped {:,}) -> {}".format(
        n_out, args.min_samples, n_drop, out_bed))


if __name__ == "__main__":
    main()
