#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
eval_truth.py -- evaluate the ATAC-Quant v4 output against the ground truth from
simulate_mini.py.

Ground-truth counting: every non-duplicate, aligned fragment is counted on the features
          that its "span on the true source copy" falls on (split evenly when it spans two
          features; falls on no feature -> background).
Reference: N_Primary = the number of fragments whose primary alignment falls on the feature
          (bowtie2 picks the primary alignment at random; this is the featureCounts view).
Checks:
  1. conservation: EM feature total + background = number of non-duplicate fragments
  2. feature level: r / mean absolute error / slope of EM vs ground truth, compared with
     N_Primary
  3. HRU group level and lineage level (called+projected) shares
  4. homeologous projection: whether a truly open element is covered by a feature on its
     copy (especially on copies where no peak was called)
  5. indistinguishable identical copies (A=B): Resolved_Frac should be low, while the group
     total should still be right
"""
import os
import sys
import json
import math
import argparse
from collections import defaultdict


def pearson(a, b):
    n = len(a)
    ma, mb = sum(a) / n, sum(b) / n
    sa = math.sqrt(sum((x - ma) ** 2 for x in a))
    sb = math.sqrt(sum((y - mb) ** 2 for y in b))
    return sum((x - ma) * (y - mb) for x, y in zip(a, b)) / (sa * sb) if sa and sb else float("nan")


def slope(y, x):
    mx, my = sum(x) / len(x), sum(y) / len(y)
    sxx = sum((a - mx) ** 2 for a in x)
    return sum((a - mx) * (b - my) for a, b in zip(x, y)) / sxx if sxx else float("nan")


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--sim", required=True)
    ap.add_argument("--out", required=True, help="atacquant_v2 output directory")
    ap.add_argument("--sample", default="SIM")
    args = ap.parse_args()

    feats = []
    with open(os.path.join(args.out, "features", "features.tsv")) as f:
        h = f.readline().rstrip("\n").split("\t")
        for line in f:
            d = dict(zip(h, line.rstrip("\n").split("\t")))
            d["Start"], d["End"] = int(d["Start"]), int(d["End"])
            feats.append(d)
    q = {}
    with open(os.path.join(args.out, "quant", args.sample, args.sample + ".quant.tsv")) as f:
        h = f.readline().rstrip("\n").split("\t")
        for line in f:
            d = dict(zip(h, line.rstrip("\n").split("\t")))
            q[d["FeatureID"]] = d
    st = json.load(open(os.path.join(args.out, "quant", args.sample, args.sample + ".quant_stats.json")))

    by_chrom = defaultdict(list)
    for i, ft in enumerate(feats):
        by_chrom[ft["Chrom"]].append((ft["Start"], ft["End"], i))
    truth = [0.0] * len(feats)
    n_frag = n_bg = 0
    with open(os.path.join(args.sim, "truth_fragments.tsv")) as f:
        f.readline()
        for line in f:
            qn, c, s, e, src, dup, al = line.rstrip("\n").split("\t")
            if dup == "1" or al == "0":
                continue
            s, e = int(s), int(e)
            n_frag += 1
            hit = [i for (a, b, i) in by_chrom.get(c, []) if a < e and b > s]
            if not hit:
                n_bg += 1
            for i in hit:
                truth[i] += 1.0 / len(hit)

    em = [float(q[ft["FeatureID"]]["EM_Count"]) for ft in feats]
    nv = [float(q[ft["FeatureID"]]["N_Primary"]) for ft in feats]
    rf = [q[ft["FeatureID"]]["Resolved_Frac"] for ft in feats]

    fails = []
    print("=" * 70)
    print("non-duplicate aligned fragments {:,} (ground-truth background {:,}) | EM feature total {:,.1f} + background {:,.1f} = {:,.1f}".format(
        n_frag, n_bg, st["em_sum_features"], st["em_background"],
        st["em_sum_features"] + st["em_background"]))
    tot = st["em_sum_features"] + st["em_background"]
    ignore_mode = st.get("params", {}).get("outside") == "ignore"
    if ignore_mode:
        print("(outside=ignore: no background component is built, so the conservation comparison against the total fragment count is skipped; only the internal EM conservation is checked)")
        if not st["conservation_check"]["ok"]:
            fails.append("internal EM conservation check failed")
    elif abs(tot - n_frag) > 1:
        fails.append("conservation: EM total {:.1f} != fragment count {}".format(tot, n_frag))
    if st.get("n_fragments_total_nondup") is not None and st["n_fragments_total_nondup"] != n_frag:
        fails.append("fragment accounting: pipeline reports {} != ground truth {}".format(st["n_fragments_total_nondup"], n_frag))

    def report(label, idx):
        if len(idx) < 3:
            return None
        t = [truth[i] for i in idx]
        e = [em[i] for i in idx]
        n = [nv[i] for i in idx]
        r_em, r_nv = pearson(e, t), pearson(n, t)
        mae_em = sum(abs(a - b) for a, b in zip(e, t)) / len(t)
        mae_nv = sum(abs(a - b) for a, b in zip(n, t)) / len(t)
        print("{:<34} n={:>4} | EM r={:.4f} MAE={:6.2f} slope={:.3f} | naive r={:.4f} MAE={:6.2f} slope={:.3f}".format(
            label, len(idx), r_em, mae_em, slope(e, t), r_nv, mae_nv, slope(n, t)))
        return r_em, r_nv, mae_em, mae_nv

    print("-" * 70 + "\nfeature level (EM vs ground truth; naive primary-alignment count vs ground truth)")
    allr = report("all features", list(range(len(feats))))
    report("HRU called", [i for i, ft in enumerate(feats) if ft["Type"] == "called"])
    report("HRU projected", [i for i, ft in enumerate(feats) if ft["Type"] == "projected"])
    report("orphan", [i for i, ft in enumerate(feats) if ft["Type"].startswith("orphan")])
    res_idx = [i for i in range(len(feats)) if rf[i] not in ("NA",) and float(rf[i]) >= 0.5]
    report("Resolved_Frac >= 0.5", res_idx)
    if allr and allr[0] < allr[1]:
        fails.append("the r of EM ({:.4f}) is worse than that of the naive count ({:.4f})".format(allr[0], allr[1]))
    if allr and allr[2] > allr[3]:
        fails.append("the MAE of EM ({:.2f}) is larger than that of the naive count ({:.2f})".format(allr[2], allr[3]))

    # group level
    gt, ge, gn = defaultdict(float), defaultdict(float), defaultdict(float)
    for i, ft in enumerate(feats):
        gt[ft["GroupID"]] += truth[i]
        ge[ft["GroupID"]] += em[i]
        gn[ft["GroupID"]] += nv[i]
    ks = [k for k in gt if k.startswith("HRU_")]
    print("-" * 70)
    print("HRU group level n={} | EM r={:.4f} | naive r={:.4f}".format(
        len(ks), pearson([ge[k] for k in ks], [gt[k] for k in ks]),
        pearson([gn[k] for k in ks], [gt[k] for k in ks])))

    # lineage shares
    lt, le, ln = defaultdict(float), defaultdict(float), defaultdict(float)
    for i, ft in enumerate(feats):
        if ft["Type"] in ("called", "projected"):
            lt[ft["Ancestry"]] += truth[i]
            le[ft["Ancestry"]] += em[i]
            ln[ft["Ancestry"]] += nv[i]
    T, E, Nn = sum(lt.values()), sum(le.values()), sum(ln.values())
    print("-" * 70 + "\nlineage shares (HRU features)   ground truth / EM / naive")
    for l in sorted(lt):
        print("  {:<6} {:6.1%} / {:6.1%} / {:6.1%}".format(l, lt[l] / T, le[l] / E, ln[l] / Nn))
        if abs(le[l] / E - lt[l] / T) > 0.03:
            fails.append("lineage {} share is off by more than 3 percentage points (EM {:.1%} vs ground truth {:.1%})".format(
                l, le[l] / E, lt[l] / T))

    # projection coverage
    peaks = defaultdict(list)
    with open(os.path.join(args.sim, "peaks.bed")) as f:
        for line in f:
            c, s, e, n = line.split()[:4]
            peaks[c].append((int(s), int(e)))
    cov_called = cov_uncalled = n_called = n_uncalled = 0
    lam_uncalled_cov = lam_uncalled = 0.0
    ab = []
    with open(os.path.join(args.sim, "truth_elements.tsv")) as f:
        f.readline()
        for line in f:
            L, e, a, lin, c, s, t, lam, iab = line.rstrip("\n").split("\t")
            s, t, lam = int(s), int(t), float(lam)
            called = any(ps < t and pe > s for ps, pe in peaks[c])
            covered = [i for (fs, fe, i) in by_chrom.get(c, []) if fs < t and fe > s]
            if called:
                n_called += 1
                cov_called += bool(covered)
            else:
                n_uncalled += 1
                cov_uncalled += bool(covered)
                lam_uncalled += lam
                if covered:
                    lam_uncalled_cov += lam
            if iab == "1":
                ab.append((L, a, covered))
    print("-" * 70)
    print("element coverage: copies with a called peak {}/{} have a feature; copies without a called peak {}/{} are covered "
          "by projection (their true openness {:.0%} is covered)".format(cov_called, n_called, cov_uncalled, n_uncalled,
                                          lam_uncalled_cov / lam_uncalled if lam_uncalled else 0))
    if n_uncalled and cov_uncalled / n_uncalled < 0.7:
        fails.append("projection coverage is too low: {:.0%}".format(cov_uncalled / n_uncalled))

    # A=B identical copies
    rfs = []
    for L, a, cov in ab:
        for i in cov:
            if rf[i] != "NA":
                rfs.append(float(rf[i]))
    if rfs:
        print("A=B fully identical copies: mean Resolved_Frac of their features {:.2f} (should be low, meaning they were correctly marked indistinguishable)".format(
            sum(rfs) / len(rfs)))
    print("=" * 70)
    if fails:
        print("FAILED:")
        for x in fails:
            print("  - " + x)
        sys.exit(1)
    print("ALL CHECKS PASSED")


if __name__ == "__main__":
    main()
