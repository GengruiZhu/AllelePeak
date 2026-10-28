#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
qc_bias_check.py -- a minimal validation of the q_c (discriminative-site density) bias correction.

[The mechanism to be validated]
    hard-evidence count  NBU_c  proportional to  m_c x q_c
    That is: a copy with many discriminative sites naturally receives more reads that
    "look unique" -> it is systematically over-estimated.
    If this holds, the relationship between the EM alpha and the ground truth should be
        alpha_c  proportional to  truth_c x pi_c        (pi_c = the "pinning rate" of that copy)
    so it can be corrected to
        m_c_hat  =  alpha_c / pi_c
    and pi_c can be predicted from [sequence alone]: pi_c ~= 1 - exp(-q_c * L)
        (q_c = discriminative-site density; L = mean fragment length; a fragment is pinned
         as soon as it covers >=1 discriminative site)

[Computing q_c from sequence]
    For every (locus, element), align the homeologous interval of the copies with each other
    (star alignment, the longest one taken as reference), then decide column by column:
    if the base of a copy in that column is unique within the whole group -> that column is
    one of its "discriminative sites".
    q_c = number of discriminative columns / number of columns where that copy is covered
    in this interval

[The four things this script does]
    A. Compute q_c from sequence (alignment + discriminative columns)
    B. Mechanism check: does the observed pinning rate pi_obs = N_Unique/ground truth agree
       with q_c (including the prediction pi=1-exp(-qL))
    C. Correction: raw EM  vs  sequence correction  vs  observed-pinning-rate correction
       (the last one is the [upper bound of the correction])
    D. Bias at the level of within-group shares (after removing the absolute scale, check
       whether the relative bias is squeezed out)

Usage:
    python analysis/qc_bias_check.py --sim <simulation directory> --out <pipeline output directory> [--sample SIM]
"""
import os
import sys
import math
import argparse
from collections import defaultdict

import numpy as np


# --------------------------------------------------------------------------- sequence
def read_fasta(path, keep=None):
    seqs, name, buf = {}, None, []
    with open(path) as f:
        for line in f:
            if line.startswith(">"):
                if name is not None:
                    seqs[name] = "".join(buf)
                name = line[1:].split()[0]
                buf = []
            else:
                buf.append(line.strip())
    if name is not None:
        seqs[name] = "".join(buf)
    if keep is not None:
        seqs = {k: v for k, v in seqs.items() if k in keep}
    return seqs


def nw_align(a, b, match=2, mism=-3, gap=-4):
    """Global alignment (Needleman-Wunsch); returns the two aligned sequences. The sequences are short (~250bp), so this is fast enough."""
    n, m = len(a), len(b)
    NEG = -10 ** 9
    S = np.full((n + 1, m + 1), NEG, dtype=np.int32)
    T = np.zeros((n + 1, m + 1), dtype=np.uint8)      # 0=diagonal 1=up 2=left
    S[0, :] = np.arange(m + 1) * gap
    S[:, 0] = np.arange(n + 1) * gap
    A = np.frombuffer(a.encode("ascii"), dtype=np.uint8)
    B = np.frombuffer(b.encode("ascii"), dtype=np.uint8)
    for i in range(1, n + 1):
        ai = A[i - 1]
        for j in range(1, m + 1):
            d = S[i - 1, j - 1] + (match if ai == B[j - 1] else mism)
            u = S[i - 1, j] + gap
            l = S[i, j - 1] + gap
            if d >= u and d >= l:
                S[i, j], T[i, j] = d, 0
            elif u >= l:
                S[i, j], T[i, j] = u, 1
            else:
                S[i, j], T[i, j] = l, 2
    # traceback
    oa, ob = [], []
    i, j = n, m
    while i > 0 or j > 0:
        t = T[i, j]
        if i > 0 and j > 0 and t == 0:
            oa.append(a[i - 1]); ob.append(b[j - 1]); i -= 1; j -= 1
        elif i > 0 and (j == 0 or t == 1):
            oa.append(a[i - 1]); ob.append("-"); i -= 1
        else:
            oa.append("-"); ob.append(b[j - 1]); j -= 1
    return "".join(reversed(oa)), "".join(reversed(ob))


def load_truth_elements(path):
    rows = []
    with open(path, encoding="utf-8") as f:
        h = f.readline().rstrip("\n").split("\t")
        for line in f:
            d = dict(zip(h, line.rstrip("\n").split("\t")))
            d["Locus"] = int(d["Locus"]); d["Start"] = int(d["Start"]); d["End"] = int(d["End"])
            d["Lambda"] = float(d["Lambda"])
            rows.append(d)
    return rows


def compute_qc(truth_rows, genome):
    """Compute q_c = discriminative-site density for every (Locus, Element, Allele)."""
    groups = defaultdict(list)
    for r in truth_rows:
        groups[(r["Locus"], r["Element"])].append(r)
    qc, ncol = {}, {}
    for key, rows in sorted(groups.items()):
        if len(rows) < 2:
            for r in rows:
                qc[(r["Locus"], r["Element"], r["Allele"])] = float("nan")
                ncol[(r["Locus"], r["Element"], r["Allele"])] = 0
            continue
        # fetch the sequences
        seqs = {}
        for r in rows:
            s = genome.get(r["Chrom"])
            if s is None:
                continue
            seqs[r["Allele"]] = s[r["Start"]:r["End"]].upper()
        if len(seqs) < 2:
            continue
        # star alignment: the longest sequence is the reference
        ref = max(seqs, key=lambda k: len(seqs[k]))
        aln = {ref: seqs[ref]}
        for k, v in seqs.items():
            if k == ref:
                continue
            _, b2 = nw_align(seqs[ref], v)
            aln[k] = b2
        L = len(aln[ref])
        cols = [aln[k] for k in aln]
        diag = {k: 0 for k in aln}
        cov = {k: 0 for k in aln}
        for c in range(L):
            bases = {}
            for k in aln:
                ch = aln[k][c]
                if ch == "-":
                    continue
                cov[k] += 1
                bases.setdefault(ch, []).append(k)
            for ch, ks in bases.items():
                if len(ks) == 1:                       # this base belongs to only one copy in this column -> discriminative site
                    diag[ks[0]] += 1
        for k in aln:
            qc[(key[0], key[1], k)] = diag[k] / cov[k] if cov[k] else 0.0
            ncol[(key[0], key[1], k)] = cov[k]
    return qc


# --------------------------------------------------------------------------- pipeline output
def load_features(path):
    feats = []
    with open(path, encoding="utf-8") as f:
        h = f.readline().rstrip("\n").split("\t")
        for line in f:
            d = dict(zip(h, line.rstrip("\n").split("\t")))
            d["Start"] = int(d["Start"]); d["End"] = int(d["End"])
            feats.append(d)
    return feats


def load_quant(path):
    q = {}
    with open(path, encoding="utf-8") as f:
        h = f.readline().rstrip("\n").split("\t")
        for line in f:
            d = dict(zip(h, line.rstrip("\n").split("\t")))
            q[d["FeatureID"]] = d
    return q


def compute_truth_counts(sim, feats):
    """Same ground-truth convention as tests/eval_truth.py: count a fragment on the features its true source span falls on, splitting evenly when it spans several."""
    by_chrom = defaultdict(list)
    for i, ft in enumerate(feats):
        by_chrom[ft["Chrom"]].append((ft["Start"], ft["End"], i))
    truth = np.zeros(len(feats))
    n_frag = 0
    with open(os.path.join(sim, "truth_fragments.tsv"), encoding="utf-8") as f:
        f.readline()
        for line in f:
            qn, c, s, e, src, dup, al = line.rstrip("\n").split("\t")
            if dup == "1" or al == "0":
                continue
            s, e = int(s), int(e)
            n_frag += 1
            hit = [i for (a, b, i) in by_chrom.get(c, []) if a < e and b > s]
            for i in hit:
                truth[i] += 1.0 / len(hit)
    return truth, n_frag


def overlap_tagged(feats, truth_rows):
    """Tag every feature with the (Locus, Element, Allele) it overlaps most."""
    byc = defaultdict(list)
    for r in truth_rows:
        byc[r["Chrom"]].append(r)
    tag = []
    for ft in feats:
        best, bov = None, 0
        for r in byc.get(ft["Chrom"], []):
            ov = min(ft["End"], r["End"]) - max(ft["Start"], r["Start"])
            if ov > bov:
                bov, best = ov, r
        tag.append(best)
    return tag


def pearson(x, y):
    x, y = np.asarray(x, float), np.asarray(y, float)
    m = np.isfinite(x) & np.isfinite(y)
    x, y = x[m], y[m]
    if len(x) < 3:
        return float("nan"), 0
    xm, ym = x - x.mean(), y - y.mean()
    d = math.sqrt((xm ** 2).sum() * (ym ** 2).sum())
    return (float((xm * ym).sum() / d) if d else float("nan")), len(x)


def describe(b, name):
    b = np.asarray(b, float)
    b = b[np.isfinite(b)]
    if not len(b):
        return "  %-26s (no data)" % name
    return ("  %-26s n=%5d | median deviation %6.2f%% | p90 %7.2f%% | max %8.2f%% | median ratio %.3f"
            % (name, len(b), 100 * np.median(np.abs(b - 1)), 100 * np.percentile(np.abs(b - 1), 90),
               100 * np.abs(b - 1).max(), np.median(b)))


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--sim", required=True)
    ap.add_argument("--out", required=True)
    ap.add_argument("--sample", default="SIM")
    ap.add_argument("--min-truth", type=float, default=10.0, help="features whose ground truth is below this value are excluded from the statistics")
    ap.add_argument("--frag-len", type=float, default=105.0, help="mean fragment length of the simulation set (sim reports ~105bp)")
    args = ap.parse_args()

    truth_rows = load_truth_elements(os.path.join(args.sim, "truth_elements.tsv"))
    chroms = {r["Chrom"] for r in truth_rows}
    genome = read_fasta(os.path.join(args.sim, "genome.fa"), keep=chroms)
    print("=" * 96)
    print("A. Computing q_c from sequence (discriminative-site density)")
    print("=" * 96)
    qc = compute_qc(truth_rows, genome)
    vals = [v for v in qc.values() if np.isfinite(v)]
    print("  computed q_c for %d (locus,element,copy) triples | median %.4f | range %.4f ~ %.4f /bp"
          % (len(vals), np.median(vals), min(vals), max(vals)))
    by_lin = defaultdict(list)
    for r in truth_rows:
        k = (r["Locus"], r["Element"], r["Allele"])
        if k in qc and np.isfinite(qc[k]):
            by_lin[r["Lineage"]].append(qc[k])
    for lin in sorted(by_lin):
        print("    median q_c (lineage %s, n=%d) = %.4f /bp" % (lin, len(by_lin[lin]), np.median(by_lin[lin])))

    # ---------------------------------------------------------------- pipeline output
    feats = load_features(os.path.join(args.out, "features", "features.tsv"))
    q = load_quant(os.path.join(args.out, "quant", args.sample, args.sample + ".quant.tsv"))
    truth, n_frag = compute_truth_counts(args.sim, feats)
    tag = overlap_tagged(feats, truth_rows)
    em = np.array([float(q[f["FeatureID"]]["EM_Count"]) for f in feats])
    nu = np.array([float(q[f["FeatureID"]]["N_Unique"]) for f in feats])
    nbu = np.array([float(q[f["FeatureID"]]["N_BestUnique"]) for f in feats])

    qcv, lam, loc, ale = [], [], [], []
    for i, f in enumerate(feats):
        r = tag[i]
        if r is None:
            qcv.append(np.nan); lam.append(np.nan); loc.append(None); ale.append(None); continue
        k = (r["Locus"], r["Element"], r["Allele"])
        qcv.append(qc.get(k, np.nan)); lam.append(r["Lambda"]); loc.append(k); ale.append(r["Allele"])
    qcv = np.array(qcv, float); lam = np.array(lam, float)
    print("\n  features %d | with q_c %d | ground-truth fragment total %.0f" % (len(feats), np.isfinite(qcv).sum(), truth.sum()))

    m = np.isfinite(qcv) & (truth >= args.min_truth) & (em > 0)
    print("  features entering the statistics (ground truth >=%.0f and EM>0): %d" % (args.min_truth, m.sum()))

    # ---------------------------------------------------------------- B. mechanism check
    print("\n" + "=" * 96)
    print("B. Mechanism check: the relationship between the observed pinning rate pi_obs = N_Unique/ground truth and q_c")
    print("   prediction (physical model): pi ~= 1 - exp(-q_c * L), L = %.0f bp" % args.frag_len)
    print("=" * 96)
    pi_obs = nu / np.maximum(truth, 1e-9)
    pi_pred = 1.0 - np.exp(-qcv * args.frag_len)
    r1, n1 = pearson(qcv[m], pi_obs[m])
    r2, n2 = pearson(pi_pred[m], pi_obs[m])
    print("  corr( q_c          , pi_obs ) = %+.3f   (n=%d)" % (r1, n1))
    print("  corr( 1-exp(-q_c*L), pi_obs ) = %+.3f   (n=%d)   <- prediction vs observed" % (r2, n2))
    print("\n  binned by q_c:")
    print("    %-16s %6s %10s %10s %10s %10s" % ("q_c bin", "n", "median q_c", "median pi_obs", "median pi_pred", "median EM/truth"))
    edges = np.array([0, 0.004, 0.008, 0.016, 0.032, 1e9])
    for a, b in zip(edges[:-1], edges[1:]):
        s = m & (qcv >= a) & (qcv < b)
        if s.sum() < 5:
            continue
        lab = "%.3f-%.3f" % (a, min(b, 0.9)) if b < 1e8 else ">%.3f" % a
        print("    %-16s %6d %10.4f %10.3f %10.3f %10.3f"
              % (lab, s.sum(), np.median(qcv[s]), np.median(pi_obs[s]), np.median(pi_pred[s]),
                 np.median(em[s] / truth[s])))

    # ---------------------------------------------------------------- C. bias and correction
    print("\n" + "=" * 96)
    print("C. Per-copy bias (EM/ground truth) and the three estimates")
    print("=" * 96)
    # correction 1: sequence-predicted pinning rate
    pc = np.where(np.isfinite(qcv), 1.0 - np.exp(-qcv * args.frag_len), np.nan)
    pc = np.maximum(pc, 1e-6)
    em_seq = em / pc                                    # m_hat = alpha / pi
    # correction 2: observed pinning rate (oracle, the upper bound)
    po = np.maximum(pi_obs, 1e-6)
    em_obs = em / po
    # normalise to the same total as raw so that absolute scales can be compared
    for arr, nm in ((em_seq, "seq"), (em_obs, "obs")):
        f = em[m].sum() / arr[m].sum()
        arr *= f

    print("\n  [all features]")
    print(describe(em[m] / truth[m], "raw EM / ground truth"))
    print(describe(em_seq[m] / truth[m], "sequence correction / ground truth"))
    print(describe(em_obs[m] / truth[m], "observed pinning-rate correction / ground truth"))
    rr_raw, _ = pearson(qcv[m], np.log(np.maximum(em[m] / truth[m], 1e-9)))
    rr_seq, _ = pearson(qcv[m], np.log(np.maximum(em_seq[m] / truth[m], 1e-9)))
    rr_obs, _ = pearson(qcv[m], np.log(np.maximum(em_obs[m] / truth[m], 1e-9)))
    print("\n  correlation of the bias with q_c (should be positive = high q is over-estimated; should approach 0 after correction)")
    print("    corr( q_c, log(raw/ground truth)                    ) = %+.3f" % rr_raw)
    print("    corr( q_c, log(sequence correction/ground truth)    ) = %+.3f" % rr_seq)
    print("    corr( q_c, log(observed pinning-rate correction/ground truth) ) = %+.3f" % rr_obs)

    # ---------------------------------------------------------------- D. within-group shares
    print("\n" + "=" * 96)
    print("D. Bias at the level of within-group shares (remove the absolute scale and look at the relative bias)")
    print("=" * 96)
    bygrp = defaultdict(list)
    for i, f in enumerate(feats):
        bygrp[f["GroupID"]].append(i)
    rows_raw, rows_seq, rows_obs, qq = [], [], [], []
    for g, idx in bygrp.items():
        idx = [i for i in idx if np.isfinite(qcv[i]) and truth[i] >= args.min_truth and em[i] > 0]
        if len(idx) < 2:
            continue
        t = truth[idx]; s = t.sum()
        if s <= 0:
            continue
        A = em[idx] / em[idx].sum()
        B = em_seq[idx] / em_seq[idx].sum()
        C = em_obs[idx] / em_obs[idx].sum()
        T = t / s
        for k, i in enumerate(idx):
            rows_raw.append(A[k] / T[k] if T[k] > 0 else np.nan)
            rows_seq.append(B[k] / T[k] if T[k] > 0 else np.nan)
            rows_obs.append(C[k] / T[k] if T[k] > 0 else np.nan)
            qq.append(qcv[i])
    if rows_raw:
        print("\n  within-group share ratio (estimated share / true share; 1.0 is perfect)")
        print(describe(np.array(rows_raw), "raw"))
        print(describe(np.array(rows_seq), "sequence correction"))
        print(describe(np.array(rows_obs), "observed pinning-rate correction"))
        r3, n3 = pearson(qq, np.log(np.maximum(rows_raw, 1e-9)))
        r4, n4 = pearson(qq, np.log(np.maximum(rows_seq, 1e-9)))
        r5, n5 = pearson(qq, np.log(np.maximum(rows_obs, 1e-9)))
        print("\n    corr( q_c, log(within-group share ratio) ): raw %+.3f -> sequence correction %+.3f -> observed correction %+.3f (n=%d)"
              % (r3, r4, r5, n3))
        # high-q half vs low-q half
        qq = np.array(qq); rr = np.array(rows_raw); ss = np.array(rows_seq)
        if len(qq) >= 20:
            cut = np.median(qq)
            hi, lo = qq >= cut, qq < cut
            print("    split at the median q_c (median within-group share ratio):")
            print("      high-q half: raw %.3f | sequence correction %.3f" % (np.median(rr[hi]), np.median(ss[hi])))
            print("      low-q half:  raw %.3f | sequence correction %.3f" % (np.median(rr[lo]), np.median(ss[lo])))
            print("      high/low ratio: raw %.3f -> sequence correction %.3f  (closer to 1 is better)"
                  % (np.median(rr[hi]) / np.median(rr[lo]), np.median(ss[hi]) / np.median(ss[lo])))

    # ---------------------------------------------------------------- E. truly closed copies
    print("\n" + "=" * 96)
    print("E. Do truly closed copies (Lambda=0) receive any reads scattered onto them")
    print("=" * 96)
    z = (lam == 0) & np.isfinite(qcv) & (em >= 0)
    nz = (lam > 0) & np.isfinite(qcv) & (em > 0)
    if z.sum() and nz.sum():
        print("  copies with Lambda=0 (n=%d): median EM %.3f | median q_c of these copies %.4f"
              % (z.sum(), np.median(em[z]), np.median(qcv[z])))
        print("  neighbouring Lambda>0 copies in the same group (n=%d): median EM %.3f | median q_c %.4f"
              % (nz.sum(), np.median(em[nz]), np.median(qcv[nz])))
        # within a group: the share of the group EM taken by closed copies
        sh = []
        for g, idx in bygrp.items():
            idx = [i for i in idx if np.isfinite(qcv[i])]
            if not idx or em[idx].sum() <= 0:
                continue
            c0 = [i for i in idx if lam[i] == 0]
            if c0:
                sh.append(sum(em[i] for i in c0) / em[idx].sum())
        if sh:
            print("  EM share taken by [truly closed copies] within a group: median %.1f%% | p90 %.1f%% | max %.1f%%"
                  % (100 * np.median(sh), 100 * np.percentile(sh, 90), 100 * max(sh)))
            print("  (truly closed copies should get 0; the share they get is exactly the amount scattered onto them)")
    print("\n" + "=" * 96)


if __name__ == "__main__":
    main()
