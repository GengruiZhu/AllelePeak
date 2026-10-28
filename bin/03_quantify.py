#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
ATAC-Quant v4 | Step 03: alignment-compatible equivalence-class EM quantification (single sample)

This is the core of v2, corresponding to fixes B1/B6/B8/B10/B11 in the LQ root-cause report
(LQ_lineage_kmer_root_cause.md),
and adapted for ATAC:

  [B1]  v1's EM never used the alignment results (W = theta*L*P). In v2 each fragment is only
        allocated among the features it is "actually aligned to",
        and the compatibility weight is softened by alignment score: w = exp((AS_pair - AS_best) / T)
        (the unfinished item in section 7 of the LQ report)
  [B6]  De-duplicate by fragment (read name); take the union of all -k alignments of a fragment,
        and count the two mates only once
  [B8]  All secondary alignments participate (v1 only looked at the "first appearing peak")
  [B10] Effective length: an ATAC fragment counts as long as it overlaps a feature,
        eff = L + mean fragment length - 1
  [B11] Cross-cluster / cross-feature multi-mapping fragments are no longer double counted:
        one EM for the whole genome, and each fragment's total weight is always 1
  [New] Genomic background component (similar to Salmon's decoy): when a fragment has an equally
        good (or better) alignment outside features,
        part of its mass goes to the background instead of all being stuffed into a peak.
        Background eff = genome length - total feature length
  [New] PCR duplicates are removed at the fragment level (judged by the 0x400 flag of the primary
        alignment read1; all alignment records are removed together);
        "orphan secondary alignments" with no primary alignment (left over after v1 used
        sambamba -r to remove duplicates) are also removed
  [New] Orphan peaks are no longer counted by "counting all read1 records" (v1 was inflated 10.8x);
        they live in the same EM as HRU features
  [New] Identifiability output: N_Unique (number of fragments with a unique candidate),
        N_BestUnique (number of fragments with a unique best), Resolved_Frac, optional bootstrap SD

Workflow:
  Pass A (parallel by chromosome): scan read1 records overlapping features ->
        (fragment hash, feature, pair score)
  Pass B (parallel by chromosome, with --outside model): scan read1 records outside features,
        keep only fragments that "also appeared inside features"
        -> (fragment hash, pair score, count); count pure-background fragments
  Reduce (parallel by hash bucket): assemble the candidate set of each fragment -> equivalence
        classes (EC)
  EM: global EC-EM (vectorized), optional bootstrap

Output (--outdir):
  <sample>.quant.tsv        FeatureID Length EffLen EM_Count EM_SD N_Unique N_BestUnique
                            N_Candidate N_Primary Resolved_Frac Density_per_kb
  <sample>.quant_stats.json fragment accounting, multi-mapping fraction, background share,
                            convergence status, mass balance check
"""
import os
import sys
import math
import time
import bisect
import shutil
import argparse
import multiprocessing as mp
from array import array
from collections import defaultdict

import numpy as np

sys.path.insert(0, os.path.join(os.path.dirname(os.path.abspath(__file__)), "..", "lib"))
from atacq_common import (log, die, ensure_hashseed, qhash, load_features, open_bam,  # noqa
                          bam_references, write_json, revcomp, Fasta)

ensure_hashseed()

QS = 1000          # weight quantization scale
G = {}             # worker global state (inherited under fork; on platforms without fork, injected by spawn + initializer)


def _spawn_init_G(d):
    """Worker initialization on spawn platforms: inject the main process's G into the child process module globals (used on Windows)."""
    G.update(d)


# =============================================================================
# Small helpers
# =============================================================================
def build_feature_index(feats):
    by = defaultdict(list)
    for i, f in enumerate(feats):
        by[f["Chrom"]].append((f["Start"], f["End"], i))
    idx = {}
    for c, lst in by.items():
        lst.sort()
        idx[c] = ([x[0] for x in lst], [x[1] for x in lst], [x[2] for x in lst],
                  max(x[1] - x[0] for x in lst))
    return idx


def overlaps(fi, fs, fe, minov):
    starts, ends, ids, maxlen = fi
    lo = bisect.bisect_left(starts, fs - maxlen)
    hi = bisect.bisect_left(starts, fe)
    out = None
    for i in range(lo, hi):
        e = ends[i]
        if e > fs:
            s = starts[i]
            ov = (e if e < fe else fe) - (s if s > fs else fs)
            if ov >= minov:
                if out is None:
                    out = [ids[i]]
                else:
                    out.append(ids[i])
    return out


def frag_span(r, flag):
    """Fragment span on the reference (0-based half-open), derived from the read1 record + TLEN; both mates give the same span."""
    tlen = r.template_length
    if (flag & 1) and tlen != 0:
        if tlen > 0:
            s = r.reference_start
            return s, s + tlen
        s = r.next_reference_start
        return s, s - tlen
    s = r.reference_start
    e = r.reference_end
    return s, (e if e is not None else s + 1)


def pair_score(r, flag):
    try:
        a = r.get_tag("AS")
    except KeyError:
        a = 0
    if flag & 1:
        try:
            y = r.get_tag("YS")
        except KeyError:
            y = a
        return a + y
    return a


def np_u64(a):
    return np.frombuffer(a, dtype=np.uint64).copy() if len(a) else np.zeros(0, np.uint64)


def np_i32(a):
    return np.frombuffer(a, dtype=np.int32).copy() if len(a) else np.zeros(0, np.int32)


def split_by_bucket(B, q, *arrs):
    """Split into buckets by q % B, returns [(q_b, arr1_b, ...)]"""
    if len(q) == 0:
        return [tuple([q] + [a for a in arrs])] * 0 + \
               [tuple([q[:0]] + [a[:0] for a in arrs]) for _ in range(B)]
    b = (q % np.uint64(B)).astype(np.int64)
    order = np.argsort(b, kind="stable")
    cnt = np.bincount(b, minlength=B)
    cuts = np.concatenate([[0], np.cumsum(cnt)])
    out = []
    qo = q[order]
    ao = [a[order] for a in arrs]
    for k in range(B):
        sl = slice(cuts[k], cuts[k + 1])
        out.append(tuple([qo[sl]] + [a[sl] for a in ao]))
    return out


# =============================================================================
# Optional: lineage k-mer
# =============================================================================
def load_lineage_kmers(specs, k):
    names, kd = [], {}
    for li, spec in enumerate(specs):
        name, path = spec.split(":", 1)
        names.append(name.upper())
        n = 0
        with open(path) as f:
            for line in f:
                if line.startswith(">"):
                    continue
                s = line.strip().upper()
                if len(s) >= k:
                    kd.setdefault(s[:k], li)
                    n += 1
        log("[03] lineage k-mer {}: {:,} entries ({})".format(name, n, path))
    # canonical detection (the same test as LQ snipe_kmer v2)
    probe = list(kd.keys())[:2000]
    both = sum(1 for km in probe if revcomp(km) in kd)
    canon = (both / max(len(probe), 1)) < 0.05
    return names, kd, canon


def count_lineage_hits(seq, kd, k, L, canon):
    h = [0] * L
    s = seq.upper()
    n = len(s) - k + 1
    if n <= 0:
        return h
    if canon:
        rc = revcomp(s)
        m = len(s)
        for i in range(n):
            km = s[i:i + k]
            j = kd.get(km)
            if j is None:
                j = kd.get(rc[m - i - k:m - i])
            if j is not None:
                h[j] += 1
    else:
        rc = revcomp(s)
        for i in range(n):
            j = kd.get(s[i:i + k])
            if j is not None:
                h[j] += 1
            j = kd.get(rc[i:i + k])
            if j is not None:
                h[j] += 1
    return h


# =============================================================================
# Pass A
# =============================================================================
def _save_part(prefix, chrom, part, arrs):
    fn = os.path.join(G["tmp"], "{}_{}_{}.npz".format(prefix, G["cidx"][chrom], part))
    np.savez(fn, **arrs)
    return fn


def pass_a(chrom):
    """Scan one chromosome; once entries exceed --chunk, spill one block to disk to cap a single worker's memory."""
    t0 = time.time()
    fi = G["fidx"].get(chrom)
    res = {"chrom": chrom, "n_rec": 0, "n_r1_in": 0, "n_prim_in": 0, "n_dup_in": 0,
           "naive_f": None, "naive_n": None, "fl": None, "files": [], "n_entries": 0}
    if fi is None:
        return res
    B, minov, chunk = G["B"], G["minov"], G["chunk"]
    lin = G["lin"]
    L = 0
    if lin is not None:
        names, kd, canon, featlin_ok = lin
        L = len(names)
    naive = defaultdict(int)
    fl = []
    buf = {}

    def reset():
        buf["qa"], buf["fa"], buf["sa"] = array("Q"), array("i"), array("i")
        buf["dq"], buf["pq"] = array("Q"), array("Q")
        buf["lq"], buf["lm"] = array("Q"), array("b")
        buf["lh"] = [array("i") for _ in range(L)]

    def flush():
        q, f, s = np_u64(buf["qa"]), np_i32(buf["fa"]), np_i32(buf["sa"])
        if len(q):
            order = np.lexsort((-s, f, q))
            q, f, s = q[order], f[order], s[order]
            keep = np.ones(len(q), bool)
            keep[1:] = (q[1:] != q[:-1]) | (f[1:] != f[:-1])
            q, f, s = q[keep], f[keep], s[keep]
        arrs = {}
        for k, (qb, fb, sb) in enumerate(split_by_bucket(B, q, f, s)):
            arrs["q%d" % k], arrs["f%d" % k], arrs["s%d" % k] = qb, fb, sb
        for k, (qb,) in enumerate(split_by_bucket(B, np_u64(buf["dq"]))):
            arrs["d%d" % k] = qb
        for k, (qb,) in enumerate(split_by_bucket(B, np_u64(buf["pq"]))):
            arrs["p%d" % k] = qb
        if lin is not None:
            lqa = np_u64(buf["lq"])
            lma = np.frombuffer(buf["lm"], dtype=np.int8).copy() if len(buf["lm"]) else np.zeros(0, np.int8)
            lha = (np.stack([np_i32(x) for x in buf["lh"]], axis=1) if len(lqa)
                   else np.zeros((0, L), np.int32))
            for k, (qb, mb, hb) in enumerate(split_by_bucket(B, lqa, lma, lha)):
                arrs["lq%d" % k], arrs["lm%d" % k], arrs["lh%d" % k] = qb, mb, hb
        res["files"].append(_save_part("A", chrom, len(res["files"]), arrs))
        res["n_entries"] += int(len(q))
        reset()

    reset()
    n_rec = n_r1 = n_prim = n_dup = 0
    bam = open_bam(G["bam"])
    for r in bam.fetch(chrom):
        n_rec += 1
        flag = r.flag
        if flag & 0x804:            # unmapped / supplementary
            continue
        paired = flag & 1
        is_r1 = (flag & 0x40) or not paired
        if not is_r1 and lin is None:
            continue
        fs, fe = frag_span(r, flag)
        hits = overlaps(fi, fs, fe, minov)
        if hits is None:
            continue
        primary = not (flag & 0x100)
        q = qhash(r.query_name)
        if is_r1:
            n_r1 += 1
            s = pair_score(r, flag)
            qa, fa, sa = buf["qa"], buf["fa"], buf["sa"]
            for f in hits:
                qa.append(q)
                fa.append(f)
                sa.append(s)
            if primary:
                if flag & 0x400:
                    buf["dq"].append(q)
                    n_dup += 1
                else:
                    buf["pq"].append(q)
                    n_prim += 1
                    for f in hits:
                        naive[f] += 1
                    if paired and len(fl) < 100000:
                        t = abs(r.template_length)
                        if 0 < t < 2000:
                            fl.append(t)
        if lin is not None and primary and not (flag & 0x400):
            seq = r.query_sequence
            if seq and any(featlin_ok[f] for f in hits):
                h = count_lineage_hits(seq, kd, G["lin_k"], L, canon)
                if any(h):
                    buf["lq"].append(q)
                    buf["lm"].append(0 if is_r1 else 1)
                    for j in range(L):
                        buf["lh"][j].append(h[j])
        if len(buf["qa"]) >= chunk:
            flush()
    bam.close()
    flush()
    res.update(n_rec=n_rec, n_r1_in=n_r1, n_prim_in=n_prim, n_dup_in=n_dup,
               naive_f=np.array(list(naive.keys()), dtype=np.int64),
               naive_n=np.array(list(naive.values()), dtype=np.int64),
               fl=np.array(fl, dtype=np.int32), sec=time.time() - t0)
    return res


# =============================================================================
# Pass B
# =============================================================================
def _members(k):
    m = G["members"][k]
    if m is None:
        m = np.load(os.path.join(G["tmp"], "M_{}.npy".format(k)), mmap_mode="r")
        G["members"][k] = m
    return m


def pass_b(chrom):
    t0 = time.time()
    fi = G["fidx"].get(chrom)
    B, minov, chunk = G["B"], G["minov"], G["chunk"]
    st = {"chrom": chrom, "n_rec": 0, "n_bg": 0, "n_bg_dup": 0, "n_prim_mem_out": 0,
          "n_out_entries": 0, "files": []}
    acc = {"q": [], "s": [], "d": [], "p": [], "n": 0}
    bq, bs, bc = array("Q"), array("i"), array("b")

    def write_part():
        q = np.concatenate(acc["q"]) if acc["q"] else np.zeros(0, np.uint64)
        s = np.concatenate(acc["s"]) if acc["s"] else np.zeros(0, np.int32)
        if len(q):
            order = np.lexsort((s, q))
            q, s = q[order], s[order]
            new = np.ones(len(q), bool)
            new[1:] = (q[1:] != q[:-1]) | (s[1:] != s[:-1])
            idx = np.nonzero(new)[0]
            n = np.diff(np.concatenate([idx, [len(q)]])).astype(np.int32)
            q, s = q[idx], s[idx]
        else:
            n = np.zeros(0, np.int32)
        arrs = {}
        for k, (qb, sb, nb) in enumerate(split_by_bucket(B, q, s, n)):
            arrs["q%d" % k], arrs["s%d" % k], arrs["n%d" % k] = qb, sb, nb
        d = np.concatenate(acc["d"]) if acc["d"] else np.zeros(0, np.uint64)
        for k, (qb,) in enumerate(split_by_bucket(B, d)):
            arrs["d%d" % k] = qb
        p = np.concatenate(acc["p"]) if acc["p"] else np.zeros(0, np.uint64)
        for k, (qb,) in enumerate(split_by_bucket(B, p)):
            arrs["p%d" % k] = qb
        st["files"].append(_save_part("B", chrom, len(st["files"]), arrs))
        st["n_out_entries"] += int(len(q))
        acc.update(q=[], s=[], d=[], p=[], n=0)

    def flush():
        if not len(bq):
            return
        q = np_u64(bq)
        s = np_i32(bs)
        c = np.frombuffer(bc, dtype=np.int8).copy()
        mem = np.zeros(len(q), bool)
        bkt = (q % np.uint64(B)).astype(np.int64)
        for k in np.unique(bkt):
            ix = np.nonzero(bkt == k)[0]
            m = _members(int(k))
            if len(m) == 0:
                continue
            pos = np.searchsorted(m, q[ix])
            pos[pos >= len(m)] = len(m) - 1
            mem[ix] = (m[pos] == q[ix])
        prim = (c & 1) > 0
        dup = (c & 2) > 0
        st["n_bg"] += int(np.sum(~mem & prim & ~dup))
        st["n_bg_dup"] += int(np.sum(~mem & prim & dup))
        st["n_prim_mem_out"] += int(np.sum(mem & prim & ~dup))
        acc["q"].append(q[mem])
        acc["s"].append(s[mem])
        acc["d"].append(q[mem & prim & dup])
        acc["p"].append(q[mem & prim & ~dup])
        acc["n"] += int(mem.sum())
        del bq[:]
        del bs[:]
        del bc[:]
        if acc["n"] >= chunk:
            write_part()

    bam = open_bam(G["bam"])
    n_rec = 0
    for r in bam.fetch(chrom):
        n_rec += 1
        flag = r.flag
        if flag & 0x804:
            continue
        if (flag & 1) and not (flag & 0x40):
            continue
        if fi is not None:
            fs, fe = frag_span(r, flag)
            if overlaps(fi, fs, fe, minov) is not None:
                continue
        bq.append(qhash(r.query_name))
        bs.append(pair_score(r, flag))
        bc.append((0 if flag & 0x100 else 1) | (2 if flag & 0x400 else 0))
        if len(bq) >= 500000:
            flush()
    bam.close()
    flush()
    write_part()
    st["n_rec"] = n_rec
    st["sec"] = time.time() - t0
    return st


# =============================================================================
# Member sets
# =============================================================================
def build_members(k):
    parts = []
    for fn in G["A_files"]:
        with np.load(fn) as z:
            parts.append(z["q%d" % k])
    m = np.unique(np.concatenate(parts)) if parts else np.zeros(0, np.uint64)
    np.save(os.path.join(G["tmp"], "M_{}.npy".format(k)), m)
    return int(len(m))


# =============================================================================
# Reduce: fragment -> equivalence class
# =============================================================================
def reduce_bucket(k):
    T, wmin, F = G["T"], G["wmin"], G["F"]
    outside, req_primary = G["outside"], G["req_primary"]
    qs, fs, ss, ds, ps = [], [], [], [], []
    lqs, lms, lhs = [], [], []
    for fn in G["A_files"]:
        with np.load(fn) as z:
            qs.append(z["q%d" % k])
            fs.append(z["f%d" % k])
            ss.append(z["s%d" % k])
            ds.append(z["d%d" % k])
            ps.append(z["p%d" % k])
            if G["lin"] is not None and ("lq%d" % k) in z:
                lqs.append(z["lq%d" % k])
                lms.append(z["lm%d" % k])
                lhs.append(z["lh%d" % k])
    oq, os_, on = [], [], []
    if outside:
        for fn in G["B_files"]:
            with np.load(fn) as z:
                oq.append(z["q%d" % k])
                os_.append(z["s%d" % k])
                on.append(z["n%d" % k])
                ds.append(z["d%d" % k])
                ps.append(z["p%d" % k])
    cat = lambda L, dt: np.concatenate(L) if L else np.zeros(0, dt)  # noqa
    q, f, s = cat(qs, np.uint64), cat(fs, np.int32), cat(ss, np.int32)
    oq, os_, on = cat(oq, np.uint64), cat(os_, np.int32), cat(on, np.int32)
    dup = np.unique(cat(ds, np.uint64))
    prim = np.unique(cat(ps, np.uint64))

    st = {"bucket": k, "n_entries": int(len(q))}
    # ---- fragment filtering: must have a non-duplicate primary read1 (outside mode); otherwise at least not a duplicate
    if req_primary:
        keep = np.isin(q, prim)
        okeep = np.isin(oq, prim)
    else:
        keep = ~np.isin(q, dup)
        okeep = ~np.isin(oq, dup)
    st["n_frag_dropped"] = int(len(np.unique(q[~keep])))
    q, f, s = q[keep], f[keep], s[keep]
    oq, os_, on = oq[okeep], os_[okeep], on[okeep]

    empty = {"single": np.zeros(F + 1, np.float64), "bestu": np.zeros(F, np.int64),
             "cand": np.zeros(F, np.int64), "ec": {}, "st": st}
    if len(q) == 0:
        st.update(n_frag=0, n_single=0, n_multi=0, n_bgonly=0, n_lin_used=0)
        return empty

    order = np.lexsort((-s, f, q))
    q, f, s = q[order], f[order], s[order]
    dd = np.ones(len(q), bool)
    dd[1:] = (q[1:] != q[:-1]) | (f[1:] != f[:-1])     # (fragment, feature) pairs that may repeat across blocks
    q, f, s = q[dd], f[dd], s[dd].astype(np.float64)
    uq, start, cnt = np.unique(q, return_index=True, return_counts=True)
    nfr = len(uq)
    row = np.repeat(np.arange(nfr), cnt)
    best_in = np.maximum.reduceat(s, start)
    s_best = best_in.copy()
    bgw = np.zeros(nfr)
    if outside and len(oq):
        pos = np.searchsorted(uq, oq)
        pos[pos >= nfr] = nfr - 1
        ok = uq[pos] == oq
        pos, os_f, on_f = pos[ok], os_[ok].astype(np.float64), on[ok].astype(np.float64)
        best_out = np.full(nfr, -np.inf)
        np.maximum.at(best_out, pos, os_f)
        s_best = np.maximum(best_in, best_out)
        bgw = np.bincount(pos, weights=on_f * np.exp((os_f - s_best[pos]) / T), minlength=nfr)
    w = np.exp((s - s_best[row]) / T)

    # ---- optional: lineage P factor
    n_lin_used = 0
    if G["lin"] is not None and lqs:
        names, _, _, featlin_ok = G["lin"]
        Ln = len(names)
        lq_ = cat(lqs, np.uint64)
        lm_ = cat(lms, np.int8)
        lh_ = np.concatenate(lhs, axis=0) if lhs else np.zeros((0, Ln), np.int32)
        if len(lq_):
            o = np.lexsort((lm_, lq_))
            lq_, lm_, lh_ = lq_[o], lm_[o], lh_[o]
            firstm = np.ones(len(lq_), bool)
            firstm[1:] = (lq_[1:] != lq_[:-1]) | (lm_[1:] != lm_[:-1])
            lq_, lh_ = lq_[firstm], lh_[firstm]
            ul, ls = np.unique(lq_, return_index=True)
            H = np.add.reduceat(lh_, ls, axis=0).astype(np.float64)
            tot = H.sum(axis=1)
            a = G["lin_alpha"]
            P = (H + a) / (tot[:, None] + Ln * a)
            top = P.max(axis=1)
            expo = np.where(tot < G["lin_min_hits"], G["pw_none"],
                            np.where(top >= G["lin_spec"], G["pw_spec"], G["pw_conf"]))
            P = np.maximum(P, G["p_floor"])
            fpos = np.searchsorted(ul, uq)
            fpos[fpos >= len(ul)] = len(ul) - 1
            has = (ul[fpos] == uq) & (expo[fpos] > 0) if len(ul) else np.zeros(nfr, bool)
            if has.any():
                n_lin_used = int(has.sum())
                flin = G["feat_lin_idx"]              # feature -> lineage column (-1 = unknown)
                ent = has[row]
                rr = fpos[row[ent]]
                li = flin[f[ent]]
                pv = np.where(li >= 0, P[rr, np.maximum(li, 0)], 1.0 / Ln)
                w[ent] = w[ent] * pv ** expo[rr]
                bgw[has] = bgw[has] * (1.0 / Ln) ** expo[fpos[has]]
                wmax = np.maximum(np.maximum.reduceat(w, start), bgw)
                w = w / wmax[row]
                bgw = bgw / wmax

    keep = w >= wmin
    ncand = np.bincount(row[keep], minlength=nfr)
    hasbg = bgw >= wmin
    wq = np.maximum(np.round(w * QS), 1).astype(np.int32)
    bgq = np.round(bgw * QS).astype(np.int64)

    single_row = (ncand == 1) & ~hasbg
    bgonly_row = (ncand == 0) & hasbg
    multi_row = ~single_row & ~bgonly_row & ((ncand + hasbg) >= 1)

    single = np.zeros(F + 1, np.float64)
    es = keep & single_row[row]
    single[:F] = np.bincount(f[es], minlength=F)
    single[F] = float(bgonly_row.sum())

    top = keep & (wq >= QS)
    ntop = np.bincount(row[top], minlength=nfr)
    bu_row = (ntop == 1) & (bgq < QS)
    bestu = np.bincount(f[top & bu_row[row]], minlength=F).astype(np.int64)
    cand = np.bincount(f[keep & (w >= 0.5)], minlength=F).astype(np.int64)

    ec = defaultdict(int)
    em = keep & multi_row[row]
    ef, ew = f[em], wq[em]
    mrows = np.nonzero(multi_row)[0]
    lens = ncand[mrows]
    offs = np.concatenate([[0], np.cumsum(lens)])
    bgm = bgq[mrows] * hasbg[mrows]
    for i in range(len(mrows)):
        a, b = offs[i], offs[i + 1]
        ec[(ef[a:b].tobytes(), ew[a:b].tobytes(), int(bgm[i]))] += 1

    st.update(n_frag=int(nfr), n_single=int(single_row.sum()), n_multi=int(multi_row.sum()),
              n_bgonly=int(bgonly_row.sum()), n_lin_used=n_lin_used,
              n_best_outside=int(np.sum(bgw >= 1.0 - 1e-9)))
    return {"single": single, "bestu": bestu, "cand": cand, "ec": dict(ec), "st": st}


# =============================================================================
# EM
# =============================================================================
class ECData(object):
    def __init__(self, F, single, ec, eff):
        self.F = F
        self.single = single                    # length F+1 (the last one = background)
        self.eff = eff
        keys = list(ec.keys())
        self.ec_count = np.array([ec[k] for k in keys], dtype=np.float64)
        idx_parts, w_parts, lens = [], [], []
        for fb, wb, bg in keys:
            fi = np.frombuffer(fb, dtype=np.int32)
            wi = np.frombuffer(wb, dtype=np.int32)
            if bg > 0:
                fi = np.append(fi, F)
                wi = np.append(wi, bg)
            idx_parts.append(fi)
            w_parts.append(wi)
            lens.append(len(fi))
        if keys:
            self.indices = np.concatenate(idx_parts).astype(np.int64)
            self.data = np.concatenate(w_parts).astype(np.float64) / QS
            lens = np.array(lens, dtype=np.int64)
        else:
            self.indices = np.zeros(0, np.int64)
            self.data = np.zeros(0)
            lens = np.zeros(0, np.int64)
        self.indptr = np.concatenate([[0], np.cumsum(lens)]).astype(np.int64)
        self.rowid = np.repeat(np.arange(len(keys)), lens)
        # ---- v4: precompute 1/eff so that each round's division becomes a multiplication
        self.eff_inv = 1.0 / np.maximum(eff, 1e-300)
        self.N = float(single.sum() + self.ec_count.sum())

    @property
    def nbytes(self):
        """Temporary data footprint (used to judge whether a spawn platform can afford parallel bootstrap)."""
        return int(self.indices.nbytes + self.data.nbytes + self.rowid.nbytes
                   + self.indptr.nbytes + self.ec_count.nbytes + self.single.nbytes)

    # ------------------------------------------------------------------ v4: single-step EM
    def _loglik(self, alpha, ec_count, single):
        """Log-likelihood of the current alpha (for convergence diagnostics / the SQUAREM monotonicity guard)."""
        rate = alpha * self.eff_inv
        ll = float(np.dot(single, np.log(np.maximum(rate, 1e-300))))
        if len(self.data):
            val = self.data * rate[self.indices]
            rs = np.add.reduceat(val, self.indptr[:-1])
            np.maximum(rs, 1e-300, out=rs)
            ll += float(np.dot(ec_count, np.log(rs)))
        return ll

    def _step(self, alpha, ec_count, single, prior):
        """One standard EM update, returns (new_alpha, log-likelihood).

        Three changes in the v4 loop body relative to v3 (semantically equivalent, only to save
        memory bandwidth):
          - rate = alpha * eff_inv      (was alpha / eff)
          - rs zero-guard via np.maximum(out=)  (was a boolean mask over 53 million elements)
          - multiply the share back into val in place  (was one more 53-million-element temporary array)
        The log-likelihood ll is also computed along the way; the SQUAREM monotonicity guard needs it.
        """
        F1 = self.F + 1
        rate = alpha * self.eff_inv
        new = single.astype(np.float64).copy()
        if len(self.data):
            val = self.data * rate[self.indices]
            rs = np.add.reduceat(val, self.indptr[:-1])
            np.maximum(rs, 1e-300, out=rs)
            val *= (ec_count / rs)[self.rowid]
            new += np.bincount(self.indices, weights=val, minlength=F1)
            ll = float(np.dot(ec_count, np.log(rs)))
        else:
            ll = 0.0
        # pinned votes also enter the likelihood (they are equivalent to fragments "with a unique candidate")
        ll += float(np.dot(single, np.log(np.maximum(rate, 1e-300))))
        new += prior
        return new, ll

    @staticmethod
    def _converged(new, alpha, tol, min_iter, it):
        """Exactly the same convergence criterion as v3 (compare relative change only on components with new>1)."""
        if it < min_iter:
            return False
        sel = new > 1.0
        if not sel.any():
            return True
        return float(np.max(np.abs(new[sel] - alpha[sel]) / new[sel])) < tol

    def em(self, alpha0, ec_count=None, single=None, tol=1e-3, max_iter=2000, min_iter=20,
           prior=1e-8, accelerate=True, trace=False, tag=""):
        """EM iteration. With accelerate=True, SQUAREM squared extrapolation is enabled.

        SQUAREM (Varadhan & Roland 2008) accelerates the fixed-point iteration x_{k+1}=F(x_k):
            one step r = F(x)-x, two steps v = F(F(x))-F(x)-r
            extrapolate along the Cauchy direction  xs = x - 2a*r + a^2*v,  a = -||r||/||v||
            take one more EM step from xs to get the candidate x3; accept it only when
            ll(x3) >= ll(F(x)),
            otherwise fall back to the ordinary EM step x1.
        This keeps the EM log-likelihood monotonically non-decreasing (it "accelerates" rather than
        "switches algorithm"),
        the endpoint is still the same MLE, only the iteration count drops sharply.
        Each round does at most 3 EM steps, so it only pays off when the round count drops below 1/3.
        """
        F1 = self.F + 1
        ec_count = self.ec_count if ec_count is None else ec_count
        single = self.single if single is None else single
        alpha = alpha0.copy()
        has_multi = len(self.data) > 0
        it = 0
        converged = False
        for it in range(1, max_iter + 1):
            x1, ll1 = self._step(alpha, ec_count, single, prior)
            if self._converged(x1, alpha, tol, min_iter, it):
                alpha = x1
                converged = True
                if trace:
                    log("[03]   EM{} converged at iteration {} ll={:.10g}".format(tag, it, ll1))
                break
            if not (accelerate and has_multi):
                alpha = x1
                continue
            x2, _ = self._step(x1, ec_count, single, prior)
            r = x1 - alpha
            v = (x2 - x1) - r
            rn = float(np.linalg.norm(r))
            vn = float(np.linalg.norm(v))
            if not (vn > 0.0) or not (rn > 0.0):
                alpha = x1
                continue
            a_sq = -rn / vn                      # Cauchy step length (negative)
            a_sq = min(-1e-6, max(a_sq, -1e3))   # clamp it so the extrapolation cannot blow up
            xs = alpha - 2.0 * a_sq * r + a_sq * a_sq * v
            np.maximum(xs, 1e-300, out=xs)
            x3, ll3 = self._step(xs, ec_count, single, prior)
            if ll3 >= ll1:
                alpha = x3
            else:
                alpha = x1
            if trace and (it <= 3 or it % 25 == 0):
                log("[03]   EM{} iteration {:4d} ll={:.10g} a={:.3f}".format(
                    tag, it, max(ll3, ll1), a_sq))
        self.last_ll = self._loglik(alpha, ec_count, single)
        self.last_converged = bool(converged)
        if not converged:
            log("[03] [WARN] EM{} reached the iteration cap of {} rounds without satisfying tol={}; "
                "the result may not be fully converged (raise --max-iter or relax --tol)".format(tag, max_iter, tol))
        return alpha, it


# =============================================================================
# v4: parallel bootstrap
# =============================================================================
_BOOT_CFG = None


def _boot_init(cfg):
    global _BOOT_CFG
    _BOOT_CFG = cfg


def _boot_worker(draw):
    ec, alpha0, em_kw = _BOOT_CFG
    F = ec.F
    a, _ = ec.em(alpha0, ec_count=draw[F + 1:], single=draw[:F + 1], **em_kw)
    return a[:F]


def _run_bootstrap(draws, ec, alpha0, em_kw, nthreads):
    """Run the bootstrap resampling. With nthreads>1, parallelize over processes (zero-copy ECData sharing under fork)."""
    if nthreads <= 1 or len(draws) <= 1:
        _boot_init((ec, alpha0, em_kw))
        return [_boot_worker(d) for d in draws]
    has_fork = "fork" in mp.get_all_start_methods()
    if not has_fork and ec.nbytes > 3e8:
        log("[03] [NOTE] this platform has no fork, so parallel bootstrap would have to copy the temporary data "
            "into every child process; "
            "the current data volume of {:.1f} GB is too large -> automatically degraded to 1 thread. "
            "Linux has fork and will parallelize automatically.".format(ec.nbytes / 1e9))
        _boot_init((ec, alpha0, em_kw))
        return [_boot_worker(d) for d in draws]
    ctx = mp.get_context("fork" if has_fork else "spawn")
    if has_fork:
        _boot_init((ec, alpha0, em_kw))          # children inherit the parent's memory, no pickling needed
        with ctx.Pool(nthreads) as pool:
            return pool.map(_boot_worker, draws)
    with ctx.Pool(nthreads, initializer=_boot_init, initargs=((ec, alpha0, em_kw),)) as pool:
        return pool.map(_boot_worker, draws)


# =============================================================================
# main
# =============================================================================
def main():
    ap = argparse.ArgumentParser(description="ATAC-Quant v4 Step 03: EC-EM quantification")
    ap.add_argument("--bam", required=True, help="BAM that keeps -k multi-alignments and carries duplicate flags (needs .bai)")
    ap.add_argument("--features", required=True, help="features.tsv written by step 02")
    ap.add_argument("--sample", required=True)
    ap.add_argument("-o", "--outdir", required=True)
    ap.add_argument("-t", "--threads", type=int, default=8)
    ap.add_argument("--tmpdir", default=None, help="parent directory for temporary files (default outdir); atacq_tmp_<sample>/ is created inside it and needs plenty of space")
    ap.add_argument("--buckets", type=int, default=128, help="number of fragment hash buckets (more buckets means less memory per bucket)")
    ap.add_argument("--chunk", type=int, default=20000000,
                    help="how many entries a worker accumulates while scanning before spilling to disk (caps memory; about 0.5 GB per 10 million entries)")
    ap.add_argument("--reduce-threads", type=int, default=0,
                    help="parallelism of the Reduce stage (the memory hog; 0 = min(threads, 16))")
    ap.add_argument("--min-overlap", type=int, default=1, help="minimum overlap in bp between a fragment and a feature")
    ap.add_argument("--as-temp", type=float, default=2.0,
                    help="compatibility weight temperature T: w = exp((AS - AS_best)/T). bowtie2 charges 6 points "
                         "for one high-quality mismatch, so at T=2 a difference of 1 mismatch gives a weight of about 0.05")
    ap.add_argument("--w-min", type=float, default=0.01, help="candidates below this weight are ignored")
    ap.add_argument("--outside", choices=["model", "ignore"], default="model",
                    help="model = scan alignments outside features, build the background component and "
                         "deduplicate by primary alignment (recommended); "
                         "ignore = only look at alignments inside features (faster, but duplicates whose "
                         "primary alignment lies outside features cannot be identified)")
    ap.add_argument("--frag-len", type=float, default=0.0, help="mean fragment length; 0 = estimate from the BAM")
    ap.add_argument("--genome-size", type=float, default=0.0, help="0 = total reference length from the BAM header")
    ap.add_argument("--tol", type=float, default=1e-6,
                    help="EM convergence criterion (relative change). v4 defaults to 1e-6 (v3 used 1e-3): "
                         "SQUAREM stops earlier under the same criterion, so only at 1e-6 does it reach precision "
                         "equal to or better than v3; "
                         "to reproduce v3 exactly use --em-accelerate off --tol 1e-3")
    ap.add_argument("--max-iter", type=int, default=2000)
    ap.add_argument("--min-iter", type=int, default=20)
    ap.add_argument("--bootstrap", type=int, default=0, help="number of bootstrap replicates (0 = skip; 20 recommended)")
    ap.add_argument("--seed", type=int, default=2026)
    # ---- new in v4: convergence acceleration and parallel bootstrap
    ap.add_argument("--em-accelerate", choices=["squarem", "off"], default="squarem",
                    help="EM convergence acceleration: squarem = SQUAREM squared extrapolation (v4 default, "
                         "sharply reduces the iteration count); "
                         "off = disabled (line-for-line identical to v3 behaviour)")
    ap.add_argument("--em-trace", action="store_true",
                    help="print the EM log-likelihood of every iteration (for checking the acceleration effect)")
    ap.add_argument("--em-prior", type=float, default=1e-8,
                    help="per-round zero-avoiding prior added to alpha (fixed at 1e-8 in v3)")
    ap.add_argument("--bootstrap-threads", type=int, default=0,
                    help="number of parallel bootstrap processes; 0 = auto min(bootstrap, 8); 1 = sequential (v3 behaviour)")
    ap.add_argument("--bootstrap-max-iter", type=int, default=500,
                    help="EM iteration cap for each bootstrap replicate (fixed at 500 in v3)")
    # lineage k-mer (optional, off by default)
    ap.add_argument("--lineage-kmers", nargs="*", default=None, metavar="NAME:PATH",
                    help="optional lineage k-mer library (e.g. SOF:B48_k15.fasta). In one evaluated project 99.76%% of fragments had no hit, "
                         "and the LQ root-cause report shows it carries no information for fragments that multi-map "
                         "between copies; off by default")
    ap.add_argument("--lineage-k", type=int, default=15)
    ap.add_argument("--lin-alpha", type=float, default=1.0)
    ap.add_argument("--lin-min-hits", type=int, default=3)
    ap.add_argument("--lin-specific", type=float, default=0.9)
    ap.add_argument("--p-weight-specific", type=float, default=1.0)
    ap.add_argument("--p-weight-conflict", type=float, default=0.0)
    ap.add_argument("--p-weight-none", type=float, default=0.0)
    ap.add_argument("--p-floor", type=float, default=1e-3)
    ap.add_argument("--chroms", default=None,
                    help="process only these chromosomes (comma separated), for a quick trial run on real data first; "
                         "features on the other chromosomes are counted as 0")
    ap.add_argument("--keep-tmp", action="store_true")
    args = ap.parse_args()

    t_start = time.time()
    os.makedirs(args.outdir, exist_ok=True)
    # only read/write/clean inside the subdirectory we create ourselves; never delete the user-supplied --tmpdir itself
    tmp = os.path.join(args.tmpdir or args.outdir, "atacq_tmp_{}".format(args.sample))
    if os.path.exists(tmp):
        shutil.rmtree(tmp)
    os.makedirs(tmp)

    feats = load_features(args.features)
    F = len(feats)
    if F == 0:
        die("features.tsv is empty")
    fidx = build_feature_index(feats)
    refs = bam_references(args.bam)
    ref_names = [r[0] for r in refs]
    missing = [c for c in fidx if c not in set(ref_names)]
    if missing:
        log("[WARN] {} feature chromosomes are not in the BAM header (e.g. {}) -- genome version mismatch?".format(
            len(missing), missing[:3]))
    genome_size = args.genome_size or float(sum(r[1] for r in refs))
    log("[03] {} | features {:,} | chromosomes {:,} | threads {} | outside={}".format(
        args.sample, F, len(refs), args.threads, args.outside))

    lin = None
    feat_lin_idx = np.full(F, -1, np.int64)
    if args.lineage_kmers:
        names, kd, canon = load_lineage_kmers(args.lineage_kmers, args.lineage_k)
        featlin_ok = bytearray(F)
        for i, ft in enumerate(feats):
            a = ft.get("Ancestry", "").upper()
            if a in names:
                feat_lin_idx[i] = names.index(a)
            if ft.get("Type") in ("called", "projected"):
                featlin_ok[i] = 1
        lin = (names, kd, canon, featlin_ok)
        log("[03] lineage k-mer: {} | canonical={}".format(names, canon))

    G.update(fidx=fidx, bam=args.bam, B=args.buckets, chunk=args.chunk, minov=args.min_overlap, tmp=tmp,
             lin=lin, lin_k=args.lineage_k, cidx={c: i for i, c in enumerate(ref_names)},
             T=args.as_temp, wmin=args.w_min, F=F, outside=(args.outside == "model"),
             req_primary=(args.outside == "model"), members=[None] * args.buckets,
             feat_lin_idx=feat_lin_idx, lin_alpha=args.lin_alpha, lin_min_hits=args.lin_min_hits,
             lin_spec=args.lin_specific, pw_spec=args.p_weight_specific,
             pw_conf=args.p_weight_conflict, pw_none=args.p_weight_none, p_floor=args.p_floor)

    # run the large chromosomes first for a more even load
    reflen = dict(refs)
    chroms_a = sorted([c for c in ref_names if c in fidx], key=lambda c: -reflen[c])
    chroms_b = sorted(ref_names, key=lambda c: -reflen[c])
    if args.chroms:
        keep = set(x.strip() for x in args.chroms.split(",") if x.strip())
        chroms_a = [c for c in chroms_a if c in keep]
        chroms_b = [c for c in chroms_b if c in keep]
        genome_size = float(sum(reflen[c] for c in chroms_b))
        log("[03] trial mode: processing only {} chromosomes ({})".format(len(chroms_b), ",".join(chroms_b[:5])))

    def pmap(func, items, label, nproc=None, ordered=False):
        out = []
        nproc = nproc or args.threads
        if nproc > 1 and len(items) > 1:
            # fork: workers inherit G (since Python 3.14 Linux no longer defaults to fork);
            # platforms without fork (Windows) use spawn + inject G into the workers
            _ctx = mp.get_context("fork") if "fork" in mp.get_all_start_methods() else mp.get_context("spawn")
            if _ctx.get_start_method() == "fork":
                pool_cm = _ctx.Pool(nproc)
            else:
                pool_cm = _ctx.Pool(nproc, initializer=_spawn_init_G, initargs=(dict(G),))
            with pool_cm as pool:
                # v4: ordered=True uses imap to preserve order. The Reduce stage must preserve order, otherwise
                # the floating-point summation order of each bucket's result varies with scheduling, which makes
                # the last few digits of p wobble; multinomial is discrete sampling, so the last digit flips the
                # sampled result -> EM_SD is not reproducible (the old v3 problem).
                it = pool.imap(func, items) if ordered else pool.imap_unordered(func, items)
                for i, r in enumerate(it):
                    out.append(r)
                    if (i + 1) % 20 == 0 or i + 1 == len(items):
                        log("[03]   {} {}/{}".format(label, i + 1, len(items)))
        else:
            for r in map(func, items):
                out.append(r)
        return out

    # ---------------- Pass A
    log("[03] Pass A: alignments inside features")
    resA = pmap(pass_a, chroms_a, "Pass A")
    G["A_files"] = [fn for r in resA for fn in r.get("files", [])]
    naive = np.zeros(F, np.int64)
    fl_all = []
    sA = defaultdict(int)
    for r in resA:
        for k in ("n_rec", "n_r1_in", "n_prim_in", "n_dup_in", "n_entries"):
            sA[k] += int(r.get(k, 0) or 0)
        if r.get("naive_f") is not None and len(r["naive_f"]):
            np.add.at(naive, r["naive_f"], r["naive_n"])
        if r.get("fl") is not None and len(r["fl"]):
            fl_all.append(r["fl"])
    fl_all = np.concatenate(fl_all) if fl_all else np.zeros(0)
    mean_fl = args.frag_len or (float(fl_all.mean()) if len(fl_all) else 0.0)
    log("[03] Pass A: read1 records inside features {:,}, entries {:,}; primary fragments {:,} (duplicates {:,}); "
        "mean fragment length {:.1f}".format(sA["n_r1_in"], sA["n_entries"], sA["n_prim_in"],
                                            sA["n_dup_in"], mean_fl))

    # ---------------- members + Pass B
    sB = defaultdict(int)
    G["B_files"] = []
    if G["outside"]:
        n_mem = sum(pmap(build_members, list(range(args.buckets)), "members"))
        log("[03] fragments seen inside features (including duplicates) {:,}".format(n_mem))
        log("[03] Pass B: alignments outside features + pure background fragments")
        resB = pmap(pass_b, chroms_b, "Pass B")
        G["B_files"] = [fn for r in resB for fn in r.get("files", [])]
        for r in resB:
            for k in ("n_rec", "n_bg", "n_bg_dup", "n_prim_mem_out", "n_out_entries"):
                sB[k] += int(r.get(k, 0) or 0)
        log("[03] Pass B: pure background fragments {:,} (plus duplicates {:,}); out-of-feature alignment entries "
            "of member fragments {:,}".format(sB["n_bg"], sB["n_bg_dup"], sB["n_out_entries"]))

    # ---------------- Reduce
    log("[03] Reduce: assembling equivalence classes")
    rthreads = args.reduce_threads or min(args.threads, 16)
    resR = pmap(reduce_bucket, list(range(args.buckets)), "reduce", nproc=rthreads,
                ordered=True)      # v4: preserve order -> reproducible results (see the note inside pmap)
    single = np.zeros(F + 1)
    bestu = np.zeros(F, np.int64)
    cand = np.zeros(F, np.int64)
    ec = defaultdict(int)
    sR = defaultdict(int)
    for r in resR:
        single += r["single"]
        bestu += r["bestu"]
        cand += r["cand"]
        for key, c in r["ec"].items():
            ec[key] += c
        for k2, v in r["st"].items():
            if k2 != "bucket":
                sR[k2] += v
    single[F] += sB["n_bg"]
    n_unique = single[:F].astype(np.int64)
    log("[03] fragments: entering features {:,} | unique candidate {:,} ({:.1%}) | multiple candidates {:,} "
        "({:.1%}) | background only {:,} | equivalence classes {:,}".format(
            sR["n_frag"], sR["n_single"], sR["n_single"] / max(sR["n_frag"], 1),
            sR["n_multi"], sR["n_multi"] / max(sR["n_frag"], 1), sR["n_bgonly"], len(ec)))

    # ---------------- EM
    lens = np.array([ft["End"] - ft["Start"] for ft in feats], dtype=np.float64)
    if mean_fl > 0:
        eff_f = lens + mean_fl - 1.0
    else:
        eff_f = lens.copy()
    proc = set(chroms_b)
    on = np.array([ft["Chrom"] in proc for ft in feats], dtype=bool)
    eff_bg = max(genome_size - float(eff_f[on].sum()), 1.0)
    eff = np.concatenate([np.maximum(eff_f, 1.0), [eff_bg]])
    data = ECData(F, single, ec, eff)
    N = data.N
    log("[03] EM: {:,.0f} fragments, {:,} features + background".format(N, F))
    alpha0 = eff * (N / eff.sum())
    # use the unique evidence as a gentle starting point (does not change the MLE, only speeds up convergence)
    alpha0 = 0.5 * alpha0 + 0.5 * np.maximum(single, 0) + 1e-6
    alpha0 *= N / alpha0.sum()
    t_em0 = time.time()
    alpha, n_it = data.em(alpha0, tol=args.tol, max_iter=args.max_iter, min_iter=args.min_iter,
                          prior=args.em_prior,
                          accelerate=(args.em_accelerate == "squarem"),
                          trace=args.em_trace, tag="(main)")
    log("[03] EM converged: {} rounds (acceleration={}); background share {:.2%}; elapsed {:.1f} s".format(
        n_it, args.em_accelerate, alpha[F] / max(N, 1), time.time() - t_em0))

    sd = np.full(F, np.nan)
    if args.bootstrap > 0:
        rng = np.random.default_rng(args.seed)
        base = np.concatenate([data.single, data.ec_count])
        p = base / base.sum()
        # ---- v4: generate all resamples first, then run them in parallel (v3 used a sequential for loop, about 76% of the total runtime)
        draws = [rng.multinomial(int(round(N)), p).astype(np.float64)
                 for _ in range(args.bootstrap)]
        nt = args.bootstrap_threads
        if nt <= 0:
            nt = min(args.bootstrap, 8)
        nt = max(1, min(nt, args.bootstrap))
        em_kw = dict(tol=args.tol, max_iter=args.bootstrap_max_iter, min_iter=5,
                     prior=args.em_prior,
                     accelerate=(args.em_accelerate == "squarem"))
        log("[03] bootstrap: {} resamples x at most {} EM rounds; {} parallel processes".format(
            args.bootstrap, args.bootstrap_max_iter, nt))
        t_b0 = time.time()
        results = _run_bootstrap(draws, data, alpha, em_kw, nt)
        mean = np.zeros(F)
        m2 = np.zeros(F)
        for b, x in enumerate(results):
            d1 = x - mean
            mean += d1 / (b + 1)
            m2 += d1 * (x - mean)
            if (b + 1) % 5 == 0 or (b + 1) == args.bootstrap:
                log("[03]   bootstrap {}/{} done (elapsed {:.1f} minutes)".format(
                    b + 1, args.bootstrap, (time.time() - t_b0) / 60.0))
        sd = np.sqrt(m2 / max(args.bootstrap - 1, 1))
        log("[03] bootstrap all done, elapsed {:.1f} minutes".format((time.time() - t_b0) / 60.0))

    # ---------------- output
    est = alpha[:F]
    qpath = os.path.join(args.outdir, "{}.quant.tsv".format(args.sample))
    with open(qpath, "w") as fo:
        fo.write("FeatureID\tLength\tEffLen\tEM_Count\tEM_SD\tN_Unique\tN_BestUnique\t"
                 "N_Candidate\tN_Primary\tResolved_Frac\tDensity_per_kb\n")
        for i, ft in enumerate(feats):
            c = est[i]
            rf = min(1.0, bestu[i] / c) if c > 0.5 else float("nan")
            fo.write("{}\t{}\t{:.0f}\t{:.3f}\t{}\t{}\t{}\t{}\t{}\t{}\t{:.4f}\n".format(
                ft["FeatureID"], int(lens[i]), eff[i], c,
                "NA" if np.isnan(sd[i]) else "{:.3f}".format(sd[i]),
                int(n_unique[i]), int(bestu[i]), int(cand[i]), int(naive[i]),
                "NA" if np.isnan(rf) else "{:.4f}".format(rf), c / eff[i] * 1000.0))

    n_frag_total = sA["n_prim_in"] + sB["n_prim_mem_out"] + sB["n_bg"] if G["outside"] else None
    stats = {
        "sample": args.sample,
        "bam": args.bam,
        "features": args.features,
        "n_features": F,
        "mean_frag_len": mean_fl,
        "genome_size": genome_size,
        "passA": dict(sA),
        "passB": dict(sB),
        "reduce": {k: int(v) for k, v in sR.items()},
        "n_equivalence_classes": len(ec),
        "n_fragments_total_nondup": n_frag_total,
        "n_fragments_in_features": int(sR["n_frag"]),
        "n_fragments_em": N,
        "em_iterations": int(n_it),
        "em_accelerate": args.em_accelerate,
        "em_seconds": round(float(time.time() - t_em0), 2),
        "em_converged": bool(getattr(data, "last_converged", False)),
        "em_loglik": float(getattr(data, "last_ll", float("nan"))),
        "em_sum_features": float(est.sum()),
        "em_background": float(alpha[F]),
        "conservation_check": {
            "sum_alpha": float(alpha.sum()),
            "n_used": N,
            "ok": bool(abs(alpha.sum() - N) <= max(1.0, 1e-4 * N) + (F + 1) * 1e-6),
        },
        "member_primary_check": {
            "fragments_in_features_after_filter": int(sR["n_frag"]),
            "primary_read1_records": int(sA["n_prim_in"] + sB["n_prim_mem_out"]) if G["outside"] else None,
        },
        "frac_multi_candidate": sR["n_multi"] / max(sR["n_frag"], 1),
        "bootstrap": args.bootstrap,
        "params": {k: (v if not isinstance(v, list) else list(v)) for k, v in vars(args).items()},
        "runtime_sec": time.time() - t_start,
    }
    write_json(stats, os.path.join(args.outdir, "{}.quant_stats.json".format(args.sample)))
    if not args.keep_tmp:
        shutil.rmtree(tmp, ignore_errors=True)
    log("[03] done: {} ({:.1f} minutes)".format(qpath, (time.time() - t_start) / 60))


if __name__ == "__main__":
    main()
