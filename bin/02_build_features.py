#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
ATAC-Quant v4 | Step 02: define the quantitative features -- Homeologous Regulatory Unit (HRU)

Why this step (the problem with v1):
  v1 ran EM only between the peaks called independently on each copy. But the primary alignments
  of homeologous reads are distributed at random among the 14 copies: a genuinely open copy may
  fail to produce a peak because its signal is diluted, and a closed copy may show a peak from
  mismapped reads. EM can only redistribute among peaks that were already called, so a missed
  copy is always 0.

What v2 does:
  1. Anchor each peak to the nearest TSS gene (TSS -3kb/+1kb, same as v1), and write the anchor
     table for every peak (v1 kept only the clusters with >=2 peaks; the rest were labelled
     UNKNOWN / Intergenic).
  2. For peaks anchored to an allele-cluster gene, project them by sequence (k-mer diagonal
     anchoring) onto the promoter window of every other copy in the same cluster and locate the
     homeologous interval -- whether or not that copy has a called peak.
     If no homeologous sequence is found (that copy lacks this element) no feature is built:
     reads would not align there in the first place.
  3. The projected intervals and the original peaks are merged by overlap; one connected set =
     one HRU group (the corresponding intervals of the same regulatory element on every copy).
     The interval on each copy = one quantitative feature, whose on/off state is decided by the
     EM in Step 03.
  4. Peaks that are not on an allele cluster (distal / non-cluster gene promoters) become
     independent features (orphan) and also enter EM (to handle multi-mapping correctly), but
     they get no homeologous projection.

Output (--outdir):
  features.tsv      every quantitative feature (the columns of EM): FeatureID Chrom Start End Length GroupID
                    ClusterID Gene Ancestry Type PeakIDs ProjCoverage
                    Type: called / projected / orphan_noncluster / orphan_distal
  features.bed
  groups.tsv        HRU group summary (copies per group, lineage composition, called/projected counts)
  peak_anchor.tsv   per-peak anchor gene, distance to TSS, cluster, lineage, class
  projections.tsv   result of every projection (number of anchors, coverage), for QC
  features_stats.json
"""
import os
import sys
import bisect
import argparse
import multiprocessing as mp
from collections import defaultdict, Counter

sys.path.insert(0, os.path.join(os.path.dirname(os.path.abspath(__file__)), "..", "lib"))
from atacq_common import (log, die, Fasta, revcomp, parse_gff_genes, load_clusters,  # noqa
                          load_ancestry, parse_lineage_map, read_peaks, write_json,
                          DEFAULT_LINEAGE_MAP, FEATURE_COLS)

# ----------------------------------------------------------------- worker state
W = {}


def init_worker(genome, genes, args):
    W["fa"] = Fasta(genome)
    W["genes"] = genes
    W["args"] = args


def tss_of(g):
    chrom, s, e, strand = g
    return s if strand == "+" else e - 1


def window_of(g, up, down, pad, chrom_len):
    chrom, s, e, strand = g
    t = tss_of(g)
    if strand == "+":
        w0, w1 = t - up - pad, t + down + pad + 1
    else:
        w0, w1 = t - down - pad, t + up + pad + 1
    return max(0, w0), min(chrom_len, w1)


def kmer_index(seq, k, cap):
    """k-mer -> list of positions; k-mers occurring more than cap times in the window (low
    complexity / tandem repeats) are discarded."""
    idx = {}
    rep = set()
    check_n = "N" in seq
    for j in range(len(seq) - k + 1):
        km = seq[j:j + k]
        if check_n and "N" in km:
            continue
        if km in rep:
            continue
        lst = idx.get(km)
        if lst is None:
            idx[km] = [j]
        else:
            lst.append(j)
            if len(lst) > cap:
                rep.add(km)
                del idx[km]
    return idx


def project(S, idx, k, band, min_anchor, min_cov):
    """
    Project the source sequence S (already in the source gene's orientation) onto the target
    window (already in the target gene's orientation).
    Main diagonal = the (j - i) shared by the most k-mers; anchors within +/-band are kept
    (tolerating small indels).
    Returns (start_w, end_w, n_anchor, coverage) or None.
    """
    L = len(S)
    if L < k:
        return None
    diag = Counter()
    hits = []
    check_n = "N" in S
    for i in range(L - k + 1):
        km = S[i:i + k]
        if check_n and "N" in km:
            continue
        js = idx.get(km)
        if not js:
            continue
        for j in js:
            diag[j - i] += 1
            hits.append((i, j))
    if not diag:
        return None
    D = diag.most_common(1)[0][0]
    best = {}
    for i, j in hits:
        dd = abs((j - i) - D)
        if dd <= band:
            if i not in best or dd < abs((best[i] - i) - D):
                best[i] = j
    if len(best) < min_anchor:
        return None
    covered = 0
    last_end = -1
    for i in sorted(best):
        s0 = max(i, last_end)
        e0 = i + k
        if e0 > s0:
            covered += e0 - s0
        last_end = max(last_end, e0)
    cov = covered / float(L)
    if cov < min_cov:
        return None
    i0 = min(best)
    i1 = max(best)
    start_w = best[i0] - i0
    end_w = best[i1] + (L - i1)
    if end_w <= start_w:
        return None
    return start_w, end_w, len(best), cov


def project_cluster(task):
    """One allele cluster: project each anchored peak of the cluster onto the other copies."""
    cid, gene_list, peaks = task
    fa, genes, a = W["fa"], W["genes"], W["args"]
    out = []
    win = {}
    for g in gene_list:
        gi = genes[g]
        clen = fa.length(gi[0]) if gi[0] in fa.index else 0
        if clen == 0:
            continue
        w0, w1 = window_of(gi, a.tss_up, a.tss_down, a.proj_pad, clen)
        seq = fa.fetch(gi[0], w0, w1)
        if gi[3] == "-":
            seq = revcomp(seq)
        win[g] = (w0, w1, kmer_index(seq, a.proj_k, a.proj_cap))
    for pname, chrom, s, e, src in peaks:
        S = fa.fetch(chrom, s, e)
        if not S:
            continue
        if genes[src][3] == "-":
            S = revcomp(S)
        for t in gene_list:
            if t == src or t not in win:
                continue
            w0, w1, idx = win[t]
            r = project(S, idx, a.proj_k, a.proj_band, a.proj_min_anchors, a.proj_min_cov)
            if r is None:
                out.append((cid, pname, src, t, None, None, None, 0, 0.0))
                continue
            sw, ew, na, cov = r
            tchrom, _, _, tstrand = genes[t]
            if tstrand == "+":
                gs, ge = w0 + sw, w0 + ew
            else:
                gs, ge = w1 - ew, w1 - sw
            gs, ge = max(gs, w0), min(ge, w1)
            if ge <= gs:
                out.append((cid, pname, src, t, None, None, None, 0, 0.0))
                continue
            out.append((cid, pname, src, t, tchrom, gs, ge, na, cov))
    return out


# ----------------------------------------------------------------- union-find
class UF(object):
    def __init__(self, n):
        self.p = list(range(n))

    def find(self, x):
        p = self.p
        while p[x] != x:
            p[x] = p[p[x]]
            x = p[x]
        return x

    def union(self, a, b):
        ra, rb = self.find(a), self.find(b)
        if ra != rb:
            self.p[max(ra, rb)] = min(ra, rb)


def main():
    ap = argparse.ArgumentParser(description="ATAC-Quant v4 Step 02: build HRU features")
    ap.add_argument("--peaks", required=True, help="consensus_peaks.bed or narrowPeak")
    ap.add_argument("--gff", required=True)
    ap.add_argument("--clusters", required=True, help="GraphAllele my_cluster.tsv")
    ap.add_argument("--ancestry", default="",
                    help="Clean_Gene_Ancestry_Mapping.tsv. v4.3: may be omitted when ploidy=auto/diploid (the lineage layer is not produced)")
    ap.add_argument("--ploidy", choices=["diploid", "auto", "allo"], default="allo",
                    help="v4.3: ploidy mode (only used for logging and guardrail hints; the projection logic is the same in all three modes)")
    ap.add_argument("--genome", required=True, help="unmasked genome FASTA (needs .fai)")
    ap.add_argument("--lineage-map", default="")   # v3: species mapping is off by default; supplied by the config / driver
    ap.add_argument("-o", "--outdir", required=True)
    ap.add_argument("-t", "--threads", type=int, default=8)
    ap.add_argument("--tss-up", type=int, default=3000)
    ap.add_argument("--tss-down", type=int, default=1000)
    ap.add_argument("--proj-pad", type=int, default=500, help="extra length by which the projection search window extends beyond the anchor window")
    ap.add_argument("--proj-k", type=int, default=15)
    ap.add_argument("--proj-cap", type=int, default=8, help="k-mers occurring more than this many times in the window count as repeats and are not used as anchors")
    ap.add_argument("--proj-band", type=int, default=30, help="diagonal tolerance (bp), allows small indels")
    ap.add_argument("--proj-min-anchors", type=int, default=4)
    ap.add_argument("--proj-min-cov", type=float, default=0.3,
                    help="minimum fraction of the source peak that the anchors must cover; below this the copy is considered to lack this element")
    ap.add_argument("--merge-gap", type=int, default=0, help="gap allowed when merging features")
    ap.add_argument("--no-projection", action="store_true",
                    help="skip homeologous projection (use only the called peaks; for ablation controls only)")
    args = ap.parse_args()

    os.makedirs(args.outdir, exist_ok=True)
    lmap = parse_lineage_map(args.lineage_map)

    # ---------------------------------------------------------- 1. Load
    genes = parse_gff_genes(args.gff)
    # v4.3: auto/diploid produces no lineage layer -> ancestry may be absent (anc stays an empty dict, the Ancestry column stays empty)
    anc = load_ancestry(args.ancestry, lmap) if args.ancestry else {}
    # v4.3: without a lineage layer (ploidy=auto/diploid) Ancestry is always written empty --
    #       otherwise it would be written as "UNKNOWN" and downstream would mistake that for
    #       "has a lineage layer" and emit a meaningless lineage table.
    #       Note: Ancestry is only used for the annotation columns (peak_anchor / features / groups);
    #       it does not affect projection or EM.
    lin_default = "UNKNOWN" if args.ancestry else ""
    c2g, g2c = load_clusters(args.clusters)
    fa = Fasta(args.genome)
    peaks = read_peaks(args.peaks)
    log("[02] genes {:,} | allele clusters {:,} | lineage table {:,} | peaks {:,} | ploidy={}".format(
        len(genes), len(c2g), len(anc), len(peaks), args.ploidy))

    # ---- v4.3 guardrail: when the allele clusters are almost all singletons, EM has nothing to do
    #      -> the output is equivalent to ordinary peak counting.
    #      The most common cause: the reference is a [haplotype consensus] while the user runs
    #      diploid/auto (there are no separable copies).
    if c2g:
        _sizes = [len(v) for v in c2g.values()]
        _n1 = sum(1 for s in _sizes if s <= 1)
        _frac1 = _n1 / len(_sizes)
        _med = sorted(_sizes)[len(_sizes) // 2]
        if _frac1 >= 0.9 and _med <= 1:
            log("[02] [WARN] {:.1%} of the allele clusters are singletons (median cluster size={}) -> [no separable copies].".format(_frac1, _med))
            log("[02]        In that case EM takes the 'pinned' branch throughout and the output is equivalent to an ordinary peak count (~featureCounts),")
            log("[02]        not homeologous copy deconvolution. If you expect deconvolution, check that the reference contains several copies of the same gene")
            log("[02]        (a diploid needs a [two-haplotype] reference; a polyploid needs every copy to be assembled).")

    n_cg = sum(len(v) for v in c2g.values())
    n_cg_found = sum(1 for g in g2c if g in genes)
    log("[02] cluster genes found in the GFF {:,}/{:,} ({:.1%})".format(
        n_cg_found, n_cg, n_cg_found / max(n_cg, 1)))
    if n_cg_found < 0.5 * n_cg:
        log("[WARN] more than half of the cluster genes are missing from the GFF; please check that the gene IDs are consistent")
    chrom_missing = {g[0] for g in genes.values() if g[0] not in fa.index}
    if chrom_missing:
        log("[WARN] {} GFF chromosomes are not present in the FASTA (e.g. {})".format(
            len(chrom_missing), ",".join(sorted(chrom_missing)[:3])))
    for c, s, e, n in peaks[:1000]:
        if c not in fa.index:
            die("peak chromosome {} is not present in the genome FASTA -- are the peaks and the genome the same version?".format(c))
    n_oob = sum(1 for c, s, e, n in peaks if c in fa.index and e > fa.length(c))
    if n_oob:
        log("[WARN] {} peaks extend past the chromosome length (they will be truncated) -- check that the peaks and the genome are the same version".format(n_oob))

    # ---------------------------------------------------------- 2. peak anchoring
    win_by_chrom = defaultdict(list)
    for g, gi in genes.items():
        t = tss_of(gi)
        if gi[3] == "+":
            ws, we = t - args.tss_up, t + args.tss_down + 1
        else:
            ws, we = t - args.tss_down, t + args.tss_up + 1
        win_by_chrom[gi[0]].append((max(0, ws), we, t, g))
    wstarts = {}
    maxw = args.tss_up + args.tss_down + 1
    for c in win_by_chrom:
        win_by_chrom[c].sort()
        wstarts[c] = [w[0] for w in win_by_chrom[c]]

    anchor = []   # (name, chrom, s, e, gene, dist, cid, lin, cls)
    for c, s, e, name in peaks:
        if c in fa.index:
            e = min(e, fa.length(c))
        best, bd, bdist = None, None, None
        wl = win_by_chrom.get(c)
        if wl:
            mid = (s + e) / 2.0
            lo = bisect.bisect_left(wstarts[c], s - maxw)
            hi = bisect.bisect_left(wstarts[c], e)
            for i in range(lo, hi):
                ws, we, t, g = wl[i]
                if we > s and ws < e:
                    d = abs(mid - t)
                    if bd is None or d < bd:
                        bd, best = d, g
                        bdist = (mid - t) if genes[g][3] == "+" else (t - mid)
        if best is None:
            cls = "distal"
            cid = lin = ""
        elif best in g2c:
            cls = "cluster_promoter"
            cid, lin = g2c[best], anc.get(best, lin_default)
        else:
            cls = "noncluster_promoter"
            cid, lin = "", anc.get(best, lin_default)
        anchor.append((name, c, s, e, best or "", bdist, cid, lin, cls))

    with open(os.path.join(args.outdir, "peak_anchor.tsv"), "w") as fo:
        fo.write("PeakID\tChrom\tStart\tEnd\tAnchorGene\tDistToTSS\tClusterID\tAncestry\tClass\n")
        for r in anchor:
            fo.write("{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\n".format(
                r[0], r[1], r[2], r[3], r[4], "" if r[5] is None else int(round(r[5])),
                r[6], r[7], r[8]))
    cls_ct = Counter(r[8] for r in anchor)
    log("[02] peak anchoring: " + ", ".join("{}={:,}".format(k, v) for k, v in cls_ct.most_common()))

    # ---------------------------------------------------------- 3. homeologous projection
    tasks = defaultdict(list)
    for r in anchor:
        if r[8] == "cluster_promoter":
            tasks[r[6]].append((r[0], r[1], r[2], r[3], r[4]))
    proj = []
    if not args.no_projection and tasks:
        task_list = []
        for cid, pk in tasks.items():
            gl = [g for g in c2g.get(cid, []) if g in genes and genes[g][0] in fa.index]
            if len(gl) >= 2:
                task_list.append((cid, gl, pk))
        log("[02] homeologous projection: {:,} clusters, {:,} anchored peaks ({} threads)".format(
            len(task_list), sum(len(t[2]) for t in task_list), args.threads))
        if args.threads > 1:
            # use fork when available (Linux servers); platforms without fork (Windows) fall back
            # to spawn (the worker arguments are self-contained)
            _mctx = "fork" if "fork" in mp.get_all_start_methods() else "spawn"
            pool = mp.get_context(_mctx).Pool(args.threads, initializer=init_worker,
                                              initargs=(args.genome, genes, args))
            for i, res in enumerate(pool.imap_unordered(project_cluster, task_list, chunksize=16)):
                proj.extend(res)
                if (i + 1) % 5000 == 0:
                    log("[02]   ... {:,}/{:,} clusters".format(i + 1, len(task_list)))
            pool.close()
            pool.join()
        else:
            init_worker(args.genome, genes, args)
            for t in task_list:
                proj.extend(project_cluster(t))
    n_try = len(proj)
    n_ok = sum(1 for p in proj if p[4] is not None)
    log("[02] projection attempts {:,}, homeologous intervals found {:,} ({:.1%})".format(
        n_try, n_ok, n_ok / max(n_try, 1)))

    with open(os.path.join(args.outdir, "projections.tsv"), "w") as fo:
        fo.write("ClusterID\tSourcePeak\tSourceGene\tTargetGene\tChrom\tStart\tEnd\tN_Anchors\tCoverage\n")
        for p in proj:
            fo.write("{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{:.3f}\n".format(
                p[0], p[1], p[2], p[3], p[4] or "", "" if p[5] is None else p[5],
                "" if p[6] is None else p[6], p[7], p[8]))

    # ---------------------------------------------------------- 4. build features + groups
    # item: (chrom, s, e, kind, gene, cid, peakname, cov)
    items = []
    peak_item = {}
    for r in anchor:
        name, c, s, e, g, _, cid, lin, cls = r
        if e <= s:
            continue
        kind = "called" if cls == "cluster_promoter" else (
            "orphan_noncluster" if cls == "noncluster_promoter" else "orphan_distal")
        peak_item[name] = len(items)
        items.append((c, s, e, kind, g, cid, name, None))
    proj_links = []
    for p in proj:
        if p[4] is None:
            continue
        cid, pname, src, t, c, s, e, na, cov = p
        proj_links.append((len(items), peak_item[pname]))
        items.append((c, s, e, "projected", t, cid, "", cov))

    uf = UF(len(items))
    for a, b in proj_links:
        uf.union(a, b)
    by_chrom = defaultdict(list)
    for i, it in enumerate(items):
        by_chrom[it[0]].append(i)
    for c, lst in by_chrom.items():
        lst.sort(key=lambda i: (items[i][1], items[i][2]))
        cur_end, cur_rep = -1, None
        for i in lst:
            s, e = items[i][1], items[i][2]
            if cur_rep is not None and s < cur_end + args.merge_gap + (1 if args.merge_gap > 0 else 0):
                uf.union(i, cur_rep)
                cur_end = max(cur_end, e)
            else:
                cur_rep, cur_end = i, e

    comp = defaultdict(list)
    for i in range(len(items)):
        comp[uf.find(i)].append(i)

    # within each group: merge overlapping intervals per chromosome -> features
    feats = []
    groups = []
    rank_kind = {"called": 0, "projected": 1, "orphan_noncluster": 2, "orphan_distal": 3}
    for root, members in comp.items():
        is_hru = any(items[i][3] in ("called", "projected") for i in members)
        cids = sorted({items[i][5] for i in members if items[i][5]})
        per_chrom = defaultdict(list)
        for i in members:
            per_chrom[items[i][0]].append(i)
        gfeats = []
        for c, lst in per_chrom.items():
            lst.sort(key=lambda i: items[i][1])
            block = []
            bend = -1
            for i in lst + [None]:
                if i is not None and block and items[i][1] < bend + max(args.merge_gap, 0) + (1 if args.merge_gap > 0 else 0):
                    block.append(i)
                    bend = max(bend, items[i][2])
                    continue
                if block:
                    gfeats.append((c, block))
                if i is not None:
                    block = [i]
                    bend = items[i][2]
        groups.append((is_hru, cids, gfeats))

    for is_hru, cids, gfeats in groups:
        recs = []
        for c, block in gfeats:
            s = min(items[i][1] for i in block)
            e = max(items[i][2] for i in block)
            kinds = sorted({items[i][3] for i in block}, key=lambda k: rank_kind[k])
            kind = kinds[0]
            if kind == "called" or ("called" in kinds):
                typ = "called"
            elif kind == "projected":
                typ = "called" if any(k.startswith("orphan") for k in kinds) else "projected"
            else:
                typ = kind
            gset = sorted({items[i][4] for i in block if items[i][4]})
            # the gene(s) of an HRU feature = the genes in the cluster (projected/called); orphan uses the anchor gene
            if is_hru:
                cg = sorted({items[i][4] for i in block if items[i][3] in ("called", "projected")})
                gset = cg or gset
            gene = ",".join(gset)
            lins = sorted({anc.get(g, lin_default) for g in gset}) if gset else []
            lin = lins[0] if len(lins) == 1 else ("MIXED" if lins else lin_default)
            pids = sorted({items[i][6] for i in block if items[i][6]})
            covs = [items[i][7] for i in block if items[i][7] is not None]
            recs.append([c, s, e, typ, gene, lin, ",".join(pids),
                         "{:.3f}".format(max(covs)) if covs else "NA"])
        feats.append((is_hru, ";".join(cids), recs))

    # sort and number
    fa_len = fa.lengths()
    chrom_rank = {c: i for i, c in enumerate(fa.order)}
    flat = []
    for gi, (is_hru, cidstr, recs) in enumerate(feats):
        for r in recs:
            flat.append((chrom_rank.get(r[0], 10 ** 9), r[1], r[2], gi, r))
    flat.sort()
    grp_id = {}
    n_hru = n_og = 0
    for _, _, _, gi, r in flat:
        if gi not in grp_id:
            if feats[gi][0]:
                n_hru += 1
                grp_id[gi] = "HRU_{:07d}".format(n_hru)
            else:
                n_og += 1
                grp_id[gi] = "OG_{:07d}".format(n_og)

    fpath = os.path.join(args.outdir, "features.tsv")
    bpath = os.path.join(args.outdir, "features.bed")
    gstats = defaultdict(lambda: {"cid": "", "n": 0, "called": 0, "projected": 0,
                                  "genes": [], "lin": Counter(), "orphan": 0})
    type_ct = Counter()
    with open(fpath, "w") as fo, open(bpath, "w") as fb:
        fo.write("\t".join(FEATURE_COLS) + "\n")
        for n, (_, _, _, gi, r) in enumerate(flat, 1):
            fid = "F{:08d}".format(n)
            c, s, e, typ, gene, lin, pids, cov = r
            if c in fa_len:
                e = min(e, fa_len[c])
            gid = grp_id[gi]
            cidstr = feats[gi][1]
            fo.write("\t".join(str(x) for x in [fid, c, s, e, e - s, gid, cidstr, gene, lin,
                                               typ, pids, cov]) + "\n")
            fb.write("{}\t{}\t{}\t{}\t0\t.\n".format(c, s, e, fid))
            type_ct[typ] += 1
            st = gstats[gid]
            st["cid"] = cidstr
            st["n"] += 1
            if typ == "called":
                st["called"] += 1
            elif typ == "projected":
                st["projected"] += 1
            else:
                st["orphan"] += 1
            if gene:
                st["genes"].append(gene)
            st["lin"][lin] += 1

    with open(os.path.join(args.outdir, "groups.tsv"), "w") as fo:
        fo.write("GroupID\tClusterID\tN_Features\tN_Called\tN_Projected\tN_Orphan\tLineages\tGenes\n")
        for gid in sorted(gstats):
            st = gstats[gid]
            fo.write("{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\n".format(
                gid, st["cid"], st["n"], st["called"], st["projected"], st["orphan"],
                ",".join("{}:{}".format(k, v) for k, v in sorted(st["lin"].items())),
                ";".join(st["genes"])))

    hru_sizes = Counter(min(gstats[g]["n"], 15) for g in gstats if g.startswith("HRU_"))
    stats = {
        "n_peaks": len(peaks),
        "peak_class": dict(cls_ct),
        "projection_attempts": n_try,
        "projection_found": n_ok,
        "n_features": len(flat),
        "feature_types": dict(type_ct),
        "n_hru_groups": n_hru,
        "n_orphan_groups": n_og,
        "hru_group_size_hist": {str(k): v for k, v in sorted(hru_sizes.items())},
        "params": vars(args),
    }
    write_json(stats, os.path.join(args.outdir, "features_stats.json"))
    log("[02] {:,} features: {} | HRU groups {:,}, orphan {:,}".format(
        len(flat), ", ".join("{}={:,}".format(k, v) for k, v in type_ct.most_common()), n_hru, n_og))
    log("[02] output: {}".format(fpath))


if __name__ == "__main__":
    main()
