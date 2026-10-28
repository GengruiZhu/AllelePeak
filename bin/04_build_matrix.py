#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
ATAC-Quant v4 | Step 04: multi-sample count matrices (feature / HRU group / group x lineage / gene)

Input: features.tsv + <sample>.quant.tsv for every sample
Output (--outdir):
  feature_annotation.tsv       feature annotation (same as features.tsv, with the Lineage column name unified)
  feature_counts.tsv           feature x sample, EM expected fragment counts rounded to integers (for DESeq2)
  feature_counts_float.tsv     not rounded
  feature_resolved_frac.tsv    Resolved_Frac per feature x sample (unique-best fragments / EM count)
  feature_sd.tsv               bootstrap standard deviation (when a bootstrap was run)
  group_counts.tsv             HRU group (regulatory element summed over all copies) x sample  -- the most robust level
  group_lineage_counts.tsv     HRU group x lineage x sample
  gene_counts.tsv              sum of all features on the promoter of each gene (copy) x sample  -- links to LQ allele expression
  lineage_summary.tsv          per sample, count totals and resolved share by lineage x feature type

Recommendation: for differential analysis prefer group_counts / group_lineage_counts; draw
copy-level (feature_counts) conclusions only on features with a high Resolved_Frac.
"""
import os
import sys
import argparse
from collections import defaultdict, OrderedDict

sys.path.insert(0, os.path.join(os.path.dirname(os.path.abspath(__file__)), "..", "lib"))
from atacq_common import log, die, load_features, FEATURE_COLS  # noqa: E402


def read_quant(path):
    d = {}
    with open(path) as f:
        h = f.readline().rstrip("\n").split("\t")
        for line in f:
            p = line.rstrip("\n").split("\t")
            r = dict(zip(h, p))
            d[r["FeatureID"]] = r
    return d


def fnum(x):
    try:
        v = float(x)
        return v
    except (TypeError, ValueError):
        return float("nan")


def main():
    ap = argparse.ArgumentParser(description="ATAC-Quant v4 Step 04: count matrices")
    ap.add_argument("--features", required=True)
    ap.add_argument("--quant", nargs="+", required=True, metavar="SAMPLE:PATH")
    ap.add_argument("-o", "--outdir", required=True)
    ap.add_argument("--resolved-cut", type=float, default=0.5,
                    help="Resolved_Frac threshold for 'resolved' in lineage_summary")
    ap.add_argument("--ploidy", choices=["diploid", "auto", "allo"], default="allo",
                    help="v4.3: auto/diploid have no lineage layer -> group_lineage_counts / lineage_summary are not written")
    args = ap.parse_args()
    os.makedirs(args.outdir, exist_ok=True)

    feats = load_features(args.features)
    # v4.3: whether there is a lineage layer (only allo has one). Also tolerates the case
    # where ancestry was supplied but is empty everywhere.
    has_lin = any((ft.get("Ancestry") or "").strip() for ft in feats)
    if not has_lin:
        log("[04] ploidy={} and features have no Ancestry -> skipping the lineage-layer output"
            " (group_lineage_counts.tsv / lineage_summary.tsv not written)".format(args.ploidy))
    samples = OrderedDict()
    for spec in args.quant:
        s, p = spec.split(":", 1)
        samples[s] = read_quant(p)
        log("[04] read {} ({:,} features)".format(s, len(samples[s])))
        if len(samples[s]) != len(feats):
            die("{} has {} features but features.tsv has {} -- did you use a different features.tsv?".format(
                s, len(samples[s]), len(feats)))
    names = list(samples.keys())

    def W(path, header, rows):
        with open(os.path.join(args.outdir, path), "w") as fo:
            fo.write("\t".join(header) + "\n")
            for r in rows:
                fo.write("\t".join(str(x) for x in r) + "\n")

    # --- feature level
    W("feature_annotation.tsv", FEATURE_COLS,
      ([ft.get(c, "") for c in FEATURE_COLS] for ft in feats))
    cnt = {s: [fnum(samples[s][ft["FeatureID"]]["EM_Count"]) for ft in feats] for s in names}
    W("feature_counts.tsv", ["FeatureID"] + names,
      ([ft["FeatureID"]] + [int(round(cnt[s][i])) for s in names] for i, ft in enumerate(feats)))
    W("feature_counts_float.tsv", ["FeatureID"] + names,
      ([ft["FeatureID"]] + ["{:.3f}".format(cnt[s][i]) for s in names] for i, ft in enumerate(feats)))
    W("feature_resolved_frac.tsv", ["FeatureID"] + names,
      ([ft["FeatureID"]] + [samples[s][ft["FeatureID"]]["Resolved_Frac"] for s in names]
       for ft in feats))
    if any(samples[s][feats[0]["FeatureID"]]["EM_SD"] != "NA" for s in names):
        W("feature_sd.tsv", ["FeatureID"] + names,
          ([ft["FeatureID"]] + [samples[s][ft["FeatureID"]]["EM_SD"] for s in names]
           for ft in feats))

    # --- group / group x lineage / gene
    g_sum = defaultdict(lambda: defaultdict(float))
    g_meta = {}
    gl_sum = defaultdict(lambda: defaultdict(float))
    gl_n = defaultdict(int)
    gene_sum = defaultdict(lambda: defaultdict(float))
    gene_meta = {}
    for i, ft in enumerate(feats):
        gid = ft["GroupID"]
        m = g_meta.setdefault(gid, {"cid": ft["ClusterID"], "n": 0, "types": defaultdict(int)})
        m["n"] += 1
        m["types"][ft["Type"]] += 1
        key = (gid, ft["Ancestry"])
        gl_n[key] += 1
        for s in names:
            g_sum[gid][s] += cnt[s][i]
            gl_sum[key][s] += cnt[s][i]
        if ft["Gene"] and ft["Type"] in ("called", "projected", "orphan_noncluster"):
            for g in ft["Gene"].split(","):
                gene_meta.setdefault(g, (ft["ClusterID"], ft["Ancestry"]))
                for s in names:
                    gene_sum[g][s] += cnt[s][i]

    W("group_counts.tsv", ["GroupID", "ClusterID", "N_Features", "Types"] + names,
      ([gid, g_meta[gid]["cid"], g_meta[gid]["n"],
        ",".join("{}:{}".format(k, v) for k, v in sorted(g_meta[gid]["types"].items()))] +
       [int(round(g_sum[gid][s])) for s in names] for gid in sorted(g_sum)))
    if has_lin:
        W("group_lineage_counts.tsv", ["GroupID", "Lineage", "N_Features"] + names,
          ([k[0], k[1], gl_n[k]] + [int(round(gl_sum[k][s])) for s in names]
           for k in sorted(gl_sum) if k[0].startswith("HRU_")))
    W("gene_counts.tsv", ["Gene", "ClusterID", "Ancestry"] + names,
      ([g, gene_meta[g][0], gene_meta[g][1]] + [int(round(gene_sum[g][s])) for s in names]
       for g in sorted(gene_sum)))

    # --- lineage summary (only written when allo has a lineage layer)
    if has_lin:
        rows = []
        for s in names:
            agg = defaultdict(lambda: [0.0, 0.0, 0])
            for i, ft in enumerate(feats):
                key = (ft["Type"], ft["Ancestry"])
                c = cnt[s][i]
                rf = fnum(samples[s][ft["FeatureID"]]["Resolved_Frac"])
                agg[key][0] += c
                if rf == rf and rf >= args.resolved_cut:
                    agg[key][1] += c
                agg[key][2] += 1
            tot = sum(v[0] for v in agg.values()) or 1.0
            for (typ, lin), (c, cr, n) in sorted(agg.items()):
                rows.append([s, typ, lin, n, "{:.1f}".format(c), "{:.4f}".format(c / tot),
                             "{:.4f}".format(cr / c if c > 0 else 0.0)])
        W("lineage_summary.tsv", ["Sample", "Type", "Lineage", "N_Features", "EM_Count",
                                  "Share_of_all", "Resolved_Share"], rows)
    log("[04] matrices written to {}".format(args.outdir))


if __name__ == "__main__":
    main()
