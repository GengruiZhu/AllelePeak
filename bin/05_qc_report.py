#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
ATAC-Quant v4 | Step 05: QC report

Summarises the key diagnostics per sample, to decide whether the results are usable
and down to which level:
  1. Fragment accounting: total / duplicate / on features / pure background, plus the
     conservation check (EM total = input fragments)
     -- directly prevents the orphan-peak count inflation of v1 (in v1 orphan peak
     counts reached 6.7x the number of sequenced fragments)
  2. Multi-mapping: among fragments on features, the fraction with more than one
     candidate (the "applicability boundary" x-axis recommended by the LQ report)
  3. Identifiability: the fraction of EM counts coming from "unique-best fragments",
     split by lineage x feature type
  4. EM vs naive counting (primary-alignment counting, i.e. the featureCounts
     convention; the LQ report points out that this is the correct comparison)
  5. Feature construction: peak classes, homeologous projection success rate,
     HRU group sizes

Output: qc_summary.tsv, qc_report.md
"""
import os
import sys
import math
import argparse
from collections import defaultdict

sys.path.insert(0, os.path.join(os.path.dirname(os.path.abspath(__file__)), "..", "lib"))
from atacq_common import log, load_features, read_json  # noqa: E402


def fnum(x):
    try:
        return float(x)
    except (TypeError, ValueError):
        return float("nan")


def pearson(a, b):
    n = len(a)
    if n < 3:
        return float("nan")
    ma, mb = sum(a) / n, sum(b) / n
    sa = math.sqrt(sum((x - ma) ** 2 for x in a))
    sb = math.sqrt(sum((y - mb) ** 2 for y in b))
    if sa == 0 or sb == 0:
        return float("nan")
    return sum((x - ma) * (y - mb) for x, y in zip(a, b)) / (sa * sb)


def main():
    ap = argparse.ArgumentParser(description="ATAC-Quant v4 Step 05: QC report")
    ap.add_argument("--features", required=True)
    ap.add_argument("--features-stats", default=None)
    ap.add_argument("--quant", nargs="+", required=True, metavar="SAMPLE:QUANT_TSV")
    ap.add_argument("-o", "--outdir", required=True)
    ap.add_argument("--resolved-cut", type=float, default=0.5)
    ap.add_argument("--ploidy", choices=["diploid", "auto", "allo"], default="allo",
                    help="v4.3: ploidy mode (decides whether the lineage layer is written/shown, and the guidance shown)")
    args = ap.parse_args()
    os.makedirs(args.outdir, exist_ok=True)

    feats = load_features(args.features)
    fmeta = {ft["FeatureID"]: ft for ft in feats}
    # v4.3: whether there is a lineage layer (only allo has one)
    has_lin = any((ft.get("Ancestry") or "").strip() for ft in feats)
    _PMODE = {"allo": "allopolyploid (has lineage layer)",
              "auto": "autopolyploid (no lineage layer)",
              "diploid": "diploid, two-haplotype reference (no lineage layer)"}
    md = ["# ATAC-Quant v4 QC report", "",
          "**Ploidy mode**: `{}` -- {}".format(args.ploidy, _PMODE.get(args.ploidy, args.ploidy)), ""]

    if args.features_stats and os.path.exists(args.features_stats):
        fs = read_json(args.features_stats)
        md += ["## 1. Feature construction", ""]
        md.append("- total peaks: {:,}".format(fs.get("n_peaks", 0)))
        for k, v in sorted(fs.get("peak_class", {}).items()):
            md.append("  - {}: {:,}".format(k, v))
        pa, pf = fs.get("projection_attempts", 0), fs.get("projection_found", 0)
        md.append("- homeologous projection: {:,} attempts, {:,} homeologous intervals found ({:.1%})".format(pa, pf, pf / max(pa, 1)))
        md.append("- features: {:,} ({})".format(fs.get("n_features", 0), ", ".join(
            "{}={:,}".format(k, v) for k, v in sorted(fs.get("feature_types", {}).items()))))
        md.append("- HRU groups {:,}; group size distribution (features: groups) {}".format(
            fs.get("n_hru_groups", 0), fs.get("hru_group_size_hist", {})))
        md.append("")

    summary_rows = []
    md += ["## 2. Per-sample", ""]
    for spec in args.quant:
        s, qpath = spec.split(":", 1)
        spath = qpath.replace(".quant.tsv", ".quant_stats.json")
        st = read_json(spath) if os.path.exists(spath) else {}
        rows = []
        with open(qpath) as f:
            h = f.readline().rstrip("\n").split("\t")
            for line in f:
                rows.append(dict(zip(h, line.rstrip("\n").split("\t"))))
        em = [fnum(r["EM_Count"]) for r in rows]
        naive = [fnum(r["N_Primary"]) for r in rows]
        bestu = [fnum(r["N_BestUnique"]) for r in rows]
        tot_em = sum(em)
        tot_bu = sum(min(b, e) for b, e in zip(bestu, em))

        red = st.get("reduce", {})
        n_in = red.get("n_frag", 0)
        n_multi = red.get("n_multi", 0)
        cons = st.get("conservation_check", {})
        mpc = st.get("member_primary_check", {})

        md.append("### {}".format(s))
        md.append("")
        md.append("| Item | Value |")
        md.append("|---|---|")
        tot = st.get("n_fragments_total_nondup")
        md.append("| Total non-duplicate fragments | {} |".format("{:,}".format(tot) if tot else "NA (outside=ignore)"))
        pa = st.get("passA", {})
        pb = st.get("passB", {})
        md.append("| Fragments whose primary alignment is a duplicate (in features / background) | {:,} / {:,} |".format(
            pa.get("n_dup_in", 0), pb.get("n_bg_dup", 0)))
        md.append("| Fragments on features | {:,} |".format(n_in))
        md.append("| of which multi-candidate (need EM allocation) | {:,} ({:.1%}) |".format(n_multi, n_multi / max(n_in, 1)))
        md.append("| Pure background fragments | {:,} |".format(pb.get("n_bg", 0)))
        md.append("| Number of equivalence classes | {:,} |".format(st.get("n_equivalence_classes", 0)))
        md.append("| EM iterations | {} |".format(st.get("em_iterations", "NA")))
        md.append("| EM total on features / background | {:,.0f} / {:,.0f} |".format(
            st.get("em_sum_features", tot_em), st.get("em_background", 0)))
        md.append("| Conservation check (EM total = input fragments) | {} |".format(
            "pass" if cons.get("ok") else "**FAILED**"))
        if mpc.get("primary_read1_records") is not None:
            ok = mpc["primary_read1_records"] == mpc["fragments_in_features_after_filter"]
            md.append("| Fragments in features = their primary-alignment records | {} ({:,} vs {:,}) |".format(
                "pass" if ok else "**mismatch**", mpc["fragments_in_features_after_filter"],
                mpc["primary_read1_records"]))
        md.append("| Share of EM counts from unique-best fragments | {:.1%} |".format(tot_bu / max(tot_em, 1)))
        md.append("| Mean fragment length | {:.1f} |".format(st.get("mean_frag_len", 0)))
        md.append("")

        # by lineage x type (v4.3: by type only when there is no lineage layer)
        agg = defaultdict(lambda: [0.0, 0.0, 0.0, 0])
        for r in rows:
            ft = fmeta[r["FeatureID"]]
            key = (ft["Type"], ft["Ancestry"]) if has_lin else (ft["Type"], "")
            e = fnum(r["EM_Count"])
            agg[key][0] += e
            agg[key][1] += min(fnum(r["N_BestUnique"]), e)
            agg[key][2] += fnum(r["N_Primary"])
            agg[key][3] += 1
        if has_lin:
            md.append("| Type | Lineage | N_Features | EM_Count | Share | BestUnique_Share | N_Primary |")
            md.append("|---|---|---|---|---|---|---|")
        else:
            md.append("| Type | N_Features | EM_Count | Share | BestUnique_Share | N_Primary |")
            md.append("|---|---|---|---|---|---|")
            md.append("")
            md.append("> This mode (`{}`) **has no lineage layer**, so the table below is not split by lineage.".format(args.ploidy))
        for (typ, lin), (e, bu, nv, n) in sorted(agg.items(), key=lambda x: (x[0][0], -x[1][0])):
            if has_lin:
                md.append("| {} | {} | {:,} | {:,.0f} | {:.1%} | {:.1%} | {:,.0f} |".format(
                    typ, lin, n, e, e / max(tot_em, 1), bu / e if e > 0 else 0, nv))
            else:
                md.append("| {} | {:,} | {:,.0f} | {:.1%} | {:.1%} | {:,.0f} |".format(
                    typ, n, e, e / max(tot_em, 1), bu / e if e > 0 else 0, nv))
            summary_rows.append([s, typ, lin, n, "{:.1f}".format(e), "{:.4f}".format(e / max(tot_em, 1)),
                                 "{:.4f}".format(bu / e if e > 0 else 0), "{:.0f}".format(nv)])
        md.append("")

        # EM vs naive
        pairs = [(math.log1p(a), math.log1p(b)) for a, b in zip(em, naive) if a + b > 0]
        r_log = pearson([p[0] for p in pairs], [p[1] for p in pairs]) if pairs else float("nan")
        big = sum(1 for a, b in zip(em, naive) if (a + 1) / (b + 1) > 2 or (b + 1) / (a + 1) > 2)
        md.append("- EM vs naive primary counts: log1p Pearson r = {:.3f}; features differing by >2x {:,} ({:.1%})".format(
            r_log, big, big / max(len(em), 1)))
        md.append("  (naive count = how many fragments have their primary alignment on this feature; for multi-mapping fragments bowtie2 picks the primary at random)")
        md.append("")

        warn = []
        if not cons.get("ok", True):
            warn.append("Conservation check failed: the EM total does not match the number of input fragments; please send quant_stats.json to the developer")
        if n_in and n_multi / n_in > 0.5:
            _lv = "HRU group / lineage totals" if has_lin else "HRU group / gene totals (this mode has no lineage layer)"
            warn.append("More than half of the in-feature fragments are multi-candidate: copy-level results rely mainly on EM, prefer {}".format(_lv))
        if tot_em and tot_bu / tot_em < 0.3:
            warn.append("Unique-best fragments account for only {:.0%} of EM counts: copy-level identifiability is low".format(tot_bu / tot_em))
        if not has_lin:
            _sizes = defaultdict(int)
            for ft in feats:
                _sizes[ft["GroupID"]] += 1
            _n1 = sum(1 for v in _sizes.values() if v <= 1)
            if _sizes and _n1 / len(_sizes) >= 0.9:
                warn.append("**Allele clusters are almost all singletons ({:.0%}) -> no separable copies**: every EM step takes the 'pinned' branch, "
                            "so this output is equivalent to ordinary peak counting (approx. featureCounts), **not homeologous copy deconvolution**. "
                            "If you expected deconvolution, check that the reference contains multiple copies of the same gene "
                            "(a diploid run needs a [two-haplotype] reference).".format(_n1 / len(_sizes)))
        if warn:
            md.append("**Notes**")
            md.append("")
            for w in warn:
                md.append("- " + w)
            md.append("")

    with open(os.path.join(args.outdir, "qc_summary.tsv"), "w") as fo:
        fo.write("Sample\tType\tLineage\tN_Features\tEM_Count\tShare\tBestUnique_Share\tN_Primary\n")
        for r in summary_rows:
            fo.write("\t".join(str(x) for x in r) + "\n")
    _fallback = "HRU group or lineage totals" if has_lin else "HRU group / gene totals (this mode has no lineage layer)"
    _MODE_TIP = {
        "allo": ["- **Recommended analysis priority** (allopolyploid): **lineage level -> element level -> copy level**.",
                 "  There are many cross-lineage discriminating sites, so the lineage level is the most stable (measured share error <1pp; a difference of >=2x is enough to assign direction); "
                 "draw copy-level conclusions only on features with a high Resolved_Frac."],
        "auto": ["- **Recommended analysis priority** (autopolyploid): **element level / gene-level totals -> copy level**.",
                 "  In an autopolyploid every pair of copies diverges by only 0.1-1% (the hardest 'within-lineage' case of an allopolyploid), "
                 "so copy-level resolution needs a **>=3-5x** difference to be reliable, and the usable range is far smaller than for an allopolyploid.",
                 "  The concept of 'subgenome dominance' **does not apply to autopolyploids** (there are no subgenomes to distinguish).",
                 "  For a genuine copy-level signal, consider **between-individual haplotype differences / population dosage designs** rather than within-individual copies."],
        "diploid": ["- **Recommended analysis priority** (diploid, two-haplotype reference): **element level / gene-level totals -> allele level**.",
                    "  Only **1** competitor has to be excluded (a 13-copy case has to exclude 12), so allele resolution is in principle easier than in a polyploid.",
                    "  **Note**: if the reference is a haplotype consensus (every gene appears only once), this pipeline **has nothing to do** -- "
                    "the output is equivalent to ordinary peak counting, not allele resolution."],
    }
    md += ["## 3. How to read these numbers", "",
           "- **BestUnique share**: how much of the EM count comes from fragments that align best only on this copy. High -> copy-level conclusions are trustworthy; "
           "low -> this copy is almost identical in sequence to its homeologous copies, the copy-level value comes mostly from proportional EM allocation, and you should use {} instead.".format(_fallback),
           "- **projected features**: no peak was called on this copy, but the corresponding sequence was found at the homeologous position. Its count can be 0 "
           "(genuinely closed), or sizeable (open copies that v1 would have missed).",
           "- **Background**: some fragments align equally well outside the features, so part of the mass is assigned to genomic background instead of being pushed into peaks.", ""]
    md += ["### Interpretation notes for this mode (`{}`)".format(args.ploidy), ""] + _MODE_TIP.get(args.ploidy, []) + [""]
    with open(os.path.join(args.outdir, "qc_report.md"), "w") as fo:
        fo.write("\n".join(md) + "\n")
    log("[05] QC report: {}".format(os.path.join(args.outdir, "qc_report.md")))


if __name__ == "__main__":
    main()
