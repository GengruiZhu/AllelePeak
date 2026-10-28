#!/usr/bin/env python3
# =============================================================================
# ATAC-Quant v2 statistics layer prototype (Step 10)
#
# Purpose: on top of the existing 03/04 outputs, turn "can a copy-level conclusion
# be drawn" into an operational standard:
#   1) per-feature tier list (Tier A-E) -- decides which features allow a
#      copy-level conclusion
#   2) HRU group-level tie-aware informative test -- on the "resolvable fragments
#      (N_BestUnique)", run a multinomial test to see whether the openness of the
#      copies in the group deviates significantly from the chance baseline
#   3) minimum detectable difference (MDD) -- for each group / copy at the current
#      fragment count, the smallest detectable between-copy fold difference
#      (80% power, alpha=0.05)
#
# Conventions and assumptions (detailed in README_copy_level_stats.md):
#   - informative evidence = N_BestUnique (unique-best fragments), whose assignment
#     is nearly certain and does not rely on EM allocation
#   - chance baseline opportunity: the N_Candidate of each copy in the group
#     (number of near-best fragments with weight>=0.5) normalised; the null
#     hypothesis is that every copy has the same "fragment opportunity to be
#     assigned" and the same degree of openness
#   - test: multinomial Monte Carlo (B draws), statistic = Pearson chi2, BH-FDR correction
#   - this test is an approximation of a conservative prototype; the strict version
#     (full-weight likelihood ratio) will be merged into the 03 engine later
#
# Dependencies: numpy, pandas (only these two)
# Usage:
#   python 10_copy_level_stats.py \
#       --quant PROJECT_out/quant/SAMPLE/SAMPLE.quant.tsv \
#       --features PROJECT_out/features/features.tsv \
#       --groups PROJECT_out/features/groups.tsv \
#       --outdir analysis/out --mc 2000 --seed 2026
# =============================================================================
import argparse
import math
import numpy as np
import pandas as pd

Z_A = 1.959963984540054  # alpha=0.05 two-sided
Z_B = 0.8416212335729143  # power 0.8


def bh_fdr(p):
    """Benjamini-Hochberg FDR; returns the array of q values (same order as input)."""
    p = np.asarray(p, dtype=float)
    n = len(p)
    if n == 0:
        return p
    order = np.argsort(p)
    q = np.empty(n, dtype=float)
    prev = 1.0
    for i in range(n - 1, -1, -1):
        idx = order[i]
        val = min(prev, p[idx] * n / (i + 1))
        q[idx] = val
        prev = val
    return np.clip(q, 0, 1)


def phi(x):
    return 0.5 * (1.0 + math.erf(x / math.sqrt(2.0)))


def mdd_fold(n, p0, target_power=0.8):
    """Detection limit: given n informative fragments and a chance baseline share p0,
    the smallest fold difference f (copy vs the rest pooled) that reaches target_power.
    Binomial normal approximation; returns None = unreachable (>cap)."""
    if n <= 0 or p0 <= 0 or p0 >= 1:
        return None
    for f in (1.2, 1.3, 1.5, 1.75, 2.0, 2.5, 3.0, 4.0, 5.0, 7.0, 10.0, 15.0, 20.0, 50.0):
        p1 = f * p0 / (f * p0 + (1.0 - p0))
        se0 = math.sqrt(p0 * (1.0 - p0) / n)
        se1 = math.sqrt(p1 * (1.0 - p1) / n)
        if se1 <= 0:
            continue
        z = (abs(p1 - p0) - Z_A * se0) / se1
        if phi(z) >= target_power:
            return f
    return None


def main():
    ap = argparse.ArgumentParser(description="ATAC-Quant v2 Step 10: copy-level statistics layer prototype")
    ap.add_argument("--quant", required=True)
    ap.add_argument("--features", required=True)
    ap.add_argument("--groups", required=True)
    ap.add_argument("--outdir", required=True)
    ap.add_argument("--min-nbu", type=int, default=12,
                    help="minimum number of informative fragments required for the group-level test (total NBU)")
    ap.add_argument("--min-cand", type=int, default=5,
                    help="minimum N_Candidate for a copy to enter the chance baseline")
    ap.add_argument("--mc", type=int, default=2000, help="number of Monte Carlo draws")
    ap.add_argument("--seed", type=int, default=2026)
    args = ap.parse_args()

    rng = np.random.default_rng(args.seed)

    # ---------------- Read input ----------------
    q = pd.read_csv(args.quant, sep="\t")
    f = pd.read_csv(args.features, sep="\t",
                    usecols=["FeatureID", "Chrom", "Start", "End",
                             "GroupID", "ClusterID", "Gene", "Ancestry", "Type", "ProjCoverage"])
    d = q.merge(f, on="FeatureID", how="left")
    print("[10] features %d | quant %d" % (len(d), len(q)))

    for c in ["EM_Count", "EM_SD", "N_Unique", "N_BestUnique", "N_Candidate",
              "N_Primary", "Resolved_Frac", "Length"]:
        d[c] = pd.to_numeric(d[c], errors="coerce")
    d["RF"] = d["Resolved_Frac"].fillna(0.0)
    d["EM"] = d["EM_Count"].fillna(0.0)
    d["NBU"] = d["N_BestUnique"].fillna(0).astype(np.int64)
    d["CAND"] = d["N_Candidate"].fillna(0).astype(np.int64)
    d["GroupID"] = d["GroupID"].fillna("")
    d["is_hru"] = d["GroupID"].str.startswith("HRU")

    # ---------------- 1. Per-feature tiering ----------------
    # Pure ladder (RF, NBU): A directly usable / B usable with test support / C detection limit only / D too little information / E abandon
    tier = np.full(len(d), "", dtype=object)
    called = (d["Type"] == "called").to_numpy()
    rf = d["RF"].to_numpy(); em = d["EM"].to_numpy(); nbu = d["NBU"].to_numpy()
    tA = called & (rf >= 0.5) & (nbu >= 10)
    tB = called & (~tA) & (rf >= 0.2) & (nbu >= 5)
    tC = called & (~tA) & (~tB) & (nbu >= 5)
    tD = called & (~tA) & (~tB) & (~tC) & (nbu >= 3)
    tE = called & (~tA) & (~tB) & (~tC) & (~tD)
    tier[tA] = "A"; tier[tB] = "B"; tier[tC] = "C"; tier[tD] = "D"; tier[tE] = "E"
    tier[~called] = d.loc[~called, "Type"].to_numpy()
    d["Tier"] = tier

    # ---------------- 2. HRU group-level test ----------------
    hru = d[d["is_hru"]].copy()
    # sort, then process group by group
    hru = hru.sort_values("GroupID")
    gids = hru["GroupID"].to_numpy()
    bounds = np.flatnonzero(np.r_[True, gids[1:] != gids[:-1], True])

    copy_rows = []      # one row per (group, copy)
    group_rows = []     # one row per group
    n_tested = 0
    B = args.mc

    for bi in range(len(bounds) - 1):
        sl = slice(bounds[bi], bounds[bi + 1])
        sub = hru.iloc[sl]
        m = (sub["Type"] == "called") & (sub["CAND"] >= args.min_cand)
        test = sub[m]
        nbu_all = int(test["NBU"].sum())
        if len(test) < 2 or nbu_all < args.min_nbu:
            group_rows.append({
                "GroupID": sub["GroupID"].iloc[0], "N_Called_Tested": int(len(test)),
                "NBU_total": nbu_all, "Tested": False,
                "p": np.nan, "q": np.nan, "chi2": np.nan, "MDD_fold": np.nan,
                "Max_ratio_copy": "", "Max_ratio": np.nan})
            continue

        obs = test["NBU"].to_numpy().astype(np.float64)
        opp = test["CAND"].to_numpy().astype(np.float64)
        exp_share = opp / opp.sum()
        em_cnt = test["EM"].to_numpy().astype(np.float64)
        em_share = em_cnt / em_cnt.sum() if em_cnt.sum() > 0 else np.full(len(test), 1.0 / len(test))

        chi2_obs = float(((obs - nbu_all * exp_share) ** 2 / (nbu_all * exp_share)).sum())
        sims = rng.multinomial(int(nbu_all), exp_share, size=B)
        chi2_sim = ((sims - nbu_all * exp_share) ** 2 / (nbu_all * exp_share)).sum(axis=1)
        p_val = (1.0 + float((chi2_sim >= chi2_obs).sum())) / (B + 1.0)

        # detection limit: compute f for each copy and take the minimum (the strictest one still detectable)
        mdd_list = [mdd_fold(nbu_all, p0) for p0 in exp_share]
        mdd_grp = min((x for x in mdd_list if x is not None), default=None)

        # copy-level rows
        inf_share = obs / nbu_all
        ratio = np.divide(inf_share, exp_share, out=np.full_like(inf_share, np.nan),
                          where=exp_share > 0)
        ratio_em = np.divide(inf_share, em_share, out=np.full_like(inf_share, np.nan),
                             where=em_share > 0)
        for j, (_, r) in enumerate(test.iterrows()):
            copy_rows.append({
                "GroupID": r["GroupID"], "FeatureID": r["FeatureID"],
                "Gene": r["Gene"], "Ancestry": r["Ancestry"],
                "NBU": int(obs[j]), "EM_Count": round(float(em_cnt[j]), 2),
                "N_Candidate": int(opp[j]),
                "Opp_Share": round(float(exp_share[j]), 4),
                "Inform_Share": round(float(inf_share[j]), 4),
                "Ratio_vs_Opp": round(float(ratio[j]), 3) if np.isfinite(ratio[j]) else np.nan,
                "Ratio_vs_EM": round(float(ratio_em[j]), 3) if np.isfinite(ratio_em[j]) else np.nan,
                "Group_NBU_total": nbu_all,
            })
        jmax = int(np.nanargmax(np.where(np.isfinite(ratio), ratio, -np.inf))) if np.isfinite(ratio).any() else -1
        group_rows.append({
            "GroupID": test["GroupID"].iloc[0], "N_Called_Tested": int(len(test)),
            "NBU_total": nbu_all, "Tested": True,
            "p": p_val, "chi2": round(chi2_obs, 3),
            "MDD_fold": mdd_grp if mdd_grp is not None else np.nan,
            "Max_ratio_copy": test["FeatureID"].iloc[jmax] if jmax >= 0 else "",
            "Max_ratio": round(float(ratio[jmax]), 3) if jmax >= 0 else np.nan})
        n_tested += 1

    gdf = pd.DataFrame(group_rows)
    cdf = pd.DataFrame(copy_rows)
    tested = gdf["Tested"].fillna(False).to_numpy()
    gdf.loc[tested, "q"] = bh_fdr(gdf.loc[tested, "p"].to_numpy())
    # merge the group-level p/q/MDD back into the copy rows
    cdf = cdf.merge(gdf[["GroupID", "p", "q", "MDD_fold"]], on="GroupID", how="left")
    cdf = cdf.rename(columns={"p": "Group_p", "q": "Group_q"})
    cdf["Group_Significant"] = cdf["Group_q"] < 0.05

    # attach group-level information to the feature level
    d2 = d.merge(gdf[["GroupID", "p", "q", "MDD_fold", "NBU_total"]],
                 on="GroupID", how="left")

    # ---------------- 3. Output ----------------
    import os
    os.makedirs(args.outdir, exist_ok=True)
    cols_t = ["FeatureID", "Chrom", "Start", "End", "Length", "GroupID", "ClusterID",
              "Gene", "Ancestry", "Type", "Tier", "EM_Count", "EM_SD", "N_BestUnique",
              "N_Candidate", "Resolved_Frac", "NBU_total", "q", "MDD_fold"]
    d2[cols_t].to_csv(os.path.join(args.outdir, "copy_level_tiers.tsv"),
                      sep="\t", index=False, float_format="%.4g")
    gdf.to_csv(os.path.join(args.outdir, "hru_group_tests.tsv"),
               sep="\t", index=False, float_format="%.4g")
    cdf.to_csv(os.path.join(args.outdir, "copy_level_informative_tests.tsv"),
               sep="\t", index=False, float_format="%.4g")

    # ---------------- 4. Summary ----------------
    print("\n===== Tier summary (called features) =====")
    for t in "ABCDE":
        m = d["Tier"] == t
        print("  Tier %s: %6d features | EM_Count %12s | NBU %10s"
              % (t, int(m.sum()), f"{d.loc[m,'EM'].sum():,.0f}", f"{d.loc[m,'NBU'].sum():,.0f}"))
    # v4.3: skip the by-lineage statistics when there is no lineage layer (ploidy=auto/diploid)
    _has_lin = bool((d["Ancestry"].astype(str).str.strip() != "").any()) if "Ancestry" in d.columns else False
    if _has_lin:
        print("\nTier feature counts by lineage:")
        tab = d[d["Type"] == "called"].pivot_table(index="Ancestry", columns="Tier",
                                                   values="FeatureID", aggfunc="count", fill_value=0)
        print(tab.to_string())
    else:
        print("\n(no lineage layer: ploidy=auto/diploid, skipping the by-lineage Tier statistics)")

    sig = int((gdf["q"] < 0.05).sum())
    strong = int(((gdf["q"] < 0.05) & (gdf["Max_ratio"] >= 3) &
                  (gdf["MDD_fold"] <= 3)).sum())
    print("\n===== Group-level tests (min_nbu=%d, MC=%d) =====" % (args.min_nbu, B))
    print("  Testable groups: %d / %d (total HRU groups)" % (n_tested, len(gdf)))
    print("  Significant groups (q<0.05): %d (%.1f%%)" % (sig, 100.0 * sig / max(n_tested, 1)))
    print("  Strong-signal groups (q<0.05 and ratio>=3 and MDD<=3): %d" % strong)
    print("  Median group-level MDD: %s" %
          (round(float(gdf.loc[tested, 'MDD_fold'].median()), 2)
           if tested.sum() else "NA"))

    if sig:
        top = cdf[cdf["Group_q"] < 0.05].copy()
        top = top[top["Group_Significant"]]
        print("\n  Top 10 by Ratio_vs_Opp among significant groups (group, copy, gene, ratio):")
        tt = top.sort_values("Ratio_vs_Opp", ascending=False).head(10)
        for _, r in tt.iterrows():
            print("   %s %s %s ratio=%.2f NBU=%d EM=%.1f q=%.2g"
                  % (r["GroupID"], r["FeatureID"], (str(r["Gene"])[:24]), r["Ratio_vs_Opp"],
                     r["NBU"], r["EM_Count"], r["Group_q"]))
    print("\n[10] done: %s" % args.outdir)


if __name__ == "__main__":
    main()
