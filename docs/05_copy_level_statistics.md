# Statistics-layer prototype v0.1 - when a copy-level call is safe (ATAC-Quant v2)

> Established 2026-10-09. Input = the stage 03/04 output of one project run (called `pipeline_out` here), zero re-running; 33 seconds in total.
> Rerun command (needs python3 + numpy + pandas):
> ```bash
> cd analysis
> python3 10_copy_level_stats.py \
>     --quant ../pipeline_out/quant/<sample>/<sample>.quant.tsv \
>     --features ../pipeline_out/features/features.tsv \
>     --groups ../pipeline_out/features/groups.tsv \
>     --outdir out --mc 2000 --seed 2026
> ```

## 1. Why a statistics layer is needed

The QC report says "copy-level identifiability is low" (unique-best fragments account for only 9.4% of EM counts), but that does not mean "copy level is unusable".
To pin down the usable range, three things are needed: a tier label for every feature, a test built on informative fragments,
and the smallest difference each group can detect at its current fragment count. This directory is the prototype implementation of those three things.

## 2. Output files

| File | Rows | Content |
|---|---|---|
| `copy_level_tiers.tsv` | 492,503 | the full feature table + `Tier` + the group's `q` / `MDD_fold` |
| `hru_group_tests.tsv` | 25,793 | HRU groups: test eligibility, chi-square, p, q, MDD, most deviating copy |
| `copy_level_informative_tests.tsv` | 46,196 | copy rows: NBU / EM / chance share / informative share / ratios against both baselines / group q |

Key columns (copy-row table): `Opp_Share` = chance-baseline share (normalised over N_Candidate within the group);
`Inform_Share` = informative share (normalised over N_BestUnique within the group); `Ratio_vs_Opp` = Inform_Share/Opp_Share;
`Ratio_vs_EM` = Inform_Share / EM share; `Group_q` = q value after group-level BH-FDR correction.

## 3. Tier definitions (for called features; numbers from the reference run)

| Tier | Condition | Features | Meaning and use |
|---|---|---|---|
| **A** | RF >= 0.5 and NBU >= 10 | 12,493 | directly usable: point estimate + test |
| **B** | RF >= 0.2 and NBU >= 5 | 14,568 | usable once the test supports it |
| **C** | NBU >= 5 and RF < 0.2 | 3,272 | informative fragments exist but the model dominates: report the detection limit only / aggregation only |
| **D** | 3 <= NBU < 5 | 3,283 | too little information |
| **E** | NBU < 3 | 38,500 | give up on the copy level |

- RF = Resolved_Frac (the resolved fraction) = N_BestUnique / EM_Count; NBU = N_BestUnique (number of unique-best fragments).
- **A+B = 27,061 features (37.5% of called), carrying 93% of the called informative fragments** - the main battleground for copy-level analysis.
- By lineage (A/B feature counts): L1 7,468/6,778, L2 2,742/4,787, L3 2,283/3,003 (example lineage labels).

## 4. Group-level tie-aware test design

- **Informative fragments only**: the assignment of an NBU is close to certain (that fragment scores highest on one copy only),
  it does not depend on EM allocation, and it is therefore a weak-assumption statistic.
- **Null hypothesis**: all copies in the group are equally open -> the expected distribution of NBUs is proportional to the chance baseline (normalised over N_Candidate within the group).
- **Statistic**: Pearson chi-square = sum of (o-e)^2/e; p from a Monte Carlo multinomial simulation (B=2000); group-level BH-FDR.
- **Eligibility**: >=2 called copies within the group (and N_Candidate >= 5), total NBU >= 12 -> 9,078 groups (35%) qualified in this run.
- Significant groups with q<0.05 total 8,177 (90.1%), **but see the confounding caveat in section 6 - they cannot all be taken at face value**.

## 5. Detection limit (MDD)

- Definition: at that group's total NBU, the **smallest fold difference between copies** that can be detected with 80% power (alpha=0.05, two-sided).
- Formula (binomial normal approximation, copy c against the rest of the group pooled): find the smallest f such that
  `Phi( (|p1-p0| - z_alpha*sqrt(p0(1-p0)/n)) / sqrt(p1(1-p1)/n) ) >= 0.8`, where `p1 = f*p0/(f*p0 + 1-p0)`.
- Use: place `MDD_fold` next to the observed `Ratio_vs_Opp` - "this group can only detect >=X-fold differences, and Y-fold was observed".
- This run: median group-level MDD 2.5x; **strong-signal groups (q<0.05 and ratio>=3 and MDD<=3) = 3,266**.

## 6. Results and the key caveat (important)

**Numbers**: 9,078 groups eligible for testing -> 8,177 significant (90.1%) -> 3,266 strong-signal groups;
3,926 copy rows exceed 2x against both baselines (chance share and EM share).

**Why the 90% cannot all be taken at face value - the confounding is not removed**:

- Case HRU_0003342: copy Chr1F has NBU=18, EM=20 (RF=0.89, informative content maxed out),
  an NBU share of 90% but an EM share of only 6%; in the same group copy Chr1G has EM=190 yet only NBU=2.
  The two statistics give opposite answers about "who is more open".
- Mechanism: NBU_c is proportional to m_c x q_c (m_c = openness, **q_c = the density of "specifically discriminating sites" on that copy, an intrinsic sequence property**);
  the EM share additionally passes through the "pool geometry" (the mutual similarity structure within the group) and through the feature's effective length. When q_c is asymmetric within a group,
  both statistics are biased towards the more diverged copy (higher q).
- Conclusion: **absolute single-sample comparisons of "who is more open" across copies must be made with care** - either correct with a sequence-level q_c (v0.2),
  or use a conditional paired design (q_c cancels between conditions). The 90% significance mainly says that "within-group heterogeneity is ubiquitous",
  but the average direction is affected by q confounding.

## 7. What you can safely say (claim rules)

1. **Lower-bound statements**: for Tier A/B features you may write "copy c carries at least NBU=... distinguishable fragments at this element" -
   this is a bound under weak assumptions, so it is safe.
2. **High-confidence candidate list**: the 3,266 groups with `q<0.05` and `ratio>=3` and `MDD<=3` as priorities for follow-up / experimental validation;
   when writing directional conclusions, state that "q confounding is uncorrected".
3. **Conditional contrasts (the recommended main use)**: when comparing the same feature x copy across conditions, q_c and the pool geometry cancel ->
   the copy-level difference test is valid. The multi-sample matrices (`feature_counts.tsv` / `group_counts.tsv`) already support this.
4. **Avoid**: writing conclusions straight from single-sample absolute comparisons; reporting point estimates for features with RF<0.2.

## 8. Roadmap for v0.2

1. **q_c estimation (the key item)**: using `genome.fa` plus the features intervals, compute for every copy of every group a "discriminating k-mer density"
   (the fraction of sliding-window 15-31-mers that do not occur in the other copies of the group; the idea can reuse the unique k-mer harvesting of v1).
   With q_c in hand the test baseline becomes m_c x q_c, which allows q-corrected inference about open-fraction shares.
2. **In-engine LRT**: run a full-weight likelihood-based heterogeneity test for each HRU group inside stage 03 (replacing this prototype's approximation).
3. **Multi-sample paired design + simulation validation** (verify the calibration of the test on mocks with known q).
4. **Fold tiers and MDD into the stage 05 QC report** so that they become standard pipeline output.

## 9. Limitations

- NBU uses a strict "unique best" definition (tied-best fragments are not counted) -> the statistic is conservative;
- the chance baseline (proportional to N_Candidate) is an approximation, not strictly "the distinguishable opportunity";
- `min_nbu=12`, `min_cand=5` are empirical starting points and can be adjusted (script parameters);
- this prototype is a "screening / description / detection limit" framework; strict calibration of directional p values awaits v0.2.
