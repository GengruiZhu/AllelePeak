# From peaks to copy fractions: a walk through the pipeline

(ATAC-Quant v3 - with measured numbers from one full production run - 2026-10-11)

> This is an "end to end" explainer: it starts from the fundamental problem, then walks through the six stages in execution order, and expands stage 3 (quantification / EM) in full. The measured numbers come from a full run on a 13-copy allopolyploid (the same BAM + peak set as the samples about to be run); the hand-worked toy examples and the measured numbers corroborate each other.
> Companions: [`00_quickstart.md`](00_quickstart.md) (how to run it); this document (why each step exists, and what it actually does).

---

## 0. The fundamental problem, and three actions

In the reference allopolyploid every gene has 13-14 homeologous copies (three subgenome lineages, labeled L1/L2/L3 here). Can the read count over a peak simply be "counted" onto the copy it sits on? **No.** The same read pair frequently aligns equally well to the positions of several copies at once - measured on our task-A dataset: **94.5% of pairs are multi-mapping**, and raising alignment sensitivity does not change that (under the most sensitive setting, 98.2% of the unique evidence stays unique).

So the whole pipeline is built around one thing: **splitting an account that cannot be told apart in an explainable, verifiable way.** There are only three actions:

1. **Nominate** - MACS2 peaks say "there is open-chromatin signal around here" (we borrow only their **position**, never their value);
2. **Spread** - project each nominated position onto every copy **by sequence**, yielding candidate intervals of the form "element x copy" (492,503 features);
3. **Allocate** - every read is split only among the features it touches, in proportion to the evidence (a global EM); once split, the account is conserved exactly.

### The six stages at a glance

| Stage | Script | What it does | Main output |
|---|---|---|---|
| 0 | `bin/00_map_call.sh` | alignment + peak calling (only for samples missing BAM/peaks) | sample BAM, narrowPeak |
| 1 | `bin/01_consensus_peaks.py` | merge peaks across samples into a consensus set | `peaks/consensus_peaks.bed` |
| 2 | `bin/02_build_features.py` | classification -> anchoring -> sequence projection -> feature definition | `features/features.tsv` |
| 3 | `bin/03_quantify.py` | scan and register -> score and fold -> global EM allocation | `quant/<sample>/<sample>.quant.tsv` |
| 4 | `bin/04_build_matrix.py` | aggregate matrices | 9 tables under `matrix/` |
| 5 | `bin/05_qc_report.py` | QC report | `qc/qc_report.md` |
| Statistics layer | `analysis/10_copy_level_stats.py` | tiers / group-level tests / MDD | `stats/<sample>/` |

In this run BAM and peaks were both supplied -> stage 0 was skipped automatically, and what actually ran was **1 -> 2 -> 3 -> 4 -> 5 -> statistics layer**.

---

## 1. Stages 0-1 (alignment and peak calling, consensus peaks)

**Stage 0 (skipped in this run)**: aligns samples that have no BAM with bowtie2 (`-k 20`, concordant pairs only, fragment window <=1000) -> sort -> mark duplicates -> ATAC filtering / Tn5 offset correction -> MACS2 peak calling (gsize=1e10, q=0.05). When future samples arrive as FASTQ only, this stage takes over automatically.

**Stage 1 (consensus peaks)**: when several samples are compared, all of them must be counted on **the same feature set**, otherwise the matrices cannot be aligned. The method: pool the narrowPeak files of all samples, merge overlapping intervals (gap <= merge_gap=0), and keep sites that are supported by at least min_samples samples.
- Single sample (this run): **passed through unchanged** -> **445,041** consensus peaks.
- Output: `consensus_peaks.bed` (coordinates + number of supporting samples), `consensus_peaks.members.tsv` (which samples contributed to each site).

---

## 2. Stage 2: feature construction = homeologous projection

This step translates "peak coordinates" into candidate intervals of the form "element x copy". **It is the step that distinguishes ATAC-Quant from an ordinary ATAC workflow.**

### 2.1 Classification: stamp every peak first

For every peak: which gene's **TSS window** (3 kb upstream to 1 kb downstream) does it fall in? And is that gene present in the **allele cluster table** (`my_cluster.tsv`)?

| Case | Class | Count | What happens next |
|---|---|---|---|
| inside the promoter window of a cluster gene | `cluster_promoter` | 68,110 | **enters projection** |
| inside the promoter window of a non-cluster gene | `noncluster_promoter` | 110,111 | standalone feature (orphan) |
| not in any gene promoter window | `distal` | 266,820 | standalone feature (orphan) |

A real row (from `peak_anchor.tsv`): `peak_2` lies at Chr10A:86,437-86,605, is anchored to gene `gene_Chr10A000006` (35 bp upstream of the TSS), cluster `Global_Cluster_036323`, lineage L2, class `cluster_promoter`.

### 2.2 What the allele cluster table looks like

Line 3 of `my_cluster.tsv`: `Global_Cluster_000002 -> gene_Chr1A000004, ..., gene_Chr1N000004` (13 copies). Looking up the cluster table = getting this "sibling list", i.e. knowing which positions to search for the corresponding interval.

### 2.3 Projection: 15-mer diagonal anchoring

Goal: find "the same stretch of sequence" on the sibling copies (an absent peak does not mean an absent sequence). Mechanism:

- take the source peak sequence (cut from the genome at its coordinates, with padding);
- take the promoter window of the target copy (500 bp padding);
- collect shared 15-mers: in a genuinely homeologous interval the shared k-mers fall on **one and the same diagonal** (consistent relative position, small indels up to 30 bp tolerated); k-mers occurring more than 8 times inside the window are discarded outright (low complexity / repeats);
- decision rule: at least **4** colinear anchors and at least **30%** of the source sequence covered -> the homeologous interval counts as found, and its coordinates are recorded (with N_Anchors and Coverage).

Why these thresholds: the random collision rate of a 15-mer is about 1/4^15, roughly one in a billion, while the genome (both strands) is about 22 billion bases - so **a single 15-mer hits about 20 places genome-wide by chance alone**, and one hit means nothing. 4 or more anchors, on the same diagonal, covering >=30%, is what essentially rules out coincidence.

A real projection record (`projections.tsv`): source peak 8577 (from `gene_Chr10D002008`) -> target `gene_Chr10F001903`, interval Chr10F:28,193,932-28,194,052, **N_Anchors=24, Coverage=0.667**.

Genome-wide: 474,370 projection attempts, 353,971 successes (**74.6%**). The 25.4% that fail are not "failures": the target copy genuinely lacks that sequence (deletion / rearrangement) - nothing is found, so no feature is built (reads would not align there either).

### 2.4 Definition: three feature types + HRU groups

Projected intervals and the original peaks are merged by overlap; one connected set = one **HRU group** (the corresponding intervals of the same regulatory element across copies). Each interval on a copy = one quantified feature, with a Type out of three:

- **called**: that copy has its own MACS2 peak (72,116);
- **projected**: that copy has no peak, but projection found a homeologous interval (49,469);
- **orphan**: outside the projection system (265,274 distal + 105,644 non-cluster promoter), kept as single-copy features (so that the read account stays closed).

**492,503 features and 25,793 HRU groups in total** (338 groups of size 13, 203 of size 14).

### 2.5 The key insight (worth highlighting)

**Projection answers "is this sequence there, and where" - not "is it open".** Example HRU_0000001 (= Global_Cluster_000002) has 13 rows: 7 called (each with its own peak) + 6 projected (Chr1B/1D/1F/1H/1I/1K, no peak, but projection coverage 100% or 99.4%). Whether they are open is decided in stage 3, by the EM, from reads: the EM_Count of a projected copy can be 0 (genuinely closed), or substantial (an open copy that v1 would have missed permanently).

---

## 3. Stage 3: quantification - splitting every read's account

There are only two inputs: **a BAM (all alignments) and features.tsv (492k intervals)**. The output is one number per feature (EM_Count). Eight steps.

### 3.1 Scan and register: every fragment registers "all the features it touched"

The genome is scanned **exactly once**. Scanning means running the feature-interval table over every record of the BAM and asking: **which intervals of the table does this record's span overlap?** (>=1 bp of overlap counts as touching). Whatever it touches is registered in that fragment's candidate list.

- one fragment = all records grouped by read name (mate1/mate2 deduplicated, the same landing site counted once);
- a fragment can touch several copies of the same cluster at once (the most common case), intervals of another cluster (repeat regions), or nothing at all (= a pure background fragment).

Example (registration of four fragments F1-F4 within one cluster):

| Fragment | Landing sites in the BAM | Registration |
|---|---|---|
| F1 | a single site, this cluster - Chr1A | {this cluster - Chr1A} |
| F2 | two sites: this cluster - Chr1C + this cluster - Chr1F | {this cluster - Chr1C, this cluster - Chr1F} |
| F3 | two sites: this cluster - Chr1A + another cluster - X | {this cluster, another cluster} |
| F4 | two sites: this cluster - Chr1A + outside any feature | {this cluster, background} |

**Wording matters: "touched" is not "belongs to".** The scan registers only what was touched; which feature it belongs to, and how much, is the proportion that the EM works out later. During counting the "cluster" is invisible - the engine sees feature IDs only, and does not know who is whose sibling.

Measured: 59,581,567 fragments inside features; 39,258,379 pure-background fragments (all landing sites outside features).

### 3.2 Scoring, folding, background

For every candidate record of every fragment:

**Scoring**: `w = exp( (score of this record - global best score of this fragment) / 2 )` (T=2 is the softening temperature of the discount; "global best" = max(best inside features, best outside features); weights are therefore calibrated so that "the fragment's best alignment = 1.0"). Score gap Delta against w: 0 -> 1.00; 2 -> 0.37; 4 -> 0.14; **9.2 -> 0.01**; 10 -> 0.0067 (at T=2 one mismatch is worth about 6 points -> w is about 0.05).

**Background candidate**: 99% of the genome lies outside features (total feature length 114.7 Mbp = 1.04% of 11.07 Gb). A read outside features very likely has an equally good landing site somewhere (repeats). Treatment: all out-of-feature records of the fragment are converted with the same formula and **summed** into a virtual feature number 492,504 - the background candidate (17 tied landing sites -> background weight 17).

**On "normalisation" (erratum)**: the default path has **no** "divide by the maximum" step - calibrating the weights to "this fragment's best = 1" is enough; the background candidate is simply the **sum** of the converted out-of-feature records (it can exceed 1). Dividing feature candidates and background by their common maximum exists only in the optional lineage k-mer branch (off by default).

**Folding**: candidates with weight < **0.01** are dropped. Because weights are calibrated to "best = 1", the 0.01 cut-off **is equivalent to an absolute score gap of 9.2** (T*ln100 ~ 9.21) - it looks only at how far a candidate is from this fragment's best, **independently of background strength**; background competitiveness is cashed out entirely inside the EM (`w x alpha/eff` competition, see 3.5). Example: a fragment with a feature candidate 4.6 points behind (w=0.1, kept) plus 17 full-score out-of-feature sites (bgw=17) -> it has a background candidate and cannot be pinned -> it enters the EM, where background competes with a 17x density advantage.

**Three-way branch** (counting candidates after folding):

| What remains | Destination | Measured share (of the 59,581,567 in-feature fragments) |
|---|---|---|
| 1 feature, background eliminated | **pinned**: counted directly, not fed to the EM | 4,809,061 (8.1%) |
| 0 features, background only | the whole fragment goes to background | 1,805,635 (3.0%) |
| >=2 candidates (>=2 features, or 1 feature + background) | **enters EM allocation** | 52,966,871 (88.9%) |

(In addition, 39,258,379 pure-background fragments go straight to background and never enter the EM.)

### 3.3 Packing: equal-length equivalence classes

Fragments that enter the EM are packed first: **fragments with an identical candidate list and identical weights are merged into "one row x count"** (their split is identical in every round, so computing it once and multiplying by the count is exactly equivalent). For example, "500 fragments all with the signature {this cluster - Chr1C:1.0, this cluster - Chr1F:0.37}" -> one row x500. Measured: 52,966,871 fragments -> **18,099,762 packs** (about 2.9 fragments per pack on average).

### 3.4 alpha and eff: the two quantities in the allocation formula

- **eff (a fixed constant, computed once)**: a feature's **effective length** = feature length + mean fragment length - 1 (211 bp -> 304); background eff = genome - total feature length, about 10.95 Gb.
- **alpha (the unknown the EM estimates)**: one number per feature, meaning "the current total read count of that feature". **Why it has to exist**: an ambiguous fragment must be divided among its candidates according to how open each of them is - and that openness is exactly what we are trying to obtain (a circularity); the EM guesses first and then iterates, and alpha is "the current guess". **alpha/eff = reads per base** (converting "who has more total" into "who is denser", so that a long interval does not win just by occupying more ground).

**Initialisation (three fixed lines)**: (1) assume equal density per bp and split the total N in proportion to eff; (2) mix half and half with the pinned votes; (3) normalise. Toy example (X: eff 300, pinned 10; Y: eff 600, pinned 5; 30 ambiguous fragments at even weights; N=45): 15/30 -> 12.5/17.5 -> **18.75 / 26.25**. When there is enough information the starting value does not change the end point; it is only a starting point.

### 3.5 EM iterations: two rules, rolled until "plugging it back changes nothing"

- **Allocation**: the share of an ambiguous fragment given to a candidate = (w x alpha/eff) / sum over candidates of (w x alpha/eff);
- **Update**: new alpha = pinned votes + sum of (share x count) (replaced wholesale every round, not accumulated; plus a tiny 1e-8 prior).

Toy example (continued): the 30 ambiguous fragments split between X and Y by density -

| Round | alpha_X | alpha_Y | Fraction of the ambiguous pool going to X |
|---|---|---|---|
| start | 18.75 | 26.25 | 58.8% |
| 1 | 27.7 | 17.4 | 76.1% |
| 2 | 32.8 | 12.2 | 84.4% |
| 3 | 35.3 | 9.7 | 87.9% |
| **stop** | **37.1** | **7.9** | **90.4%** |

**Convergence = plugging the current alpha into the allocation formula reproduces the same alpha** (check: 37.1/300=0.1237, 7.9/600=0.0132, share = 90.4% -> new alpha_X = 10 + 30x0.904 = 37.1, consistent). The number it stops at is **not** "the ratio of pinned votes" (that would be 30:15): **the end point is set self-consistently by three inputs together - pinned votes, weight w, and density alpha/eff**; here Y is twice as long, and the density term moves the end point away from the pinned ratio to 37.1/7.9. A perfectly symmetric, uninformative account (twin copies) stops exactly at half (inside HRU_0000001 this is Chr1I/Chr1J at 30.153 each).

### 3.6 One global ledger: nothing counted twice, nothing lost

All fragments and all features are allocated inside **one and the same EM**; each fragment is split only among its own candidates, its **total weight is always 1**, and it is split exactly once. Where F1-F4 finally go:

- F1 (single candidate) is pinned and counted straight into Chr1A;
- F2 is split between Chr1C and Chr1F (two copies of the same cluster);
- F3 is split between "this cluster" and "another cluster" - viewed cluster by cluster the account **does not duplicate**: the same sum exists only once, the "this cluster" view sees most of it and the "another cluster" view sees the remainder; reviewing it from both viewpoints does not turn one sum into two;
- F4 is split between "this cluster" and background.

**Mass balance is the sentinel**: total over features + background = all fragments. Measured: 57,308,082.7 + 41,531,863.3 = **98,839,946.00** (closing exactly on the 98,839,946 fragments, leaving only a floating-point tail). Any double counting would blow the identity up - v1's counts inflated 10.8x precisely because it had no such global ledger.

### 3.7 Output after convergence

Each feature's final alpha (= pinned votes + accumulated allocations) is written to `quant.tsv`: `EM_Count` (the primary result), `N_Primary` (naive primary-alignment count, a comparison column), `N_BestUnique` (number of unique-best pieces of evidence), `Resolved_Frac` (the resolved fraction: share of hard evidence), `EM_SD` (bootstrap jitter: resample [pinned; equal-length equivalence classes] 20 times and rerun from the convergence point; the spread is reported). Measured: features 57,308,082.7 / background 41,531,863.3 (background is 42% of the whole); **run time 65,735 s, about 18.3 h**; the EM ran the full 2000-round cap (the tail of a large ledger converges slowly, a known item; mass balance and the main conclusions are unaffected).

The 13 final numbers of HRU_0000001 (in the order 1A -> 1N):

`8.956 | 0.000 | 0.000 | 4.701 | 38.799 | 28.067 | 0.000 | 30.153 | 30.153 | 0.000 | 0.000 | 19.329 | 19.329`

How to read it: Chr1C has a peak yet ends at 0.000 (its pile-up was borrowed from elsewhere); Chr1F has no peak yet ends at 28.067 (an open copy rescued by the allocation); Chr1I/Chr1J are exactly half and half (indistinguishable).

### 3.8 Boundaries (the honesty written into the result)

Perfectly indistinguishable copies can only be split evenly - that is an information-theoretic limit, and no algorithm can invent a difference. When pinned-vote ratios are unequal, the split carries a directional bias. A calibration experiment measured the usable range: **across lineages, >=2x is called in the right direction 94% of the time; within a lineage 3-5x is needed; <=1.5x is always flagged "unreliable"** (the flagging is carried by Tier/MDD).

---

## 4. Stages 4-5: matrices and QC

**Stage 4 (matrices)**: 9 tables under `matrix/`: `feature_counts.tsv` / `feature_counts_float.tsv` (count / EM convention), `feature_sd.tsv` (bootstrap fluctuation), `feature_resolved_frac.tsv`, `feature_annotation.tsv`, `gene_counts.tsv`, `group_counts.tsv` (aggregated by HRU group), `group_lineage_counts.tsv`, `lineage_summary.tsv`. One sample is one column; multiple samples can be compared side by side directly.

**Stage 5 (QC)**: `qc/qc_report.md` + `qc/qc_summary.tsv` - peak counts and classes, projection rate, feature composition, HRU group size distribution, per-sample fragment accounting, EM vs naive comparison, and notes on how to read the tables. **Read this first when a run finishes.**

Acceptance check (a homeologous input should reproduce these): 492,503 features; 72,116 called; projection 74.6%; multi-candidate 88.9%; unique-best share 9.4%; mass balance closing exactly.

## 5. Statistics layer

Turning quantification into quotable conclusions: evidence strength is bucketed into Tiers A-E; HRU groups get tie-aware tests; MDD (minimum detectable difference). Measured: of the 72,116 called features, Tier A 12,493 / B 14,568 / C 3,272 / D 3,283 / E 38,500 (A+B = 27,061 carry 93% of the called hard evidence); 9,078 groups qualified for testing -> 8,177 significant (90.1%), 3,266 strong-signal, median MDD 2.5x. The standard reading example HRU_0000006: 750 candidate fragments -> 13 unique pieces of evidence -> EM=26.18, chi-square=50.6, q=6.2e-4, MDD=5x.

## 6. Methodology and credibility

EM allocation of multi-mapping reads is a standard tool in the quantification field (it is the core mechanism of RSEM/kallisto/Salmon; both the "equivalence class" and the "decoy background" come from that lineage). Three layers of validation:

1. **truth simulation** (13 isostructural copies, known truth): r=0.9957, slope=1.005, exact mass balance, projection 764/764;
2. **calibration experiment** (injected indistinguishable degeneracy): yields the detection limits of section 3.8;
3. **real-data diagnostics** (tasks A/B): the alignment-layer error is bounded - 1.8% false-unique, 1.8% misattributed, 0.02% missed true source.

## 7. Cheat sheet (memorise these five)

1. Type (called/projected/orphan) answers only "how did this homeologous copy get into the table", **never "is it open"**;
2. "how open" is `EM_Count`; "how hard is the evidence" is `N_BestUnique` / `Resolved_Frac` / `EM_SD`;
3. by default rely on **HRU-group / lineage aggregates**; use copy-level numbers only from Tier A/B, or from the subset with MDD >=3-5x;
4. when a run finishes, read `qc_report.md` first and compare against the acceptance numbers of section 4;
5. if in doubt, go back to: three actions (nominate / spread / allocate) -> two formulas (`w=exp(Delta/2)`; share = ratio of `w x alpha/eff` across candidates) -> three classes -> two boundaries (symmetry: an even split is the ceiling; tiny differences are not decidable).

## Appendix: glossary

| Term | In one line |
|---|---|
| consensus peak | the unified site set obtained by merging peaks across samples |
| allele cluster table | the full copy list of each homeologous gene cluster (my_cluster.tsv) |
| projection | spreading an element from a source copy onto its sibling copies by 15-mer diagonal anchoring |
| called / projected / orphan | has its own peak / no peak but projected onto / cannot enter the projection system |
| HRU group | the set of corresponding intervals of one regulatory element across copies |
| pinned (vote) | a fragment left with a single candidate after folding, counted directly |
| equal-length equivalence class / packing | fragments with identical candidates and weights merged into "one row x count" |
| alpha | the current total read count of each feature (the unknown the EM estimates) |
| eff | effective length: feature length + mean fragment length - 1; background = genome - features |
| weight w | `exp((score - best)/2)`; calibrated to "best = 1"; folded below 0.01 (~ a 9.2 point gap) |
| background candidate | the virtual candidate number F+1 formed by summing the converted out-of-feature sites |
| mass balance | total over features + background = all fragments (closing exactly) |
| Resolved_Frac / N_BestUnique / EM_SD | hard-evidence share / unique-best evidence count / bootstrap jitter |
| Tier / MDD | evidence tier / minimum detectable difference |
| T=2 / 0.01 folding line | score discount temperature; folding is equivalent to an absolute gap of 9.2 points |
