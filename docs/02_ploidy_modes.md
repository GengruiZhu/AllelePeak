# v4.3 - Ploidy modes (`[inputs] ploidy`)

**In one sentence**: one line in the configuration, `ploidy = diploid | auto | allo`, has the user declare the ploidy; the three modes **change no counts at all** - they only decide **whether the lineage layer is produced, which inputs are required, and how results should be read**.

- Version: v4.3 (2026-10-13)
- Naming and semantics aligned with the sister project **LQ (LineageQuant)** and its `--ploidy {auto, allo}`; `diploid` is an extension of this pipeline

---

## 1. First, dispel a misunderstanding: **there are two algorithmic paths, not three**

ATAC-Quant is about "assigning reads to homeologous copies". The three ploidies are **completely identical in the core EM (stage 03)** - lineage never took part in the allocation: stage `03` reads `Ancestry` only in the **optional lineage k-mer channel** (off by default).

| | **with the lineage layer** | **without the lineage layer** |
|---|---|---|
| who takes it | `allo` | `diploid` / `auto` |
| difference | one extra **aggregation dimension**: split by lineage first, then split into copies within a lineage | split into copies directly |

**So `diploid` and `auto` are algorithmically identical**; they differ in exactly four places:

| # | Difference | Explanation |
|---|---|---|
| 1 | whether `ancestry` is required | required for `allo`; not needed for the other two (if supplied it is ignored, with a notice) |
| 2 | interpretation rules (written into the QC report) | recommended analysis level, expected detection limits, caveats |
| 3 | guardrail messages | an explicit warning when the allele clusters are almost all singletons |
| 4 | biological interpretation | in an autopolyploid the concept of "subgenome dominance" does not apply |

**This matters**: once it is clear that "the difference is not in the algorithm", the implementation stays bounded and low-risk.

---

## 2. Behaviour of the three modes

| Dimension | `diploid` | `auto` | `allo` (default) |
|---|---|---|---|
| applies to | diploid, **two-haplotype reference** (2 copies per cluster) | autopolyploid (subgenomes undifferentiated) | allopolyploid (subgenomes differentiated) |
| `clusters` | required (2 copies per cluster) | required | required |
| `ancestry` | **not needed** | **not needed** | **required** |
| **core EM** | same | same | same |
| lineage-layer output<br>(`group_lineage_counts.tsv` / `lineage_summary.tsv` / per-lineage tiers) | skipped | skipped | produced |
| lineage k-mer channel | no | no | optional |
| expected minimum detectable difference | must exclude **1** competitor (theoretically the best case) | all pairs >=3-5x | >=2x across lineages / >=3-5x within a lineage |
| recommended analysis level | element/gene totals -> allele level | **element/gene totals** (copy level is hard) | **lineage level** -> element level -> copy level |

**A hard truth** (both `auto` and `diploid` have to accept it): **the fewer the copies and the lower the divergence, the harder copy-level resolution becomes.**
In an autopolyploid every pair of copies is only 0.1-1% diverged - **which amounts to applying the hardest allopolyploid case ("within a lineage") to every pair of copies**.
The pipeline's own QC report has always said what to do in that situation: **"use HRU group or gene totals instead"**.

---

## 3. The two `diploid` situations - why only a two-haplotype reference is supported

| | **haplotype reference** (the vast majority of diploid species) | **two-haplotype reference** (phased assembly / pangenome) |
|---|---|---|
| each gene in the reference | **appears once** | appears twice (hap1/hap2) |
| `clusters` | all **singletons** | **two** genes per cluster |
| this pipeline | **nothing to do** (each peak has only one candidate; everything is pinned) | **fully applicable** |
| what you actually get | **equivalent to ordinary peak counting** (~featureCounts) | true allele-resolved accessibility |

**In v4.3, `diploid` targets the two-haplotype reference only.** With a haplotype consensus reference the pipeline **will not silently hand you a deconvolution result that merely looks normal** -
it warns explicitly both in stage 02 and in the QC report (see section 4).

---

## 4. Guardrail: an explicit warning when the allele clusters are almost all singletons

**Trigger**: `>=90%` of the allele clusters are singletons (median cluster size <=1).

**The stage 02 log**:

```
[02] [WARNING] 100.0% of allele clusters are singletons (median cluster size=1) -> [no separable copies].
[02]         In that case the EM takes the 'pinned' branch throughout, the output is equivalent to ordinary
[02]         peak counting (~featureCounts) and is NOT homeologous-copy deconvolution. If you expected
[02]         deconvolution, check that the same gene exists in several copies in your reference
[02]         (a diploid needs a [two-haplotype] reference; a polyploid needs every copy assembled).
```

**The "notes" section of the QC report says the same thing** (that is where the user will actually read it).

**Why this guardrail is mandatory**: without it, a user with a haplotype consensus reference would get a result that **looks perfectly normal but is not deconvolution at all** - the most dangerous failure mode there is.

---

## 5. Consistency with LQ

LQ takes **`--ploidy {auto, allo}`** (`auto` = autopolyploid, copy level only; `allo` = allopolyploid, lineage assignment; default `allo`) and has **no `diploid`** - because LQ is an RNA quantifier and simply quantifies diploids as ordinary data.

| | LQ | ATAC-Quant v4.3 |
|---|---|---|
| values | `auto` / `allo` | `auto` / `allo` / **`diploid` (extension)** |
| `auto` semantics | subgenomes undifferentiated, copy level only | **same** |
| `allo` semantics | assign lineage first, then split copies within a lineage | **same** |
| option name | `--ploidy` | `[inputs] ploidy` |

> Note that under `auto` LQ **still accepts `--ancestry`** (it simply does not run the lineage step). ATAC-Quant is stricter: under `auto`/`diploid` **a supplied ancestry is ignored, with a notice**, and the `Ancestry` column of the feature table is written empty - so that nothing downstream can believe a lineage layer is available.

---

## 6. How to use it

```ini
[inputs]
genome = /path/genome.fasta
gff = /path/annotation.gff3
clusters = /path/my_cluster.tsv
ploidy = allo                 # diploid | auto | allo
ancestry = /path/ancestry.tsv # required for allo only; for auto/diploid the whole line can be deleted
```

```bash
python run_atacquant.py --config config.ini --check   # inspect the plan first (the ploidy mode is printed)
python run_atacquant.py --config config.ini
```

The `--check` output:

```
  lineage table : (unused; ploidy=auto produces no lineage layer)
  ploidy mode   : auto (autopolyploid -> no lineage layer: copies split directly; copy-level resolution is hard)
```

**What differs in the output** (only the two lineage-layer tables):

| File | `allo` | `auto` / `diploid` |
|---|---|---|
| `matrix/group_lineage_counts.tsv` | produced | not produced |
| `matrix/lineage_summary.tsv` | produced | not produced |
| `matrix/group_counts.tsv` / `gene_counts.tsv` / `feature_counts.tsv` / `quant.tsv` | produced | produced, **content completely identical** |
| `Ancestry` column of `features.tsv` | L1/L2/... | **empty** |

### Runnable examples (simulated data, shipped with the package)

| Configuration | Mode | Description |
|---|---|---|
| `_selftest/config.test.ini` | `allo` | full example (with ancestry, produces the lineage layer) |
| `_selftest/config.auto.ini` | `auto` | autopolyploid example (no ancestry, lineage layer skipped) |
| `_selftest/config.diploid.ini` | `diploid` | diploid example (no ancestry, lineage layer skipped) |

```bash
# _selftest/sim must exist first (generated by: bash tests/run_test.sh ./_selftest)
python run_atacquant.py --config _selftest/config.auto.ini --check
python run_atacquant.py --config _selftest/config.auto.ini
```

---

## 7. Validation (full chain on local simulated data)

| Item | Result |
|---|---|
| **`allo` regression** | **19 of 23 files byte-identical** to before the change; the 4 differences are all metadata (`params.ploidy`, the timings `em_seconds`/`runtime_sec`, log text, the new QC report section) |
| **counts unchanged** | `quant.tsv` / `feature_counts.tsv` / `feature_counts_float.tsv` / `group_counts.tsv` / `feature_resolved_frac.tsv` / `feature_sd.tsv` are **byte-identical across all three modes** (verified) |
| **`features.tsv` difference** | **the `Ancestry` column only** (empty vs L1/L2/L3) (verified) |
| **lineage layer skipped** | the matrix directory of `auto`/`diploid` holds 7 files (allo holds 9); neither lineage table is produced (verified) |
| **guardrail** | an all-singleton cluster table (164 singletons) -> stage 02 **and** the QC report both warn (verified) |
| **error handling** | `allo` without ancestry -> explicit hint to switch to auto/diploid (verified); the illegal value `triploid` -> rejected (verified); `auto` with ancestry -> notice that it will be ignored (verified) |
| **QC framing** | the report header shows the ploidy mode; with no lineage the tables are no longer split by lineage and an explanation is added; reading notes are given per mode (verified) |

### A third old bug fixed along the way

In `atacquant_core.py`:

```python
for k in ("genome", "gff", "clusters", "ancestry", "samples"):
    setattr(args, k, os.path.abspath(getattr(args, k)))
```

**`os.path.abspath("")` returns the current working directory** - so "ancestry not supplied" was silently turned into "a directory path", the downstream check `if args.ancestry:` concluded it had been supplied, and the CWD was passed to stage 02 as the lineage table (`PermissionError`).
v3 never hit this because `--ancestry` was mandatory; making ancestry optional is what exposed it. **Fixed** (empty values are no longer passed through abspath).

---

## 8. Honest boundaries

| # | Boundary |
|---|---|
| 1 | **This change was validated end to end on simulated data only** (7 copies, known truth). **No real data** has been run in any of the three modes |
| 2 | **No real two-haplotype reference exists to test `diploid` with** - its semantics and implementation follow the "2 copies per cluster" design, which is equivalent to `auto` with a copy number of 2 |
| 3 | **Haplotype references are explicitly out of scope** (by user decision). The pipeline warns, but **does not block the run** - blocking would get in the way of people who want to try it first |
| 4 | **The actual detection limit of `auto` is unmeasured**: it is expected to be worse than `allo` across lineages and on a par with ">=3-5x within a lineage", but it has never been calibrated separately |
| 5 | The `diploid` detection limit (only one competitor to exclude -> theoretically easier) **is unmeasured**; it is a mechanistic inference only |
| 6 | Without a lineage layer, `gene_counts.tsv` still carries an `Ancestry` column (written empty) - this is **deliberate**: it keeps the matrix schema identical across the three modes, so downstream scripts need no change |
