# v4.2 - ChIP-seq support

**In one sentence**: from v4.2 on, one extra configuration line, `[map] assay = atac | chip`, is all it takes to run ChIP-seq.
**Not a single downstream line changes (stages 01-05 plus the statistics layer)** - because the pipeline's core proposition, "assign reads to homeologous copies", does not depend on where the reads came from.

- Version: v4.2 (2026-10-13)
- Files involved: `config.ini` / `config.template.ini` (new `[map] assay` and `macs2_extra`),
  `bin/00_map_call.sh`, `bin/00b_map_call_sensitive.sh`, `atacquant_core.py`, `run_atacquant.py`

---

## 1. Why only the upstream part needs changing

ATAC-Quant really has only three inputs: **a BAM with multi-mapping enabled + a set of intervals (peaks -> features) + a cluster table / lineage table**.
It neither knows nor cares whether those reads were cut out by Tn5 or enriched by an antibody. The physical model behind the EM allocation is
"every fragment is divided among its candidate intervals in proportion to *evidence x density*" - which has nothing to do with the assay.

| Stage | ATAC-specific? |
|---|---|
| 00 alignment + peak calling | **two ATAC-specific spots** (see section 2) |
| 01 consensus peak merging | no - generic |
| 02 HRU feature construction (15-mer projection) | no - generic |
| 03 EC-EM quantification | no - fully generic |
| 04 matrices / 05 QC / statistics layer | no - generic |

---

## 2. The two ATAC-specific spots, now parameterised in v4.2

### (1) Tn5 offset correction (controlled by `assay`)

```bash
# assay = atac (default): the Tn5 transposase nicks are staggered (+4/-5 on the two strands);
# shift back to the true cut site first
alignmentSieve -b "$PRIM" -o ... --ATACshift -p "$THREADS"

# assay = chip: antibody enrichment has no such physical step -> skip it
# (hand primary.bam straight to MACS2)
```

Keeping it -> every peak is shifted by 4-5 bp overall. The numerical effect is small, but **it is conceptually wrong**, and it will be questioned in the methods.

### (2) The primary-alignment subset must be paired correctly (v4.2: **automatic detection**)

```bash
# old: samtools view -b -f 2 -F 0xD04 ...     <- -f 2 selects nothing on single-end data -> no peaks come out
# new: sample 20,000 records automatically to decide single- vs paired-end
#      paired-end -> -f 2 plus MACS2 -f BAMPE
#      single-end -> no -f 2, MACS2 -f BAM (and a warning is printed)
```

---

## 3. Hard requirement: **paired-end only**

`frag_span()` in `bin/03_quantify.py`:

```python
def frag_span(r, flag):
    tlen = r.template_length
    if (flag & 1) and tlen != 0:        # paired-end: TLEN gives the true fragment span
        ...
    s = r.reference_start               # single-end: degrades to the read's own span (~read length)
    e = r.reference_end
    return s, (e if e is not None else s + 1)
```

| Data type | Fragment span | Consequence |
|---|---|---|
| **paired-end ChIP** | TLEN -> the true fragment (150-300 bp) | all the logic holds |
| **single-end ChIP** | the read itself (50-75 bp) | (1) `EffLen = L + mean fragment length - 1` is computed wrongly (2) reads at the peak edges may not overlap at all |

**Single-end data will "run through", but the fragment semantics are wrong, so the conclusions are not trustworthy.** Stage 00 detects this automatically and prints a warning.
(Really supporting single-end data would require changing `frag_span` to estimate the fragment length by cross-correlation and extend by a fixed amount - not done.)

---

## 4. ChIP types ranked by fitness

| Fitness | Type | Why |
|---|---|---|
| **5/5** | **histone modifications** (H3K27ac / H3K4me3 / H3K9ac / H3K27me3) | **broad peaks** (hundreds of bp to kb), the peak shape closest to ATAC -> high 15-mer projection success; and **H3K27ac is closer to "regulatory activity" than ATAC** |
| 4/5 | broad-peak TFs (CTCF, cohesin) | the peaks are wide enough |
| 2/5 | **narrow-peak TFs** (most transcription factors) | the peaks are only 100-300 bp -> projection requires `--proj-min-anchors 4` + `--proj-min-cov 0.3`, so **the projection rate drops markedly** (relaxing the thresholds introduces spurious projections instead - a trade-off) |

**Histone-modification ChIP is the ideal second modality.**

---

## 5. How to use it

```ini
[map]
assay = chip          # atac = Tn5 transposase / chip = antibody enrichment
qvalue = 0.01         # broad histone peaks are often tightened to 0.01; use 0.05 for TFs
gsize = 1.0e10        # set this to your effective genome size
macs2_extra = --broad # (optional) --broad for broad marks; the product is _peaks.broadPeak, recognised automatically
```

```bash
python run_atacquant.py --config config_chip.ini --check    # inspect the plan first
python run_atacquant.py --config config_chip.ini
```

`--check` prints the assay explicitly:

```
  steps    : map
  assay    : chip (antibody enrichment -> Tn5 offset correction skipped); extra MACS2 args: --broad
```

### Other points

| # | Item | Explanation |
|---|---|---|
| 1 | **paired-end** | mandatory (section 3). When starting from FASTQ keep `-1/-2` in `[samples]` |
| 2 | existing BAM without peaks | write `-` for peaks in `[samples]`; if `steps` includes `map`, peaks are called automatically |
| 3 | `macs2_extra` | appended verbatim to the MACS2 command line (e.g. `--nomodel --extsize 200` for TFs) |
| 4 | fragment length | **do not set it** - the pipeline estimates it from the BAM |
| 5 | `[quant]` section | unchanged (`as_temp` / `bootstrap` / `lineage_map` all stay the same) |
| 6 | `alignmentSieve` | **chip mode no longer requires deeptools to be installed** (validation is now assay-aware) |

---

## 6. Two old bugs fixed along the way in v4.2 (both of which blocked ChIP users)

| # | Bug | Symptom | Fix |
|---|---|---|---|
| 1 | **`-` in the inline `[samples]` form was not normalised to an empty value** | the documented "FASTQ only, write `-` for bam/peaks" was rejected by validation with the error `BAM not found: -` | the inline path now applies the same rule as an external sample table (`-`/`NA`/`None`/empty -> empty) |
| 2 | **single-token `*_extra` arguments failed in argparse** | `macs2_extra = --broad` raised `expected one argument` (argparse only treats the value as a positional argument when it **contains a space**) | values are now passed as `--opt=value`; `quant_extra` / `feature_extra` were fixed as well |

---

## 7. Honest boundaries

| # | Boundary |
|---|---|
| 1 | **This change was not tested on real ChIP data.** The logic passed bash syntax checks, argument validation, full-chain `--check` / `--dry-run` verification and a unit test of the paired-end detection; **but no real ChIP BAM has ever been run** (the development machine had no samtools/bowtie2/MACS2) |
| 2 | **Single-end ChIP is not supported** (section 3); only a warning is printed, not a block - blocking would get in the way of people who want to try it first |
| 3 | **The projection rate for narrow-peak TFs is unmeasured**: it is expected to drop, but by how much there is no data |
| 4 | The broadPeak output of `--broad` is confirmed to be parseable by `read_peaks` (only the first 4 columns are used), but **broad mode has never been tested end to end** |
| 5 | The effect of peak-width distributions on the stage-02 projection thresholds (`--proj-min-anchors` / `--proj-min-cov`) is unquantified |
