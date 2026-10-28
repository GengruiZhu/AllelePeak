# ATAC-Quant v4

**Copy-resolved quantification of ATAC-seq and ChIP-seq in polyploids.**

A peak-based workflow that stops treating a polyploid genome as if every regulatory
element existed once. For every consensus peak it builds a *Homeologous Regulatory
Unit* (HRU) - the set of homeologous copies of that element - projects the peak onto
each copy, and then allocates the ambiguous fragments to those copies by EM. The
result is a count matrix indexed by **element x copy x sample** instead of
**element x sample**, plus an explicit measure of how much of each count is
unambiguous evidence.

The pipeline is species-agnostic. Every threshold and every grouping comes from the
files you supply; nothing about a particular organism is hard-coded.

---

## What problem it solves

In a diploid, reads that align to several positions are a nuisance. In a polyploid
they are the signal: reads from copy A of a gene look like reads from copy B, so plain
counting collapses 10-14 real copies into one number and any statement about
"subgenome dominance", "copy-specific accessibility" or "homeolog expression bias"
rests on a number that never existed.

ATAC-Quant keeps them apart without pretending the ambiguity is gone:

- fragments that hit exactly one copy are counted as such (`N_Unique`);
- fragments compatible with several copies are split by an EM over the equivalence
  class, weighted by effective length;
- every feature carries `Resolved_Frac`, the fraction of its count that is
  undisputed - so a copy-level claim can be tiered by how much evidence is behind it.

EM results are mass-conserving: the sum over the copies of a group equals the number
of fragments that entered it. That property is checked in the QC report, not assumed.

---

## Requirements

```bash
conda env create -f environment.yml      # python>=3.7 + numpy + pysam
# mapping from FASTQ additionally needs: bowtie2, samtools, deeptools (alignmentSieve), MACS2
```

Python 3.7+ with `numpy` and `pysam` is the hard requirement. The mapping tools are
needed only when you start from FASTQ or ask the pipeline to call peaks.

**Paired-end is mandatory for ChIP-seq** and strongly recommended for ATAC-seq: step 03
infers fragment span from the TLEN field, which degenerates to read length for
single-end data and silently corrupts every effective length.

Self-test (1-3 minutes, generates its own simulated data):

```bash
bash tests/run_test.sh ./_selftest      # last line should be: ALL CHECKS PASSED
```

---

## Three steps to a result

### 1. Fill in one configuration file

```bash
cp config.template.ini config.ini
# edit the five inputs under [inputs]; declare samples under [samples] at the end
```

The five inputs are the only things the pipeline needs to know about your data:

| Key | What it is | Format | Where it comes from |
|---|---|---|---|
| `genome` | Genome FASTA | one file, any species | your assembly |
| `gff` | GFF3 annotation | at least `gene` lines with `ID=` | your annotation |
| `clusters` | Homeologous-copy clusters | `ClusterID <TAB> gene1,gene2,...` | your copy-grouping tool (e.g. PolyAllele) |
| `ancestry` | Gene to lineage assignment | a header with an `Ancestry`/`Lineage` column, or two columns (gene, lineage) | subgenome annotation |
| `samples` | Sample sheet | `sample r1 r2 bam peaks`, missing values written as `-` | your data |

Set `ploidy` to match the biology: `allo` (allopolyploid, lineage layer produced),
`auto` (autopolyploid, no lineage layer), or `diploid` (two-haplotype reference,
2 copies per cluster). See [docs/02_ploidy_modes.md](docs/02_ploidy_modes.md).

### 2. Validate, do not guess

```bash
python run_atacquant.py --config config.ini --check
```

`--check` validates every path, every required tool and the bowtie2 index, then prints
the exact plan and the exact command line it would run. It returns before touching the
filesystem, so it is safe to run on a login node. `--dry-run` prints the commands
without executing the steps.

### 3. Run

```bash
python run_atacquant.py --config config.ini          # checkpoint resume is automatic
```

Interrupted? Run the same command again: completed steps are recorded in
`.checkpoints/` and skipped. On a scheduler, edit `run_atacquant.pbs` and `qsub` it.

---

## What you get (under `outdir/`)

| Location | Content |
|---|---|
| `qc/qc_report.md` | **read this first**: mass balance, identifiability overview |
| `features/features.tsv`, `groups.tsv` | The HRU feature table (element x copy) and the group table |
| `quant/<sample>/<sample>.quant.tsv` | **the core result**: per-copy quantification (EM counts and evidence columns) |
| `matrix/` | Four count matrices: feature, group, group x lineage, gene |
| `stats/<sample>/` | Statistics layer: tier stratification, group-level tests, detection limits (MDD) |
| `run_atacquant.log` | Master log, with a pointer to each step's own log |

Per-sample convergence diagnostics live in `quant/<sample>/<sample>.quant_stats.json`
(`em_iterations`, `em_converged`, `em_seconds`). Check `em_converged` before quoting
bootstrap standard deviations.

---

## Reading the result honestly

Copy-level numbers are not all equally trustworthy, and the pipeline reports enough to
tell which is which:

- a feature's `Resolved_Frac` says how much of its count is undisputed evidence;
- the tier stratification puts each feature in a band by identifiability;
- the detection limit (MDD) says how large a difference must be before it is
  distinguishable from allocation noise.

Practical rule of thumb: **across lineage >= 2x, within lineage >= 3-5x** before
believing a copy-level difference. Below that, use the lineage-level summary, which is
far more robust than any individual copy. Details in
[docs/05_copy_level_statistics.md](docs/05_copy_level_statistics.md).

---

## Notes on this distribution

- English-only and species-agnostic: no organism-specific values, paths or thresholds
  are baked in, and no Chinese text remains anywhere in the package.
- It carries the v4.4 duplicate-marking fix: duplicate marking always goes through
  `samtools markdup` and the `sambamba` branch was removed. `sambamba markdup` is
  tuned for ~10-20M read pairs and degrades badly past that; one 145M-pair library
  ran for 62 hours without progress and had to be killed. `sambamba` is therefore no
  longer a dependency.

---

## Documentation

| Document | What it covers |
|---|---|
| [docs/00_quickstart.md](docs/00_quickstart.md) | The short version: setup, inputs, run, FAQ |
| [docs/01_from_peaks_to_copy_fractions.md](docs/01_from_peaks_to_copy_fractions.md) | Step by step from a peak to a copy fraction: the projection rule, the EM, the evidence columns |
| [docs/02_ploidy_modes.md](docs/02_ploidy_modes.md) | `allo` / `auto` / `diploid`, what actually changes (less than you think) and what does not |
| [docs/03_chip_seq_support.md](docs/03_chip_seq_support.md) | Running ChIP-seq instead of ATAC-seq: what changes, what to avoid |
| [docs/04_performance_v4.md](docs/04_performance_v4.md) | What v4 changed, measured speedups, how to verify them |
| [docs/05_copy_level_statistics.md](docs/05_copy_level_statistics.md) | The statistics layer: tiers, tie-aware group tests, MDD, and when not to use them |

---

## Command-line reference

```
python run_atacquant.py --config CONFIG [--check] [--dry-run] [--force]
                        [--steps map,peaks,features,quant,matrix,qc]
                        [--no-stats]
```

| Flag | Effect |
|---|---|
| `--check` | Validate inputs, tools and index; print the plan; run nothing and write nothing |
| `--dry-run` | Print the commands that would be executed, do not execute them |
| `--force` | Ignore checkpoints and redo the requested steps |
| `--steps` | Run only the listed steps, in pipeline order |
| `--no-stats` | Skip the statistics layer even if `[stats] run = true` |

---

## Layout of this package

```
config.template.ini     the only input definition (copy to config.ini)
samples.template.tsv    sample-sheet format, if you prefer a separate file
environment.yml         conda environment
run_atacquant.py        entry point: validation, planning, orchestration
atacquant_core.py       step dispatcher
run_atacquant.pbs       PBS job template
docs/                   the manual
bin/                    00 map+call peaks, 01 consensus peaks, 02 features,
                        03 quantify (EM), 04 matrices, 05 QC report
lib/                    shared helpers (FASTA, GFF, clusters, ancestry, mock BAM)
analysis/               statistics layer and bias diagnostics
tests/                  self-test, simulator, ground-truth evaluation, benchmark
```

---

## Honest boundaries

Stated up front, because they decide which conclusions the data can support:

1. **Copy-level resolution has a floor.** Copies that differ by less than the
   detection limit are, for practical purposes, one number. The MDD column tells you
   where that floor is for each feature; it is not optional reading.
2. **`projected` features are inferred, not observed.** They have no peak of their
   own; they exist because a homeologous copy is expected there. They are marked, and
   their detection rate is much lower than for `called` features.
3. **Autopolyploids are much harder than allopolyploids.** Copy divergence of 0.1-1%
   pushes the detection limit to 3-5x, and "subgenome dominance" is not even a
   meaningful concept in that case. The pipeline says so rather than producing a
   plausible-looking deconvolution.
4. **A consensus/haploid reference makes deconvolution impossible.** If every gene
   appears once in the reference, there is nothing to split. The pipeline detects this
   (nearly all clusters are singletons) and warns explicitly instead of returning a
   result that looks like deconvolution but is just peak counting.
