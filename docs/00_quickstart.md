# ATAC-Quant v4 (species-agnostic) - Start Here

**In one sentence**: a polyploid ATAC-seq pipeline that quantifies homeologous copies and resolves lineage, for any species. One configuration file defines every input; one command takes you from FASTQ (or an existing BAM) to the final matrices and statistics tables.

> **v4 = v3 + performance work** (faster EM convergence + parallel bootstrap + one reproducibility fix).
> Algorithm, output formats, decision rules and directory layout are **all unchanged**; results reproduce v3 bit for bit.
> **For details see [`04_performance_v4.md`](04_performance_v4.md)**: what changed, measured speed-ups, expected gains at production scale, how to verify it on your own server, and seven honest boundaries.

## Relationship to v3 / v2

- `ATAC-Quant-pipeline-v2/` was a single-project working tree: the engine plus that project's results, diagnostics and calibration output, all mixed together.
- `ATAC-Quant-pipeline-v3/` = the same engine in a general-purpose wrapper; **results produced with v3 remain valid** (v4 reproduces them bit for bit with `em_accelerate = off` + `em_tol = 1e-3` + `bootstrap_threads = 1`).
- **v4 = the v3 code + three performance / reproducibility changes**; apart from `bin/03_quantify.py`, `atacquant_core.py`, `run_atacquant.py` and the two config files, every file is byte-identical to v3.
- **When to use v4**: use v4 for any new run (the default configuration is both faster and more accurate). Switch back to the v3-equivalent settings only when you want to reproduce an old v3 result.

## New configuration options in v4 (`[quant]`, all with v3-compatible defaults)

| Option | Default | Meaning |
|---|---|---|
| `em_accelerate` | `squarem` | EM convergence acceleration; `off` = line-for-line identical to v3 |
| `em_tol` | `1e-6` | main EM convergence criterion (v3 used 1e-3; v4 is stricter -> **more accurate**) |
| `bootstrap_threads` | `0` | number of worker processes for bootstrap; `0` = automatic `min(bootstrap, 8)` |
| `bootstrap_max_iter` | `500` | EM iteration cap for each bootstrap replicate |
| `em_max_iter` / `em_prior` / `em_trace` | `2000` / `1e-8` / `false` | same semantics as v3 |

**Three fields are all you need to read the speed-up** (in `quant/<sample>/<sample>.quant_stats.json`):
`em_iterations` (number of rounds), `em_converged` (**whether it really converged**), `em_seconds` (EM wall clock).

## New in v4.3: ploidy modes (`[inputs] ploidy`)

```ini
[inputs]
ploidy = allo                 # diploid | auto | allo (names aligned with the sister project LQ's --ploidy)
ancestry = /path/ancestry.tsv # required for allo only; for auto / diploid the whole line can be deleted
```

| Mode | Applies to | ancestry | Lineage layer |
|---|---|---|---|
| `allo` (default) | allopolyploid (subgenomes differentiated) | **required** | produced |
| `auto` | autopolyploid (subgenomes undifferentiated, "the notion of ancestry does not apply") | not needed | skipped |
| `diploid` | diploid, **two-haplotype reference** (2 copies per cluster) | not needed | skipped |

**Key insight**: **there are only two algorithmic paths (with and without the lineage layer), not three** - lineage never takes part in the EM allocation.
So under all three modes the count matrices such as `quant.tsv` / `feature_counts.tsv` / `group_counts.tsv` are **byte-identical**; the differences are confined to "the two lineage-layer tables + required inputs + interpretation rules + guardrails".

- **The hard truth about `auto`**: in an autopolyploid every pair of copies is only 0.1-1% diverged -> minimum detectable difference >=3-5x (much harder than the >=2x needed across lineages in an allopolyploid);
  the concept of "subgenome dominance" **does not apply**; the real copy-level signal may live in **haplotype differences between individuals**, not in copies within an individual.
- **`diploid` supports a two-haplotype reference only**. With a haplotype consensus (each gene appears once) the pipeline **has nothing to do** -
  it warns explicitly that "the output is equivalent to ordinary peak counting, not deconvolution", and **it will not hand you a deconvolution result that looks normal but is wrong**.
- Details, the three-mode comparison, guardrails, the differences from LQ and the honest boundaries -> see [`02_ploidy_modes.md`](02_ploidy_modes.md)

## New in v4.2: ChIP-seq support (`[map] assay`)

```ini
[map]
assay = atac          # atac = Tn5 transposase (apply +4/-5 Tn5 offset correction before peak calling)
                      # chip = antibody enrichment (skip that correction; ChIP has no staggered Tn5 nicks)
macs2_extra =         # (optional) appended verbatim to MACS2, e.g. --broad for broad histone marks
```

- **Not a single downstream line changes** - stages 01-05 and the statistics layer are fully generic (the core proposition, "assign reads to homeologous copies", does not depend on where the reads came from)
- **Hard requirement: paired-end only**. `frag_span()` infers the fragment span from TLEN; on single-end data it degrades to the read span -> the effective length is computed wrongly
- **The best-fit ChIP assays are histone modifications** (H3K27ac / H3K4me3): broad peaks, high projection success rate, and semantically closer to regulatory activity than ATAC
- Details, fitness ranking, two old bugs fixed along the way, honest boundaries -> see [`03_chip_seq_support.md`](03_chip_seq_support.md)

## Three steps to a run

### 1. Environment

```bash
conda env create -f environment.yml     # or reuse an existing environment
# requires: python(>=3.7) + numpy + pysam
# the full FASTQ-to-results path additionally needs: bowtie2, samtools, deeptools(alignmentSieve), MACS2
bash tests/run_test.sh ./_selftest      # self-test (1-3 minutes); the last line should report that all checks passed
```

### 2. Fill in the configuration (the single place where inputs are defined)

```bash
cp config.template.ini config.ini
# edit config.ini: fill in the five inputs under [inputs]; write the samples under [samples] at the end of the file
```

What the five inputs are:

| Option | What it is | Format | Where it comes from |
|---|---|---|---|
| genome | genome FASTA | any species, single file | your assembly |
| gff | GFF3 annotation | at least the `gene` rows (with an `ID=` attribute) | your annotation |
| clusters | cluster table | one line per cluster: `clusterID <TAB> gene1,gene2,...` | a homeolog-grouping tool such as GraphAllele |
| ancestry | gene -> lineage table | either a header with `Ancestry`/`Lineage` columns, or two columns (gene, lineage) | subgenome annotation |
| samples | sample table | `sample r1 r2 bam peaks` (BAM or FASTQ, write `-` for whichever is absent) | your own data |

### 3. Run

```bash
python run_atacquant.py --config config.ini --check   # validate first and print the execution plan
python run_atacquant.py --config config.ini           # real run (checkpoint resume is on by default)
# on a cluster with qsub: set CONFIG= in run_atacquant.pbs, then: qsub run_atacquant.pbs
```

## Where the output goes (under outdir/)

| Location | Content |
|---|---|
| `qc/qc_report.md` | **read this first**: mass balance, identifiability overview |
| `features/features.tsv`, `groups.tsv` | the HRU feature table (element x copy) and the group table |
| `quant/<sample>/<sample>.quant.tsv` | **the core result**: copy-level quantification (EM counts and friends) |
| `matrix/` | four count matrices: feature / group / group x lineage / gene |
| `stats/<sample>/` | the three statistics-layer tables: tier assignment, group-level tests, MDD detection limits |
| `run_atacquant.log` | the master log (with a pointer to each subtask's log) |

## Common operations

- Validate only: `--check`; print the commands only: `--dry-run`; force a re-run: `--force`; run a subset of stages: `--steps features,quant`
- Resume after an interruption: run the same command again (each completed stage writes to `.checkpoints/` and is skipped automatically)
- Two ways to write the sample table, pick one: inline under `[samples]` (recommended), or point `[inputs] samples` at a separate sample-table file

## FAQ

**Q1 I only have FASTQ?** Put r1/r2 under `[samples]` (write `-` for bam/peaks), keep `steps = all`; the pipeline aligns (bowtie2 `-k 20`) and calls peaks (MACS2) for you.
**Q2 I have a BAM but no peaks?** Set bam to its path and write `-` for peaks, keep `map` in `steps`; peaks are then called from the BAM (no re-alignment).
**Q3 Do lineage names need standardising?** Only when you want to abbreviate long names (e.g. `LongName=L1,...` in `lineage_map`); if you leave it out, names are used as-is in upper case.
**Q4 Can the pipeline be used on another species?** Yes - every decision rule in the pipeline comes from the cluster table / lineage table / annotation you supply, and no species assumption is hard-coded; just remember to set `gsize` to your genome size.
**Q5 How can copy-level results be used reliably?** Check copy-level conclusions against the tiers and detection limits in [`05_copy_level_statistics.md`](05_copy_level_statistics.md) (reliable at >=2x across lineages and >=3-5x within a lineage); for low-confidence features, use lineage-level aggregates only.

## Change history

See `CHANGELOG.md` in the project root (version numbers v3.x).
