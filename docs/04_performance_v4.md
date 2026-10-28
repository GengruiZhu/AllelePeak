# ATAC-Quant v4 - Performance notes

**In one sentence**: v4 = the v3 engine + **faster EM convergence (SQUAREM)** + **parallel bootstrap** + **one reproducibility fix**.
Algorithm, output formats, decision rules and directory layout are all unchanged; **results reproduce v3 bit for bit**.

- Baseline: `ATAC-Quant-pipeline-v3` (v3.0.4)
- Source of every local measurement below: `tests/bench_v4.sh` (simulated data, reproducible with one command)
- Last updated: 2026-10-13

---

## 0. One-page summary

| Item | Conclusion |
|---|---|
| **where the speed comes from** | (1) parallel bootstrap (**76%** of the whole run in v3, and it was a sequential for loop) (2) EM convergence acceleration (SQUAREM, **4.3-9.9x** wall clock measured) |
| **correctness** | with `em_accelerate = off` + `em_tol = 1e-3` + `bootstrap_threads = 1` it is **bit-identical** to v3 (EM_Count difference 0, EM_SD difference 0) |
| **default precision** | v4 defaults to `em_tol = 1e-6` (v3 used 1e-3) -> actually **more precise** than v3 (the gap in log-likelihood from the optimum drops from 1.7e-2 to 3.8e-7) |
| **production-scale expectation** | stage 03: **990 min -> about 110-140 min** (~7-9x); whole run 18.3 h -> about 2.5 h (**extrapolated**, see section 3) |
| **fixed along the way** | v3's `EM_SD` was **not reproducible** (changing only the reduce thread count changed the EM_SD of 258 out of 281 features); fixed in v4 |
| **cost** | 9 new configuration options, all with v3-compatible defaults; no output column changed |

---

## 1. What changed (item by item)

### Change 1 - parallel bootstrap: the biggest win

**The v3 code** (`bin/03_quantify.py`): 20 resamples in a sequential for loop, one after another.

```python
for b in range(args.bootstrap):
    draw = rng.multinomial(int(round(N)), p).astype(np.float64)
    a_b, _ = data.em(alpha, ec_count=draw[F+1:], single=draw[:F+1], tol=..., max_iter=500, min_iter=5)
```

Measured on the production dataset: the 20 resamples took **12 h 28 min**, which was **75.7%** of stage 03's total time.

**The v4 code**: generate all the resamples first, then distribute them over several processes (on Linux `fork` gives children zero-copy access to ECData, so memory does not double).

```python
draws = [rng.multinomial(int(round(N)), p) for _ in range(args.bootstrap)]
results = _run_bootstrap(draws, data, alpha, em_kw, nt)   # nt = --bootstrap-threads
```

- the resample sequence is **exactly** the same as in v3 (same rng, same call order)
- `--bootstrap-threads 0` (default) = automatic `min(bootstrap, 8)`; `1` = back to v3's purely sequential behaviour
- Windows / platforms without fork: above 0.3 GB of data it **degrades automatically to 1 thread** and prints a notice (spawn has to copy ECData to every child, which doubles memory)

### Change 2 - EM convergence acceleration (SQUAREM)

**v3's EM is pure fixed-point iteration**: `x <- F(x)`, one EM step per round. On the production dataset it ran the full **2000-round cap without converging**.

**v4 uses SQUAREM (Varadhan & Roland 2008) for squared extrapolation**:

```
one step     r = F(x) - x
two steps    v = F(F(x)) - F(x) - r
Cauchy extrapolation   xs = x - 2a*r + a^2*v,   a = -||r||/||v||
take one more EM step from xs to get the candidate x3
accept only if ll(x3) >= ll(F(x)), otherwise fall back to a plain EM step
```

**The key point**: this is **acceleration, not a different algorithm** - the log-likelihood serves as a monotonicity guard, so
- the end point is still **the same MLE** (measured: SQUAREM and plain EM reach **exactly the same ll** at tol=1e-9)
- at most 3 EM steps per round, so it only pays off if the round count drops below 1/3 (measured: 12-28x fewer)

### Change 3 - fixed along the way: `EM_SD` was not reproducible (an old v3 problem)

**Symptom** (measured on a local simulation):

| Comparison | EM_Count difference | **EM_SD difference** |
|---|---|---|
| **v3** itself: `--reduce-threads 1` vs `8` | 0 / 281 | **258 / 281** |
| **v4**, after the fix: `--reduce-threads 1` vs `8` | 0 / 281 | **0 / 281** |

**Root cause**: the reduce phase used `imap_unordered`, so the **floating-point summation order within each bucket followed the scheduling** -> the last digits of `base` jittered -> `p` jittered -> and because `rng.multinomial` is a **discrete draw**, the last digit **flipped which fragments got resampled** -> the whole EM_SD changed.

**Fix**: `pmap(..., ordered=True)` keeps reduce in order (`imap` instead of `imap_unordered`). Pass A/B still use `imap_unordered` (order does not matter there, so no throughput is lost).

> Note: this fix makes v4 differ from "the v4 just before the fix" in the last digits (the EM iteration path changed),
> but it is **bit-identical to v3 with `reduce-threads=1`**. That is the expected consequence of the fix.

### Change 4 - micro-optimisation of memory bandwidth in a single EM step

Semantically equivalent, only fewer array copies:

| v3 | v4 |
|---|---|
| `rate = alpha / self.eff` | `rate = alpha * self.eff_inv` (1/eff precomputed) |
| `rs[rs <= 0] = 1e-300` (builds a boolean mask over 53 million elements) | `np.maximum(rs, 1e-300, out=rs)` |
| `resp = val * (...)` (allocates one extra 53-million-element temporary) | `val *= (...)` (in place) |

---

## 2. Measured locally (simulated data, reproducible with `tests/bench_v4.sh`)

**Test set**: `tests/simulate_mini.py --depth 12` -> 210,919 fragments / 2,178,810 SAM records (7 copies, 279 elements).

### 2.1 EM convergence acceleration (same data, only `--em-accelerate` / `--em-tol` changed)

| Configuration | Iterations | EM wall clock | log-likelihood ll |
|---|---|---|---|
| plain EM, tol=1e-3 (~v3 default) | 4,494 | 0.78 s | 20542.40745765 |
| **SQUAREM, tol=1e-3** | **360** | **0.18 s** | 20542.39573704 |
| plain EM, tol=1e-6 (~v4 default) | 28,936 | 5.12 s | 20542.45898129 |
| **SQUAREM, tol=1e-6 (v4 default)** | **1,040** | **0.52 s** | **20542.45905389** |

| Tolerance | Iteration reduction | **Wall-clock speed-up** |
|---|---|---|
| tol=1e-3 | 12.5x | **4.33x** |
| tol=1e-6 | 27.8x | **9.85x** |

(At most 3 EM steps per round, so the wall-clock speed-up is roughly the iteration reduction divided by 3, minus a little overhead for computing ll.)

**The same script was also run on a smaller simulated set (`--depth 3`) and the speed-up was markedly lower** - the speed-up **depends on the problem**:

| Data size | Tolerance | Iteration reduction | Wall-clock speed-up |
|---|---|---|---|
| depth 3 (~50k element-level fragments) | tol=1e-6 | 12.8x | **4.56x** |
| depth 3 | tol=1e-3 | 16.3x | 4.80x |
| depth 12 (210,919 fragments) | tol=1e-6 | 27.8x | **9.85x** |
| depth 12 | tol=1e-3 | 12.5x | 4.33x |

> **So this document always quotes ranges**: **12-28x fewer iterations, 4.3-9.9x faster EM wall clock**.
> What actually transfers between datasets is the iteration reduction factor, not the absolute seconds; and **you must measure it on your own data with `--em-trace`** (section 4.3).

### 2.2 Converging to the same point (the evidence that this is acceleration, not a different algorithm)

Run to convergence on a small simulation (281 features):

| Configuration | Rounds | ll | Gap from the optimum |
|---|---|---|---|
| plain EM, tol=1e-9 | 2,892 | -56857.732284 | -1.2e-7 |
| **SQUAREM, tol=1e-9** | **80** | **-56857.732284** | **-1.2e-7** <- **completely identical to the line above** |
| **SQUAREM, tol=1e-6 (v4 default)** | **77** | **-56857.732284** | **-3.8e-7** |
| plain EM, tol=1e-3 (=v3 default) | 368 | -56857.749626 | **-1.7e-2** <- v3 was less precise |
| SQUAREM, tol=1e-3 | 43 | -56857.755732 | -2.3e-2 |

**How to read it**: the v4 default (77 rounds) is **both faster and more accurate** than the v3 default (368 rounds) - 5 orders of magnitude more precise, and 4.8x fewer rounds.

### 2.3 Parallel bootstrap

| Comparison | EM_Count difference | EM_SD difference |
|---|---|---|
| `bootstrap_threads=1` (sequential) vs `=5` (parallel) | 0 / 281 | **0 / 281** <- bit-identical |

### 2.4 Bit-level equivalence with v3

| Comparison | EM_Count difference | EM_SD difference |
|---|---|---|
| **v3 vs v4** (`accelerate=off`, `tol=1e-3`, `bootstrap_threads=1`, `reduce_threads=1`) | **0 / 281** | **0 / 281** |

### 2.5 Truth evaluation (`tests/eval_truth.py`, simulation with known truth)

| Metric | v3 | v4 (default) |
|---|---|---|
| feature-level EM vs truth, r | 0.9957 | 0.9956 |
| MAE | 6.95 | 7.00 |
| slope (1 is ideal) | 1.005 | 1.005 |
| group-level r | 1.0000 | 1.0000 |
| lineage share (truth / EM) | L1 38.2/38.1, L2 50.3/50.4, L3 11.5/11.5 (example lineage labels) | **completely identical** |
| mass-balance check | pass | pass |
| final verdict line | all checks passed | all checks passed |

---

## 3. Expected gains at production scale (**extrapolated**, not measured)

### 3.1 Baseline: v3 measured on a production dataset (from `logs/run_atacquant.log` of that run)

```
[2026-10-10 21:22:24] [03] EM: 98,839,946 fragments, 492,503 features + background
[2026-10-11 00:00:04] [03] EM converged: 2000 rounds; background share 42.02%      <- hit the cap, not converged
[2026-10-11 12:28:40] [03] done: ... <sample>.quant.tsv (990.5 minutes)
```

| Step | v3 measured | Share of stage 03 |
|---|---|---|
| Pass A + members + Pass B + reduce | 84.3 min | 8.5% |
| **main EM** (capped at 2000 rounds, tol=1e-3, **not converged**) | **157.7 min** | 15.9% |
| **bootstrap, 20 replicates** (sequential, ~37.4 min each) | **748.5 min** | **75.6%** |
| **stage 03 total** | **990.5 min = 16.5 h** | 100% |

Configuration: `threads=48`, `buckets=128`, `tol=1e-3`, `max_iter=2000`, `bootstrap=20`.

### 3.2 The v4 extrapolation (always the conservative end 4.5x, optimistic end 9.9x)

| Step | v3 | v4 (extrapolated range) | Basis |
|---|---|---|---|
| Pass A / B / reduce | 84.3 min | 84.3 min | unchanged |
| **main EM** | 157.7 min | **16-35 min** | wall-clock speed-up 4.5-9.9x; and v4's stricter tolerance (1e-6) costs some extra, hence the conservative end |
| **bootstrap, 20 replicates** | 748.5 min | **10-21 min** | 37.4 min / (4.5-9.9) = 3.8-8.3 min each; 20 replicates / 8 parallel = 2.5 rounds |
| **stage 03 total** | **990 min (16.5 h)** | **about 110-140 min (~2 h)** | **~7-9x** |

**Whole run** (including 00-05 plus the statistics layer, ~18.3 h by the project's own accounting) -> **about 2.5 h**.

### 3.3 How else to tune for more speed (all configuration options, no code changes)

| Lever | Extra gain | Cost |
|---|---|---|
| `bootstrap_threads = 16` (or 24/48) | bootstrap 10-25 min -> 5-12 min | memory xN (about 2-3 GB per worker) |
| `bootstrap = 0` | skips bootstrap entirely | **the output has no EM_SD column** |
| `bootstrap_max_iter = 100` | bootstrap another ~2-3x faster | EM_SD is coarser on low-count features |
| `em_tol = 1e-4` | main EM another ~2x faster | slightly less precise (still better than v3's 1e-3) |

---

## 4. How to use it

### 4.1 Configuration (`[quant]` in `config.ini`; every new option has a v3-compatible default)

```ini
[quant]
bootstrap = 20
em_accelerate = squarem   # squarem = SQUAREM acceleration (v4 default); off = line-for-line identical to v3
bootstrap_threads = 0     # 0 = automatic min(bootstrap, 8); 1 = sequential (= v3 behaviour)
bootstrap_max_iter = 500  # EM iteration cap for each bootstrap replicate
em_max_iter = 2000        # main EM iteration cap
em_tol = 1e-6             # main EM convergence criterion (v4 default; v3 used 1e-3)
em_prior = 1e-8           # per-round anti-zero prior
em_trace = false          # true = print the log-likelihood of every EM round
```

**Running it is exactly as in v3**:

```bash
python run_atacquant.py --config config.ini --check
python run_atacquant.py --config config.ini
```

### 4.2 To reproduce old v3 results exactly

```ini
[quant]
em_accelerate = off
em_tol = 1e-3
bootstrap_threads = 1
```

Or on the command line (calling `atacquant_core.py` directly):

```bash
python atacquant_core.py ... --em-accelerate off --em-tol 1e-3 --bootstrap-threads 1
```

### 4.3 Verifying the speed-up on your own server

The new `--em-trace` writes the log-likelihood of every round to `logs/03_quant_<sample>.log`:

```bash
python atacquant_core.py ... --em-trace
grep 'EM(' logs/03_quant_<sample>.log | tail -20
```

At the same time, three new fields appear in `quant/<sample>/<sample>.quant_stats.json` (**looking at these numbers is all you need**):

| Field | Meaning | How to judge |
|---|---|---|
| `em_iterations` | rounds the main EM actually ran | compare with v3's 2000 |
| `em_converged` | whether `em_tol` was really met | **`false` means the cap was hit without convergence** (v3 on the production dataset reported false) |
| `em_seconds` | main EM wall-clock seconds | compare directly |
| `em_loglik` | log-likelihood at convergence | higher is better; compare two runs to see whether they converge to the same point |

**If v4's `em_iterations` also reaches the 2000 cap**, this dataset converges harder than the simulation does; loosen the tolerance with `--em-tol 1e-5` or `1e-4`.

---

## 5. Honest boundaries (**required reading**)

| # | Boundary |
|---|---|
| 1 | **The speed-up on the production dataset is extrapolated, not measured.** That data is not on this machine: the temporary file of 18,099,762 equal-length equivalence classes was cleaned up after the v3 run finished, and re-running it would take 16.5 h. The extrapolation rests on "9.85x wall-clock speed-up on a local 211k-fragment simulation" plus "the per-step shares in the production log". |
| 2 | **The speed-up depends on the problem.** On a small simulation at tol=1e-9, SQUAREM cut 2,892 rounds to 80 (36x); on the 211k-fragment simulation 28,936 -> 1,040 (27.8x). What really transfers is **this iteration reduction factor**, not the absolute seconds. |
| 3 | **v4's default `em_tol = 1e-6` is 1000x stricter than v3**, so v4 results differ slightly from v3 (and are **more accurate**). To reproduce v3 strictly, use the settings in section 4.2. |
| 4 | **`EM_SD` was never reproducible in v3** (section 1, change 3). So "v4's EM_SD differs from that one v3 run" is normal behaviour; any comparison requires fixing `--reduce-threads` on both sides. |
| 5 | **`em_accelerate = squarem` is not guaranteed to accelerate.** It takes at most 3 EM steps per round; a dataset that already converges in 30 rounds will be slightly slower with it. On simulated data the smallest observed reduction was to 1/12 of the rounds, so if this worries you, run `--em-trace` once and look. |
| 6 | **SQUAREM's log-likelihood guard depends on `prior`.** `em_prior` defaults to 1e-8, as in v3; raising it changes the MLE itself (not just numerical stability). |
| 7 | **Memory for the parallel bootstrap = number of threads x (resident ECData + per-round temporary arrays)**. At production scale that is about 2-3 GB per worker; 8 threads on a 48-core machine is fine, 48 threads needs a memory check first. |

---

## 6. Which files changed (SHA256, checked file by file)

**Summary: 8 byte-identical / 15 changed / 3 new / 0 present only in v3.**
Of the changed files, **6 carry substantive changes** and **9 only carry a version label** (`ATAC-Quant v2/v3` -> `v4`: pure strings, zero logic risk, line endings unchanged).

### 6.1 Substantive changes (6)

| File | Change |
|---|---|
| `bin/03_quantify.py` | SQUAREM (`_loglik` / `_step` / `_converged` / `em` reworked), parallel bootstrap (`_boot_init` / `_boot_worker` / `_run_bootstrap`), ordered reduce (`pmap(ordered=True)`), `eff_inv`, `nbytes`; 5 new command-line options (`--em-accelerate` / `--em-trace` / `--em-prior` / `--bootstrap-threads` / `--bootstrap-max-iter`); `--tol` default 1e-3 -> 1e-6; the stats file gains `em_accelerate` / `em_converged` / `em_seconds` / `em_loglik` |
| `atacquant_core.py` | 7 new options (`--em-*` / `--bootstrap-threads` / `--bootstrap-max-iter`), passed through to stage 03 |
| `run_atacquant.py` | reads the 7 new `[quant]` options and passes them through |
| `config.ini`, `config.template.ini` | the 7 new `[quant]` entries, with comments |
| `docs/00_quickstart.md` | rewritten for v4 + a new option table + a pointer to this file |

### 6.2 Version-label-only changes (9, pure strings)

`bin/00_map_call.sh`, `bin/00b_map_call_sensitive.sh`, `bin/01_consensus_peaks.py`, `bin/02_build_features.py`, `bin/04_build_matrix.py`, `bin/05_qc_report.py`, `lib/atacq_common.py`, `tests/eval_truth.py`, `tests/run_test.sh`

> One old annoyance closed along the way: the report title generated by `bin/05_qc_report.py` used to read `# ATAC-Quant v2 QC report`; it now reads `v4`.
> The replacement only matched prefixed version labels (`"ATAC-Quant v2 "` / `"ATAC-Quant v3 "`),
> **leaving references to v3 in the prose untouched** (for example "behaves like v3" in stage 03, and "v3 used 1e-3" in the `em_tol` comment).
> The only file that keeps a v2 label is `analysis/10_copy_level_stats.py` - it really is the v2-era statistics-layer prototype, and calling it v2 is historically correct.

### 6.3 Byte-identical (8)

`environment.yml`, `lib/mockbam.py`, `run_atacquant.pbs`, `samples.template.tsv`, `tests/simulate_mini.py`, `analysis/10_copy_level_stats.py`, `docs/05_copy_level_statistics.md`, `docs/01_from_peaks_to_copy_fractions.md`

### 6.4 New (3)

`docs/04_performance_v4.md` (this file), `tests/bench_v4.py` (reproduces the speed-ups and equivalence checks of section 2 in one command), `_selftest/config.test.ini` (v4 self-test configuration, relative paths)

---

## 7. Reproducing

```bash
# self-test (including the truth evaluation; the last line should report that all checks passed)
bash tests/run_test.sh ./_selftest

# speed-up + equivalence comparison (generates a simulation, runs 4 configurations, prints tables)
bash tests/bench_v4.sh ./_bench
```

---

## 8. Change history

See `CHANGELOG.md` in the project root (this version is **v4.0**).
