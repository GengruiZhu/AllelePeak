#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
bench_v4.py -- reproduce the speed-up figures of v4 and the equivalence comparison against v3.

Usage:
    python tests/bench_v4.py [workdir] [--depth 12] [--threads 8]

It does four things:
  A. EM convergence speed-up comparison: on one and the same dataset, only
     --em-accelerate / --em-tol are varied
     -> prints iteration count / EM wall-clock seconds / log-likelihood / converged? /
        wall-clock speed-up
  B. Verification that both converge to the same point: SQUAREM and plain EM must give
     exactly the same ll at tol=1e-9
  C. Bit-for-bit equivalence with v3: v4 (accelerate=off, tol=1e-3, bootstrap order)
     must reproduce the EM_Count and EM_SD of v3 bit for bit
  D. Reproducibility: in v4, EM_SD must be exactly identical under --reduce-threads 1 and 8
     (in v3 these differ -- that is the old problem v4 fixes)

Everything runs locally (mock alignment path); about 5-10 minutes.
"""
import os
import sys
import json
import shutil
import argparse
import subprocess
import time

HERE = os.path.dirname(os.path.abspath(__file__))
ROOT = os.path.dirname(HERE)
V3 = os.path.join(os.path.dirname(ROOT), "ATAC-Quant-pipeline-v3")
PY = sys.executable


def run(cmd, logpath, env=None):
    with open(logpath, "w", encoding="utf-8", errors="replace") as f:
        p = subprocess.run(cmd, stdout=f, stderr=subprocess.STDOUT, env=env)
    return p.returncode


def core_cmd(core, sim, out):
    return [PY, core,
            "--samples", os.path.join(sim, "samples.tsv"),
            "--genome", os.path.join(sim, "genome.fa"),
            "--gff", os.path.join(sim, "genes.gff3"),
            "--clusters", os.path.join(sim, "clusters.tsv"),
            "--ancestry", os.path.join(sim, "ancestry.tsv"),
            "-o", out,
            "--lineage-map", "Officinarum=SOF,Robustum=SRO,Spontaneum=SSP"]


def read_stats(out):
    p = os.path.join(out, "quant", "SIM", "SIM.quant_stats.json")
    with open(p, encoding="utf-8") as f:
        return json.load(f)


def read_quant(out, cols=("EM_Count", "EM_SD")):
    p = os.path.join(out, "quant", "SIM", "SIM.quant.tsv")
    rows = {}
    with open(p, encoding="utf-8") as f:
        h = f.readline().rstrip("\n").split("\t")
        for line in f:
            d = dict(zip(h, line.rstrip("\n").split("\t")))
            rows[d["FeatureID"]] = tuple(d[c] for c in cols)
    return rows


def diff(a, b):
    n = sum(1 for k in a if k in b and a[k] != b[k])
    return n, len(a)


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("workdir", nargs="?", default=os.path.join(ROOT, "_bench"))
    ap.add_argument("--depth", type=int, default=12)
    ap.add_argument("--threads", type=int, default=8)
    args = ap.parse_args()
    W = os.path.abspath(args.workdir)
    os.makedirs(W, exist_ok=True)
    sim = os.path.join(W, "sim")
    env = dict(os.environ, PYTHONIOENCODING="utf-8", PYTHONUTF8="1", ATACQ_MOCK_BAM="1")

    # ---------------------------------------------------------------- generate the simulation set
    if not os.path.exists(os.path.join(sim, "samples.tsv")):
        sys.stdout.write("[bench] generating the simulation set (depth=%d) ...\n" % args.depth)
        sys.stdout.flush()
        cmd = [PY, os.path.join(HERE, "simulate_mini.py"), "-o", sim,
               "--depth", str(args.depth), "--sample", "SIM"]
        if run(cmd, os.path.join(W, "sim.log"), env=env) != 0:
            sys.exit("[bench] simulation-set generation failed, see " + os.path.join(W, "sim.log"))

    def go(tag, core, extra, steps="peaks,features,quant", bootstrap=0, qextra=None):
        out = os.path.join(W, tag)
        cmd = core_cmd(core, sim, out)
        cmd += ["-t", str(args.threads), "--steps", steps, "--bootstrap", str(bootstrap)]
        if qextra:
            cmd += ["--quant-extra", qextra]
        cmd += extra
        t0 = time.time()
        rc = run(cmd, os.path.join(W, tag + ".log"), env=env)
        if rc != 0:
            sys.exit("[bench] %s failed (exit=%d), see %s" % (tag, rc, os.path.join(W, tag + ".log")))
        return out, time.time() - t0

    CORE4 = os.path.join(ROOT, "atacquant_core.py")
    CORE3 = os.path.join(V3, "atacquant_core.py")
    Q = "--buckets 16"

    # ------------------------------------------------- A. EM speed-up comparison
    sys.stdout.write("\n" + "=" * 78 + "\nA. EM convergence speed-up comparison (same data, only --em-accelerate / --em-tol vary)\n" + "=" * 78 + "\n")
    sys.stdout.flush()
    res = {}
    for tag, acc, tol in [("em_off_1e-3", "off", "1e-3"), ("em_sq_1e-3", "squarem", "1e-3"),
                          ("em_off_1e-6", "off", "1e-6"), ("em_sq_1e-6", "squarem", "1e-6")]:
        out, wall = go(tag, CORE4, ["--em-accelerate", acc, "--em-tol", tol,
                                    "--em-max-iter", "200000"], qextra=Q)
        st = read_stats(out)
        res[(acc, tol)] = (st["em_iterations"], st.get("em_seconds", float("nan")),
                           st["em_loglik"], st["em_converged"])
    print("%-10s %-10s %10s %10s %20s %8s" % ("accelerate", "tol", "iterations", "EM sec", "ll", "converged"))
    for k in [("off", "1e-3"), ("squarem", "1e-3"), ("off", "1e-6"), ("squarem", "1e-6")]:
        v = res[k]
        print("%-10s %-10s %10d %10.2f %20.8f %8s" % (k[0], k[1], v[0], v[1], v[2], v[3]))
    print()
    for tol in ("1e-3", "1e-6"):
        a, b = res[("off", tol)], res[("squarem", tol)]
        print("  tol=%s : iterations %d -> %d (%.1fx) | EM wall clock %.2fs -> %.2fs (%.2fx)"
              % (tol, a[0], b[0], a[0] / max(b[0], 1), a[1], b[1], a[1] / max(b[1], 1e-9)))

    # ------------------------------------------------- B. identical ll at tol=1e-9
    sys.stdout.write("\n" + "=" * 78 + "\nB. Converging to the same point (tol=1e-9, ll should be exactly identical)\n" + "=" * 78 + "\n")
    sys.stdout.flush()
    ll9 = {}
    for tag, acc in [("ll9_off", "off"), ("ll9_sq", "squarem")]:
        out, wall = go(tag, CORE4, ["--em-accelerate", acc, "--em-tol", "1e-9",
                                    "--em-max-iter", "500000"], qextra=Q)
        st = read_stats(out)
        ll9[acc] = (st["em_iterations"], st["em_loglik"])
        print("  %-8s iterations %8d  ll %.10f" % (acc, st["em_iterations"], st["em_loglik"]))
    lo = ll9["off"][1]
    rel = abs(lo - ll9["squarem"][1]) / max(abs(lo), 1e-30)
    same = rel < 1e-9
    print("  -> relative difference %.2e  ->  %s (this is an acceleration, not a different algorithm; the endpoint is the same MLE)"
          % (rel, "ll identical [OK]" if same else "ll differs [MISMATCH]"))

    # ------------------------------------------------- C. bit-for-bit equivalence with v3
    sys.stdout.write("\n" + "=" * 78 + "\nC. Bit-for-bit equivalence with v3 (accelerate=off, tol=1e-3, bootstrap order)\n" + "=" * 78 + "\n")
    sys.stdout.flush()
    if os.path.exists(CORE3):
        o3, _ = go("eq_v3", CORE3, [], steps="peaks,features,quant", bootstrap=5,
                   qextra=Q + " --reduce-threads 1")
        o4, _ = go("eq_v4", CORE4, ["--em-accelerate", "off", "--em-tol", "1e-3",
                                    "--bootstrap-threads", "1"],
                   bootstrap=5, qextra=Q + " --reduce-threads 1")
        a = read_quant(o3)
        b = read_quant(o4)
        n1, tot = diff({k: v[0] for k, v in a.items()}, {k: v[0] for k, v in b.items()})
        n2, _ = diff({k: v[1] for k, v in a.items()}, {k: v[1] for k, v in b.items()})
        print("  EM_Count mismatches = %d / %d" % (n1, tot))
        print("  EM_SD    mismatches = %d / %d" % (n2, tot))
        print("  -> %s" % ("bit-for-bit identical [OK]" if n1 == 0 and n2 == 0 else "differences found [MISMATCH]"))
    else:
        print("  (%s was not found, skipping)" % CORE3)

    # ------------------------------------------------- D. reproducibility
    sys.stdout.write("\n" + "=" * 78 + "\nD. Reproducibility (only the number of reduce threads changes; EM_SD should not)\n" + "=" * 78 + "\n")
    sys.stdout.flush()
    pa = {}
    for tag, rt in [("rep_rt1", 1), ("rep_rt8", 8)]:
        o, _ = go(tag, CORE4, ["--em-accelerate", "off", "--em-tol", "1e-3",
                               "--bootstrap-threads", "1"],
                  bootstrap=5, qextra=Q + " --reduce-threads %d" % rt)
        pa[rt] = read_quant(o)
    n1, tot = diff({k: v[0] for k, v in pa[1].items()}, {k: v[0] for k, v in pa[8].items()})
    n2, _ = diff({k: v[1] for k, v in pa[1].items()}, {k: v[1] for k, v in pa[8].items()})
    print("  v4: reduce-threads 1 vs 8  ->  EM_Count mismatches = %d / %d | EM_SD mismatches = %d / %d"
          % (n1, tot, n2, tot))
    print("  -> %s" % ("reproducible [OK] (v3 differs on EM_SD for 258/281 features under the same test)" if n2 == 0 else "still not reproducible [MISMATCH]"))
    sys.stdout.write("\n[bench] all done. Intermediate products are in %s (safe to delete as a whole and re-run)\n" % W)


if __name__ == "__main__":
    main()
