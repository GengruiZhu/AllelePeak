#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
ATAC-Quant v4 engine master - homeologous-copy / lineage-resolved quantification of polyploid ATAC-seq (any species)

  00  mapping + mark duplicates + peak calling     (per sample; skipped automatically when its BAM/peaks already exist)
  01  merge the per-sample peaks -> one consensus peak set
  02  homeologous regulatory unit (HRU): peak anchoring + homeologous copy sequence projection + grouping
  03  mapping-compatibility equivalence-class EM quantification          (per sample)
  04  multi-sample count matrix (feature / HRU group / groupxlineage / gene)
  05  QC report

Recommended usage (config-file driven, inputs defined once):
  python run_atacquant.py --config config.ini

Direct usage (skip the config file and pass every argument by hand; same as v2):
  python atacquant_core.py --samples samples.tsv --genome genome.fasta --gff ann.gff3 \\
      --clusters my_cluster.tsv --ancestry ancestry.tsv \\
      --outdir OUT --threads 48

samples.tsv (tab-separated, with a header; write - for an empty value):
  sample   r1                 r2                 bam                    peaks
  S1       -                  -                  /path/S1.dedup.bam     /path/S1_peaks.narrowPeak
  S2       /path/S2_R1.fq.gz  /path/S2_R2.fq.gz  -                      -

Checkpoint resume: each finished step writes <outdir>/.checkpoints/*.done; --force reruns from scratch; --steps runs only the named steps.
"""
import os
import sys
import time
import shlex
import shutil
import argparse
import subprocess

HERE = os.path.dirname(os.path.abspath(__file__))
BIN = os.path.join(HERE, "bin")
sys.path.insert(0, os.path.join(HERE, "lib"))
from atacq_common import DEFAULT_LINEAGE_MAP  # noqa: E402

STEPS = ["map", "peaks", "features", "quant", "matrix", "qc"]


def now():
    return time.strftime("%Y-%m-%d %H:%M:%S")


def log(msg):
    print("[{}] {}".format(now(), msg), flush=True)


def run(cmd, name, logdir, done, env_extra=None, dry=False):
    if done and os.path.exists(done):
        log("[SKIP] {} (already done)".format(name))
        return
    log("[RUN ] {}".format(name))
    print("       $ " + " ".join(shlex.quote(str(c)) for c in cmd), flush=True)
    if dry:
        return
    env = dict(os.environ)
    env["PYTHONHASHSEED"] = "0"
    env.setdefault("PYTHONUNBUFFERED", "1")
    if env_extra:
        env.update(env_extra)
    os.makedirs(logdir, exist_ok=True)
    lp = os.path.join(logdir, name.replace(" ", "_").replace("/", "_") + ".log")
    t0 = time.time()
    with open(lp, "w") as lf:
        p = subprocess.Popen([str(c) for c in cmd], stdout=subprocess.PIPE,
                             stderr=subprocess.STDOUT, env=env, universal_newlines=True, bufsize=1)
        for line in p.stdout:
            sys.stdout.write("       " + line)
            lf.write(line)
        rc = p.wait()
    if rc != 0:
        log("[FAIL] {} exited with code {}; log: {}".format(name, rc, lp))
        sys.exit(rc)
    if done:
        with open(done, "w") as f:
            f.write("done {}\n".format(now()))
    log("[ OK ] {} ({:.1f} minutes)".format(name, (time.time() - t0) / 60))


def read_samples(path):
    out = []
    with open(path, encoding="utf-8") as f:
        header = None
        for line in f:
            if not line.strip() or line.startswith("#"):
                continue
            p = [x.strip() for x in line.rstrip("\n").split("\t")]
            if header is None:
                header = [x.lower() for x in p]
                if "sample" not in header:
                    sys.exit("[FATAL] samples.tsv needs a header: sample r1 r2 bam peaks")
                continue
            d = dict(zip(header, p + [""] * (len(header) - len(p))))
            for k in ("r1", "r2", "bam", "peaks"):
                v = d.get(k, "")
                d[k] = "" if v in ("", "-", "NA", "na", "None") else os.path.abspath(v)
            out.append(d)
    if not out:
        sys.exit("[FATAL] no samples in samples.tsv")
    names = [d["sample"] for d in out]
    if len(set(names)) != len(names):
        sys.exit("[FATAL] duplicate sample names")
    return out


def preflight(args, samples, steps):
    ok = True
    try:
        import numpy  # noqa
    except ImportError:
        print("[FATAL] numpy is missing")
        ok = False
    if os.environ.get("ATACQ_MOCK_BAM") != "1":
        try:
            import pysam  # noqa
        except ImportError:
            print("[FATAL] pysam is missing (conda install -c bioconda pysam)")
            ok = False
    need_map = "map" in steps and any(not s["bam"] for s in samples)
    need_peaks = "map" in steps and any(not s["peaks"] for s in samples)
    tools = []
    if need_map:
        tools += ["bowtie2", "samtools"]
    if need_peaks:
        tools += ["samtools", "alignmentSieve"]
    for t in sorted(set(tools)):
        if shutil.which(t) is None:
            print("[FATAL] cannot find {}".format(t))
            ok = False
    for p in (args.genome, args.gff, args.clusters, args.ancestry):
        if p and not os.path.exists(p):        # v4.3: ancestry may be empty (auto/diploid)
            print("[FATAL] file does not exist: {}".format(p))
            ok = False
    for s in samples:
        for k in ("r1", "r2", "bam", "peaks"):
            if s[k] and not os.path.exists(s[k]):
                print("[FATAL] {}: {} does not exist: {}".format(s["sample"], k, s[k]))
                ok = False
        if not s["bam"] and not (s["r1"] and s["r2"]):
            print("[FATAL] {}: has neither bam nor r1/r2".format(s["sample"]))
            ok = False
    if not ok:
        sys.exit(1)


def main():
    ap = argparse.ArgumentParser(description="ATAC-Quant v4 master", formatter_class=argparse.RawDescriptionHelpFormatter,
                                 epilog=__doc__)
    g = ap.add_argument_group("Inputs")
    g.add_argument("--samples", required=True)
    g.add_argument("--genome", required=True, help="unmasked genome FASTA (v2 no longer needs .masked)")
    g.add_argument("--gff", required=True)
    g.add_argument("--clusters", required=True, help="GraphAllele my_cluster.tsv")
    g.add_argument("--ancestry", default="",
                   help="gene->lineage table. Required when ploidy=allo; unused for auto/diploid (no lineage layer is produced)")
    # ---- v4.3: ploidy mode (aligned with --ploidy in the sister project LQ; diploid is an extension of this pipeline)
    g.add_argument("--ploidy", choices=["diploid", "auto", "allo"], default="allo",
                   help="allo=allopolyploid (lineage layer, needs ancestry); auto=autopolyploid (no lineage layer); "
                        "diploid=diploid two-haplotype reference (no lineage layer, 2 copies per cluster)")
    g.add_argument("--bt2-index", default=None, help="bowtie2 index prefix (default = genome with its extension stripped)")
    g.add_argument("--lineage-map", default="")   # v3: empty by default; the species mapping arrives via --lineage-map/config
    g.add_argument("--lineage-kmers", nargs="*", default=None, metavar="NAME:PATH",
                   help="optional: lineage k-mers (off by default, rationale in the README)")
    o = ap.add_argument_group("Run")
    o.add_argument("-o", "--outdir", required=True)
    o.add_argument("-t", "--threads", type=int, default=16)
    o.add_argument("--steps", default="all", help="comma-separated: " + ",".join(STEPS))
    o.add_argument("--force", action="store_true", help="ignore checkpoints and rerun the selected steps")
    o.add_argument("--dry-run", action="store_true")
    o.add_argument("--tmpdir", default=None, help="temporary directory for 03 (needs a fair amount of space)")
    m = ap.add_argument_group("00 mapping/peak calling")
    m.add_argument("--bt2-k", type=int, default=20)
    m.add_argument("--maxins", type=int, default=1000)
    m.add_argument("--gsize", default="1.0e10")
    m.add_argument("--qvalue", default="0.05")
    m.add_argument("--macs2-env", default="", help="name of the conda environment to activate before peak calling (v1 used ATAC)")
    m.add_argument("--macs2", default="macs2")
    # ---- v4.2: assay switch (controls the Tn5 offset correction in 00 and the MACS2 peak-calling details) ----
    m.add_argument("--assay", choices=["atac", "chip"], default="atac",
                   help="atac = Tn5 transposase (+4/-5 offset correction before peak calling); "
                        "chip = antibody enrichment (that correction is skipped)")
    m.add_argument("--macs2-extra", default="",
                   help="extra arguments appended verbatim to MACS2 callpeak (e.g. \"--broad\" for broad peaks)")
    q = ap.add_argument_group("03 quantification")
    q.add_argument("--outside", choices=["model", "ignore"], default="model")
    q.add_argument("--as-temp", type=float, default=2.0)
    q.add_argument("--bootstrap", type=int, default=0)
    q.add_argument("--quant-extra", default="", help="extra arguments passed verbatim to 03_quantify.py")
    q.add_argument("--feature-extra", default="", help="extra arguments passed verbatim to 02_build_features.py")
    # ---- new in v4: EM convergence acceleration and bootstrap parallelism (forwarded to 03_quantify.py)
    q.add_argument("--em-accelerate", choices=["squarem", "off"], default="squarem",
                   help="EM acceleration: squarem = SQUAREM squared extrapolation (v4 default); off = same as v3")
    q.add_argument("--em-trace", action="store_true", help="print the log-likelihood of every EM iteration")
    q.add_argument("--em-prior", default="1e-8", help="per-iteration zero-avoiding prior")
    q.add_argument("--em-max-iter", type=int, default=2000, help="main EM iteration cap")
    q.add_argument("--em-tol", default="1e-6", help="main EM convergence tolerance (v4 default 1e-6; v3 used 1e-3)")
    q.add_argument("--bootstrap-threads", type=int, default=0,
                   help="number of parallel bootstrap processes; 0 = auto; 1 = sequential")
    q.add_argument("--bootstrap-max-iter", type=int, default=500,
                   help="EM iteration cap for each bootstrap replicate")
    q.add_argument("--min-samples", type=int, default=1, help="minimum number of supporting samples when merging peaks in 01")
    args = ap.parse_args()

    steps = STEPS if args.steps == "all" else [s.strip() for s in args.steps.split(",")]
    for s in steps:
        if s not in STEPS:
            sys.exit("[FATAL] unknown step {}".format(s))
    for k in ("genome", "gff", "clusters", "ancestry", "samples"):
        # v4.3 fix: os.path.abspath("") returns the CURRENT WORKING DIRECTORY, turning "not provided" into
        # "a directory path", so the downstream `if args.ancestry:` wrongly treated it as provided
        # (v3 never hit this because ancestry was mandatory).
        _v = getattr(args, k)
        if _v:
            setattr(args, k, os.path.abspath(_v))
    out = os.path.abspath(args.outdir)
    ck = os.path.join(out, ".checkpoints")
    logs = os.path.join(out, "logs")
    os.makedirs(ck, exist_ok=True)
    if args.force:
        for s in steps:
            for fn in os.listdir(ck):
                if fn.startswith(s):
                    os.remove(os.path.join(ck, fn))
    samples = read_samples(args.samples)
    preflight(args, samples, steps)
    py = sys.executable
    bt2 = args.bt2_index or os.path.splitext(args.genome)[0]
    log("ATAC-Quant v4 | {} sample(s) | output {}".format(len(samples), out))

    # ------------------------------------------------------------ 00
    for s in samples:
        sd = os.path.join(out, "samples", s["sample"])
        if not s["bam"]:
            s["bam"] = os.path.join(sd, s["sample"] + ".marked.bam")
            need = True
        else:
            need = False
        if not s["peaks"]:
            # v4.2: MACS2 --broad writes _peaks.broadPeak; accept both (read_peaks only takes the first 4 columns)
            for _ext in ("narrowPeak", "broadPeak"):
                _c = os.path.join(sd, s["sample"] + "_peaks." + _ext)
                if os.path.exists(_c):
                    s["peaks"] = _c
                    break
            if not s["peaks"]:
                s["peaks"] = os.path.join(sd, s["sample"] + "_peaks.narrowPeak")
            need = True
        if "map" in steps and need:
            cmd = ["bash", os.path.join(BIN, "00_map_call.sh"), "--outdir", sd,
                   "--prefix", s["sample"], "--threads", args.threads,
                   "--bt2-k", args.bt2_k, "--maxins", args.maxins, "--gsize", args.gsize,
                   "--qvalue", args.qvalue, "--macs2", args.macs2,
                   "--assay", args.assay]          # v4.2: atac(Tn5 offset correction) / chip(no correction)
            if args.macs2_extra:
                cmd += ["--macs2-extra", args.macs2_extra]
            if args.macs2_env:
                cmd += ["--macs2-env", args.macs2_env]
            if s["bam"].startswith(sd) and not os.path.exists(s["bam"] + ".bai"):
                cmd += ["--r1", s["r1"], "--r2", s["r2"], "--index", bt2]
            else:
                cmd += ["--bam", s["bam"]]
            run(cmd, "00_map_{}".format(s["sample"]), logs,
                os.path.join(ck, "map_{}.done".format(s["sample"])), dry=args.dry_run)

    # ------------------------------------------------------------ 01
    pk_dir = os.path.join(out, "peaks")
    cons = os.path.join(pk_dir, "consensus_peaks.bed")
    if "peaks" in steps:
        run([py, os.path.join(BIN, "01_consensus_peaks.py"), "-o", pk_dir,
             "--min-samples", args.min_samples, "--peaks"] +
            ["{}:{}".format(s["sample"], s["peaks"]) for s in samples],
            "01_consensus_peaks", logs, os.path.join(ck, "peaks.done"), dry=args.dry_run)

    # ------------------------------------------------------------ 02
    ft_dir = os.path.join(out, "features")
    feats = os.path.join(ft_dir, "features.tsv")
    if "features" in steps:
        c02 = [py, os.path.join(BIN, "02_build_features.py"), "--peaks", cons, "--gff", args.gff,
               "--clusters", args.clusters, "--genome", args.genome,
               "--lineage-map", args.lineage_map, "--ploidy", args.ploidy,
               "-o", ft_dir, "-t", args.threads] + shlex.split(args.feature_extra)
        if args.ancestry:                    # v4.3: ancestry is not passed for auto/diploid
            c02 += ["--ancestry", args.ancestry]
        run(c02, "02_build_features", logs, os.path.join(ck, "features.done"), dry=args.dry_run)

    # ------------------------------------------------------------ 03
    qfiles = []
    for s in samples:
        qd = os.path.join(out, "quant", s["sample"])
        qfiles.append((s["sample"], os.path.join(qd, s["sample"] + ".quant.tsv")))
        if "quant" in steps:
            cmd = [py, os.path.join(BIN, "03_quantify.py"), "--bam", s["bam"], "--features", feats,
                   "--sample", s["sample"], "-o", qd, "-t", args.threads,
                   "--outside", args.outside, "--as-temp", args.as_temp,
                   "--bootstrap", args.bootstrap,
                   # ---- v4: EM acceleration + bootstrap parallelism
                   "--em-accelerate", args.em_accelerate,
                   "--em-prior", args.em_prior,
                   "--max-iter", str(args.em_max_iter),
                   "--tol", args.em_tol,
                   "--bootstrap-threads", str(args.bootstrap_threads),
                   "--bootstrap-max-iter", str(args.bootstrap_max_iter)]
            if args.em_trace:
                cmd += ["--em-trace"]
            if args.tmpdir:
                cmd += ["--tmpdir", os.path.abspath(args.tmpdir)]
            if args.lineage_kmers:
                cmd += ["--lineage-kmers"] + args.lineage_kmers
            cmd += shlex.split(args.quant_extra)
            run(cmd, "03_quant_{}".format(s["sample"]), logs,
                os.path.join(ck, "quant_{}.done".format(s["sample"])), dry=args.dry_run)

    # ------------------------------------------------------------ 04 / 05
    mx = os.path.join(out, "matrix")
    if "matrix" in steps:
        run([py, os.path.join(BIN, "04_build_matrix.py"), "--features", feats, "-o", mx,
             "--ploidy", args.ploidy, "--quant"] +
            ["{}:{}".format(a, b) for a, b in qfiles],
            "04_matrix", logs, os.path.join(ck, "matrix.done"), dry=args.dry_run)
    if "qc" in steps:
        run([py, os.path.join(BIN, "05_qc_report.py"), "--features", feats,
             "--features-stats", os.path.join(ft_dir, "features_stats.json"),
             "--ploidy", args.ploidy,
             "-o", os.path.join(out, "qc"), "--quant"] + ["{}:{}".format(a, b) for a, b in qfiles],
            "05_qc", logs, os.path.join(ck, "qc.done"), dry=args.dry_run)
    log("all done. Start with {}".format(os.path.join(out, "qc", "qc_report.md")))


if __name__ == "__main__":
    main()
