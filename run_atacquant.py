#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
run_atacquant.py - ATAC-Quant v4 unified entry point (config-file driven, works for any species)

Usage:
  python run_atacquant.py --config config.ini            read config -> validate -> run the whole pipeline
  python run_atacquant.py --config config.ini --check    validate and print the execution plan only, no run
  python run_atacquant.py --config config.ini --dry-run  print every command but do not run it
  python run_atacquant.py --config config.ini --force    ignore checkpoints, rerun the selected steps
  python run_atacquant.py --config config.ini --steps features,quant   override the steps in the config
  python run_atacquant.py --config config.ini --no-stats skip the statistics layer at the end

All inputs are defined in the config file (template: config.template.ini). This script does only three things:
  1) parse and validate the config (missing files and missing tools are reported before anything runs)
  2) assemble and call the engine atacquant_core.py (same origin as v2)
  3) once the main pipeline has finished, (optionally) run the statistics layer 10_copy_level_stats.py per sample

Path rules: "relative paths" in the config file and in the sample sheet always resolve against
"the directory containing the config file".
"""
import argparse
import configparser
import os
import shlex
import shutil
import subprocess
import sys
import time

HERE = os.path.dirname(os.path.abspath(__file__))
CORE = os.path.join(HERE, "atacquant_core.py")
STATS = os.path.join(HERE, "analysis", "10_copy_level_stats.py")

ALL_STEPS = ["map", "peaks", "features", "quant", "matrix", "qc"]


def now():
    return time.strftime("%Y-%m-%d %H:%M:%S")


def log(msg):
    print("[{}] {}".format(now(), msg), flush=True)


def die(msg):
    print("[{}] [FATAL] {}".format(now(), msg), flush=True)
    sys.exit(2)


def gopt(cp, sec, key, default=None):
    if cp.has_section(sec):
        for k, v in cp.items(sec):     # Key names are case-insensitive (needed once optionxform=str preserves sample-name case)
            if k.lower() == key.lower():
                v = v.strip()
                return v if v != "" else default
    return default


def gint(cp, sec, key, default):
    v = gopt(cp, sec, key)
    if v is None:
        return default
    try:
        return int(v)
    except ValueError:
        die("[{}.{}] expects an integer, got: {}".format(sec, key, v))


def resolve(base, p):
    if p is None or p.strip() in ("", "-"):
        return ""
    p = os.path.expanduser(p.strip())
    if os.path.isabs(p):
        return os.path.normpath(p)
    if p.startswith("/"):
        # Compatibility: an absolute path written for Linux ("/...") is not absolute under
        # Windows isabs; joining it directly would produce a drive path (e.g. D:\data\genome.fa) -> keep it as is
        return p
    return os.path.normpath(os.path.join(base, p))


def fmt_cmd(cmd):
    return " ".join(shlex.quote(str(c)) for c in cmd)


# Canonical spelling of a missing column (inline paths and external sample-sheet paths must agree, see the v4.2 fix)
MISSING = ("", "-", "NA", "na", "None", "none")


def load_samples_inline(cp, base):
    """[samples] inline: key = sample name, value = 'r1 r2 bam peaks' (write - for anything missing)."""
    rows = []
    if cp.has_section("samples"):
        for name, val in cp.items("samples"):
            val = val.strip()
            if not val:
                continue
            f = val.split()
            if len(f) != 4:
                die("config [samples] entry '{}' should have 4 fields (r1 r2 bam peaks), got {}: {}".format(
                    name, len(f), val))
            # v4.2 fix: `-` must be normalised to an empty value. The old code kept "-" verbatim, so the
            # "FASTQ only, write - for bam/peaks" usage - which the docs explicitly recommend - was blocked
            # by validation ("BAM not found: -"). The external sample-sheet path already did this; now they agree.
            rows.append([name.strip()] +
                        ["" if x in MISSING else resolve(base, x) for x in f])
    return rows


def write_sample_sheet(rows, path):
    with open(path, "w", encoding="utf-8", newline="\n") as f:
        f.write("sample\tr1\tr2\tbam\tpeaks\n")
        for r in rows:
            f.write("\t".join(r) + "\n")


def read_sheet(path):
    rows = []
    header = None
    with open(path, encoding="utf-8") as f:
        for line in f:
            if not line.strip() or line.startswith("#"):
                continue
            p = [x.strip() for x in line.rstrip("\n").split("\t")]
            if header is None:
                header = [x.lower() for x in p]
                if "sample" not in header:
                    die("sample sheet is missing its header (need: sample r1 r2 bam peaks): " + path)
                continue
            d = dict(zip(header, p + [""] * (len(header) - len(p))))
            for k in ("r1", "r2", "bam", "peaks"):
                v = d.get(k, "")
                d[k] = "" if v in ("", "-", "NA", "na", "None") else v
            rows.append(d)
    return rows


def run_stream(cmd, log_file, env):
    log("$ " + fmt_cmd(cmd))
    t0 = time.time()
    with open(log_file, "a", encoding="utf-8") as lf:
        lf.write("\n[{}]\n$ {}\n".format(now(), fmt_cmd(cmd)))
        lf.flush()
        p = subprocess.Popen([str(c) for c in cmd], stdout=subprocess.PIPE,
                             stderr=subprocess.STDOUT, env=env,
                             universal_newlines=True, encoding="utf-8",
                             errors="replace", bufsize=1)
        for line in p.stdout:
            sys.stdout.write(line)
            sys.stdout.flush()
            lf.write(line)
        rc = p.wait()
    log("[{}] took {:.1f} minutes".format("OK" if rc == 0 else "FAIL rc={}".format(rc),
                                          (time.time() - t0) / 60.0))
    return rc


def main():
    ap = argparse.ArgumentParser(description="ATAC-Quant v4 config-driven entry point")
    ap.add_argument("--config", required=True, help="config file (a copy of config.template.ini)")
    ap.add_argument("--check", action="store_true", help="validate and print the execution plan only")
    ap.add_argument("--dry-run", action="store_true", help="print the commands but do not run them")
    ap.add_argument("--force", action="store_true", help="ignore checkpoints and rerun")
    ap.add_argument("--steps", default=None, help="override the steps in the config")
    ap.add_argument("--stats", dest="stats", action="store_true", default=None, help="force the statistics layer to run")
    ap.add_argument("--no-stats", dest="stats", action="store_false", help="skip the statistics layer")
    args = ap.parse_args()

    cfg_path = os.path.abspath(args.config)
    if not os.path.isfile(cfg_path):
        die("cannot find the config file: " + cfg_path)
    base = os.path.dirname(cfg_path)

    cp = configparser.ConfigParser(interpolation=None,
                                   inline_comment_prefixes=("#", ";"))   # allow end-of-line comments
    cp.optionxform = str    # Keep key case (sample names such as Sample1 must not be lower-cased by the configparser default)
    try:
        cp.read(cfg_path, encoding="utf-8")
    except configparser.Error as e:
        die("failed to parse the config: {}".format(e))

    name = gopt(cp, "project", "name", "ATACQuant")
    outdir = resolve(base, gopt(cp, "project", "outdir", "./{}_out".format(name)))
    threads = gint(cp, "project", "threads", 16)
    tmpdir = resolve(base, gopt(cp, "project", "tmpdir", ""))

    genome = resolve(base, gopt(cp, "inputs", "genome", ""))
    gff = resolve(base, gopt(cp, "inputs", "gff", ""))
    clusters = resolve(base, gopt(cp, "inputs", "clusters", ""))
    ancestry = resolve(base, gopt(cp, "inputs", "ancestry", ""))
    # ---- v4.3: ploidy mode (aligned with --ploidy in the sister project LQ; diploid is an ATAC-Quant extension)
    #   allo    = allopolyploid (subgenomes are differentiated) -> ancestry required, produces a lineage layer
    #   auto    = autopolyploid (subgenomes not differentiated, so no "lineage" concept) -> no ancestry, lineage layer skipped
    #   diploid = diploid (two-haplotype reference, 2 copies per cluster) -> no ancestry, lineage layer skipped
    ploidy = gopt(cp, "inputs", "ploidy", "allo").strip().lower()
    if ploidy not in ("diploid", "auto", "allo"):
        die("[inputs] ploidy must be diploid / auto / allo, got: {}".format(ploidy))
    if not (genome and gff and clusters):
        die("config [inputs] requires: genome / gff / clusters (see config.template.ini)")
    if ploidy == "allo" and not ancestry:
        die("[inputs] ploidy=allo (allopolyploid) requires ancestry (a lineage table).\n"
            "       If the subgenomes are not differentiated (autopolyploid), set ploidy = auto;\n"
            "       if this is a diploid two-haplotype reference, set ploidy = diploid.")
    if ploidy != "allo" and ancestry:
        warns_early = ("[inputs] with ploidy={} the ancestry table is IGNORED (no lineage layer is produced)".format(ploidy))
    else:
        warns_early = None
    bt2 = resolve(base, gopt(cp, "inputs", "bt2_index", ""))
    sheet_opt = resolve(base, gopt(cp, "inputs", "samples", ""))

    inline = load_samples_inline(cp, base)
    sheet_generated = None
    if inline:
        sheet = os.path.join(outdir, "samples.tsv")
        sheet_generated = sheet
        samples = [{"sample": r[0], "r1": r[1], "r2": r[2], "bam": r[3], "peaks": r[4]}
                   for r in inline]
    elif sheet_opt:
        if not os.path.isfile(sheet_opt):
            die("cannot find the sample sheet: " + sheet_opt)
        sheet = sheet_opt
        samples = read_sheet(sheet)
    else:
        die("either point [inputs] samples at a sample sheet, or list the samples inline under [samples]")

    steps = args.steps or gopt(cp, "run", "steps", "all")
    steps_list = ALL_STEPS[:] if steps == "all" else [s.strip() for s in steps.split(",") if s.strip()]
    for s in steps_list:
        if s not in ALL_STEPS:
            die("unknown step: {} (available: {})".format(s, ",".join(ALL_STEPS)))

    bt2_k = gopt(cp, "map", "bt2_k", "20")
    maxins = gopt(cp, "map", "maxins", "1000")
    gsize = gopt(cp, "map", "gsize", "1.0e10")
    qvalue = gopt(cp, "map", "qvalue", "0.05")
    macs2_env = gopt(cp, "map", "macs2_env", "")
    macs2 = gopt(cp, "map", "macs2", "macs2")
    # ---- v4.2: assay switch (atac = Tn5 transposase / chip = antibody enrichment)
    #      controls the Tn5 (+4/-5) offset correction in step 00 and the extra arguments passed on to MACS2
    assay = gopt(cp, "map", "assay", "atac").strip().lower()
    macs2_extra = gopt(cp, "map", "macs2_extra", "")
    if assay not in ("atac", "chip"):
        die("[map] assay must be atac or chip, got: {}".format(assay))

    as_temp = gopt(cp, "quant", "as_temp", "2.0")
    bootstrap = gint(cp, "quant", "bootstrap", 0)
    outside = gopt(cp, "quant", "outside", "model")
    lineage_map = gopt(cp, "quant", "lineage_map", "")
    lineage_kmers = gopt(cp, "quant", "lineage_kmers", "")
    # ---- new in v4: EM convergence acceleration and bootstrap parallelism
    em_accelerate = gopt(cp, "quant", "em_accelerate", "squarem")   # squarem / off
    em_trace = gopt(cp, "quant", "em_trace", "false").lower() in ("1", "true", "yes", "on")
    em_prior = gopt(cp, "quant", "em_prior", "1e-8")
    em_max_iter = gint(cp, "quant", "em_max_iter", 2000)
    em_tol = gopt(cp, "quant", "em_tol", "1e-6")
    bootstrap_threads = gint(cp, "quant", "bootstrap_threads", 0)   # 0 = auto
    bootstrap_max_iter = gint(cp, "quant", "bootstrap_max_iter", 500)
    if em_accelerate not in ("squarem", "off"):
        die("[quant] em_accelerate must be squarem or off, got: {}".format(em_accelerate))
    kmers_list = []
    if lineage_kmers:
        for item in lineage_kmers.split(","):
            item = item.strip()
            if not item:
                continue
            if ":" not in item:
                die("lineage_kmers should be formatted as NAME:PATH, comma-separated")
            nm, pp = item.split(":", 1)
            kmers_list.append(nm.strip() + ":" + resolve(base, pp))

    min_samples = gopt(cp, "run", "min_samples", "1")
    quant_extra = gopt(cp, "quant", "quant_extra", "")
    feature_extra = gopt(cp, "run", "feature_extra", "")

    stats_run = gopt(cp, "stats", "run", "true").lower() in ("1", "true", "yes", "on")
    if args.stats is not None:
        stats_run = args.stats
    mc = gint(cp, "stats", "mc", 2000)
    seed = gint(cp, "stats", "seed", 2026)

    # ------------------------------------------------------------ validation
    errors, warns = [], []
    if warns_early:
        warns.append(warns_early)

    def chk_file(p, label):
        if not p or not os.path.isfile(p):
            errors.append("missing {}: {}".format(label, p or "(not set)"))

    chk_file(genome, "genome FASTA")
    chk_file(gff, "GFF3 annotation")
    chk_file(clusters, "homeologous cluster table")
    if ploidy == "allo":
        chk_file(ancestry, "lineage table")
    if not samples:
        errors.append("the sample sheet contains no samples")

    for s in samples:
        nm = s["sample"]
        bam, peaks = s.get("bam", ""), s.get("peaks", "")
        r1, r2 = s.get("r1", ""), s.get("r2", "")
        if bam:
            if not os.path.isfile(bam):
                errors.append("[{}] BAM does not exist: {}".format(nm, bam))
        else:
            if not (r1 and r2):
                errors.append("[{}] has neither a BAM nor R1/R2".format(nm))
            for lab, p in (("R1", r1), ("R2", r2)):
                if p and not os.path.isfile(p):
                    errors.append("[{}] {} does not exist: {}".format(nm, lab, p))
            if "map" not in steps_list:
                errors.append("[{}] only FASTQ given -> steps must include map".format(nm))
        if peaks:
            if not os.path.isfile(peaks):
                errors.append("[{}] peak file does not exist: {}".format(nm, peaks))
        elif "map" not in steps_list:
            errors.append("[{}] peak file missing and steps does not include map (cannot build it from the BAM)".format(nm))

    # Mapping / peak-calling checks: the corresponding tools and index are only required when some sample really needs them
    need_align = "map" in steps_list and any(not s.get("bam") for s in samples)
    need_peakgen = "map" in steps_list and any((not s.get("bam")) or (not s.get("peaks")) for s in samples)
    if need_peakgen:
        # v4.2: alignmentSieve only serves ATAC Tn5 (+4/-5) offset correction;
        #       assay=chip skips that step, so deeptools must not be required
        peak_tools = ["samtools"] + (["alignmentSieve"] if assay == "atac" else [])
        for t in peak_tools:
            if shutil.which(t) is None:
                errors.append("a sample needs mapping/peak calling (map), but the tool is not on the PATH: " + t)
        if not macs2_env and shutil.which(macs2) is None:
            warns.append("cannot find {} on the PATH (use [map] macs2_env to name a conda environment that has MACS2)".format(macs2))
    if need_align:
        if shutil.which("bowtie2") is None:
            errors.append("a sample needs mapping (map), but the tool is not on the PATH: bowtie2")
        idx = bt2 or os.path.splitext(genome)[0]
        if not (os.path.isfile(idx + ".1.bt2") or os.path.isfile(idx + ".1.bt2l")):
            errors.append("cannot find the bowtie2 index: {} (run bowtie2-build first, or set [inputs] bt2_index to its prefix)".format(idx))

    # ------------------------------------------------------------ plan and check output
    print("")
    log("ATAC-Quant v4 | project: {} | config: {}".format(name, cfg_path))
    print("  Genome       : {}".format(genome))
    print("  Annotation   : {}".format(gff))
    print("  Cluster table: {}".format(clusters))
    print("  Ancestry     : {}".format(ancestry if ancestry else "(unused; ploidy={} produces no lineage layer)".format(ploidy)))
    _pm = {"allo": "allo (allopolyploid -> lineage layer present: split lineages first, then copies within a lineage)",
           "auto": "auto (autopolyploid -> no lineage layer: split copies directly; copy-level resolution is harder)",
           "diploid": "diploid (diploid, two-haplotype reference -> no lineage layer; 2 copies per cluster)"}[ploidy]
    print("  Ploidy mode  : {}".format(_pm))
    print("  Output       : {}".format(outdir))
    print("  Steps        : {}".format(",".join(steps_list)))
    if "map" in steps_list:
        _a = {"atac": "atac (T5/Tn5 transposase -> +4/-5 offset correction before peak calling)",
              "chip": "chip (antibody enrichment -> that correction is skipped)"}[assay]
        _mx = "; extra MACS2 arguments: {}".format(macs2_extra) if macs2_extra else ""
        print("  Assay        : {}{}".format(_a, _mx))
    n_map = sum(1 for s in samples if not s.get("bam"))
    print("  Samples ({})  : {}".format(
        len(samples),
        "; ".join("{}{}".format(s["sample"], " [FASTQ->auto-align]" if not s.get("bam") else " [BAM]")
                  for s in samples[:8]) + (" ..." if len(samples) > 8 else "")))
    for w in warns:
        print("  [WARN] " + w)
    if errors:
        for e in errors:
            print("  [ERROR] " + e)
        die("validation failed ({} problem(s)); fix the config and try again".format(len(errors)))

    # ------------------------------------------------------------ build the command
    cmd = [sys.executable, CORE, "--samples", sheet, "--genome", genome, "--gff", gff,
           "--clusters", clusters, "-o", outdir,
           "-t", str(threads), "--steps", ",".join(steps_list),
           "--ploidy", ploidy]                     # v4.3: diploid / auto / allo
    if ancestry:
        cmd += ["--ancestry", ancestry]
    if tmpdir:
        cmd += ["--tmpdir", tmpdir]
    if bt2:
        cmd += ["--bt2-index", bt2]
    if "map" in steps_list:
        cmd += ["--bt2-k", bt2_k, "--maxins", maxins, "--gsize", gsize,
                "--qvalue", qvalue, "--macs2", macs2,
                "--assay", assay]                   # v4.2: atac(Tn5) / chip(antibody enrichment)
        if macs2_extra:
            cmd += ["--macs2-extra=" + macs2_extra]   # use the = form, see the comment below
    if macs2_env:
        cmd += ["--macs2-env", macs2_env]
    cmd += ["--as-temp", as_temp, "--bootstrap", str(bootstrap), "--outside", outside]
    # ---- v4: forward the EM acceleration and bootstrap parallelism settings
    cmd += ["--em-accelerate", em_accelerate, "--em-prior", em_prior,
            "--em-max-iter", str(em_max_iter), "--em-tol", em_tol,
            "--bootstrap-threads", str(bootstrap_threads),
            "--bootstrap-max-iter", str(bootstrap_max_iter)]
    if em_trace:
        cmd += ["--em-trace"]
    if lineage_map:
        cmd += ["--lineage-map", lineage_map]
    if kmers_list:
        cmd += ["--lineage-kmers"] + kmers_list
    cmd += ["--min-samples", min_samples]
    # v4.2 fix: pass the "extra arguments" string as "--opt=value".
    #   The old form ["--quant-extra", "--buckets 8"] only worked when the value contained a space, so argparse
    #   treated it as a positional; a single-token value (e.g. macs2_extra = "--broad") was mistaken for an
    #   option -> "expected one argument".
    if quant_extra:
        cmd += ["--quant-extra=" + quant_extra]
    if feature_extra:
        cmd += ["--feature-extra=" + feature_extra]
    if args.force:
        cmd += ["--force"]
    if args.dry_run:
        cmd += ["--dry-run"]

    if stats_run and "quant" in steps_list:
        print("  Statistics   : on (per sample -> {}/stats/<sample>/, mc={}, seed={})".format(outdir, mc, seed))
    else:
        print("  Statistics   : off")
    print("  Plan command : {}".format(fmt_cmd(cmd)))

    if args.check:
        log("check-only mode (--check): the lines above are the execution plan; nothing was run")
        return

    # ------------------------------------------------------------ run
    os.makedirs(outdir, exist_ok=True)
    if sheet_generated:
        write_sample_sheet(inline, sheet_generated)
        log("generated the sample sheet from config [samples]: " + sheet_generated)
    log_file = os.path.join(outdir, "run_atacquant.log")

    os.chdir(base)   # so relative paths in the sample sheet also resolve against "the directory containing the config file"
    env = os.environ.copy()
    env.setdefault("PYTHONHASHSEED", "0")
    env.setdefault("PYTHONUNBUFFERED", "1")
    env.setdefault("PYTHONUTF8", "1")            # uniform UTF-8: reliable across platforms and locales
    env.setdefault("PYTHONIOENCODING", "utf-8")

    rc = run_stream(cmd, log_file, env)
    if rc != 0:
        die("the main pipeline failed (exit code {}). Finished steps have checkpoints, so after fixing the problem just rerun the same command to resume. Log: {}".format(rc, log_file))

    # ------------------------------------------------------------ statistics layer (optional)
    if stats_run and not args.dry_run and "quant" in steps_list:
        feat = os.path.join(outdir, "features", "features.tsv")
        grp = os.path.join(outdir, "features", "groups.tsv")
        n_fail = 0
        for s in samples:
            nm = s["sample"]
            q = os.path.join(outdir, "quant", nm, nm + ".quant.tsv")
            if not os.path.isfile(q):
                log("[SKIP] statistics {}: cannot find {}".format(nm, q))
                continue
            scmd = [sys.executable, STATS, "--quant", q, "--features", feat, "--groups", grp,
                    "--outdir", os.path.join(outdir, "stats", nm), "--mc", str(mc), "--seed", str(seed)]
            rc2 = run_stream(scmd, log_file, env)
            if rc2 != 0:
                n_fail += 1
                log("[WARN] statistics failed for {} (rc={}); continuing with the remaining samples".format(nm, rc2))
        if n_fail:
            die("the main pipeline is done, but the statistics layer failed for {} sample(s) (numpy+pandas required)".format(n_fail))

    log("all done. Start with {}/qc/qc_report.md".format(outdir))


if __name__ == "__main__":
    main()
