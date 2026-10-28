# -*- coding: utf-8 -*-
"""
atacq_common.py -- ATAC-Quant v4 shared utilities

Contains:
  - logging / timing
  - deterministic hashing (PYTHONHASHSEED=0 guard; read-name hashing that is consistent across processes)
  - FASTA random access (.fai; when it is missing, build the index with samtools faidx or in pure Python)
  - parsers for GFF3 genes, allele clusters and the lineage table (lineage names normalized automatically, Robustum->SRO etc.)
  - BAM opening (pysam; can switch to mockbam for testing)

Compatible with Python >= 3.7.
"""
import os
import sys
import time
import subprocess
from collections import defaultdict

# --------------------------------------------------------------------------
# Logging
# --------------------------------------------------------------------------
def now():
    return time.strftime("%Y-%m-%d %H:%M:%S")


def log(msg):
    print("[{}] {}".format(now(), msg), flush=True)


def die(msg, code=1):
    print("[{}] [FATAL] {}".format(now(), msg), file=sys.stderr, flush=True)
    sys.exit(code)


# --------------------------------------------------------------------------
# Deterministic hashing
# --------------------------------------------------------------------------
MASK64 = (1 << 64) - 1


def ensure_hashseed():
    """
    Read names must hash identically across processes and across stages. Python randomizes str
    hashing by default, so this guarantees PYTHONHASHSEED=0; if the current interpreter was not
    started with it, re-exec itself with the same arguments.
    """
    if os.environ.get("PYTHONHASHSEED") != "0":
        env = dict(os.environ)
        env["PYTHONHASHSEED"] = "0"
        os.execve(sys.executable, [sys.executable] + sys.argv, env)


def qhash(name):
    """read name -> 64-bit unsigned integer (requires ensure_hashseed first)."""
    return hash(name) & MASK64


# --------------------------------------------------------------------------
# Sequence utilities
# --------------------------------------------------------------------------
_RC = str.maketrans("ACGTNacgtn", "TGCANtgcan")


def revcomp(s):
    return s.translate(_RC)[::-1]


class Fasta(object):
    """FASTA random access on top of .fai; returns uppercase sequence."""

    def __init__(self, path):
        self.path = path
        fai = path + ".fai"
        if not os.path.exists(fai):
            self._build_fai(fai)
        self.index = {}
        self.order = []
        with open(fai) as f:
            for line in f:
                p = line.rstrip("\n").split("\t")
                if len(p) < 5:
                    continue
                name, ln, off, lb, lw = p[0], int(p[1]), int(p[2]), int(p[3]), int(p[4])
                self.index[name] = (ln, off, lb, lw)
                self.order.append(name)
        self._fh = None
        self._verify_index()

    def _verify_index(self):
        """Check that the .fai and the actual byte layout of the FASTA agree (Windows CRLF /
        newline-difference guard). On mismatch rebuild the index from the real bytes -- otherwise
        fetching sequence by offset drifts cumulatively with the coordinates."""
        if not self.order:
            return
        c0 = self.order[0]
        ln, off, lb, lw = self.index[c0]
        if lb <= 0 or lw <= 0:
            return
        try:
            with open(self.path, "rb") as fh:
                fh.seek(off)
                probe = fh.read(lb + 4)
        except OSError:
            return
        bad = len(probe) < lb or any(b in (13, 10) for b in probe[:lb])
        if not bad:
            w = 0
            while lb + w < len(probe) and probe[lb + w] in (13, 10):   # \r=13 \n=10
                w += 1
            if w > 0 and (lw - lb) != w:
                bad = True
        if bad:
            log("[WARN] genome FASTA and .fai disagree (line width / offset; typical of Windows CRLF) -- index rebuilt: "
                + self.path + ".fai")
            self._build_fai(self.path + ".fai")
            self.index = {}
            self.order = []
            with open(self.path + ".fai") as f:
                for line in f:
                    p = line.rstrip("\r\n").split("\t")
                    if len(p) < 5:
                        continue
                    self.index[p[0]] = (int(p[1]), int(p[2]), int(p[3]), int(p[4]))
                    self.order.append(p[0])

    def _build_fai(self, fai):
        try:
            subprocess.run(["samtools", "faidx", self.path], check=True,
                           stdout=subprocess.DEVNULL, stderr=subprocess.DEVNULL)
            if os.path.exists(fai):
                return
        except Exception:
            pass
        log("[INFO] samtools not found, building .fai for {} in Python (slow on large genomes)".format(self.path))
        entries = []
        with open(self.path, "rb") as f:
            name = None
            seqlen = off = lb = lw = 0
            pos = 0
            for raw in f:
                if raw.startswith(b">"):
                    if name is not None:
                        entries.append((name, seqlen, off, lb, lw))
                    name = raw[1:].split()[0].decode()
                    seqlen = 0
                    off = pos + len(raw)
                    lb = lw = 0
                else:
                    if lb == 0:
                        lb = len(raw.rstrip(b"\r\n"))
                        lw = len(raw)
                    seqlen += len(raw.rstrip(b"\r\n"))
                pos += len(raw)
            if name is not None:
                entries.append((name, seqlen, off, lb, lw))
        with open(fai, "w") as fo:
            for e in entries:
                fo.write("{}\t{}\t{}\t{}\t{}\n".format(*e))

    def length(self, chrom):
        return self.index[chrom][0]

    def lengths(self):
        return {c: self.index[c][0] for c in self.order}

    def fetch(self, chrom, start, end):
        """0-based half-open interval [start, end). Out-of-range coordinates are truncated."""
        if chrom not in self.index:
            return ""
        ln, off, lb, lw = self.index[chrom]
        start = max(0, start)
        end = min(ln, end)
        if end <= start:
            return ""
        if self._fh is None:
            self._fh = open(self.path, "rb")
        b0 = off + (start // lb) * lw + (start % lb)
        b1 = off + ((end - 1) // lb) * lw + ((end - 1) % lb) + 1
        self._fh.seek(b0)
        raw = self._fh.read(b1 - b0)
        return raw.replace(b"\n", b"").replace(b"\r", b"").decode().upper()


# --------------------------------------------------------------------------
# Lineage name normalization
# --------------------------------------------------------------------------
# Kept for backward compatibility only: the engine enables NO default mapping. Supply yours via
# --lineage-map or config [quant] lineage_map; when it is empty the original lineage name is used
# uppercased. Syntax example: "LongName=L1,OtherName=L2".
DEFAULT_LINEAGE_MAP = ""


def parse_lineage_map(spec):
    m = {}
    if not spec:
        return m
    for item in spec.split(","):
        item = item.strip()
        if not item or "=" not in item:
            continue
        k, v = item.split("=", 1)
        m[k.strip().upper()] = v.strip().upper()
    return m


def normalize_lineage(name, lmap):
    u = (name or "").strip().upper()
    if not u:
        return "UNKNOWN"
    return lmap.get(u, u)


# --------------------------------------------------------------------------
# Input parsing
# --------------------------------------------------------------------------
def strip_version(gid):
    return gid.rsplit(".", 1)[0] if "." in gid else gid


def load_ancestry(path, lmap):
    """
    Two formats are supported:
      1) with a header that has an Ancestry column (e.g. Clean_Gene_Ancestry_Mapping.tsv:
         Gene_ID Chromosome Subgenome Ancestry)
      2) two columns without a header: GeneID  Lineage
    Returns {gene: normalized lineage}; IDs with the version suffix stripped are registered too.
    """
    anc = {}
    with open(path, encoding="utf-8") as f:
        first = f.readline()
        cols = first.rstrip("\n").split("\t")
        if len(cols) == 1:
            cols = first.split()
        low = [c.strip().lower() for c in cols]
        if "ancestry" in low:
            gi = 0
            ai = low.index("ancestry")
            header = True
        elif "lineage" in low:
            gi = 0
            ai = low.index("lineage")
            header = True
        else:
            gi, ai, header = 0, 1, False
        lines = f if header else [first] + list(f)
        for line in lines:
            if not line.strip() or line.startswith("#"):
                continue
            p = line.rstrip("\n").split("\t")
            if len(p) <= ai:
                p = line.split()
            if len(p) <= ai:
                continue
            g = p[gi].strip()
            lin = normalize_lineage(p[ai], lmap)
            anc[g] = lin
            anc.setdefault(strip_version(g), lin)
    return anc


def load_clusters(path):
    """cid \\t g1,g2,... -> (cid->[genes], gene->cid)"""
    c2g, g2c = {}, {}
    with open(path, encoding="utf-8") as f:
        for line in f:
            if not line.strip() or line.startswith("#"):
                continue
            p = line.rstrip("\n").split("\t")
            if len(p) < 2:
                continue
            cid = p[0].strip()
            genes = [x.strip() for x in p[1].split(",") if x.strip() and x.strip() != "-"]
            if cid.lower() in ("clusterid", "cluster_id", "cluster"):
                continue
            c2g[cid] = genes
            for g in genes:
                g2c[g] = cid
    return c2g, g2c


def parse_gff_genes(path, feature_types=("gene",)):
    """Returns {gene_id: (chrom, start0, end, strand)} (0-based half-open)."""
    genes = {}
    with open(path, encoding="utf-8") as f:
        for line in f:
            if not line.strip() or line.startswith("#"):
                continue
            p = line.rstrip("\n").split("\t")
            if len(p) < 9 or p[2] not in feature_types:
                continue
            gid = None
            for item in p[8].split(";"):
                item = item.strip()
                if item.startswith("ID="):
                    gid = item[3:].strip()
                    break
            if not gid:
                continue
            genes[gid] = (p[0], int(p[3]) - 1, int(p[4]), p[6] if p[6] in "+-" else "+")
    return genes


def read_peaks(path):
    """BED / narrowPeak -> list of (chrom, start, end, name)"""
    out = []
    with open(path) as f:
        n = 0
        for line in f:
            if not line.strip() or line.startswith(("#", "track", "browser")):
                continue
            p = line.rstrip("\n").split("\t")
            n += 1
            name = p[3] if len(p) > 3 and p[3] not in ("", ".") else "peak_{}".format(n)
            out.append((p[0], int(p[1]), int(p[2]), name))
    return out


# --------------------------------------------------------------------------
# Feature table
# --------------------------------------------------------------------------
FEATURE_COLS = ["FeatureID", "Chrom", "Start", "End", "Length", "GroupID", "ClusterID",
                "Gene", "Ancestry", "Type", "PeakIDs", "ProjCoverage"]


def load_features(path):
    """Read the features.tsv written by 02_build_features -> list of dict (in file order, i.e. the feature index)"""
    feats = []
    with open(path, encoding="utf-8") as f:
        header = f.readline().rstrip("\n").split("\t")
        for line in f:
            if not line.strip():
                continue
            p = line.rstrip("\n").split("\t")
            d = dict(zip(header, p))
            d["Start"] = int(d["Start"])
            d["End"] = int(d["End"])
            feats.append(d)
    return feats


# --------------------------------------------------------------------------
# BAM
# --------------------------------------------------------------------------
def open_bam(path):
    """
    Use pysam for production runs. When the environment variable ATACQ_MOCK_BAM=1 is set (or the
    input is a .sam and pysam is not available), use the mockbam written for the tests (a pure
    Python SAM reader that implements only the interfaces this pipeline needs).
    """
    use_mock = os.environ.get("ATACQ_MOCK_BAM") == "1"
    if not use_mock:
        try:
            import pysam  # noqa
        except ImportError:
            if path.endswith(".sam"):
                use_mock = True
            else:
                die("pysam is required: conda install -c bioconda pysam")
    if use_mock:
        from mockbam import MockAlignmentFile
        return MockAlignmentFile(path)
    import pysam
    return pysam.AlignmentFile(path, "rb")


def bam_references(path):
    b = open_bam(path)
    refs = list(zip(b.references, b.lengths))
    b.close()
    return refs


def write_json(obj, path):
    import json
    with open(path, "w") as f:
        json.dump(obj, f, indent=2, ensure_ascii=False, sort_keys=True)


def read_json(path):
    import json
    with open(path) as f:
        return json.load(f)
