#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
simulate_mini.py -- generate a small "sugarcane-style" polyploid ATAC test dataset (with
ground truth) for end-to-end validation.

Genome: a parameterised set of homeologous copies (7 by default: chromosomes Chr1A..Chr1G)
  SOF: A B C D   (0.3% divergence between copies; at 1/4 of the loci A and B are exactly
                  identical -> in principle indistinguishable)
  SRO: E F       (about 2% divergence from SOF)
  SSP: G         (about 3-4% divergence from SOF/SRO)
  --n-sof/--n-sro/--n-ssp set the copy numbers (e.g. 5+6+2=13 to mimic the FJDY
  configuration);
  --within-sof/--within-sro/--within-ssp set the between-copy divergence rate within each
  lineage.
Each locus: spacer + 3 kb upstream promoter + TSS + gene body; about half of the loci are on
  the minus strand;
  open elements E1 (core promoter, TSS-200..+50) and E2 (upstream enhancer,
  TSS-1600..-1400);
  some copies carry a small indel between the elements; some loci are missing on some
  copies; transposon (TE) copies are scattered through the spacers;
  plus 2 non-clustered genes that exist only on the last copy, and 1 distal open site of TE
  origin.
Ground truth: the openness of every (locus, element, copy) (some copies are closed; the
  second copy that is identical to the first SOF copy is deliberately set closed).
Sequencing: NFR/mononucleosome mixture of fragment lengths
  (--nfr-frac/--nfr-mean/--mono-mean); paired-end --read-len (2x75 by default); 0.3%
  sequencing error; 5% PCR duplicates (0x400 on the primary alignment); 30% background
  fragments.
Alignment: a pure-Python simple aligner mimicking bowtie2 -k 20 --no-mixed --no-discordant
  (AS = -6 x mismatch, YS = mate AS, primary = highest score, ties broken at random).
Output (--outdir): genome.fa(.fai) genes.gff3 clusters.tsv ancestry.tsv peaks.bed
  sim.sam sim_R1.fq sim_R2.fq (for reproducing with a real aligner; disable with
  --no-fastq)
  truth_fragments.tsv truth_elements.tsv lineage_{SOF,SRO,SSP}_k15.fasta samples.tsv
"""
import os
import sys
import math
import random
import argparse
from collections import defaultdict

ANC_NAME = {"SOF": "Officinarum", "SRO": "Robustum", "SSP": "Spontaneum"}
SUBG = {"SOF": "Noble", "SRO": "Wild", "SSP": "Wild"}
BASES = "ACGT"
RC = str.maketrans("ACGT", "TGCA")


def rc(s):
    return s.translate(RC)[::-1]


def rand_seq(n, rng):
    return "".join(rng.choice("AACCGGTTAT") for _ in range(n))


def mutate(seq, rate, rng):
    s = list(seq)
    for i in range(len(s)):
        if rng.random() < rate:
            s[i] = rng.choice([b for b in BASES if b != s[i]])
    return "".join(s)


# block coordinates (gene orientation): [0,3000) spacer; [3000,6000) upstream; TSS=6000; [6000,8000) gene body; [8000,9000) spacer
BLOCK = 9000
TSS = 6000
E1 = (5800, 6050)
E2 = (4400, 4600)
INDEL_ZONE = (4700, 5700)      # between E2 and E1


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("-o", "--outdir", required=True)
    ap.add_argument("--loci", type=int, default=24)
    ap.add_argument("--depth", type=float, default=1.0)
    ap.add_argument("--seed", type=int, default=7, help="random seed for the genome")
    ap.add_argument("--frag-seed", type=int, default=None,
                    help="random seed for openness/fragments/sequencing (defaults to --seed); use a different value for each sample sharing one genome")
    ap.add_argument("--sample", default="SIM")
    ap.add_argument("--within-sof", type=float, default=0.003, help="between-copy divergence rate within SOF")
    ap.add_argument("--n-sof", type=int, default=4, help="number of SOF copies (default 4)")
    ap.add_argument("--n-sro", type=int, default=2, help="number of SRO copies (default 2)")
    ap.add_argument("--n-ssp", type=int, default=1, help="number of SSP copies (default 1)")
    ap.add_argument("--within-sro", type=float, default=0.005, help="between-copy divergence rate within SRO")
    ap.add_argument("--within-ssp", type=float, default=0.002, help="between-copy divergence rate within SSP")
    ap.add_argument("--nfr-frac", type=float, default=0.6, help="fraction of NFR (short) fragments")
    ap.add_argument("--nfr-mean", type=float, default=100.0, help="mean NFR fragment length")
    ap.add_argument("--mono-mean", type=float, default=200.0, help="mean mononucleosome fragment length")
    ap.add_argument("--read-len", type=int, default=75, help="read length (one read for each end)")
    ap.add_argument("--no-fastq", action="store_true", help="do not write FASTQ (written by default, for use with a real aligner)")
    args = ap.parse_args()
    rng = random.Random(args.seed)
    os.makedirs(args.outdir, exist_ok=True)
    O = lambda n: os.path.join(args.outdir, n)  # noqa

    te = rand_seq(400, rng)

    # ------------------------------------------------------------ copy layout (7 copies by default: 4 SOF + 2 SRO + 1 SSP)
    n_tot = args.n_sof + args.n_sro + args.n_ssp
    if not (1 <= n_tot <= 26):
        raise SystemExit("the total number of copies must be within 1..26")
    ALLELES_L = list("ABCDEFGHIJKLMNOPQRSTUVWXYZ"[:n_tot])
    sof_l = ALLELES_L[:args.n_sof]
    sro_l = ALLELES_L[args.n_sof:args.n_sof + args.n_sro]
    ssp_l = ALLELES_L[args.n_sof + args.n_sro:]
    LIN_L = {}
    LIN_L.update({a: "SOF" for a in sof_l})
    LIN_L.update({a: "SRO" for a in sro_l})
    LIN_L.update({a: "SSP" for a in ssp_l})
    idcopy2 = sof_l[1] if len(sof_l) >= 2 else None   # the copy that is exactly identical to the first SOF copy

    def pick(lst, i):
        return {lst[i]} if 0 <= i < len(lst) else set()

    # ------------------------------------------------------------ build the locus blocks of every copy
    chrom_seq = {a: [] for a in ALLELES_L}
    chrom_len = {a: 0 for a in ALLELES_L}
    genes = []          # (gid, allele, locus, chrom, start0, end, strand)
    elements = []       # (locus, elem, allele, chrom, start, end)
    missing = {5: pick(sof_l, 3) | pick(sro_l, 1), 11: pick(sof_l, 2), 17: pick(sro_l, 0)}
    identicalAB = {i for i in range(args.loci) if i % 4 == 1}
    te_loci = set(rng.sample(range(args.loci), min(10, args.loci)))
    for L in range(args.loci):
        anc = rand_seq(BLOCK, rng)
        sofsro = mutate(anc, 0.02, rng)
        cons = {"SOF": mutate(sofsro, 0.01, rng), "SRO": mutate(sofsro, 0.01, rng),
                "SSP": mutate(anc, 0.03, rng)}
        strand = "+" if L % 2 == 0 else "-"
        base_A = mutate(cons["SOF"], args.within_sof, rng)
        for a in ALLELES_L:
            if a in missing.get(L, set()):
                blk = rand_seq(BLOCK, rng)
                present = False
            else:
                lin = LIN_L[a]
                if a == sof_l[0]:
                    blk = base_A
                elif a == idcopy2 and L in identicalAB:
                    blk = base_A
                else:
                    blk = mutate(cons[lin], {"SOF": args.within_sof, "SRO": args.within_sro,
                                             "SSP": args.within_ssp}[lin], rng)
                present = True
            e1, e2 = list(E1), list(E2)
            if present and rng.random() < 0.4 and not (a == idcopy2 and L in identicalAB):
                p = rng.randint(*INDEL_ZONE)
                d = rng.choice([-8, -5, -3, -1, 1, 2, 4, 9])
                if d > 0:
                    blk = blk[:p] + rand_seq(d, rng) + blk[p:]
                else:
                    blk = blk[:p] + blk[p - d:]
                # p lies after E2 and before E1: only E1, the TSS and everything after shift
                e1 = [e1[0] + d, e1[1] + d]
                tss = TSS + d
            else:
                tss = TSS
            if L in te_loci and rng.random() < 0.6:
                p = rng.randint(200, 2400)
                blk = blk[:p] + mutate(te, 0.02, rng) + blk[p + 400:]
            n = len(blk)
            off = chrom_len[a]
            if strand == "+":
                seq = blk
                g0, g1 = off + tss, off + tss + 2000
                m = lambda x: off + x  # noqa
                E1g = (m(e1[0]), m(e1[1]))
                E2g = (m(e2[0]), m(e2[1]))
            else:
                seq = rc(blk)
                # gene-orientation coordinate x -> genome off + n - 1 - x
                g0, g1 = off + n - (tss + 2000), off + n - tss
                E1g = (off + n - e1[1], off + n - e1[0])
                E2g = (off + n - e2[1], off + n - e2[0])
            chrom = "Chr1" + a
            chrom_seq[a].append(seq)
            chrom_len[a] += n
            if present:
                gid = "SIM_Chr1{}{:06d}".format(a, L + 1)
                genes.append((gid, a, L, chrom, g0, g1, strand))
                elements.append((L, "E1", a, chrom, E1g[0], E1g[1]))
                if L % 3 != 2:
                    elements.append((L, "E2", a, chrom, E2g[0], E2g[1]))
    # non-clustered gene (only on the last copy) + distal TE open site (on the first copy)
    last_l = ALLELES_L[-1]
    first_l = ALLELES_L[0]
    extra_elems = []
    for j in range(2):
        blk = rand_seq(BLOCK, rng)
        off = chrom_len[last_l]
        chrom_seq[last_l].append(blk)
        chrom_len[last_l] += BLOCK
        gid = "SIM_Chr1{}{:06d}".format(last_l, 900 + j)
        genes.append((gid, last_l, 900 + j, "Chr1" + last_l, off + TSS, off + TSS + 2000, "+"))
        extra_elems.append((900 + j, "E1", last_l, "Chr1" + last_l, off + E1[0], off + E1[1]))
    blk = rand_seq(3000, rng)
    blk = blk[:1300] + mutate(te, 0.01, rng) + blk[1700:]
    off = chrom_len[first_l]
    chrom_seq[first_l].append(blk)
    chrom_len[first_l] += 3000
    extra_elems.append((999, "TE", first_l, "Chr1" + first_l, off + 1300, off + 1700))
    genome = {"Chr1" + a: "".join(chrom_seq[a]) for a in ALLELES_L}

    # ------------------------------------------------------------ ground-truth openness
    if args.frag_seed is not None:
        rng = random.Random(args.frag_seed)
    truth_el = []
    for (L, e, a, c, s, t) in elements + extra_elems:
        lin = LIN_L[a]
        if L >= 900:
            lam = 60.0
        else:
            base = rng.gammavariate(2.0, 30.0) * (1.0 if e == "E1" else 0.6)
            eff = {"SOF": 0.7, "SRO": 2.0, "SSP": 1.0}[lin]
            lam = base * eff * rng.lognormvariate(0, 0.5)
            if rng.random() < 0.25:
                lam = 0.0
            if a == idcopy2 and L in identicalAB:
                lam = 0.0               # identical to A but closed -> indistinguishable
        truth_el.append((L, e, a, c, s, t, lam * args.depth))

    # ------------------------------------------------------------ fragments
    def frag_len():
        if rng.random() < args.nfr_frac:
            return int(min(180, max(55, rng.gauss(args.nfr_mean, 25))))
        return int(min(300, max(150, rng.gauss(args.mono_mean, 30))))

    frags = []   # (chrom, start, end, source)
    for (L, e, a, c, s, t, lam) in truth_el:
        n = 0
        # Poisson
        x, p = 0, 1.0
        lim = math.exp(-lam) if lam < 500 else 0
        if lam > 0:
            if lam < 500:
                while True:
                    p *= rng.random()
                    if p < lim:
                        break
                    x += 1
                n = x
            else:
                n = int(rng.gauss(lam, math.sqrt(lam)))
        for _ in range(n):
            ln = frag_len()
            ctr = rng.randint(s - 20, t + 20)
            st = max(0, ctr - ln // 2)
            frags.append((c, st, min(st + ln, len(genome[c])), "{}:{}:{}".format(L, e, a)))
    n_el = len(frags)
    n_bg = int(n_el * 0.3 / 0.7)
    chroms = sorted(genome)
    tot_len = sum(len(genome[c]) for c in chroms)
    for _ in range(n_bg):
        r = rng.randrange(tot_len)
        for c in chroms:
            if r < len(genome[c]):
                break
            r -= len(genome[c])
        ln = frag_len()
        st = max(0, min(r, len(genome[c]) - ln))
        frags.append((c, st, st + ln, "background"))
    # PCR duplicates
    dups = []
    for fr in frags:
        if rng.random() < 0.05:
            dups.append(fr)
    frags = [(f, False) for f in frags] + [(f, True) for f in dups]
    rng.shuffle(frags)

    # ------------------------------------------------------------ simple aligner
    K = 20
    index = defaultdict(list)
    for c in chroms:
        g = genome[c]
        for i in range(0, len(g) - K + 1):
            index[g[i:i + K]].append((c, i))
    # high-copy k-mer (TE) restriction
    RL = args.read_len

    def align_read(r):
        """Return [(chrom, pos, strand, mm)]; strand '+' means the read aligns to the forward strand"""
        out = {}
        for strand, seq in (("+", r), ("-", rc(r))):
            L_ = len(seq)
            offs = sorted({0, 18, 36, max(0, L_ - K)})
            for o in offs:
                if o + K > L_:
                    continue
                for (c, i) in index.get(seq[o:o + K], ()):
                    p = i - o
                    if p < 0 or p + L_ > len(genome[c]):
                        continue
                    key = (c, p, strand)
                    if key in out:
                        continue
                    g = genome[c]
                    mm = 0
                    for x in range(L_):
                        if g[p + x] != seq[x]:
                            mm += 1
                            if mm > 4:
                                break
                    if mm <= 4:
                        out[key] = mm
        return [(c, p, s, mm) for (c, p, s), mm in out.items()]

    def seqerr(s):
        s = list(s)
        for i in range(len(s)):
            if rng.random() < 0.003:
                s[i] = rng.choice([b for b in BASES if b != s[i]])
        return "".join(s)

    recs = []
    truth_rows = []
    n_unal = 0
    fq1 = fq2 = None
    if not args.no_fastq:
        fq1 = open(O("sim_R1.fq"), "w")
        fq2 = open(O("sim_R2.fq"), "w")
    for fi, ((c, st, en, src), isdup) in enumerate(frags):
        qn = "SIM:{:07d}".format(fi)
        fseq = genome[c][st:en]
        L_ = min(RL, len(fseq))
        left = seqerr(fseq[:L_])
        right = seqerr(rc(fseq[-L_:]))
        flip = rng.random() < 0.5
        r1, r2 = (right, left) if flip else (left, right)
        if fq1 is not None:
            fq1.write("@{}\n{}\n+\n{}\n".format(qn, r1, "I" * len(r1)))
            fq2.write("@{}\n{}\n+\n{}\n".format(qn, r2, "I" * len(r2)))
        truth_rows.append((qn, c, st, en, src, int(isdup)))
        a1 = align_read(r1)
        a2 = align_read(r2)
        pairs = []
        for (c1, p1, s1, m1) in a1:
            for (c2, p2, s2, m2) in a2:
                if c1 != c2 or s1 == s2:
                    continue
                if s1 == "+":
                    fs, fe = p1, p2 + len(r2)
                else:
                    fs, fe = p2, p1 + len(r1)
                if fe - fs < max(len(r1), len(r2)) or fe - fs > 1000:
                    continue
                if s1 == "+" and p2 < p1:
                    continue
                if s1 == "-" and p1 < p2:
                    continue
                pairs.append((-6 * m1 - 6 * m2, c1, p1, s1, m1, p2, s2, m2, fs, fe))
        if not pairs:
            n_unal += 1
            truth_rows[-1] = truth_rows[-1] + (0,)
            continue
        truth_rows[-1] = truth_rows[-1] + (1,)
        rng.shuffle(pairs)
        pairs.sort(key=lambda x: -x[0])
        pairs = pairs[:20]
        best = pairs[0][0]
        nbest = sum(1 for p in pairs if p[0] == best)
        mapq = 42 if nbest == 1 and len(pairs) == 1 else (1 if nbest > 1 else 23)
        for k, (sc, cc, p1, s1, m1, p2, s2, m2, fs, fe) in enumerate(pairs):
            tl = fe - fs
            prim = (k == 0)
            xs1 = [p for p in pairs if p is not pairs[k]]
            for mate in (1, 2):
                if mate == 1:
                    pos, s, mm, seq, mpos, ms, mmm = p1, s1, m1, r1, p2, s2, m2
                else:
                    pos, s, mm, seq, mpos, ms, mmm = p2, s2, m2, r2, p1, s1, m1
                flag = 0x1 | 0x2 | (0x40 if mate == 1 else 0x80)
                if s == "-":
                    flag |= 0x10
                if ms == "-":
                    flag |= 0x20
                if not prim:
                    flag |= 0x100
                if prim and isdup:
                    flag |= 0x400
                sseq = seq if s == "+" else rc(seq)
                leftmost = pos < mpos or (pos == mpos and mate == 1)
                tlen = tl if leftmost else -tl
                tags = ["AS:i:{}".format(-6 * mm), "YS:i:{}".format(-6 * mmm), "NM:i:{}".format(mm)]
                if xs1:
                    tags.append("XS:i:{}".format(max(-6 * (p[4] if mate == 1 else p[7]) for p in xs1)))
                recs.append((cc, pos, "\t".join([qn, str(flag), cc, str(pos + 1), str(mapq if prim else 0),
                                                  "{}M".format(len(seq)), "=", str(mpos + 1), str(tlen),
                                                  sseq, "I" * len(seq)] + tags)))
    if fq1 is not None:
        fq1.close()
        fq2.close()
    rank = {c: i for i, c in enumerate(chroms)}
    recs.sort(key=lambda x: (rank[x[0]], x[1]))

    # ------------------------------------------------------------ write output files
    with open(O("genome.fa"), "w", newline="\n") as f:   # explicit LF: keeps the .fai line-width assumption consistent (cross-platform)
        for c in chroms:
            f.write(">{}\n".format(c))
            g = genome[c]
            for i in range(0, len(g), 60):
                f.write(g[i:i + 60] + "\n")
    with open(O("genome.fa.fai"), "w") as f:
        off = 0
        for c in chroms:
            off += len(">{}\n".format(c))
            n = len(genome[c])
            f.write("{}\t{}\t{}\t60\t61\n".format(c, n, off))
            off += n + (n + 59) // 60
    with open(O("sim.sam"), "w") as f:
        f.write("@HD\tVN:1.6\tSO:coordinate\n")
        for c in chroms:
            f.write("@SQ\tSN:{}\tLN:{}\n".format(c, len(genome[c])))
        for _, _, line in recs:
            f.write(line + "\n")
    with open(O("genes.gff3"), "w") as f:
        f.write("##gff-version 3\n")
        for gid, a, L, c, s, e, strand in genes:
            f.write("{}\tsim\tgene\t{}\t{}\t.\t{}\t.\tID={}\n".format(c, s + 1, e, strand, gid))
    by_locus = defaultdict(list)
    for gid, a, L, c, s, e, strand in genes:
        by_locus[L].append(gid)
    with open(O("clusters.tsv"), "w") as f:
        for L in sorted(by_locus):
            if len(by_locus[L]) >= 2:
                f.write("Global_Cluster_{:06d}\t{}\n".format(L, ",".join(by_locus[L])))
    with open(O("ancestry.tsv"), "w") as f:
        f.write("Gene_ID\tChromosome\tSubgenome\tAncestry\n")
        for gid, a, L, c, s, e, strand in genes:
            f.write("{}\t{}\t{}\t{}\n".format(gid, c[3:], SUBG[LIN_L[a]], ANC_NAME[LIN_L[a]]))
    # called peaks: elements with high enough signal (with a random offset) + TE site + non-clustered gene + 2 false positives
    with open(O("peaks.bed"), "w") as f:
        n = 0
        for (L, e, a, c, s, t, lam) in truth_el:
            if lam >= 25 * args.depth:
                n += 1
                j1, j2 = rng.randint(-30, 30), rng.randint(-30, 30)
                f.write("{}\t{}\t{}\tsim_peak_{}\n".format(c, max(0, s + j1), t + j2, n))
        for c in ("Chr1" + ALLELES_L[min(1, n_tot - 1)],
                  "Chr1" + ALLELES_L[min(args.n_sof, n_tot - 1)]):
            n += 1
            p = rng.randint(1000, 2000)
            f.write("{}\t{}\t{}\tsim_peak_{}\n".format(c, p, p + 200, n))
    with open(O("truth_fragments.tsv"), "w") as f:
        f.write("QName\tChrom\tStart\tEnd\tSource\tIsDup\tAligned\n")
        for r in truth_rows:
            f.write("\t".join(str(x) for x in r) + "\n")
    with open(O("truth_elements.tsv"), "w") as f:
        f.write("Locus\tElement\tAllele\tLineage\tChrom\tStart\tEnd\tLambda\tIdenticalAB\n")
        for (L, e, a, c, s, t, lam) in truth_el:
            f.write("{}\t{}\t{}\t{}\t{}\t{}\t{}\t{:.2f}\t{}\n".format(
                L, e, a, LIN_L[a], c, s, t, lam,
                int(L in identicalAB and (a == sof_l[0] or a == idcopy2))))
    # lineage-specific 15-mers (present only on copies of that lineage)
    k = 15
    kset = {l: set() for l in ("SOF", "SRO", "SSP")}
    for c in chroms:
        g = genome[c]
        l = LIN_L[c[-1]]
        for i in range(0, len(g) - k + 1, 3):
            kset[l].add(g[i:i + k])
    for l in kset:
        others = set().union(*[kset[o] for o in kset if o != l])
        spec = sorted(x for x in kset[l] if x not in others and rc(x) not in others)
        rng.shuffle(spec)
        with open(O("lineage_{}_k15.fasta".format(l)), "w") as f:
            for i, x in enumerate(spec[:4000]):
                f.write(">0.57|{}_k15_0.55-0.60\n{}\n".format(l, x))
    with open(O("samples.tsv"), "w") as f:
        f.write("sample\tr1\tr2\tbam\tpeaks\n")
        f.write("{}\t-\t-\t{}\t{}\n".format(args.sample, os.path.abspath(O("sim.sam")),
                                              os.path.abspath(O("peaks.bed"))))
    print("[sim] genome {:.2f} Mb | genes {} | elements {} | fragments {} (duplicates {}) | unaligned {} | SAM records {}".format(
        tot_len / 1e6, len(genes), len(truth_el), len(frags), len(dups), n_unal, len(recs)))


if __name__ == "__main__":
    main()
