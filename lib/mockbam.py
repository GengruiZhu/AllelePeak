# -*- coding: utf-8 -*-
"""
mockbam.py -- for testing only: reads SAM text and emulates the pysam.AlignmentFile interface
used by this pipeline.

It is never used in production (it is enabled only when ATACQ_MOCK_BAM=1 is set).
"""
import re
from collections import defaultdict

_CIG = re.compile(r"(\d+)([MIDNSHP=X])")


class MockRead(object):
    __slots__ = ("query_name", "flag", "reference_name", "reference_start",
                 "reference_end", "mapping_quality", "template_length",
                 "next_reference_start", "query_sequence", "_tags")

    def __init__(self, p):
        self.query_name = p[0]
        self.flag = int(p[1])
        self.reference_name = p[2]
        self.reference_start = int(p[3]) - 1
        self.mapping_quality = int(p[4])
        span = 0
        if p[5] != "*":
            for n, op in _CIG.findall(p[5]):
                if op in "MDN=X":
                    span += int(n)
        self.reference_end = self.reference_start + span if span else None
        self.next_reference_start = int(p[7]) - 1
        self.template_length = int(p[8])
        self.query_sequence = None if p[9] == "*" else p[9]
        self._tags = {}
        for t in p[11:]:
            k, typ, v = t.split(":", 2)
            self._tags[k] = int(v) if typ == "i" else v

    # flag properties
    is_paired = property(lambda s: bool(s.flag & 0x1))
    is_proper_pair = property(lambda s: bool(s.flag & 0x2))
    is_unmapped = property(lambda s: bool(s.flag & 0x4))
    is_reverse = property(lambda s: bool(s.flag & 0x10))
    is_read1 = property(lambda s: bool(s.flag & 0x40))
    is_read2 = property(lambda s: bool(s.flag & 0x80))
    is_secondary = property(lambda s: bool(s.flag & 0x100))
    is_duplicate = property(lambda s: bool(s.flag & 0x400))
    is_supplementary = property(lambda s: bool(s.flag & 0x800))

    def get_tag(self, k):
        return self._tags[k]

    def has_tag(self, k):
        return k in self._tags


class MockAlignmentFile(object):
    def __init__(self, path, mode="r"):
        self.references = []
        self.lengths = []
        self._by_ref = defaultdict(list)
        self._all = []
        with open(path) as f:
            for line in f:
                if line.startswith("@"):
                    if line.startswith("@SQ"):
                        d = dict(x.split(":", 1) for x in line.rstrip("\n").split("\t")[1:])
                        self.references.append(d["SN"])
                        self.lengths.append(int(d["LN"]))
                    continue
                p = line.rstrip("\n").split("\t")
                r = MockRead(p)
                self._all.append(r)
                if p[2] != "*":
                    self._by_ref[p[2]].append(r)
        for k in self._by_ref:
            self._by_ref[k].sort(key=lambda r: r.reference_start)

    def fetch(self, contig=None, start=None, end=None, until_eof=False):
        if contig is None:
            for r in self._all:
                yield r
            return
        for r in self._by_ref.get(contig, []):
            if start is not None and (r.reference_end or r.reference_start + 1) <= start:
                continue
            if end is not None and r.reference_start >= end:
                continue
            yield r

    def close(self):
        pass

    def __enter__(self):
        return self

    def __exit__(self, *a):
        self.close()
