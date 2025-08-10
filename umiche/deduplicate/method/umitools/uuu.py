#!/usr/bin/env python3
# -*- coding: utf-8 -*-

import os
import sys
from multiprocessing import Pool
from typing import List, Tuple

import pysam as ps


def _umi_from_qname(qname: str) -> str:
    parts = qname.split("_")
    if len(parts) < 2:
        raise ValueError(f"Read name lacks UMI (expected qname_UMI_...): {qname}")
    return parts[1]


def _mark_one_reference(args: Tuple[str, str]) -> Tuple[str, int]:
    """
    For a single reference, execute a two-pass process:
        First pass: On the leftmost mate (not is_reverse) only, collect duplicate QNAMEs using the original key.
        Second pass: Write the partitioned BAM and mark duplicate reads with 0x400.
        Return (chrom, dup_count_on_chrom).
    """
    infile, chrom = args
    bam = ps.AlignmentFile(infile, "rb")
    out = ps.AlignmentFile(f"{infile}.{chrom}.bam", "wb", template=bam)

    r1_seen, r2_seen = {}, {}
    r1_dup, r2_dup = {}, {}

    for read in bam.fetch(reference=chrom):
        qn = read.query_name
        print(qn)
        umi = _umi_from_qname(qn)
        read5 = read.reference_start
        tlen = read.template_length

        if (not read.is_reverse) and read.is_read1:
            key = (read5, umi, tlen)
            if key in r1_seen:
                r1_dup[qn] = 1
            else:
                r1_seen[key] = 0

        elif (not read.is_reverse) and read.is_read2:
            key = str(read5) + umi + str(tlen)
            if key in r2_seen:
                r2_dup[qn] = 1
            else:
                r2_seen[key] = 0

    bam.close()

    # Write the partitioned BAM and flag duplicates with 0x400.
    bam = ps.AlignmentFile(infile, "rb")
    for read in bam.fetch(reference=chrom):
        if (read.query_name in r1_dup) or (read.query_name in r2_dup):
            read.flag = read.flag | 0x400
        out.write(read)
    bam.close()
    out.close()

    return chrom, (len(r1_dup) + len(r2_dup))


def _merge_parts(infile: str, refs: List[str]) -> str:
    """
    Merge the split files into the final BAM and build the index.
        For input xxx.sorted.bam, the output will be xxx.sorted.dedup.sorted.bam.
    """
    prefix = infile[:-4] if infile.endswith(".bam") else infile
    out_bam = f"{prefix}.deumi.sorted.bam"

    tmpl = ps.AlignmentFile(infile, "rb")
    out = ps.AlignmentFile(out_bam, "wb", template=tmpl)
    tmpl.close()

    for chrom in sorted(refs):
        part = f"{infile}.{chrom}.bam"
        ps.index(part)
        bam = ps.AlignmentFile(part, "rb")
        for r in bam.fetch(until_eof=True):
            out.write(r)
        bam.close()
        try:
            os.remove(part)
        except FileNotFoundError:
            pass
        try:
            os.remove(part + ".bai")
        except FileNotFoundError:
            pass

    out.close()
    ps.index(out_bam)
    return out_bam


def _count_total_and_nondub(bam_path: str) -> Tuple[int, int]:
    """samtools view -c -F 0x400"""
    total = nondub = 0
    bam = ps.AlignmentFile(bam_path, "rb")
    for r in bam.fetch(until_eof=True):
        total += 1
        if (r.flag & 0x400) == 0:
            nondub += 1
    bam.close()
    return total, nondub


class UMItoolsReplica:
    """
    Replicates Weng's PCR duplicate marking implementation.
        out_bam, total, nondup, dup = UMItoolsReplica(processes=8).run("input.sorted.bam")
    """

    def __init__(self, processes: int = 8):
        self.processes = int(processes)

    def run(self, infile: str) -> Tuple[str, int, int, int]:
        # indexinig
        bam = ps.AlignmentFile(infile, "rb")
        if not bam.has_index():
            ps.index(infile)
        refs = list(bam.references)
        bam.close()

        # Stage 1: Per-reference two-pass (collect → flag) → write splits
        if self.processes > 1:
            with Pool(self.processes) as pool:
                pool.map(_mark_one_reference, [(infile, r) for r in refs])
        else:
            for r in refs:
                _mark_one_reference((infile, r))

        # merging
        out_bam = _merge_parts(infile, refs)

        # stats
        total, nondup = _count_total_and_nondub(out_bam)
        dup = total - nondup
        return out_bam, total, nondup, dup


if __name__ == "__main__":
    IN_BAM = "/mnt/d/Document/Programming/Python/umiche/umiche/data/r1/umitools/umitools.test.RNA-seq.sorted.bam"
    tool = UMItoolsReplica(processes=8)
    OUT_BAM, TOTAL, NONDUP, DUP = tool.run(IN_BAM)
    print(f"[DONE] wrote: {OUT_BAM}")
    print(f"[DONE] Total={TOTAL}, NonDup={NONDUP}, PCR_Duplicates={DUP}")
