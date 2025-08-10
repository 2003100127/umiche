"""
The original umitools implementation from the Weng lab parses UMI from read names.
Modern sequencing/alignment tools now typically store UMI in the MB (Mapped-Barcode) tag
(or tags like RX/UB/UR). This version exclusively uses the MB tag, abandoning read-name parsing.

Raises ValueError if MB is missing, preventing silent fallback to empty strings as UMI.

PCR-duplicate criteria remain unchanged:
chrom + read5 + template_length + UMI (+ PE mate read).

Deduplication results now use a custom FU:i tag instead of flag 0x400:

FU:i:1 → duplicate

FU:i:0 → unique.
"""

import collections
import os
import sys
import time
from multiprocessing import Pool
from pathlib import Path
from typing import Dict, List, Tuple, Set

import pysam

__all__ = ["UmiDeduplicator", "deduplicate"]


class RsqRunStats(dict):
    def __init__(self):
        super().__init__(n_read=0, n_duplicate=0, n_non_duplicate=0)

class UmiDeduplicator:
    """Mark PCR duplicates using UMI stored in **`MB` tag** (string)."""

    def __init__(
        self,
        bam_fpn: str | Path,
        *,
        n_proc: int = 8,
        debug: bool = False,
        count_locus: bool = False,
    ) -> None:
        self.bam_fpn = str(bam_fpn)
        self.n_proc = max(1, n_proc)
        self.debug = debug
        self.count_locus = count_locus

        self.stats = RsqRunStats()
        self.clusters: Dict[Tuple[str, int, int, str], List[str]] = collections.defaultdict(list)
        self.assigned: Dict[str, str] = {}
        self.dedup_count = 0
        self.output_bam: str | None = None

    # ------------------------------------------------------------------
    def _log(self, *msg):
        if self.debug:
            print(time.strftime("[%H:%M:%S]"), *msg, file=sys.stderr)

    # ------------------------------------------------------------------
    @staticmethod
    def _get_umi(read: pysam.AlignedSegment) -> str:
        """Return UMI from `MB` tag; raise if missing."""
        try:
            return read.get_tag("MB")
        except KeyError:
            raise ValueError("Read is missing MB tag: " + read.query_name)

    def _first_pass(self, bam: pysam.AlignmentFile, chrom: str):
        r1_seen: Set[Tuple[int, str, int]] = set()
        r2_seen: Set[Tuple[int, str, int]] = set()
        dup_qnames: Set[str] = set()
        counts_loc: Dict[str, int] = collections.Counter()

        for read in bam.fetch(reference=chrom):
            umi = self._get_umi(read)
            read5 = read.reference_start
            tlen = read.template_length if read.template_length >= 0 else -read.template_length
            key = (read5, umi, tlen)
            if not read.is_reverse and read.is_read1:
                if key in r1_seen:
                    dup_qnames.add(read.query_name)
                else:
                    r1_seen.add(key)
            elif not read.is_reverse and read.is_read2:
                if key in r2_seen:
                    dup_qnames.add(read.query_name)
                else:
                    r2_seen.add(key)

            self.clusters[(chrom, read5, tlen, umi)].append(read.query_name)
            if self.count_locus:
                counts_loc[f"{read5},{tlen}"] += 1

        return dup_qnames, counts_loc

    def _process_chrom(self, chrom: str):
        self._log(f"chromosome {chrom} … start")
        bam = pysam.AlignmentFile(self.bam_fpn, "rb")
        dup_qnames, counts_loc = self._first_pass(bam, chrom)
        bam.close()

        bam = pysam.AlignmentFile(self.bam_fpn, "rb")
        tmp_bam_fpn = f"{self.bam_fpn}.{chrom}.bam"
        out = pysam.AlignmentFile(tmp_bam_fpn, "wb", template=bam)
        dup_count = 0
        for read in bam.fetch(reference=chrom):
            is_dup = read.query_name in dup_qnames
            read.set_tag("FU", 1 if is_dup else 0, value_type="i", replace=True)
            if is_dup:
                dup_count += 1
            out.write(read)
        out.close(); bam.close()

        if self.count_locus:
            with open(f"{self.bam_fpn}.{chrom}.loc_count", "w") as fh:
                for loc, cnt in counts_loc.items():
                    fh.write(f"{chrom}\t{loc}\t{cnt}\n")

        self._log(f"chromosome {chrom} … done")
        return tmp_bam_fpn, dup_count

    def run(self):
        t0 = time.time()
        if not Path(f"{self.bam_fpn}.bai").exists():
            pysam.index(self.bam_fpn)

        chromosomes = pysam.AlignmentFile(self.bam_fpn, "rb").references
        with Pool(self.n_proc) as pool:
            results = pool.map(self._process_chrom, chromosomes)

        merged_fpn = f"{Path(self.bam_fpn).with_suffix('')}.fu_tagged.bam"
        template_bam = pysam.AlignmentFile(self.bam_fpn, "rb")
        merged = pysam.AlignmentFile(merged_fpn, "wb", template=template_bam)
        template_bam.close()

        for tmp_bam_fpn, dup_count in results:
            self.dedup_count += dup_count
            tmp_bam = pysam.AlignmentFile(tmp_bam_fpn, "rb")
            for read in tmp_bam.fetch(until_eof=True):
                merged.write(read)
            tmp_bam.close()
            os.remove(tmp_bam_fpn)
            if Path(f"{tmp_bam_fpn}.bai").exists():
                os.remove(f"{tmp_bam_fpn}.bai")

        merged.close(); pysam.index(merged_fpn)
        self.output_bam = merged_fpn
        self.stats["n_duplicate"] = self.dedup_count
        self.stats["n_read"] = sum(len(v) for v in self.clusters.values())
        self.stats["n_non_duplicate"] = self.stats["n_read"] - self.dedup_count
        self._log(f"finished in {time.time() - t0:.1f}s – duplicates: {self.dedup_count}")
        return self

def deduplicate(bam_fpn: str | Path, **kwargs):
    """Run deduplication; return `(clusters, dup_count, assigned)`.

    Output BAM is `<input>.fu_tagged.bam` with `FU` duplicate tags.
    """
    d = UmiDeduplicator(bam_fpn, **kwargs).run()
    return d.clusters, d.dedup_count, d.assigned


if __name__ == "__main__":
    from umiche.path import to
    # samtools view /mnt/d/Document/Programming/Python/umiche/umiche/data/r1/umitools/umitools.test.RNA-seq.sorted.bam | head -n 5

    clusters, dedup_num, assigned = deduplicate(
        bam_fpn='/mnt/d/Document/Programming/Python/umiche/umiche/data/r1/umitools/umitools.test.RNA-seq.sorted.tagged.sorted.bam',
        n_proc=12,
        debug=True,
    )

    print("PCR duplicates:", dedup_num)
    print("c:", clusters)
    print(len(assigned))
    print("", {k: assigned[k] for k in list(assigned)[:5]})