__version__ = "v1.3"  # add cell+gene mode: dedup by UMI within (chrom, cell, gene) only
__copyright__ = "Copyright 2025"
__license__ = "GPL-3.0"
__developer__ = "Jianfeng Sun"
__maintainer__ = "Jianfeng Sun"
__email__="jianfeng.sunmt@gmail.com"

from typing import Dict, List, Tuple, Set, Optional

import os
import pandas as pd
from umiche.util.Console import Console


class UMItoolsReplica:
    """
    Weng-style duplicate detection with PD tag write-back.

    Modes:
      1) If BOTH `cell_col` and `gene_col` are provided:
         - For each (chrom, cell, gene) group, mark duplicates based on UMI ONLY
           (i.e., reads sharing the same UMI in the group → keep first, others PD=1).
         - This matches single-cell practice where read names are unique and
           position/tlen jitter should not block dedup within a gene.
      2) Otherwise (no cell/gene grouping):
         - Classic UMItools behavior on leftmost mate:
           * read1 duplicates: duplicated over (pos, UMI, tlen)
           * read2 duplicates: duplicated over str(pos)+UMI+str(tlen)

    Write-back:
      - Each read gets SAM tag PD:i:{0|1}, no use of 0x400.
      - Writing is chromosome-wise and then merged, consistent with previous behavior.
    """

    def __init__(
        self,
        bam_fpn: str,
        umi_col: str = "MB",
        cell_col: Optional[str] = None,
        gene_col: Optional[str] = None,
        verbose: bool = False,
    ):
        import pysam
        self.pysam = pysam
        self.bam_fpn = bam_fpn
        self.umi_col = umi_col
        self.cell_col = cell_col
        self.gene_col = gene_col

        self.console = Console()
        self.console.verbose = verbose

    # === Calculate duplicate QNAME sets per chromosome ===
    def compute_duplicates_by_chrom(self, df: pd.DataFrame) -> Dict[str, Set[str]]:
        # Required columns
        req = {"qname", "chrom", "pos", "tlen", "is_read1", "is_read2", "is_reverse", self.umi_col}
        if self.cell_col: req.add(self.cell_col)
        if self.gene_col: req.add(self.gene_col)
        miss = req - set(df.columns)
        if miss:
            raise ValueError(f"DataFrame missing required columns: {sorted(miss)}")

        # Same as fetch(reference=...): only keep aligned records
        df = df[df["chrom"].notna()].copy()

        dup_by_chrom: Dict[str, Set[str]] = {}

        # --- Single-cell gene-level mode: BOTH cell_col & gene_col are set ---
        if self.cell_col and self.gene_col:
            group_keys = ["chrom", self.cell_col, self.gene_col]
            for key_vals, sub in self.console._tqdm(
                df.groupby(group_keys, sort=False),
                desc="[group chrom/cell/gene by UMI]",
                unit="grp",
                position=0,
                leave=True,
                dynamic_ncols=False,
            ):
                # key_vals -> (chrom, cell, gene); dedup by UMI only within the group
                chrom = key_vals[0] if isinstance(key_vals, tuple) else key_vals
                dups: Set[str] = set()

                # Mark second+ occurrences of the same UMI as duplicates
                dup_mask = sub.duplicated([self.umi_col], keep="first")
                if dup_mask.any():
                    dups.update(sub.loc[dup_mask, "qname"].tolist())

                if dups:
                    dup_by_chrom.setdefault(str(chrom), set()).update(dups)
            return dup_by_chrom

        # --- Classic mode (no cell/gene grouping): UMItools leftmost logic ---
        group_keys = ["chrom"]
        tag_title = "group chrom"
        for key_vals, sub in self.console._tqdm(
            df.groupby(group_keys, sort=False),
            desc=f"[{tag_title}]",
            unit="grp",
            position=0,
            leave=True,
            dynamic_ncols=False,
        ):
            chrom = key_vals if isinstance(key_vals, str) else key_vals[0]
            dups: Set[str] = set()

            # R1: tuple key
            r1 = sub[(~sub["is_reverse"]) & (sub["is_read1"])]
            if not r1.empty:
                dup_mask = r1.duplicated(["pos", self.umi_col, "tlen"], keep="first")
                if dup_mask.any():
                    dups.update(r1.loc[dup_mask, "qname"].tolist())

            # R2: string key
            r2 = sub[(~sub["is_reverse"]) & (sub["is_read2"])]
            if not r2.empty:
                k = r2["pos"].astype(str) + r2[self.umi_col] + r2["tlen"].astype(str)
                dup_mask = k.duplicated(keep="first")
                if dup_mask.any():
                    dups.update(r2.loc[dup_mask, "qname"].tolist())

            if dups:
                dup_by_chrom.setdefault(str(chrom), set()).update(dups)

        return dup_by_chrom

    # === Merge helper ===
    def _merge_parts(self, infile: str, refs: List[str]) -> str:
        prefix = infile[:-4] if infile.endswith(".bam") else infile
        out_bam = f"{prefix}.deumi.sorted.bam"

        tmpl = self.pysam.AlignmentFile(infile, "rb")
        out = self.pysam.AlignmentFile(out_bam, "wb", template=tmpl)
        tmpl.close()

        for chrom in sorted(refs):
            part = f"{infile}.{chrom}.bam"
            self.pysam.index(part)
            bam = self.pysam.AlignmentFile(part, "rb")
            for r in bam.fetch(until_eof=True):
                out.write(r)
            bam.close()
            for p in (part, part + ".bai"):
                try: os.remove(p)
                except FileNotFoundError: pass

        out.close()
        self.pysam.index(out_bam)
        return out_bam

    def _count_total_and_nondub_by_pd(self, bam_fpn: str) -> Tuple[int, int]:
        total = nondup = 0
        bam = self.pysam.AlignmentFile(bam_fpn, "rb")
        for r in bam.fetch(until_eof=True):
            total += 1
            pd_val = r.get_tag("PD") if r.has_tag("PD") else 0
            if int(pd_val) == 0:
                nondup += 1
        bam.close()
        return total, nondup

    # === Write back: per chromosome, mark PD=1 for duplicate QNAMEs in that chromosome, else PD=0 ===
    @Console.vignette()
    def mark_and_merge(self, bam_fpn: str, dup_by_chrom: Dict[str, Set[str]]) -> Tuple[str, int, int, int]:
        bam = self.pysam.AlignmentFile(bam_fpn, "rb")
        if not bam.has_index():
            self.pysam.index(bam_fpn)
        refs = list(bam.references)
        bam.close()

        for chrom in self.console._tqdm(
            refs, desc="[chrom merge]", unit="chr", position=0, leave=True, dynamic_ncols=False
        ):
            dup_q = dup_by_chrom.get(chrom, set())
            inb = self.pysam.AlignmentFile(bam_fpn, "rb")
            out = self.pysam.AlignmentFile(f"{bam_fpn}.{chrom}.bam", "wb", template=inb)
            for read in inb.fetch(reference=chrom):
                is_dup = 1 if (read.query_name in dup_q) else 0
                read.set_tag("PD", is_dup, value_type="i")
                out.write(read)
            inb.close()
            out.close()

        out_bam = self._merge_parts(bam_fpn, refs)
        total, nondup = self._count_total_and_nondub_by_pd(out_bam)
        return out_bam, total, nondup, total - nondup

    # === End-to-end entry point ===
    def run(self, df_bam: pd.DataFrame):
        dup_by_chrom = self.compute_duplicates_by_chrom(df_bam)
        return self.mark_and_merge(self.bam_fpn, dup_by_chrom)


if __name__ == "__main__":
    # Example usage (Reader is in your project; pass column names explicitly there)
    from umiche.bam.Reader import ReaderChunk

    # bam_fpn = "/mnt/d/Document/Programming/Python/umiche/umiche/data/r1/umitools/umitools.test.RNA-seq.sorted.tagged.bam"
    bam_fpn = "/mnt/d/Document/Programming/Python/umiche/umiche/data/r1/trumicount/10xn9k_10c/10xn9k_10c_tagged.sorted.filtered.sorted.bam"

    # df_bam must contain: qname, chrom, pos, tlen, is_read1, is_read2, is_reverse, and UMI column (default MB)
    # For single-cell gene-level dedup, set both cell_col and gene_col in the tool below.
    df_bam = ReaderChunk(
        bam_fpn=bam_fpn,
        bam_fields=None,
        # tag_whitelist=["MB"],  # Using MB as UMI
        tag_whitelist=["MB", "CB", "XF"],  # Using MB as UMI
        categorize=["chrom"],
        verbose=True,
    ).todf(chunk_size=2_000_000)

    tool = UMItoolsReplica(
        bam_fpn=bam_fpn,
        umi_col="MB",
        cell_col="CB",
        gene_col="XF",
        # cell_col=None,
        # gene_col=None,
        verbose=True,
    )
    OUT_BAM, TOTAL, NONDUP, DUP = tool.run(df_bam)
    print(f"[DONE] wrote: {OUT_BAM}")
    print(f"[DONE] Total={TOTAL}, NonDup={NONDUP}, PCR_Duplicates={DUP}")