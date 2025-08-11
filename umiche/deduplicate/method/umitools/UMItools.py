__version__ = "v1.0"
__copyright__ = "Copyright 2025"
__license__ = "GPL-3.0"
__developer__ = "Jianfeng Sun"
__maintainer__ = "Jianfeng Sun"
__email__="jianfeng.sunmt@gmail.com"

from typing import Iterable, Dict, List, Tuple, Set

import os
import pandas as pd
import pysam as ps


class UMItoolsReplica:
    """
    Weng：
      - compute_duplicates_from_df(df) -> Set[qname]
      - compute_duplicates_from_chunks(chunk_iter) -> Set[qname]（内存友好）
      - mark_and_merge(in_bam, dup_qnames) -> <prefix>.deumi.sorted.bam
      - run_from_dataframe(df, in_bam) -> (out_bam, total, nondup, dup)
      - run_from_chunks(chunk_iter, in_bam) -> (out_bam, total, nondup, dup)  ← 大文件推荐

    rule kept the same：
      - Only check for duplicates on the leftmost mate (not is_reverse)
      - read1： (pos, UMI, tlen)
      - read2： str(pos) + UMI + str(tlen)
      - The original implementation processes each chromosome individually;
      - here the calculation can be done all at once/in blocks,
      but the write back is still done by chromosome and merged to keep the product consistent
    """

    # @@ /*** ---- calc PCR duplicates using df ---- ***/
    @staticmethod
    def compute_duplicates_by_chrom_from_df(df: pd.DataFrame) -> Dict[str, Set[str]]:
        req = {"qname", "MB", "chrom", "pos", "tlen", "is_read1", "is_read2", "is_reverse"}
        miss = req - set(df.columns)
        if miss:
            raise ValueError(f"DataFrame missing required columns: {miss}")

        # Only keep records that match the reference (same behavior as fetch(reference=...)
        df = df[df["chrom"].notna()].copy()

        dup_by_chrom: Dict[str, Set[str]] = {}
        for chrom, sub in df.groupby("chrom", sort=False):
            dups: Set[str] = set()
            # R1: key = (pos, UMI, tlen)
            r1 = sub[(~sub["is_reverse"]) & (sub["is_read1"])]
            if not r1.empty:
                dup_mask = r1.duplicated(["pos", "MB", "tlen"], keep="first")
                dups.update(r1.loc[dup_mask, "qname"].tolist())
            # R2: key = str(pos)+UMI+str(tlen)
            r2 = sub[(~sub["is_reverse"]) & (sub["is_read2"])]
            if not r2.empty:
                k = (r2["pos"].astype(str) + r2["MB"] + r2["tlen"].astype(str))
                dup_mask = k.duplicated(keep="first")
                dups.update(r2.loc[dup_mask, "qname"].tolist())
            dup_by_chrom[chrom] = dups
        return dup_by_chrom

    # @@ /*** ---- calc PCR duplicates using chunks ---- ***/
    @staticmethod
    def compute_duplicates_by_chrom_from_chunks(chunks: Iterable[pd.DataFrame]) -> Dict[str, Set[str]]:
        seen_r1: Dict[str, Set[Tuple[int, str, int]]] = {}  # chrom -> set of (pos, UMI, tlen)
        seen_r2: Dict[str, Set[str]] = {}  # chrom -> set of f"{pos}{UMI}{tlen}"
        dup_by_chrom: Dict[str, Set[str]] = {}

        for df in chunks:
            req = {"qname", "MB", "chrom", "pos", "tlen", "is_read1", "is_read2", "is_reverse"}
            miss = req - set(df.columns)
            if miss:
                raise ValueError(f"Chunk missing required columns: {miss}")

            # Only keep the matched records
            df = df[df["chrom"].notna()]
            if df.empty:
                continue

            # R1
            r1 = df[(~df["is_reverse"]) & (df["is_read1"])]
            for qn, chrom, pos, umi, tlen in zip(
                    r1["qname"], r1["chrom"].astype(str), r1["pos"], r1["MB"], r1["tlen"]
            ):
                S1 = seen_r1.setdefault(chrom, set())
                D = dup_by_chrom.setdefault(chrom, set())
                key = (int(pos), str(umi), int(tlen))
                if key in S1:
                    D.add(qn)
                else:
                    S1.add(key)

            # R2
            r2 = df[(~df["is_reverse"]) & (df["is_read2"])]
            for qn, chrom, pos, umi, tlen in zip(
                    r2["qname"], r2["chrom"].astype(str), r2["pos"], r2["MB"], r2["tlen"]
            ):
                S2 = seen_r2.setdefault(chrom, set())
                D = dup_by_chrom.setdefault(chrom, set())
                key = f"{int(pos)}{str(umi)}{int(tlen)}"
                if key in S2:
                    D.add(qn)
                else:
                    S2.add(key)

        # Make sure all occurrences of chrom are collected
        return {chrom: dup_by_chrom.get(chrom, set()) for chrom in set(list(seen_r1.keys()) + list(seen_r2.keys()) + list(dup_by_chrom.keys()))}

    # --------- Mark and merge: write back using the original BAM (keep format/header information) ----------
    @staticmethod
    def _merge_parts(infile: str, refs: List[str]) -> str:
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
            for p in (part, part + ".bai"):
                try: os.remove(p)
                except FileNotFoundError: pass

        out.close()
        ps.index(out_bam)
        return out_bam

    @staticmethod
    def _count_total_and_nondub(bam_path: str) -> Tuple[int, int]:
        total = nondub = 0
        bam = ps.AlignmentFile(bam_path, "rb")
        for r in bam.fetch(until_eof=True):
            total += 1
            if (r.flag & 0x400) == 0:
                nondub += 1
        bam.close()
        return total, nondub

    @staticmethod
    def mark_and_merge(in_bam: str, dup_by_chrom: Dict[str, Set[str]]) -> Tuple[str, int, int, int]:
        bam = ps.AlignmentFile(in_bam, "rb")
        if not bam.has_index():
            ps.index(in_bam)
        refs = list(bam.references)
        bam.close()

        # Chromosome-by-chromosome: only mark the duplicate qnames of this chromosome
        for chrom in refs:
            dup_q = dup_by_chrom.get(chrom, set())
            inb = ps.AlignmentFile(in_bam, "rb")
            out = ps.AlignmentFile(f"{in_bam}.{chrom}.bam", "wb", template=inb)
            for read in inb.fetch(reference=chrom):
                if read.query_name in dup_q:
                    read.flag = read.flag | 0x400
                out.write(read)
            inb.close()
            out.close()

        out_bam = UMItoolsReplica._merge_parts(in_bam, refs)
        total, nondup = UMItoolsReplica._count_total_and_nondub(out_bam)
        return out_bam, total, nondup, total - nondup

    # ---- End-to-end entry (changed to flow by chromosome results) ----
    def run_from_dataframe(self, df_bam: pd.DataFrame, in_bam_path: str):
        dup_by_chrom = self.compute_duplicates_by_chrom_from_df(df_bam)
        return self.mark_and_merge(in_bam_path, dup_by_chrom)

    def run_from_chunks(self, chunk_iter: Iterable[pd.DataFrame], in_bam_path: str):
        dup_by_chrom = self.compute_duplicates_by_chrom_from_chunks(chunk_iter)
        return self.mark_and_merge(in_bam_path, dup_by_chrom)


if __name__ == "__main__":
    from umiche.bam.Reader import ReaderChunk

    bam_fpn = "/mnt/d/Document/Programming/Python/umiche/umiche/data/r1/umitools/umitools.test.RNA-seq.sorted.tagged.bam"

    # 1) 大文件推荐: 分块DataFrame → 去重 → 写回
    # chunk_iter = ReadChunk.iter_chunks(
    #     bam_fpn,
    #     chunk_size=2_000_000,
    #     tag_whitelist=None,            # 例：["CB","MB","XF"]
    #     project_cols=None,             # 默认 ReadChunk.DEFAULT_PROJECT_COLS
    #     categorize=["chrom"]           # 可选
    # )
    # deduper = UMItoolsReplica()
    # OUT_BAM, TOTAL, NONDUP, DUP = deduper.run_from_chunks(chunk_iter, IN_BAM)
    # print(f"[DONE] wrote: {OUT_BAM}")
    # print(f"[DONE] Total={TOTAL}, NonDup={NONDUP}, PCR_Duplicates={DUP}")

    # 2) 小/中等数据：一次性DataFrame作输入
    df_bam = ReaderChunk(
        bam_fpn=bam_fpn,
        bam_fields=None,
        tag_whitelist=['MB'],
        categorize=["chrom"],
    ).todf(chunk_size=2_000_000)
    print(df_bam)
    OUT_BAM, TOTAL, NONDUP, DUP = UMItoolsReplica().run_from_dataframe(df_bam, bam_fpn)
    print(f"[DONE] wrote: {OUT_BAM}")
    print(f"[DONE] Total={TOTAL}, NonDup={NONDUP}, PCR_Duplicates={DUP}")