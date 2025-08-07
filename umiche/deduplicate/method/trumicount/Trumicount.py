#!/usr/bin/env python3
# -*- coding: utf-8 -*-

import sys, math, statistics, gzip, bz2
from collections import defaultdict
from typing import Dict, List, Tuple, Iterable, Any, Optional

import click

try:
    import pysam
except ImportError as e:
    raise SystemExit("This script requires pysam. Install with `pip install pysam` or conda.") from e

VERSION = "0.3.0-no-umitools"

# ------------------------ I/O helpers ------------------------
def open_byext(path: str, mode: str = "rt"):
    if path.endswith((".gz", ".gzip")):
        return gzip.open(path, mode)
    if path.endswith((".bz2", ".bzip2")):
        return bz2.open(path, mode)
    return open(path, mode)

def write_tsv(path: str, header: List[str], rows: Iterable[Iterable[Any]]):
    fh = open_byext(path, "wt")
    fh.write("\t".join(header) + "\n")
    for r in rows:
        fh.write("\t".join(str(x) for x in r) + "\n")
    fh.close()

# ------------------------ TRUmiCount core ------------------------
class TRUmiCountPy:
    """
    TRUmiCount pipeline without umi_tools:
      - Read CB/UB (and optionally GX) tags directly from BAM
      - Group by (cell,gene) or (cell,contig) and UMI
      - Aggregate read counts per UMI (exact UMI only; no UMI error-correction merge here)
      - Estimate global/group models (hooks unchanged)
      - Write counts + assigned (raw_umi -> final_umi identity)

    Notes:
      * If --paired is set, we count only read1 to avoid double-counting mates.
      * Secondary/supplementary alignments are skipped.
      * If --per-contig is used (or gene tag missing), 'gene' column will hold reference name.
    """
    def __init__(
        self,
        molecules: int = 1,
        threshold: Optional[int] = None,
        threshold_quantile: float = 0.05,
        genewise_min_umis: int = 3,
        mapping_quality: int = 0,
        paired: bool = False,
        gene_tag: Optional[str] = "GX",
        per_contig: bool = False,
        cell_tag: Optional[str] = "CB",
        umi_tag: Optional[str] = "UB",
        verbose: bool = False,
    ):
        self.molecules = molecules
        self.threshold = threshold
        self.threshold_quantile = threshold_quantile
        self.genewise_min_umis = genewise_min_umis
        self.mapping_quality = mapping_quality
        self.paired = paired
        self.gene_tag = gene_tag
        self.per_contig = per_contig
        self.cell_tag = cell_tag
        self.umi_tag = umi_tag
        self.verbose = verbose

    # ---------- iterate BAM → group-like records ----------
    def _iter_bam_tag_records(self, bam_path: str) -> Iterable[Dict[str, Any]]:
        """
        Yield records mimicking umi_tools `group --group-out` rows:
          { "cell_id": <CB>, "gene": <GX or contig>, "final_umi": <UB>, "umi": <UB>, "umi_count": 1 }

        Filtering:
          - skip unmapped, secondary, supplementary
          - mapping_quality >= self.mapping_quality
          - if paired=True: count only read1 (is_read1), skip read2 (to avoid double counting)

        If per_contig=True or gene_tag is missing/unset, gene field is filled with reference_name.
        """
        af = pysam.AlignmentFile(bam_path, "rb")
        fetch_iter = af.fetch(until_eof=True)  # works for sorted BAM as well

        missing = 0
        produced = 0
        for aln in fetch_iter:
            if aln.is_unmapped or aln.is_secondary or aln.is_supplementary:
                continue
            if aln.mapping_quality < self.mapping_quality:
                continue
            if self.paired:
                # count only read1 to avoid double counting mates
                if not aln.is_read1:
                    continue

            # tags
            cell = None
            umi  = None
            try:
                if self.cell_tag:
                    cell = aln.get_tag(self.cell_tag)
                if self.umi_tag:
                    umi = aln.get_tag(self.umi_tag)
            except KeyError:
                pass

            if not cell or not umi:
                missing += 1
                continue

            # gene or contig
            gene_val = None
            if not self.per_contig and self.gene_tag:
                try:
                    gene_val = aln.get_tag(self.gene_tag)
                except KeyError:
                    gene_val = None
            if not gene_val:
                # fallback to reference name (per-contig mode)
                gene_val = af.get_reference_name(aln.reference_id)

            yield {
                "cell_id": cell,
                "gene": gene_val,
                "final_umi": umi,   # no correction here
                "umi": umi,
                "umi_count": 1      # one read contributes 1 line; aggregator will sum
            }
            produced += 1

        if self.verbose:
            sys.stderr.write(f"[trumicount] iter_bam: produced={produced}, skipped_missing_tags={missing}\n")

    # ---------- aggregate to UMI-level + build assigned ----------
    @staticmethod
    def _aggregate_to_umis(group_recs):
        """
        输入：_iter_bam_tag_records 的逐 read 记录（每条 umi_count=1）。
        逻辑与之前一致：同 (cell,gene,final_umi) 的 reads = max(umi_count最大值, 行数)
        同时构建 raw_umi→final_umi 的 assigned（此处为恒等映射）。
        """
        per_final_max = {}                 # (cell,gene,final) -> max umi_count
        per_final_sum = defaultdict(int)   # (cell,gene,final) -> #lines
        per_mapping   = defaultdict(int)   # (cell,gene,raw,final) -> #lines

        for rec in group_recs:
            cell  = rec.get("cell_id", "")
            gene  = rec.get("gene", "NA")
            final = rec.get("final_umi", "")
            raw   = rec.get("umi", final)
            cnt   = int(rec.get("umi_count", 1))

            kf = (cell, gene, final)
            per_final_sum[kf] += 1
            if kf not in per_final_max or cnt > per_final_max[kf]:
                per_final_max[kf] = cnt

            per_mapping[(cell, gene, raw, final)] += 1

        umi_rows = []
        all_keys = set(per_final_sum) | set(per_final_max)
        for k in all_keys:
            reads = max(per_final_sum.get(k, 0), per_final_max.get(k, 0))
            umi_rows.append({"cell": k[0], "gene": k[1], "final_umi": k[2], "reads": reads})

        assigned_rows = [(k[0], k[1], k[2], k[3], v) for k, v in per_mapping.items()]
        return umi_rows, assigned_rows

    # ---------- threshold auto ----------
    @staticmethod
    def _auto_threshold(counts: List[int], q: float) -> int:
        if not counts: return 0
        s = sorted(counts)
        k = max(0, min(len(s) - 1, int(q * (len(s) - 1))))
        return int(s[k])

    # ---------- model: Poisson-depth × efficiency (approx) ----------
    @staticmethod
    def _loss_prob(threshold: int, eff: float, depth: float) -> float:
        import math
        lam = max(1e-9, depth * eff)
        p = 0.0
        for k in range(0, threshold):
            p += math.exp(-lam) * (lam ** k) / math.factorial(k)
        return min(max(p, 0.0), 1.0)

    def estimate_global_model(self, read_counts: List[int], threshold: int) -> Dict[str, float]:
        if not read_counts:
            return dict(efficiency=0.5, depth=1.0, loss=0.5)
        m = statistics.mean(read_counts)
        depth = max(1e-6, m)
        frac_below = sum(c < threshold for c in read_counts) / max(1, len(read_counts))
        best_e, best_diff = 0.5, 1e9
        for e in [i / 100 for i in range(5, 100)]:
            d = abs(self._loss_prob(threshold, e, depth) - frac_below)
            if d < best_diff:
                best_e, best_diff = e, d
        return dict(
            efficiency=best_e,
            depth=depth,
            loss=self._loss_prob(threshold, best_e, depth)
        )

    def estimate_group_models(
        self,
        group_counts: Dict[Tuple[str,str], List[int]],
        threshold: int,
        global_model: Dict[str, float],
    ) -> Dict[Tuple[str,str], Dict[str, Any]]:
        out = {}
        for gid, counts in group_counts.items():
            n_obs = len(counts)
            if n_obs == 0:
                eff, dep, loss = global_model["efficiency"], global_model["depth"], global_model["loss"]
            else:
                m = self.estimate_global_model(counts, threshold)
                eff, dep, loss = m["efficiency"], m["depth"], m["loss"]
            n_umis = n_obs
            n_tot = n_umis / max(1e-8, 1.0 - loss)  # bias correction
            out[gid] = dict(
                n_umis=n_umis, n_tot=n_tot,
                efficiency=eff, depth=dep, loss=loss, n_obs=n_obs
            )
        return out

    # ---------- main run ----------
    def run_from_bam(
        self,
        bam_path: str,
        output_counts: Optional[str] = None,
        output_assigned: Optional[str] = None,
    ):
        # 1) 直接从 BAM tags 读取记录
        recs = self._iter_bam_tag_records(bam_path)

        # 2) 聚合到 UMI 粒度 & 构建 assigned
        umi_rows, assigned_rows = self._aggregate_to_umis(recs)

        # 3) 阈值
        if self.threshold is None:
            thr = self._auto_threshold([r["reads"] for r in umi_rows], self.threshold_quantile)
            self.threshold = thr
            if self.verbose:
                print(f"[trumicount] autothreshold={self.threshold}", file=sys.stderr)

        # 4) 过滤低 reads/UMI
        if self.threshold > 0:
            umi_rows = [r for r in umi_rows if r["reads"] >= self.threshold]

        # 5) 分组到 (cell,gene)
        group_counts: Dict[Tuple[str,str], List[int]] = defaultdict(list)
        for r in umi_rows:
            key = (r["cell"], r["gene"])
            group_counts[key].append(int(r["reads"]))

        # 6) 全局模型 + 分组模型
        global_model = self.estimate_global_model([r for counts in group_counts.values() for r in counts],
                                                  self.threshold)
        group_models = self.estimate_group_models(group_counts, self.threshold, global_model)

        # 7) 输出 counts
        rows = []
        for (cell, gene), cnts in group_counts.items():
            m = group_models[(cell, gene)]
            rows.append([
                cell, gene,
                int(m["n_umis"]),
                round(m["n_tot"], 3),
                round(m["efficiency"], 3),
                round(m["depth"], 3),
                round(m["loss"], 3),
                int(m["n_obs"]),
            ])
        header = ["cell", "gene", "n.umis", "n.tot", "efficiency", "depth", "loss", "n.obs"]
        if output_counts:
            write_tsv(output_counts, header, rows)
        else:
            write_tsv("/dev/stdout", header, rows)

        # 8) 输出 assigned（仅保留阈值过滤后仍保留的 final_umi）
        if output_assigned:
            kept = {(r["cell"], r["gene"], r["final_umi"]) for r in umi_rows}
            hdr = ["cell", "gene", "raw_umi", "final_umi", "reads"]
            outrows = []
            for cell, gene, raw, final, nr in assigned_rows:
                if (cell, gene, final) in kept:
                    outrows.append([cell, gene, raw, final, nr])
            write_tsv(output_assigned, hdr, outrows)

        return rows

# ------------------------ CLI ------------------------
@click.command()
@click.option("--input-bam", required=True, type=click.Path(exists=True), help="Input BAM with CB/UB (and optionally GX) tags.")
@click.option("--molecules", default=1, show_default=True, type=int, help="Assumed molecules per UMI (compatibility only).")
@click.option("--threshold", default=None, type=int, help="Reads/UMI threshold. If omitted, use --threshold-quantile.")
@click.option("--threshold-quantile", default=0.05, show_default=True, type=float, help="Auto threshold quantile when --threshold not set.")
@click.option("--genewise-min-umis", default=3, show_default=True, type=int, help="Minimum UMIs per (cell,gene) for group-wise model (kept for compatibility).")
@click.option("--output-counts", required=True, type=click.Path(), help="Output counts table.")
@click.option("--output-assigned", default=None, type=click.Path(), help="Output mapping raw_umi→final_umi (optional).")
@click.option("--mapping-quality", default=0, show_default=True, type=int, help="Minimum MAPQ to keep a read.")
@click.option("--paired", is_flag=True, help="If set, only count read1 per template to avoid double-counting.")
@click.option("--per-contig", is_flag=True, help="Use reference contig as 'gene' when gene_tag missing/unset.")
@click.option("--gene-tag", default="GX", show_default=True, help="BAM tag for gene id (e.g. GX). Ignored when --per-contig.")
@click.option("--cell-tag", default="CB", show_default=True, help="BAM tag for cell barcode (e.g. CB).")
@click.option("--umi-tag",  default="UB", show_default=True, help="BAM tag for UMI (e.g. UB).")
@click.option("--verbose", is_flag=True, help="Verbose logs.")
def main(
    input_bam, molecules, threshold, threshold_quantile, genewise_min_umis,
    output_counts, output_assigned, mapping_quality, paired, per_contig, gene_tag, cell_tag, umi_tag, verbose
):
    runner = TRUmiCountPy(
        molecules=molecules,
        threshold=threshold,
        threshold_quantile=threshold_quantile,
        genewise_min_umis=genewise_min_umis,
        mapping_quality=mapping_quality,
        paired=paired,
        gene_tag=gene_tag,
        per_contig=per_contig,
        cell_tag=cell_tag,
        umi_tag=umi_tag,
        verbose=verbose,
    )
    runner.run_from_bam(
        bam_path=input_bam,
        output_counts=output_counts,
        output_assigned=output_assigned,
    )

if __name__ == "__main__":
    main()
