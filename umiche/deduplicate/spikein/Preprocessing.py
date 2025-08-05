__version__ = "v1.0"
__copyright__ = "Copyright 2025"
__license__ = "GPL-3.0"
__developer__ = "Jianfeng Sun"
__maintainer__ = "Jianfeng Sun"
__email__="jianfeng.sunmt@gmail.com"

from typing import Optional, List, Dict, Any

import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
from collections import Counter, defaultdict

import pysam


def hamming_distance(a: str, b: str) -> int:
    if len(a) != len(b):
        return max(len(a), len(b))  # 不同长时给个大值，避免被合并
    return sum(ch1 != ch2 for ch1, ch2 in zip(a, b))

def _correct_block(umis: List[str], editham: int) -> Dict[str, str]:
    """
    基于 '邻接（adjacency）' 的纠错：频数降序依次作为中心，把 HD<=editham 的未指派 UMI 并到该中心。
    返回 mapping: raw_umi -> representative_umi
    """
    counts = Counter(umis)
    # 按频数降序，频数相同按字典序稳定化
    order = sorted(counts.keys(), key=lambda u: (-counts[u], u))
    assigned = {}
    for u in order:
        if u in assigned:
            continue
        # u 成为自己簇的代表
        assigned[u] = u
        for v in order:
            if v in assigned:
                continue
            if len(v) != len(u):
                continue
            if hamming_distance(u, v) <= editham:
                assigned[v] = u
    return assigned

def return_corrected_umi(
    umi_series: pd.Series,
    editham: int = 1,
    ngram_split: Optional[int] = None
) -> pd.Series:
    """
    等价接口：对一组 UMI（同一 BC 内）做 HD<=editham 纠错，返回同序列长度的代表 UMI。
    - ngram_split: 若提供，则在该 BC 内按 UMI 的前缀 umi[:ngram_split] 分桶后分别纠错（加速）。
      （注意：这与 UMIcountR 的 ngram 阻塞思想一致，但实现细节可能与其完全一致实现有细微差异。）
    """
    vals = umi_series.fillna("").astype(str).tolist()
    if ngram_split is None or ngram_split <= 0:
        mapping = _correct_block(vals, editham)
        return pd.Series([mapping.get(u, u) for u in vals], index=umi_series.index)
    else:
        # 以前缀分桶
        buckets = defaultdict(list)
        for idx, u in zip(umi_series.index, vals):
            key = u[:ngram_split] if len(u) >= ngram_split else u
            buckets[key].append((idx, u))
        out = pd.Series(index=umi_series.index, dtype=object)
        for key, items in buckets.items():
            idxs = [i for i, _ in items]
            umis = [u for _, u in items]
            mapping = _correct_block(umis, editham)
            out.loc[idxs] = [mapping.get(u, u) for u in umis]
        return out

class Preprocessing:

    def __init__(self, ):
        pass

    def extract_spike_dat(
            self,
            bam_path: str,
            spikecontig: str = "diySpike",
            spikename: str = "g_diySpike4",
            match_seq_before_UMI: Optional[str] = "GAGCCTGGGGGAACAGGTAGG",
            match_seq_after_UMI: Optional[str] = "CTCGGAGGAGAAA",
            spikeUMI_start: Optional[int] = None,
            spikeUMI_end: Optional[int] = None,
            fixed_start_pos: Optional[int] = None,
            spikeUMI_length: Optional[int] = None,
    ) -> pd.DataFrame:
        """

        Parameters
        ----------
        bam_path
        spikecontig
            read reference==spikecontig
        spikename
            GE==spikename
        match_seq_before_UMI
        match_seq_after_UMI
        spikeUMI_start
        spikeUMI_end
        fixed_start_pos
        spikeUMI_length

        Returns
        -------
            contig:
                {rname,pos,cigar,seq}
            pos:
                {rname,pos,cigar,seq}
            CIGAR:
                {rname,pos,cigar,seq}
            seq:
                {rname,pos,cigar,seq}
            BC:
                {BC,QU,UX,UB}
            QU:
                {BC,QU,UX,UB}
            UX:
                {BC,QU,UX,UB}
                UX==""（Smart-seq3 built-in stuff）
            UB:
                {BC,QU,UX,UB}
            TSSseq:
            spikeUMI:
            seqAfterUMI:
            spikeUMI_hd1:
                grouped by BC
            spikeUMI_hd2:
                grouped by BC
        """
        # ## @@ /*** obtain spikecontig ***/
        with pysam.AlignmentFile(bam_path, "rb") as bam:
            references = bam.header.references
            lengths = bam.header.lengths
            ref2len = dict(zip(references, lengths))
            if spikecontig not in ref2len:
                raise ValueError(f"Contig '{spikecontig}' not found in BAM header.")
            reflen = ref2len[spikecontig]

        print("Reading in data from bam file...")
        rows = {
            "contig": [],
            "pos": [],
            "CIGAR": [],
            "seq": [],
            "BC": [],
            "QU": [],
            "UX": [],
            "UB": [],
        }
        with pysam.AlignmentFile(bam_path, "rb") as bam:
            for aln in bam.fetch(spikecontig, 0, reflen):
                try:
                    ge = aln.get_tag("GE")
                except KeyError:
                    continue
                if ge != spikename:
                    continue

                def _get_tag_safe(a, tag):
                    try:
                        v = a.get_tag(tag)
                        return str(v)
                    except KeyError:
                        return ""

                rows["contig"].append(aln.reference_name)
                # pos 为 1-based，pysam reference_start 为 0-based
                rows["pos"].append(aln.reference_start + 1)
                rows["CIGAR"].append(aln.cigarstring or "")
                rows["seq"].append(aln.query_sequence or "")

                rows["BC"].append(_get_tag_safe(aln, "BC"))
                rows["QU"].append(_get_tag_safe(aln, "QU"))
                rows["UX"].append(_get_tag_safe(aln, "UX"))
                rows["UB"].append(_get_tag_safe(aln, "UB"))

        dat = pd.DataFrame(rows)
        dat = dat[dat["UX"] != ""].copy()
        if fixed_start_pos is not None:
            dat = dat[dat["pos"] == int(fixed_start_pos)].copy()

        if (spikeUMI_start is None) and (match_seq_before_UMI is None):
            raise ValueError("give either spikeUMI_start or match_seq_before_UMI")

        # ## @@ /*** obtain spUMI ***/
        if (match_seq_before_UMI is not None) and (match_seq_after_UMI is not None):
            # ## @@ /*** mode 1: chars ***/
            s = dat["seq"].astype(str)
            parts1 = s.str.split(match_seq_before_UMI, n=1, expand=True)
            if parts1.shape[1] == 1:
                parts1 = pd.concat([parts1, pd.Series([np.nan] * len(parts1), index=parts1.index)], axis=1)
            dat["TSSseq"] = parts1[0]
            tmpseq = parts1[1]
            mask = tmpseq.notna()
            dat.loc[~mask, "spikeUMI"] = np.nan
            dat.loc[~mask, "seqAfterUMI"] = np.nan
            if mask.any():
                parts2 = tmpseq[mask].astype(str).str.split(match_seq_after_UMI, n=1, expand=True)
                if parts2.shape[1] == 1:
                    parts2 = pd.concat([parts2, pd.Series([np.nan] * len(parts2), index=parts2.index)], axis=1)
                dat.loc[mask, "spikeUMI"] = parts2[0]
                dat.loc[mask, "seqAfterUMI"] = parts2[1]
            # 补回被切掉的锚定序列
            dat["TSSseq"] = dat["TSSseq"].fillna("") + (match_seq_before_UMI or "")
            dat["seqAfterUMI"] = (match_seq_after_UMI or "") + dat["seqAfterUMI"].fillna("")

        else:
            # ## @@ /*** mode 2: pos ***/
            if (fixed_start_pos is None) and (spikeUMI_start is not None):
                print("Warning: using spike extraction by position in read without fixed read start position!")
            s = dat["seq"].astype(str)
            start = int(spikeUMI_start) if spikeUMI_start is not None else 1
            end = int(spikeUMI_end) if spikeUMI_end is not None else start - 1
            def substr1(x: str, i: int, j: int) -> str:
                n = len(x)
                i = max(1, i);
                j = max(1, j)
                if i > n or i > j:
                    return ""
                j = min(j, n)
                return x[i - 1:j]
            dat["TSSseq"] = [substr1(x, 1, start - 1) for x in s]
            dat["spikeUMI"] = [substr1(x, start, end) for x in s]
            dat["seqAfterUMI"] = [substr1(x, end + 1, 48) for x in s]

        # R: dat <- dat[!is.na(spikeUMI)]
        dat = dat[dat["spikeUMI"].notna()].copy()

        if spikeUMI_length is not None:
            spikeUMI_length = int(spikeUMI_length)
            dat = dat[dat["spikeUMI"].astype(str).str.len() == spikeUMI_length].copy()
            ngram_split = spikeUMI_length // 2
        else:
            ngram_split = None

        print("Hamming correct spikeUMIs...")
        # 按 BC 分组，分别生成 hd1 / hd2
        def apply_group(df_g: pd.DataFrame) -> pd.DataFrame:
            g = df_g.copy()
            g["spikeUMI_hd1"] = return_corrected_umi(g["spikeUMI"], editham=1, ngram_split=ngram_split)
            g["spikeUMI_hd2"] = return_corrected_umi(g["spikeUMI"], editham=2, ngram_split=ngram_split)
            return g
        if len(dat) > 0:
            dat = dat.groupby("BC", group_keys=False, sort=False).apply(apply_group)
        cols = [
            "contig", "pos", "CIGAR", "seq", "BC", "QU", "UX", "UB",
            "TSSseq", "spikeUMI", "seqAfterUMI", "spikeUMI_hd1", "spikeUMI_hd2"
        ]
        for c in cols:
            if c not in dat.columns:
                dat[c] = np.nan
        dat = dat[cols].reset_index(drop=True)
        return dat

    def get_overrepresented_spikes(
            self,
            dat: pd.DataFrame,
            readcutoff: int = 100,
            nbccutoff: int = 5,
            id_col: str = "spikeUMI_hd2",
            bc_col: str = "BC",
    ) -> Dict[str, Any]:
        """
        Python re-implementation of the R function `get_overrepresented_spikes`.

        Parameters
        ----------
        dat : pd.DataFrame
            输入表，必须至少包含 [id_col, bc_col] 两列；每一行代表一个原始 read（与 R: .N 语义一致）。
        readcutoff : int, default 100
            单个 spUMI 出现的原始 reads 上限阈值；超过者视为 overrepresented（与 R: N > readcutoff 一致）。
        nbccutoff : int, default 5
            单个 spUMI 出现的细胞条形码(BC)个数上限阈值；超过者视为 overrepresented（与 R: nBCs > nbccutoff 一致）。
        id_col : str, default "spikeUMI_hd2"
            spUMI 标识列（R 里按 spikeUMI_hd2 分组）。
        bc_col : str, default "BC"
            细胞条形码列（R 里 length(unique(BC))）。

        Returns
        -------
        result : dict
            {
              "plots": [fig1, fig2],               # 两张 matplotlib Figure：与 R 的 p1/p2 对应
              "over_readcutoff": List[str],        # N > readcutoff 的 spUMI 标识
              "over_nbcs": List[str],              # nBCs > nbccutoff 的 spUMI 标识
              "spikeoccurance": pd.DataFrame       # （可选）返回中间表，便于调试/二次作图
            }
        """
        required = {id_col, bc_col}
        missing = required - set(dat.columns)
        if missing:
            raise ValueError(f"Input dat is missing required columns: {missing}")

        # --- 统计：每个 spUMI 的读段数 N 与出现的 BC 个数 nBCs（与 R: .N, length(unique(BC)) 一致）
        grp = dat.groupby(id_col, sort=False)
        spikeoccurance = pd.DataFrame({
            id_col: grp.size().index,
            "N": grp.size().values,
            "nBCs": grp[bc_col].nunique().values
        })

        # ---------- 图1：按 N 升序的累计 reads ----------
        so1 = spikeoccurance.sort_values("N", kind="mergesort").reset_index(drop=True)  # 稳定排序，贴近 data.table::setorder
        so1["cs"] = so1["N"].cumsum()

        fig1, ax1 = plt.subplots()
        ax1.scatter(so1["N"], so1["cs"], s=6, rasterized=True)
        ax1.axvline(readcutoff, linestyle="--")
        ax1.set_xlabel("Reads per spike UMI")
        ax1.set_ylabel("Cumulative Reads")
        fig1.tight_layout()

        # ---------- 图2：按 nBCs 升序的累计 spikes ----------
        so2 = spikeoccurance.sort_values("nBCs", kind="mergesort").reset_index(drop=True)
        so2["cs2"] = np.arange(1, len(so2) + 1)

        fig2, ax2 = plt.subplots()
        ax2.scatter(so2["nBCs"], so2["cs2"], s=6, rasterized=True)
        ax2.axvline(nbccutoff, linestyle="--")
        ax2.set_xlabel("Spike in n BCs")
        ax2.set_ylabel("Cumulative Spikes")
        fig2.tight_layout()

        # ---------- 阈值筛选（与 R: N > readcutoff, nBCs > nbccutoff 一致） ----------
        over_read = spikeoccurance.loc[spikeoccurance["N"] > readcutoff, id_col].tolist()
        over_nbcs = spikeoccurance.loc[spikeoccurance["nBCs"] > nbccutoff, id_col].tolist()

        # 结果结构与 R 对齐：plots 列表 + 两个 ID 列表
        return {
            "plots": [fig1, fig2],
            "over_readcutoff": over_read,
            "over_nbcs": over_nbcs,
            "spikeoccurance": spikeoccurance,  # 方便你复用数据；不需要可删除
        }


if __name__ == "__main__":
    p = Preprocessing()
    df = p.extract_spike_dat(
        bam_path="/mnt/d/Document/Programming/Python/mcverse/mcverse/data/umiche/spike-in/Smartseq3.TTACCTGCCAGATTCG.bam",
        spikename="g_diySpike4",
        spikecontig="diySpike",
        match_seq_before_UMI="GAGCCTGGGGGAACAGGTAGG",
        match_seq_after_UMI="CTCGGAGGAGAAA"
    )
    print(df)
    print(df.shape)
    print(df.head())

    res = p.get_overrepresented_spikes(df, readcutoff=75, nbccutoff=5)
    fig1, fig2 = res["plots"]
    fig1.savefig("1.png", dpi=300, bbox_inches='tight')
    fig2.savefig("2.png", dpi=300, bbox_inches='tight')
    fig1.show()
    fig2.show()

    overrep_read = res["over_readcutoff"]
    overrep_bc = res["over_nbcs"]