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

from umiche.util.Console import Console


def hamming_distance(a: str, b: str) -> int:
    if len(a) != len(b):
        return max(len(a), len(b))
    return sum(ch1 != ch2 for ch1, ch2 in zip(a, b))

def _correct_block(umis: List[str], editham: int) -> Dict[str, str]:
    counts = Counter(umis)
    order = sorted(counts.keys(), key=lambda u: (-counts[u], u))
    assigned = {}
    for u in order:
        if u in assigned:
            continue
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
    vals = umi_series.fillna("").astype(str).tolist()
    if ngram_split is None or ngram_split <= 0:
        mapping = _correct_block(vals, editham)
        return pd.Series([mapping.get(u, u) for u in vals], index=umi_series.index)
    else:
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

    def __init__(
            self,
            bam_fpn: str,
            verbose=True,
    ):
        self.bam_fpn = bam_fpn
        import pysam
        self.pysam = pysam

        self.console = Console()
        self.console.verbose = verbose

    # @Console.vignette()
    def extract_spike_dat(
            self,
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
        with self.pysam.AlignmentFile(self.bam_fpn, "rb") as bam:
            references = bam.header.references
            lengths = bam.header.lengths
            ref2len = dict(zip(references, lengths))
            if spikecontig not in ref2len:
                raise ValueError(f"Contig '{spikecontig}' not found in BAM header.")
            reflen = ref2len[spikecontig]

        self.console.print("===>Read data from BAM {}".format(self.bam_fpn))
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
        with self.pysam.AlignmentFile(self.bam_fpn, "rb") as bam:
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

        dat = dat[dat["spikeUMI"].notna()].copy()

        if spikeUMI_length is not None:
            spikeUMI_length = int(spikeUMI_length)
            dat = dat[dat["spikeUMI"].astype(str).str.len() == spikeUMI_length].copy()
            ngram_split = spikeUMI_length // 2
        else:
            ngram_split = None

        self.console.print("===>Hamming correct spikeUMIs...")

        def apply_group(df_g: pd.DataFrame) -> pd.DataFrame:
            g = df_g.copy()
            g["spikeUMI_hd1"] = return_corrected_umi(g["spikeUMI"], editham=1, ngram_split=ngram_split)
            g["spikeUMI_hd2"] = return_corrected_umi(g["spikeUMI"], editham=2, ngram_split=ngram_split)
            return g
        if len(dat) > 0:
            dat = dat.groupby("BC", group_keys=False, sort=False).apply(apply_group)
        cols = [
            "contig",
            "pos",
            "CIGAR",
            "seq",
            "BC",
            "QU",
            "UX",
            "UB",
            "TSSseq",
            "spikeUMI",
            "seqAfterUMI",
            "spikeUMI_hd1",
            "spikeUMI_hd2",
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
        required = {id_col, bc_col}
        missing = required - set(dat.columns)
        if missing:
            raise ValueError(f"Input dat is missing required columns: {missing}")

        grp = dat.groupby(id_col, sort=False)
        spikeoccurance = pd.DataFrame({
            id_col: grp.size().index,
            "N": grp.size().values,
            "nBCs": grp[bc_col].nunique().values
        })

        so1 = spikeoccurance.sort_values("N", kind="mergesort").reset_index(drop=True)
        so1["cs"] = so1["N"].cumsum()

        fig1, ax1 = plt.subplots()
        ax1.scatter(so1["N"], so1["cs"], s=6, rasterized=True)
        ax1.axvline(readcutoff, linestyle="--")
        ax1.set_xlabel("Reads per spike UMI")
        ax1.set_ylabel("Cumulative Reads")
        fig1.tight_layout()

        so2 = spikeoccurance.sort_values("nBCs", kind="mergesort").reset_index(drop=True)
        so2["cs2"] = np.arange(1, len(so2) + 1)

        fig2, ax2 = plt.subplots()
        ax2.scatter(so2["nBCs"], so2["cs2"], s=6, rasterized=True)
        ax2.axvline(nbccutoff, linestyle="--")
        ax2.set_xlabel("Spike in n BCs")
        ax2.set_ylabel("Cumulative Spikes")
        fig2.tight_layout()

        over_read = spikeoccurance.loc[spikeoccurance["N"] > readcutoff, id_col].tolist()
        over_nbcs = spikeoccurance.loc[spikeoccurance["nBCs"] > nbccutoff, id_col].tolist()

        return {
            "plots": [fig1, fig2],
            "over_readcutoff": over_read,
            "over_nbcs": over_nbcs,
            "spikeoccurance": spikeoccurance,
        }


if __name__ == "__main__":
    bam_fpn = "/mnt/d/Document/Programming/Python/umiche/umiche/data/r1/umicountr/Smartseq3.TTACCTGCCAGATTCG.bam"

    p = Preprocessing(
        bam_fpn=bam_fpn
    )
    df = p.extract_spike_dat(
        spikename="g_diySpike4",
        spikecontig="diySpike",
        match_seq_before_UMI="GAGCCTGGGGGAACAGGTAGG",
        match_seq_after_UMI="CTCGGAGGAGAAA"
    )
    print(df)
    print(df.columns)
    print(df.shape)
    print(df.head())

    res = p.get_overrepresented_spikes(df, readcutoff=75, nbccutoff=5)
    # print(res)

    # fig1, fig2 = res["plots"]
    # fig1.savefig("/mnt/d/Document/Programming/Python/umiche/umiche/data/r1/umicountr/1.png", dpi=300, bbox_inches='tight')
    # fig2.savefig("/mnt/d/Document/Programming/Python/umiche/umiche/data/r1/umicountr/2.png", dpi=300, bbox_inches='tight')
    # fig1.show()
    # fig2.show()
    spikeoccurance = res['spikeoccurance']
    print(spikeoccurance)

    overrep_read = res["over_readcutoff"]
    print(overrep_read)
    overrep_bc = res["over_nbcs"]
    print(overrep_bc)