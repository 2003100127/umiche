__version__ = "v1.0"
__copyright__ = "Copyright 2025"
__license__ = "GPL-3.0"
__developer__ = "Jianfeng Sun"
__maintainer__ = "Jianfeng Sun"
__email__="jianfeng.sunmt@gmail.com"

from dataclasses import dataclass
from typing import List, Optional
import math

import pysam
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt

from umiche.deduplicate.method.ReformKit import ReformKit

@dataclass
class OverRepResult:
    plots: List[plt.Figure]
    over_readcutoff: List[str]
    over_nbcs: List[str]
    stats_df: pd.DataFrame  # for transparency / debugging


class UMIcountPy:
    """
    - extract_spike_dat(...)
    - get_overrepresented_spikes(...)
    - subsample_recompute(...)

    """

    # @@ /*** -------------------- ***/
    # UMI collapse / correction
    # @@ /*** -------------------- ***/

    # @@ /*** -------------------- ***/
    # Extract spike reads & spUMIs
    # @@ /*** -------------------- ***/
    @staticmethod
    def extract_spike_dat(
        df: pd.DataFrame,
        spikeUMI_start: Optional[int] = None,
        spikeUMI_end: Optional[int] = None,
        fixed_start_pos: Optional[int] = None,
        match_seq_before_UMI: Optional[str] = None,
        match_seq_after_UMI: Optional[str] = None,
        spikeUMI_length: Optional[int] = None,
    ) -> pd.DataFrame:
        """
        Load reads on spike contig & gene, parse spUMI by flanking motifs or by fixed positions,
        and compute HD1/HD2-corrected spUMIs (grouped within BC) — columns match R output:
        [contig, pos, CIGAR, seq, BC, QU, UX, UB, TSSseq, spikeUMI, seqAfterUMI, spikeUMI_hd1, spikeUMI_hd2].

        """
        # Extract spUMI
        TSSseq = []
        spikeUMI = []
        seqAfter = []
        if match_seq_before_UMI and match_seq_after_UMI:
            pre = match_seq_before_UMI
            aft = match_seq_after_UMI
            for s in df["seq"].tolist():
                i = s.find(pre)
                if i < 0:
                    TSSseq.append(None)
                    spikeUMI.append(None)
                    seqAfter.append(None)
                    continue
                left = s[:i]
                rem = s[i + len(pre) :]
                j = rem.find(aft)
                if j < 0:
                    TSSseq.append(None)
                    spikeUMI.append(None)
                    seqAfter.append(None)
                    continue
                umi = rem[:j]
                right = rem[j + len(aft) :]
                # include the matched motifs into TSSseq & seqAfterUMI
                TSSseq.append(left + pre)
                seqAfter.append(aft + right)
                spikeUMI.append(umi)
        else:
            # fixed positions mode
            if spikeUMI_start is None:
                raise ValueError("give either spikeUMI_start or match_seq_before_UMI")
            if spikeUMI_end is None:
                raise ValueError("spikeUMI_end must be provided with spikeUMI_start")
            for s in df["seq"].tolist():
                # R uses 1-based inclusive ranges; Python slice is 0-based, end-excl.
                tss = s[: spikeUMI_start - 1]
                umi = s[spikeUMI_start - 1 : spikeUMI_end]
                # seqAfterUMI
                seq_after = s[spikeUMI_end : 48] if len(s) >= 48 else s[spikeUMI_end :]
                TSSseq.append(tss)
                spikeUMI.append(umi)
                seqAfter.append(seq_after)

        df["TSSseq"] = TSSseq
        df["spikeUMI"] = spikeUMI
        df["seqAfterUMI"] = seqAfter

        # Drop rows with no spikeUMI
        df = df[df["spikeUMI"].notna()].copy()

        # Optional strict length filter & ngram_split parity
        if spikeUMI_length is not None:
            df = df[df["spikeUMI"].str.len() == int(spikeUMI_length)].copy()

        # Compute HD1 / HD2 corrected spUMIs grouped by BC
        def _corr(group, ed):
            return ReformKit().return_corrected_umi(group["spikeUMI"].tolist(), ed_thres=ed, collapse_mode="adjacency")

        df["spikeUMI_hd1"] = (
            df.groupby("BC", group_keys=False)["spikeUMI"]
            .apply(lambda s: pd.Series(ReformKit().return_corrected_umi(s.tolist(), ed_thres=1, collapse_mode="adjacency")))
            .values
        )
        df["spikeUMI_hd2"] = (
            df.groupby("BC", group_keys=False)["spikeUMI"]
            .apply(lambda s: pd.Series(ReformKit().return_corrected_umi(s.tolist(), ed_thres=2, collapse_mode="adjacency")))
            .values
        )
        return df.reset_index(drop=True)

    # @@ /*** -------------------- ***/
    # Over-represented spUMIs
    # @@ /*** -------------------- ***/
    @staticmethod
    def get_overrepresented_spikes(
        dat: pd.DataFrame, readcutoff: int = 100, nbccutoff: int = 5
    ) -> OverRepResult:
        """
        Find overrepresented spUMIs by read count and #BCs; also return plots:
        plots[0] == cumulative reads vs reads-per-spUMI (with dashed vline),
        plots[1] == cumulative spikes vs #BCs (with dashed vline).
        Mirrors UMIcountR::get_overrepresented_spikes.
        """
        # data.table: .N and unique(BC) by spikeUMI_hd2
        grp = dat.groupby("spikeUMI_hd2")
        N = grp.size().rename("N")
        nBCs = grp["BC"].nunique().rename("nBCs")
        stats = pd.concat([N, nBCs], axis=1).reset_index().rename(columns={"index": "spikeUMI_hd2"})

        # Plot 1: cumulative reads sorted by N
        s1 = stats.sort_values("N").copy()
        s1["cs"] = s1["N"].cumsum()

        fig1, ax1 = plt.subplots(figsize=(6, 4))
        ax1.plot(s1["N"], s1["cs"], ".", markersize=2)
        ax1.axvline(readcutoff, linestyle="--")
        ax1.set_xlabel("Reads per spike UMI")
        ax1.set_ylabel("Cumulative Reads")
        ax1.set_title("Cumulative reads by spUMI abundance")

        # Plot 2: cumulative spikes over #BCs
        s2 = stats.sort_values("nBCs").copy()
        s2["cs2"] = np.arange(1, len(s2) + 1)

        fig2, ax2 = plt.subplots(figsize=(6, 4))
        ax2.plot(s2["nBCs"], s2["cs2"], ".", markersize=2)
        ax2.axvline(nbccutoff, linestyle="--")
        ax2.set_xlabel("Spike in n BCs")
        ax2.set_ylabel("Cumulative Spikes")
        ax2.set_title("Cumulative spikes by #BCs")

        overrep1 = stats.loc[stats["N"] > readcutoff, "spikeUMI_hd2"].tolist()
        overrep2 = stats.loc[stats["nBCs"] > nbccutoff, "spikeUMI_hd2"].tolist()

        return OverRepResult(
            plots=[fig1, fig2],
            over_readcutoff=overrep1,
            over_nbcs=overrep2,
            stats_df=stats,
        )

    # ------------------------------
    # 4) Subsample & recompute UMIs
    # ------------------------------
    @staticmethod
    def subsample_recompute(
            dat: pd.DataFrame,
            mu_nSpikeUMI: int,
    ) -> pd.DataFrame:
        """
        For each BC, subsample the number of unique spikeUMI_hd2 to ~N(mu, sd=sqrt(mu)),
        then recompute UB (library UMI) HD1/HD2 ReformKit() return_corrected_umi within each BC.
        Mirrors UMIcountR::subsample_recompute.
        """
        # Keep barcodes with enough spikes (>= 0.8 * mu)
        nspikes_per_bc = (
            dat.groupby("BC")["spikeUMI_hd2"]
            .nunique()
            .rename("nspikes")
            .reset_index()
        )
        ok_bcs = nspikes_per_bc.loc[nspikes_per_bc["nspikes"] >= 0.8 * mu_nSpikeUMI, "BC"].tolist()
        df = dat[dat["BC"].isin(ok_bcs)][["BC", "spikeUMI_hd2", "UX"]].copy()

        out_rows = []
        for bc, sub in df.groupby("BC"):
            uniq_spikes = sub["spikeUMI_hd2"].unique().tolist()
            # draw target n from |N(mu, sqrt(mu))|, at least 1
            target = int(abs(round(np.random.normal(loc=mu_nSpikeUMI, scale=math.sqrt(mu_nSpikeUMI)))))
            target = max(1, target)

            if target < len(uniq_spikes):
                pick = set(np.random.choice(uniq_spikes, size=target, replace=False).tolist())
                sub2 = sub[sub["spikeUMI_hd2"].isin(pick)].copy()
            else:
                sub2 = sub.copy()

            # Recompute UB correction per BC
            ub_hd1 = ReformKit().return_corrected_umi(sub2["UX"].tolist(), ed_thres=1, collapse_mode="adjacency")
            ub_hd2 = ReformKit().return_corrected_umi(sub2["UX"].tolist(), ed_thres=2, collapse_mode="adjacency")
            sub2["UB_hd1"] = ub_hd1
            sub2["UB_hd2"] = ub_hd2
            out_rows.append(sub2)

        if not out_rows:
            return pd.DataFrame(columns=["BC", "spikeUMI_hd2", "UX", "UB_hd1", "UB_hd2", "mean_nSpikeUMI"])

        out = pd.concat(out_rows, axis=0, ignore_index=True)
        out["mean_nSpikeUMI"] = mu_nSpikeUMI
        return out


if __name__ == "__main__":
    bam_fpn = "/mnt/d/Document/Programming/Python/umiche/umiche/data/r1/umicountr/Smartseq3.TTACCTGCCAGATTCG.bam"

    from umiche.bam.Reader import ReaderChunk
    df = ReaderChunk(
        bam_fpn=bam_fpn,
        bam_fields=['contig', 'pos', 'CIGAR', 'seq'],
        tag_whitelist=['BC', 'QU', 'UX', 'UB', 'GE'],
        verbose=True,
    ).todf(chunk_size=1_000_000)
    print(df)
    spikecontig = "diySpike"
    spikename = "g_diySpike4"
    df = df[df["UX"] != ""]
    df = df[df["GE"] == spikename]
    df = df.reset_index(drop=True)
    df = df.drop("GE", axis=1)
    print(df)

    # 1) extract spike data from BAM
    df_sp = UMIcountPy.extract_spike_dat(
        df=df,
        match_seq_before_UMI="GAGCCTGGGGGAACAGGTAGG",
        match_seq_after_UMI="CTCGGAGGAGAAA",
    )
    # print(df_sp)
    # print(df_sp.columns)
    # contig, pos, CIGAR, seq, BC, QU, UX, UB, TSSseq, spikeUMI, seqAfterUMI, spikeUMI_hd1, spikeUMI_hd2。:contentReference[oaicite:5]{index=5}

    # 2) overrepr spUMI filter plot
    # overrep = UMIcountPy.get_overrepresented_spikes(spikedat, readcutoff=75)
    # fig = overrep.plots[0]  # cumulative reads 曲线图
    # fig.savefig("/mnt/d/Document/Programming/Python/umiche/umiche/data/r1/umicountr/cumulative.png", dpi=300, bbox_inches='tight')
    # fig.show()


    # 3) directional-adjacency for Smart-seq3 UMI
    df_sp["UB_directional"] = (
        df_sp.groupby("BC")["UX"]
        .transform(lambda UX: pd.Series(
            ReformKit().return_corrected_umi(
                UX.tolist(), ed_thres=1, collapse_mode="adjacency_directional"
            ),
            index=UX.index
        ))
    )
    # print(df_sp)
    # # check if dedup
    # bc = df_sp["BC"].iloc[0]
    # before = df_sp.loc[df_sp["BC"] == bc, "UX"].nunique()
    # after = df_sp.loc[df_sp["BC"] == bc, "UB_directional"].nunique()
    # print("Unique(UX) vs Unique(UB_directional) in BC =", bc, "->", before, after)

    # 4) subsample
    # sub = UMIcountPy.subsample_recompute(df_sp, mu_nSpikeUMI=100)
    # print(sub)
