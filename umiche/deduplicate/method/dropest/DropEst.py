__version__ = "v1.0"
__copyright__ = "Copyright 2025"
__license__ = "GPL-3.0"
__developer__ = "Jianfeng Sun"
__maintainer__ = "Jianfeng Sun"
__email__="jianfeng.sunmt@gmail.com"

from typing import List, Optional, Tuple

from dataclasses import dataclass
import numpy as np
import pandas as pd


def _hamming_distance(a: str, b: str) -> int:
    """Compute the Hamming distance between two equal-length strings."""
    if len(a) != len(b):
        return max(len(a), len(b))  # large so it's never considered adjacent
    # Fast path: vectorized compare via numpy from bytes
    return sum(ch1 != ch2 for ch1, ch2 in zip(a, b))


def _neighbors_within_hamming1(umis: List[str]) -> List[Tuple[str, str]]:
    """Return all unordered pairs (u, v) of UMIs with Hamming distance == 1.

    Complexity: O(n^2 * L) worst-case; for large n consider bucketing by 1-mismatch
    patterns (not implemented here to keep code compact).
    """
    pairs = []
    n = len(umis)
    for i in range(n):
        ui = umis[i]
        for j in range(i + 1, n):
            uj = umis[j]
            if _hamming_distance(ui, uj) == 1:
                pairs.append((ui, uj))
    return pairs


def _quantize(values: np.ndarray, quant_borders: np.ndarray) -> np.ndarray:
    """Quantize values by given borders (right-open bins). Returns integer bins.
    Example: borders for 10 quantiles -> indices 0..9.
    """
    # np.digitize returns 1..len(borders); shift to 0-based
    return np.digitize(values, quant_borders, right=False)


def _average_pairwise_hamming(umis: List[str]) -> float:
    """Average pairwise Hamming distance among a list of UMIs.
    Returns -1.0 when fewer than 2 UMIs are present (undefined).
    """
    n = len(umis)
    if n < 2:
        return -1.0
    tot = 0
    cnt = 0
    for i in range(n):
        for j in range(i + 1, n):
            tot += _hamming_distance(umis[i], umis[j])
            cnt += 1
    return float(tot) / float(cnt) if cnt else -1.0


@dataclass
class NBParams:
    # Likelihood tables for discrete features per class
    # Each is a 2D array: shape (n_classes=2, n_bins)
    quality_lh: np.ndarray  # P(Q=q | class)
    ratio_lh: np.ndarray    # P(R=r | class)
    cbase_lh: np.ndarray    # P(Cb=b | class)
    prior: np.ndarray       # P(class)
    # Bin borders (right-open): used to quantize new data
    quality_borders: np.ndarray
    ratio_borders: np.ndarray
    cbase_borders: np.ndarray


class DiscreteNaiveBayes:

    """Naive Bayes model
    A tiny discrete Naive Bayes over binned features with Laplace smoothing.

    Classes: 0 = negative (keep separate UMIs), 1 = positive (merge base->target)
    """

    def __init__(self, n_bins_quality=10, n_bins_ratio=10, n_bins_cbase=10, laplace=1.0):
        self.n_bins_quality = n_bins_quality
        self.n_bins_ratio = n_bins_ratio
        self.n_bins_cbase = n_bins_cbase
        self.laplace = laplace
        self.params: Optional[NBParams] = None

    def _compute_borders(self, x: np.ndarray, n_bins: int) -> np.ndarray:
        # Use empirical quantiles as borders (exclude 0 and 1)
        qs = np.linspace(0, 1, n_bins + 1)[1:-1]
        borders = np.unique(np.quantile(x, qs))
        if borders.size == 0:
            # Fallback to a single mid border
            borders = np.array([np.median(x)])
        return borders

    def fit(self, quality: np.ndarray, ratio: np.ndarray, cbase: np.ndarray, y: np.ndarray) -> None:
        # Determine borders from data
        q_b = self._compute_borders(quality, self.n_bins_quality)
        r_b = self._compute_borders(ratio, self.n_bins_ratio)
        c_b = self._compute_borders(cbase, self.n_bins_cbase)

        qz = _quantize(quality, q_b)
        rz = _quantize(ratio, r_b)
        cz = _quantize(cbase, c_b)

        classes = np.array([0, 1], dtype=int)
        prior = np.array([np.mean(y == 0), np.mean(y == 1)], dtype=float)
        # Regularize priors if degenerate
        if prior.sum() == 0 or np.any(prior == 0):
            prior = np.array([0.5, 0.5], dtype=float)

        def lh_from_bins(binned: np.ndarray, n_bins: int) -> np.ndarray:
            tab = np.zeros((2, n_bins)) + self.laplace
            for cls in classes:
                vals, counts = np.unique(binned[y == cls], return_counts=True)
                tab[cls, vals] += counts
                tab[cls, :] /= tab[cls, :].sum()
            return tab

        q_lh = lh_from_bins(qz, max(q_b.size + 1, 1))
        r_lh = lh_from_bins(rz, max(r_b.size + 1, 1))
        c_lh = lh_from_bins(cz, max(c_b.size + 1, 1))

        self.params = NBParams(
            quality_lh=q_lh,
            ratio_lh=r_lh,
            cbase_lh=c_lh,
            prior=prior,
            quality_borders=q_b,
            ratio_borders=r_b,
            cbase_borders=c_b,
        )

    def predict_proba(self, quality: np.ndarray, ratio: np.ndarray, cbase: np.ndarray) -> np.ndarray:
        assert self.params is not None, "Model not fitted"
        qz = _quantize(quality, self.params.quality_borders)
        rz = _quantize(ratio, self.params.ratio_borders)
        cz = _quantize(cbase, self.params.cbase_borders)

        # Log-space for numerical stability
        logp0 = (
            np.log(self.params.prior[0])
            + np.log(self.params.quality_lh[0, qz])
            + np.log(self.params.ratio_lh[0, rz])
            + np.log(self.params.cbase_lh[0, cz])
        )
        logp1 = (
            np.log(self.params.prior[1])
            + np.log(self.params.quality_lh[1, qz])
            + np.log(self.params.ratio_lh[1, rz])
            + np.log(self.params.cbase_lh[1, cz])
        )
        # Normalize
        m = np.maximum(logp0, logp1)
        p0 = np.exp(logp0 - m)
        p1 = np.exp(logp1 - m)
        denom = p0 + p1
        return np.vstack([p0 / denom, p1 / denom]).T


@dataclass
class DropEstUMICorrector:
    """Bayesian UMI correction inspired by dropEst, wrapped for in-Python use.

    Parameters
    ----------
    method : {'bayesian', 'classic'}
        'bayesian' uses a trained NB classifier over discretized features.
        'classic' performs directional adjacency merging with a multiplier rule.
    mult : float
        Only used for 'classic': merge u->v if count(u) <= mult * count(v).
    quality_quants : int
        Number of quantile bins for quality features in NB.
    ratio_quants : int
        Number of quantile bins for the count ratio feature in NB.
    cbase_quants : int
        Number of quantile bins for base UMI counts in NB.
    merge_threshold : float
        Posterior threshold P(merge | features) to merge base->target.
    max_iter : int
        Max iterations of the merge loop within a (cell, gene) group.
    adjust_collisions : bool
        Whether to apply a simple uniform-UMI-space occupancy correction.
    umi_space : Optional[int]
        Total distinct UMI codes (e.g., 4**L). Required if adjust_collisions=True.
    cell_col, gene_col, umi_col, qual_col : str
        Column names in the input DataFrame.
    """

    method: str = 'bayesian'
    mult: float = 1.0
    quality_quants: int = 10
    ratio_quants: int = 10
    cbase_quants: int = 10
    merge_threshold: float = 0.90
    max_iter: int = 100

    adjust_collisions: bool = False
    umi_space: Optional[int] = None

    cell_col: str = 'cell'
    gene_col: str = 'gene'
    umi_col: str = 'umi'
    qual_col: Optional[str] = 'umi_qual'

    # learned model
    nb: Optional[DiscreteNaiveBayes] = None

    def fit_transform(self, df: pd.DataFrame, return_type: str = 'counts') -> pd.DataFrame | pd.Series | pd.DataFrame:
        """Train (if needed) and perform UMI correction.

        Parameters
        ----------
        df : pandas.DataFrame
            Input per-read table.
        return_type : {'counts', 'matrix', 'umis', 'summary'}
            - 'counts': per-(cell,gene,umi) corrected counts (default)
            - 'matrix': a cell x gene count matrix
            - 'umis':   per-(cell,gene) UMI totals after correction (Series)
            - 'summary': MultiIndex (cell,gene) DataFrame with columns
                ['dedup_cnt','ave_eds','num_uniq_umis','num_diff_dedup_uniq_umis','num_diff_dedup_reads']
        """
        df = self._normalize_columns(df)

        # Aggregate per (cell,gene,umi)
        agg = self._aggregate_umis(df)

        if self.method.lower() == 'bayesian':
            if self.nb is None:
                self.nb = self._train_nb_from_dataset(agg)
            corrected = self._correct_groups_bayesian(agg)
        elif self.method.lower() == 'classic':
            corrected = self._correct_groups_classic(agg)
        else:
            raise ValueError(f"Unknown method: {self.method}")

        # Prepare outputs
        if return_type == 'summary':
            # Pre-correction unique UMI count per (cell,gene)
            pre_uniqs = (
                agg.groupby([self.cell_col, self.gene_col])[self.umi_col]
                   .nunique()
                   .rename('num_uniq_umis')
            )
            # Post-correction unique UMI count (dedup count)
            post_uniqs = (
                corrected.groupby([self.cell_col, self.gene_col])[self.umi_col]
                         .nunique()
                         .rename('dedup_cnt')
            )
            # Average pairwise Hamming distance among retained UMIs
            aveds = (
                corrected.groupby([self.cell_col, self.gene_col])[self.umi_col]
                         .apply(lambda s: _average_pairwise_hamming(list(pd.unique(s))))
                         .rename('ave_eds')
            )
            summary = pd.concat([post_uniqs, aveds, pre_uniqs], axis=1).fillna(0)
            # Differences
            summary['num_diff_dedup_uniq_umis'] = (summary['num_uniq_umis'] - summary['dedup_cnt']).astype(int)
            # In this aggregated-UMI pipeline, "reads" here refers to removed UMI representatives
            summary['num_diff_dedup_reads'] = summary['num_diff_dedup_uniq_umis'].astype(int)
            # Dtypes / order
            summary['dedup_cnt'] = summary['dedup_cnt'].astype(int)
            summary['num_uniq_umis'] = summary['num_uniq_umis'].astype(int)
            summary = summary[['dedup_cnt','ave_eds','num_uniq_umis','num_diff_dedup_uniq_umis','num_diff_dedup_reads']]
            summary.sort_index(inplace=True)

            summary.index = summary.index.map(lambda t: f"('{t[0]}', '{t[1]}')")
            summary.index.name = None

            return summary
        if return_type == 'umis':
            s = corrected.groupby([self.cell_col, self.gene_col])['umi'].nunique()
            if self.adjust_collisions:
                s = self._adjust_collisions_series(s)
            return s

        # Per-(cell,gene,umi) counts
        counts = corrected.groupby([self.cell_col, self.gene_col, 'umi'], as_index=False)['read_count'].sum()

        if return_type == 'matrix':
            mat = counts.pivot_table(index=self.cell_col, columns=self.gene_col, values='read_count', fill_value=0, aggfunc='sum')
            if self.adjust_collisions:
                # Optional: convert to UMI counts (presence) and adjust — conservative
                umis_per_gene = (counts.assign(val=1)
                                 .pivot_table(index=self.cell_col, columns=self.gene_col, values='val', fill_value=0, aggfunc='sum'))
                if self.umi_space is not None:
                    mat = self._adjust_collisions_matrix(umis_per_gene)
            return mat

        return counts

    def _normalize_columns(self, df: pd.DataFrame) -> pd.DataFrame:
        # Internals
        # Ensure required columns exist
        for col in [self.cell_col, self.gene_col, self.umi_col]:
            if col not in df.columns:
                raise KeyError(f"Missing required column: '{col}'")
        return df

    def _aggregate_umis(self, df: pd.DataFrame) -> pd.DataFrame:
        # Aggregate per (cell,gene,umi): read_count and mean quality
        group_cols = [self.cell_col, self.gene_col, self.umi_col]
        base = df.groupby(group_cols, as_index=False).size().rename(columns={'size': 'read_count'})
        if self.qual_col and self.qual_col in df.columns:
            qual_df = (
                df.groupby(group_cols, as_index=False)[self.qual_col]
                  .mean()
                  .rename(columns={self.qual_col: 'umi_quality'})
            )
            agg_df = base.merge(qual_df, on=group_cols, how='left')
        else:
            agg_df = base.copy()
            agg_df['umi_quality'] = 30.0  # a reasonable default Phred-like score
        return agg_df

    def _build_training_pairs(self, agg: pd.DataFrame) -> pd.DataFrame:
        """
        Training (Bayesian)
        Construct candidate (base -> target) pairs for NB training across all groups.

        For each (cell,gene) group, connect UMIs whose Hamming distance == 1.
        We form a *directed* pair from lower-count UMI (base) to higher-count UMI (target).
        """
        pairs = []
        for (cell, gene), sub in agg.groupby([self.cell_col, self.gene_col]):
            umis = sub[self.umi_col].tolist()
            counts = dict(zip(sub[self.umi_col], sub['read_count']))
            quals = dict(zip(sub[self.umi_col], sub['umi_quality']))
            edges = _neighbors_within_hamming1(umis)
            for u, v in edges:
                cu, cv = counts[u], counts[v]
                if cu == cv:
                    # Choose direction by quality if available, else lexicographic
                    bu, bv = quals.get(u, 30.0), quals.get(v, 30.0)
                    base, target = (u, v) if bu <= bv else (v, u)
                else:
                    base, target = (u, v) if cu < cv else (v, u)
                cb = counts[base]
                ct = counts[target]
                qb = float(quals.get(base, 30.0))
                ratio = cb / max(ct, 1)
                pairs.append((cb, ct, qb, ratio))
        if not pairs:
            return pd.DataFrame(columns=['c_base', 'c_tgt', 'q_base', 'ratio', 'y'])

        dfp = pd.DataFrame(pairs, columns=['c_base', 'c_tgt', 'q_base', 'ratio'])
        # Pseudo-labels (weak supervision):
        # Positive (merge) if base is much smaller than target and/or low quality
        y = ((dfp['ratio'] <= 0.2) | ((dfp['c_base'] <= 2) & (dfp['c_tgt'] >= 4))).astype(int)
        # Hard negatives if similar counts
        y[(dfp['ratio'] >= 0.5) & (dfp['c_base'] >= 3) & (dfp['c_tgt'] >= 3)] = 0
        dfp['y'] = y.values
        return dfp

    def _train_nb_from_dataset(self, agg: pd.DataFrame) -> DiscreteNaiveBayes:
        pairs = self._build_training_pairs(agg)
        nb = DiscreteNaiveBayes(
            n_bins_quality=self.quality_quants,
            n_bins_ratio=self.ratio_quants,
            n_bins_cbase=self.cbase_quants,
            laplace=1.0,
        )
        if len(pairs) == 0:
            # No neighbors; create a degenerate NB that never merges
            nb.fit(
                np.array([30.0]),
                np.array([1.0]),
                np.array([1.0]),
                np.array([0]),
            )
            return nb
        nb.fit(
            pairs['q_base'].to_numpy(float),
            pairs['ratio'].to_numpy(float),
            pairs['c_base'].to_numpy(float),
            pairs['y'].to_numpy(int),
        )
        return nb

    def _correct_groups_bayesian(self, agg: pd.DataFrame) -> pd.DataFrame:
        """(Bayesian)"""
        assert self.nb is not None, "NB model must be trained"
        out_rows = []
        for (cell, gene), sub in agg.groupby([self.cell_col, self.gene_col]):
            corrected = self._correct_one_group_bayes(sub)
            corrected[self.cell_col] = cell
            corrected[self.gene_col] = gene
            out_rows.append(corrected)
        return pd.concat(out_rows, ignore_index=True) if out_rows else agg

    def _correct_one_group_bayes(self, sub: pd.DataFrame) -> pd.DataFrame:
        # Work on a mutable copy
        sub = sub.copy()
        # Map umi -> (count, quality)
        counts = dict(zip(sub[self.umi_col], sub['read_count']))
        quals = dict(zip(sub[self.umi_col], sub['umi_quality']))

        # Iterative merging
        umi_set = set(counts.keys())
        for step in range(self.max_iter):
            if len(umi_set) <= 1:
                break
            umis = list(umi_set)
            # Build candidate directed edges base->target for current iteration
            cands = []
            for u, v in _neighbors_within_hamming1(umis):
                cu, cv = counts[u], counts[v]
                if cu == cv:
                    bu, bv = quals.get(u, 30.0), quals.get(v, 30.0)
                    base, target = (u, v) if bu <= bv else (v, u)
                else:
                    base, target = (u, v) if cu < cv else (v, u)
                cb, ct = counts[base], counts[target]
                qb = float(quals.get(base, 30.0))
                ratio = cb / max(ct, 1)
                cands.append((base, target, cb, ct, qb, ratio))
            if not cands:
                break

            cand_df = pd.DataFrame(cands, columns=['base', 'target', 'c_base', 'c_tgt', 'q_base', 'ratio'])
            # Predict posterior P(merge)
            proba = self.nb.predict_proba(cand_df['q_base'].to_numpy(float),
                                          cand_df['ratio'].to_numpy(float),
                                          cand_df['c_base'].to_numpy(float))
            cand_df['p_merge'] = proba[:, 1]

            # Dynamic threshold: slightly stricter for larger graphs
            dyn_thr = min(self.merge_threshold + 0.05 * np.log1p(len(umi_set)), 0.99)
            cand_df = cand_df[cand_df['p_merge'] >= dyn_thr]
            if cand_df.empty:
                break

            # Resolve conflicts: prefer merges into larger targets, then lower ratio, then lower quality
            cand_df.sort_values(by=['c_tgt', 'ratio', 'q_base'], ascending=[False, True, True], inplace=True)

            merged_any = False
            used_bases = set()
            used_targets = set()
            for row in cand_df.itertuples(index=False):
                base, target = row.base, row.target
                if base in used_bases:
                    continue
                if base not in umi_set or target not in umi_set:
                    continue
                # Apply merge: move counts from base to target; drop base
                counts[target] += counts[base]
                used_bases.add(base)
                used_targets.add(target)
                umi_set.remove(base)
                merged_any = True
            if not merged_any:
                break

        # Emit corrected per-UMI counts for this group
        rows = []
        for u in sorted(umi_set):
            rows.append({self.umi_col: u, 'read_count': counts[u], 'umi_quality': quals.get(u, 30.0)})
        return pd.DataFrame(rows)

    def _correct_groups_classic(self, agg: pd.DataFrame) -> pd.DataFrame:
        """(Classic directional)"""
        out_rows = []
        for (cell, gene), sub in agg.groupby([self.cell_col, self.gene_col]):
            umis = sub[self.umi_col].tolist()
            counts = dict(zip(sub[self.umi_col], sub['read_count']))
            quals = dict(zip(sub[self.umi_col], sub['umi_quality']))
            umi_set = set(umis)
            for u, v in _neighbors_within_hamming1(umis):
                cu, cv = counts[u], counts[v]
                if cu <= self.mult * cv:
                    loser, winner = u, v
                elif cv <= self.mult * cu:
                    loser, winner = v, u
                else:
                    continue
                if loser in umi_set and winner in umi_set:
                    counts[winner] += counts[loser]
                    umi_set.remove(loser)
            rows = [{self.umi_col: u, 'read_count': counts[u], 'umi_quality': quals.get(u, 30.0)} for u in sorted(umi_set)]
            dfc = pd.DataFrame(rows)
            dfc[self.cell_col] = cell
            dfc[self.gene_col] = gene
            out_rows.append(dfc)
        return pd.concat(out_rows, ignore_index=True) if out_rows else agg

    # ----------------- Collision adjustment (optional) ----------------- #

    def _adjust_collisions_series(self, s: pd.Series) -> pd.Series:
        """Apply a simple occupancy adjustment to per-(cell,gene) UMI counts.

        Assumes a uniform UMI space with size `umi_space`. This is a *rough* heuristic
        and differs from dropEst's DP-based adjustment; use with caution.
        """
        if not self.adjust_collisions:
            return s
        if self.umi_space is None:
            raise ValueError("umi_space must be provided when adjust_collisions=True")
        # Heuristic: invert expected unique occupancy under Poissonization
        # E[unique] = m * (1 - exp(-N/m))  => estimate N given observed unique k
        m = float(self.umi_space)
        k = s.astype(float)
        # Newton iterations to solve for N: f(N) = m*(1-exp(-N/m)) - k = 0
        def invert_k(kk: float) -> float:
            if kk <= 1e-9:
                return 0.0
            # initial guess
            N = kk
            for _ in range(20):
                f = m * (1.0 - np.exp(-N / m)) - kk
                fp = np.exp(-N / m)
                if abs(fp) < 1e-12:
                    break
                N_new = N - f / fp
                if N_new < 0:
                    N_new = N / 2
                if abs(N_new - N) < 1e-6:
                    N = N_new
                    break
                N = N_new
            return float(N)
        adjusted = s.apply(invert_k)
        return adjusted

    def _adjust_collisions_matrix(self, umis_matrix: pd.DataFrame) -> pd.DataFrame:
        s = umis_matrix.stack()
        adj = self._adjust_collisions_series(s)
        return adj.unstack(umis_matrix.columns)


# ---------------------------- Convenience function ---------------------------- #

def dropest_correct(
    df: pd.DataFrame,
    method: str = 'bayesian',
    mult: float = 1.0,
    quality_quants: int = 10,
    ratio_quants: int = 10,
    cbase_quants: int = 10,
    merge_threshold: float = 0.90,
    max_iter: int = 100,
    adjust_collisions: bool = False,
    umi_space: Optional[int] = None,
    cell_col: str = 'cell',
    gene_col: str = 'gene',
    umi_col: str = 'umi',
    qual_col: Optional[str] = 'umi_qual',
    return_type: str = 'counts',
) -> pd.DataFrame | pd.Series | pd.DataFrame:
    """One-shot convenience wrapper around DropEstUMICorrector."""
    corr = DropEstUMICorrector(
        method=method,
        mult=mult,
        quality_quants=quality_quants,
        ratio_quants=ratio_quants,
        cbase_quants=cbase_quants,
        merge_threshold=merge_threshold,
        max_iter=max_iter,
        adjust_collisions=adjust_collisions,
        umi_space=umi_space,
        cell_col=cell_col,
        gene_col=gene_col,
        umi_col=umi_col,
        qual_col=qual_col,
    )
    return corr.fit_transform(df, return_type=return_type)


if __name__ == "__main__":
    from umiche.bam.Reader import ReaderChunk

    # bam_fpn = "/mnt/d/Document/Programming/Python/umiche/umiche/data/r1/trumicount/10xn9k_10c/10xn9k_10c_tagged.sorted.filtered.sorted.bam"
    # df = ReaderChunk(
    #     bam_fpn=bam_fpn,
    #     bam_fields=None,
    #     tag_whitelist=['MB', 'CB', 'XF'],
    #     categorize=["chrom"],
    #     verbose=True,
    # ).todf(chunk_size=2_000_000)
    # print(df)

    # df.rename(columns={'CB': 'cell', 'XF': 'gene', 'MB': 'umi'}, inplace=True)

    from umiche.bam.Reader import ReaderChunk
    from umiche.util.Console import Console
    bam_fpn = "/mnt/d/Document/Programming/Python/umiche/umiche/data/r1/umicountr/Smartseq3.TTACCTGCCAGATTCG.bam"
    df_bam = ReaderChunk(
        bam_fpn=bam_fpn,
        bam_fields=['contig', 'pos', 'CIGAR', 'seq', 'read'],
        tag_whitelist=['BC', 'QU', 'UX', 'UB', 'GE'],
        verbose=True,
    ).todf(chunk_size=1_000_000)
    # print(df)
    spikecontig = "diySpike"
    spikename = "g_diySpike4"
    df_bam = df_bam[df_bam["UX"] != ""]
    df_bam = df_bam[df_bam["GE"] == spikename]
    df_bam = df_bam.reset_index(drop=True)
    df_bam = df_bam.drop("GE", axis=1)
    # print(df_bam)
    from umiche.deduplicate.spikein.SpikeUMI import SpikeUMI

    df_bam = SpikeUMI(verbose=True).extract_spike_dat(
        df=df_bam,
        match_seq_before_UMI="GAGCCTGGGGGAACAGGTAGG",
        match_seq_after_UMI="CTCGGAGGAGAAA",
    )
    Console(True).df_column_summary(df_bam)
    from umiche.deduplicate.method.trimer.Expand import Expand

    df_bam["htUMI"] = df_bam["UX"].apply(Expand().homotrimer_umi)
    mut_rate = 0.3
    from umiche.deduplicate.method.trimer.Error import Error

    df_bam["htUMI_" + str(mut_rate)] = df_bam["htUMI"].apply(lambda umi: Error().mutated(umi, mut_rate=mut_rate, mode="normal"))


    summary = dropest_correct(
        df_bam,
        method='bayesian',
        cell_col='BC', gene_col='spikeUMI', umi_col="htUMI_" + str(mut_rate),
        qual_col='umi_qual',
        return_type='summary'
    )

    print(summary.head())

    summary.to_csv(
        '/mnt/d/Document/Programming/Python/umiche/umiche/data/r1/trumicount/10xn9k_10c/df.txt',
        header=True,
        sep="\t",
        index=True,
    )

