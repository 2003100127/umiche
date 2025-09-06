import pandas as pd
from collections import defaultdict
import pysam  # added for BAM writing

class UMIS:
    """
    Minimal UMI-collapsing (umis/tagcount style) — DataFrame version.
    - Input: a DataFrame converted from BAM (one row per alignment).
    - Required columns: cell/umi/gene. Optional: NH (multi-mapping count), pos (positional dedup), is_unmapped/ref_name fallback.
    - Evidence = sum of per-alignment weights; default weight = 1/NH; if evidence >= min_evidence, the (cell,gene,umi[,pos]) counts as 1 molecule.
    - Output: a gene-by-cell matrix of unique UMI counts (integers).
    """
    def __init__(
            self,
             tag_map=None,
             min_evidence=1.0,
             weighted=True,
             positional=False,
    ):
        # Column name mapping (adjust to your DataFrame columns if needed)
        self.tag_map = {
            'cell': 'CB',          # cell barcode column
            'umi': 'MB',           # UMI column
            'gene': 'XF',          # gene tag column (fallback to ref_name if empty)
            'ref_name': 'chrom',   # reference name column (transcript/contig), used as gene fallback
            'nh': 'NH',            # multi-mapping count column
            'pos': 'pos',          # alignment start position (used only when positional=True)
            'is_unmapped': 'is_unmapped',  # unmapped flag (optional)
        }
        if tag_map:
            self.tag_map.update(tag_map)

        self.min_evidence = float(min_evidence)
        self.weighted = bool(weighted)       # True => use 1/NH weighting
        self.positional = bool(positional)

    def _resolve_gene(self, df: pd.DataFrame) -> pd.Series:
        gene_col = self.tag_map['gene']
        ref_col = self.tag_map['ref_name']
        s_gene = None
        if gene_col in df.columns:
            # XF/GX may look like "gene1,gene2", take the first entry
            s_gene = df[gene_col].astype(str).str.split(',').str[0]
            s_gene = s_gene.replace({'None': pd.NA, 'nan': pd.NA})
        if (s_gene is None or s_gene.isna().all()) and ref_col in df.columns:
            s_gene = df[ref_col].astype(str)
        if s_gene is None:
            raise ValueError("No gene-like column found; provide tag_map['gene'] or tag_map['ref_name'].")
        return s_gene

    def count_df(
            self,
            df: pd.DataFrame,
            whitelist=None,
            drop_zeros=True,
    ) -> pd.DataFrame:
        # 1) Drop unmapped reads if the flag exists
        unmapped_col = self.tag_map['is_unmapped']
        if unmapped_col in df.columns:
            df = df[~df[unmapped_col].astype(bool)].copy()

        # 2) Parse required columns
        cell_col = self.tag_map['cell']
        umi_col  = self.tag_map['umi']
        pos_col  = self.tag_map['pos']
        nh_col   = self.tag_map['nh']

        if cell_col not in df.columns or umi_col not in df.columns:
            raise ValueError(f"DataFrame is missing required columns: '{cell_col}' or '{umi_col}'. Please map them via tag_map.")

        s_cell = df[cell_col].astype(str)
        s_umi  = df[umi_col].astype(str)
        s_gene = self._resolve_gene(df)

        # 3) Optional whitelist filtering
        if whitelist is not None:
            s_mask = s_cell.isin(set(whitelist))
            s_cell, s_umi, s_gene = s_cell[s_mask], s_umi[s_mask], s_gene[s_mask]
            df = df.loc[s_mask]

        # 4) Weights (default 1/NH; fall back to 1 if NH is absent)
        if self.weighted and (nh_col in df.columns):
            nh = pd.to_numeric(df[nh_col], errors='coerce').fillna(1.0)
            nh = nh.mask(nh <= 0, 1.0)  # guard against 0
            weight = 1.0 / nh
        else:
            weight = pd.Series(1.0, index=df.index)

        # 5) Positional key (optional)
        if self.positional:
            if pos_col not in df.columns:
                raise ValueError(f"positional=True but missing position column '{pos_col}'. Please specify via tag_map.")
            s_pos = pd.to_numeric(df[pos_col], errors='coerce').fillna(-1).astype(int)
            keys_df = pd.DataFrame({'cell': s_cell, 'gene': s_gene, 'umi': s_umi, 'pos': s_pos})
            group_keys = ['cell', 'gene', 'umi', 'pos']
        else:
            keys_df = pd.DataFrame({'cell': s_cell, 'gene': s_gene, 'umi': s_umi})
            group_keys = ['cell', 'gene', 'umi']

        # Drop rows with missing key fields
        keys_df = keys_df.dropna(subset=['cell', 'gene', 'umi']).copy()
        weight  = weight.loc[keys_df.index]

        # 6) Aggregate evidence per (cell,gene,umi[,pos])
        #    Build a tuple key for grouping to reuse in multiple places
        tmp = keys_df.copy()
        tmp['_w'] = weight.values
        agg = tmp.groupby(group_keys, sort=False)['_w'].sum().reset_index(name='evidence')

        # 7) Keep UMI-molecule keys whose evidence passes threshold
        passed = agg[agg['evidence'] >= self.min_evidence].reset_index(drop=True)

        # 8) Summarize to (gene, cell): each passed UMI contributes exactly 1
        collapsed = passed.groupby(['gene', 'cell'], sort=False).size()

        # 9) Pivot to matrix (rows=gene, cols=cell), fill NA with 0, cast to int
        counts = collapsed.unstack('cell', fill_value=0)
        counts = counts.reindex(sorted(counts.index), axis=0)
        counts = counts.reindex(sorted(counts.columns), axis=1)

        summary_pairs = self.counts_to_cell_gene_index(counts, drop_zeros=drop_zeros)

        return summary_pairs

    def counts_to_cell_gene_index(self, counts: pd.DataFrame, drop_zeros: bool = True) -> pd.DataFrame:
        """
        Convert a gene-by-cell count matrix to a 1-column DataFrame whose index is
        the unique (cell, gene) pair formatted as '(cell, gene)' with no index name.

        Parameters
        ----------
        counts : pd.DataFrame
            Matrix with rows=genes and columns=cells (integers).
        drop_zeros : bool, default True
            If True, keep only entries with dedup_cnt > 0.

        Returns
        -------
        pd.DataFrame
            Index: string like '(CELL, GENE)' (no name)
            Column: 'dedup_cnt' (int)
        """
        # Stack to long format: index becomes (gene, cell)
        s = counts.stack(dropna=False).astype(int)

        if drop_zeros:
            s = s[s > 0]

        # Reorder to (cell, gene) and format as '(cell, gene)'
        idx = pd.Index([f"('{cell}', '{gene}')" for (gene, cell) in s.index], name=None)

        out = pd.DataFrame({'dedup_cnt': s.values}, index=idx)
        out.index.name = None  # ensure the index has no name
        return out


def write_dedup_bam(
    df: pd.DataFrame,
    out_bam_fpn: str,
    template_bam_path: str,
    tag_map: dict = None,
    min_evidence: float = 1.0,
    weighted: bool = True,
    positional: bool = False,
    mode: str = "filter",
    duplicate_tag: str = "PD",
):
    """
    Write a deduplicated BAM from a DataFrame (no change to your counting logic).

    Parameters
    ----------
    df : pd.DataFrame
        DataFrame converted from BAM. Must contain a 'read' column with pysam.AlignedSegment objects.
        Should also contain the columns specified by tag_map (cell, umi, gene, optional NH/pos/is_unmapped).
    out_bam_fpn : str
        Output BAM file path for the deduplicated BAM.
    template_bam_path : str
        Path to a source BAM used only to provide the output header/template (references, header text).
        This does NOT re-read alignments; we still write from df['read'].
    tag_map : dict, optional
        Column name mapping like in UMIS:
            {'cell':..., 'umi':..., 'gene':..., 'ref_name':..., 'nh':..., 'pos':..., 'is_unmapped':...}
        Defaults: {'cell':'CB','umi':'MB','gene':'XF','ref_name':'chrom','nh':'NH','pos':'pos','is_unmapped':'is_unmapped'}
    min_evidence : float, default 1.0
        Minimum evidence a (cell,gene,umi[,pos]) needs to be counted (>= threshold).
    weighted : bool, default True
        If True, weight each alignment by 1/NH (fallback to 1 if NH is absent).
    positional : bool, default False
        If True, treat (cell,gene,umi,pos) as the dedup key; otherwise (cell,gene,umi).
    mode : {"filter","mark"}, default "filter"
        "filter": write only one representative read per passed UMI-molecule key.
        "mark":   write all reads, set duplicate_tag=0 for representatives and 1 for duplicates.
    duplicate_tag : str, default "PD"
        SAM tag used when mode="mark" to indicate duplicate (1) vs kept (0).

    Notes
    -----
    - This function does not modify your counting logic and does not alter df columns.
    - It only adds a file-emitting step using the same pass/fail criterion (evidence >= threshold).
    - We require a template BAM for header; the 'read' objects in df will be written as-is.
    """
    # Default mapping consistent with UMIS
    _tag_map = {
        'cell': 'CB',
        'umi': 'MB',
        'gene': 'XF',
        'ref_name': 'chrom',
        'nh': 'NH',
        'pos': 'pos',
        'is_unmapped': 'is_unmapped',
    }
    if tag_map:
        _tag_map.update(tag_map)

    if 'read' not in df.columns:
        raise ValueError("DataFrame must contain a 'read' column with pysam.AlignedSegment objects.")

    # Drop unmapped if flag exists
    unmapped_col = _tag_map['is_unmapped']
    if unmapped_col in df.columns:
        df = df[~df[unmapped_col].astype(bool)].copy()

    cell_col = _tag_map['cell']
    umi_col  = _tag_map['umi']
    pos_col  = _tag_map['pos']
    nh_col   = _tag_map['nh']

    if cell_col not in df.columns or umi_col not in df.columns:
        raise ValueError(f"DataFrame is missing required columns: '{cell_col}' or '{umi_col}'. Map them via tag_map.")

    # Build cell/umi/gene series
    s_cell = df[cell_col].astype(str)
    s_umi  = df[umi_col].astype(str)

    # Resolve gene with fallback to ref_name
    gene_col = _tag_map['gene']
    ref_col  = _tag_map['ref_name']
    s_gene = None
    if gene_col in df.columns:
        s_gene = df[gene_col].astype(str).str.split(',').str[0]
        s_gene = s_gene.replace({'None': pd.NA, 'nan': pd.NA})
    if (s_gene is None or s_gene.isna().all()) and ref_col in df.columns:
        s_gene = df[ref_col].astype(str)
    if s_gene is None:
        raise ValueError("No gene-like column found; provide tag_map['gene'] or tag_map['ref_name'].")

    # Weights
    if weighted and (nh_col in df.columns):
        nh = pd.to_numeric(df[nh_col], errors='coerce').fillna(1.0)
        nh = nh.mask(nh <= 0, 1.0)
        weight = 1.0 / nh
    else:
        weight = pd.Series(1.0, index=df.index)

    # Key DataFrame
    if positional:
        if pos_col not in df.columns:
            raise ValueError(f"positional=True but missing position column '{pos_col}'. Specify via tag_map.")
        s_pos = pd.to_numeric(df[pos_col], errors='coerce').fillna(-1).astype(int)
        keys_df = pd.DataFrame({'cell': s_cell, 'gene': s_gene, 'umi': s_umi, 'pos': s_pos})
        group_keys = ['cell', 'gene', 'umi', 'pos']
    else:
        keys_df = pd.DataFrame({'cell': s_cell, 'gene': s_gene, 'umi': s_umi})
        group_keys = ['cell', 'gene', 'umi']

    # Keep only rows with complete keys
    keys_df = keys_df.dropna(subset=['cell', 'gene', 'umi']).copy()
    weight  = weight.loc[keys_df.index]

    # Evidence per key
    tmp = keys_df.copy()
    tmp['_w'] = weight.values
    agg = tmp.groupby(group_keys, sort=False)['_w'].sum().reset_index(name='evidence')
    passed = agg[agg['evidence'] >= float(min_evidence)].reset_index(drop=True)

    # Select one representative row index per passed key
    passed_keys = set(passed[group_keys].apply(tuple, axis=1))
    # First occurrence per (cell,gene,umi[,pos]) group (deterministic)
    rep_df = keys_df.groupby(group_keys, sort=False).head(1)
    rep_tuples = rep_df[group_keys].apply(tuple, axis=1)
    rep_mask = rep_tuples.isin(passed_keys)
    keep_idx = set(rep_df.index[rep_mask].tolist())

    # Open template for header and write out
    with pysam.AlignmentFile(template_bam_path, 'rb') as src:
        with pysam.AlignmentFile(out_bam_fpn, 'wb', template=src) as outbam:
            if mode == "filter":
                # Write only representatives
                for i in sorted(keep_idx):
                    read = df.at[i, 'read']
                    outbam.write(read)
            elif mode == "mark":
                # Write all reads, tagging duplicates
                for i, read in zip(df.index, df['read']):
                    # Keep=0, Duplicate=1
                    dup_val = 0 if i in keep_idx else 1
                    try:
                        read.set_tag(duplicate_tag, int(dup_val), value_type='i')
                    except TypeError:
                        # Older pysam may not require value_type
                        read.set_tag(duplicate_tag, int(dup_val))
                    outbam.write(read)
            else:
                raise ValueError("mode must be one of {'filter','mark'}")

    return out_bam_fpn


if __name__ == "__main__":
    # Assume df is produced by your ReaderChunk and preprocessed
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
    mut_rate = 0.1
    from umiche.deduplicate.method.trimer.Error import Error

    df_bam["htUMI_" + str(mut_rate)] = df_bam["htUMI"].apply(
        lambda umi: Error().mutated(umi, mut_rate=mut_rate, mode="normal"))

    # Example: map your column names CB/MB/XF/NH/pos/chrom to the semantic fields
    tag_map = {
        'cell':'BC',
        'umi':"htUMI_" + str(mut_rate),
        'gene':'spikeUMI',
        'ref_name':'chrom',
        'nh':'NH',
        'pos':'pos',
        'is_unmapped':'is_unmapped',
    }

    counter = UMIS(tag_map=tag_map, min_evidence=1.0, weighted=True, positional=False)
    df_sum = counter.count_df(df_bam, drop_zeros=True)
    print(df_sum.head())
    print(df_sum)


    write_dedup_bam(
        df=df_bam,
        out_bam_fpn="/mnt/d/Document/Programming/Python/umiche/umiche/data/r1/umicountr/sss.bam",
        template_bam_path=bam_fpn,
        tag_map=tag_map,
        min_evidence=1.0,
        weighted=True,
        positional=False,
        mode="filter",
        duplicate_tag="PD",
    )