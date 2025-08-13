__version__ = "v1.0"
__copyright__ = "Copyright 2025"
__license__ = "GPL-3.0"
__developer__ = "Jianfeng Sun"
__maintainer__ = "Jianfeng Sun"
__email__="jianfeng.sunmt@gmail.com"

from typing import Optional, Iterable, Dict, Any, List, Generator

import time
import pandas as pd
import pysam

from umiche.util.Console import Console


class Reader:

    def __init__(
            self,
            bam_fpn,
            verbose=False,
    ):
        """
        # pos_last = 0
        # chr_last = 0
        # if read.pos <= (pos_last + 1000) and read.reference_id == chr_last:
        # print(read.pos, pos_last, read.reference_id)
        # pos_last = read.pos
        # chr_last = read.reference_id
        Parameters
        ----------
        bam_fpn
        verbose
        """
        import pysam
        self.pysam = pysam

        # import bamnostic as bs

        self.console = Console()
        self.console.verbose = verbose
        self.console.print('===>reading the bam file... {}'.format(bam_fpn))
        read_bam_stime = time.time()
        self.pysam_bam = pysam.AlignmentFile(bam_fpn, "rb")
        # self.pysam_bam = bs.AlignmentFile(bam_fpn, "rb")
        self.console.print('===>reading BAM time: {:.2f}s'.format(time.time() - read_bam_stime))

    def bycol(self, col='sname'):
        """

        Parameters
        ----------
        col

        Returns
        -------

        """
        t = []
        if col == 'sname':
            for id, read in enumerate(self.pysam_bam):
                # print(read)
                # print(read.get_tags())
                if read.reference_id != -1:
                    tt = read.query_name
                else:
                    continue
                t.append(tt)
            return pd.DataFrame(t)

    def todf(
            self,
            tags=[],
            filter_by=None,
    ):
        """

        Notes
        -----
        11 columns from alignments deciphered by Pysam
            # read.query_name
            # read.flag
            # read.reference_id
            # read.reference_start
            # read.mapping_quality
            # read.cigar
            # read.query_sequence
            # read.next_reference_id
            # read.next_reference_start
            # read.template_length
            # read.query_qualities

        See Also
        --------
        https://pysam.readthedocs.io/en/latest/usage.html#creating-bam-cram-sam-files-from-scratch
        https://pysam.readthedocs.io/_/downloads/en/v0.12.0/pdf/

        Parameters
        ----------
        tags

        Returns
        -------

        """
        l = []
        self.console.print('=========>start converting bam to df...')
        stime = time.time()
        for id, read in enumerate(self.pysam_bam):
            # print(id)

            read_tags = read.get_tags()
            # print(read_tags)
            rt_dict = {k: v for k, v in read_tags}
            rt_keys = [*rt_dict.keys()]
            ### @@ read_tags | per read
            ### @@ rt_dict | per read
            ### @@ rt_keys | per read
            # 1st read
            # [('XA', 2), ('MD', '44T2T5'), ('NM', 2)]
            # {'XA': 2, 'MD': '44T2T5', 'NM': 2}
            # ['XA', 'MD', 'NM']
            # 2nd read
            # [('XA', 2), ('MD', '12A0C6'), ('NM', 2)]
            # {'XA': 2, 'MD': '12A0C6', 'NM': 2}
            # ['XA', 'MD', 'NM']
            # ...
            # the last read
            # [('XA', 2), ('MD', '1T4A13'), ('NM', 2)]
            # {'XA': 2, 'MD': '1T4A13', 'NM': 2}
            # ['XA', 'MD', 'NM']
            tag_keys = [rt_dict[k] if k in rt_keys else 'None' for k in tags]
            ### @@ tag_keys
            # [2, '44T2T5', 2]
            # [2, '12A0C6', 2]
            # ...
            # [2, '1T4A13', 2]
            vignette = [
                id,
                read.query_name,
                read.flag,
                read.reference_id,
                read.pos,
                read.mapping_quality,
                read.cigar,
                read.query_sequence,
                read.next_reference_id,
                read.next_reference_start,
                read.query_qualities,
                # read.template_length,
                read,
            ] + tag_keys
            l.append(vignette)
        df = pd.DataFrame(
            l,
            columns=[
                'id',
                'query_name',
                'flag',
                'reference_id',
                'genome_pos',
                'mapping_quality',
                'cigar',
                'query_sequence',
                'next_reference_id',
                'next_reference_start',
                'query_qualities',
                # 'template_length',
                'read',
            ] + tags,
        )
        if filter_by == 'sequence':
            df = df.loc[~df['query_sequence'].isin([None])] # .reset_index(drop=True)
            # df['query_sequence'] = df['query_sequence'].apply(lambda x: x[:20])
            return df
        else:
            # print(df['XA'].loc[df['reference_id'] != -1].shape)
            # print(df['MD'].loc[df['MD'] != 'None'].shape)
            # print(df['NM'].loc[df['NM'] != 'None'].shape)
            # print(df['XS'].loc[df['XS'] != 'None'].shape)
            # print(df['XT'].loc[df['XT'] != 'None'].shape)
            self.console.print('=========>time to df: {:.3f}s'.format(time.time() - stime))
            # df = df.sample(
            #     n=50,
            #     replace=True,
            # )
            return df

    def todf_trust4(
            self,
            tags=[],
            filter_by=None,
    ):
        l = []
        self.console.print('=========>start converting bam to df...')
        stime = time.time()
        for id, read in enumerate(self.pysam_bam):
            # print(id)
            # print(read.get_tag('CB'))
            # arr_query_name = read.query_name.split('_')
            # BC = arr_query_name[1]
            # UMI = arr_query_name[2]
            # read.set_tag('CB', BC)
            # read.set_tag('UB', UMI)
            # read.query_name = arr_query_name[0]
            # print(read.query_name)

            # if read.query_sequence != None:
            # read.query_sequence = None
            # len(read.query_sequence)

            # print(len(read.query_sequence))
            # print(read.query_qualities[:20])
            # t = read.query_qualities[:20]
            # u = read.query_sequence[:20]
            # print(u)
            # read.query_sequence = u
            # read.query_qualities = None

            read_tags = read.get_tags()
            # print(read_tags)
            rt_dict = {k: v for k, v in read_tags}
            rt_keys = [*rt_dict.keys()]
            ### @@ read_tags | per read
            ### @@ rt_dict | per read
            ### @@ rt_keys | per read
            # 1st read
            # [('XA', 2), ('MD', '44T2T5'), ('NM', 2)]
            # {'XA': 2, 'MD': '44T2T5', 'NM': 2}
            # ['XA', 'MD', 'NM']
            # 2nd read
            # [('XA', 2), ('MD', '12A0C6'), ('NM', 2)]
            # {'XA': 2, 'MD': '12A0C6', 'NM': 2}
            # ['XA', 'MD', 'NM']
            # ...
            # the last read
            # [('XA', 2), ('MD', '1T4A13'), ('NM', 2)]
            # {'XA': 2, 'MD': '1T4A13', 'NM': 2}
            # ['XA', 'MD', 'NM']
            tag_keys = [rt_dict[k] if k in rt_keys else 'None' for k in tags]
            ### @@ tag_keys
            # [2, '44T2T5', 2]
            # [2, '12A0C6', 2]
            # ...
            # [2, '1T4A13', 2]
            vignette = [
                id,
                read.query_name,
                read.flag,
                read.reference_id,
                read.pos,
                read.mapping_quality,
                read.cigar,
                read.query_sequence,
                read.next_reference_id,
                read.next_reference_start,
                read.query_qualities,
                # read.template_length,
                read,
            ] + tag_keys
            l.append(vignette)
        df = pd.DataFrame(
            l,
            columns=[
                'id',
                'query_name',
                'flag',
                'reference_id',
                'genome_pos',
                'mapping_quality',
                'cigar',
                'query_sequence',
                'next_reference_id',
                'next_reference_start',
                'query_qualities',
                # 'template_length',
                'read',
            ] + tags,
        )
        if filter_by == 'sequence':
            df = df.loc[~df['query_sequence'].isin([None])] # .reset_index(drop=True)
            # df['query_sequence'] = df['query_sequence'].apply(lambda x: x[:20])
            return df
        else:
            # print(df['XA'].loc[df['reference_id'] != -1].shape)
            # print(df['MD'].loc[df['MD'] != 'None'].shape)
            # print(df['NM'].loc[df['NM'] != 'None'].shape)
            # print(df['XS'].loc[df['XS'] != 'None'].shape)
            # print(df['XT'].loc[df['XT'] != 'None'].shape)
            self.console.print('=========>time to df: {:.3f}s'.format(time.time() - stime))
            # df = df.sample(
            #     n=50,
            #     replace=True,
            # )
            return df

    def todf11_depr(self, ):
        """
        Note
        ----
        Deprecated.

        Returns
        -------
        Dataframe of a bam file

        """
        l = []
        self.console.print('=========>start converting bam to df')
        import time
        stime = time.time()
        for id, read in enumerate(self.pysam_bam):
            read_piece = {
                'id': id,
                'query_name': read.query_name,
                'flag': read.flag,
                'reference_id': read.reference_id,
                'reference_start': read.reference_start,
                'mapping_quality': read.mapping_quality,
                'cigar': read.cigar,
                'query_sequence': read.query_sequence,
                'next_reference_id': read.next_reference_id,
                'next_reference_start': read.next_reference_start,
                'template_length': read.template_length,
                'query_qualities': read.query_qualities,
            }
            l.append(read_piece)
        df = pd.DataFrame.from_dict(l)
        self.console.print('=========>time to df: {:.3f}s'.format(time.time() - stime))
        return df


def _umi_from_qname(qname: str) -> str:
    parts = qname.split("_")
    if len(parts) < 2:
        raise ValueError(f"Read name lacks UMI (expected qname_UMI_...): {qname}")
    return parts[1]

class ReaderChunk:

    # @@ /*** -------- based on chuncks -------- ***/
    def __init__(
            self,
            bam_fpn: str,
            tag_whitelist: Optional[Iterable[str]],
            # sv_bam_fpn: str = None,
            bam_fields: Optional[List[str]] = None,
            categorize: Optional[List[str]] = None,
            verbose=False,
    ):
        import pysam
        self.pysam = pysam
        self.bam_fpn = bam_fpn
        # self.sv_bam_fpn = sv_bam_fpn
        self.bam_fields = bam_fields
        self.tag_whitelist = tag_whitelist
        self.categorize = categorize

        self.console = Console()
        self.console.verbose = verbose

        # Minimum necessary columns for deduplication
        self._DEFAULT_BAM_FIELDS = [
            "qname",
            "UMI",
            "chrom",
            "pos",
            "tlen",
            "is_read1",
            "is_read2",
            "is_reverse",
            "flag",
            "mapq",
            "read",
            "read_id",
        ]

    def _iter_rows(
            self,
            bam,
            tag_whitelist: Optional[Iterable[str]],
            bam_fields: List[str],
    ) -> Iterable[Dict[str, Any]]:
        bam_field_set = set(bam_fields)
        tag_running_title = 'iter rows'
        for r in self.console._tqdm(
                bam.fetch(until_eof=True),
                desc=f"[{tag_running_title}]",
                unit="read",
                position=0,
                leave=True,
                dynamic_ncols=False,
        ):
            row: Dict[str, Any] = {}
            if "qname" in bam_field_set:       row["qname"] = r.query_name
            # if "UMI" in bam_field_set:         row["UMI"] = _umi_from_qname(r.query_name)
            if "chrom" in bam_field_set:       row["chrom"] = r.reference_name
            if "pos" in bam_field_set:         row["pos"] = r.reference_start
            if "tlen" in bam_field_set:        row["tlen"] = r.template_length
            if "is_read1" in bam_field_set:    row["is_read1"] = r.is_read1
            if "is_read2" in bam_field_set:    row["is_read2"] = r.is_read2
            if "is_reverse" in bam_field_set:  row["is_reverse"] = r.is_reverse
            if "flag" in bam_field_set:        row["flag"] = r.flag
            if "mapq" in bam_field_set:        row["mapq"] = r.mapping_quality
            if "cigar" in bam_field_set:       row["cigar"] = r.cigarstring
            if "mate_chrom" in bam_field_set:  row["mate_chrom"] = r.next_reference_name
            if "read" in bam_field_set:    row["read"] = r
            if "read_id" in bam_field_set:    row["read_id"] = id(r)
            if tag_whitelist is not None:
                for tag, val in r.get_tags():
                    if tag in tag_whitelist:
                        row[tag] = val
            yield row

    def to_dataframe_chunks(
            self,
            bam_fpn: str,
            chunk_size: int = 1_000_000,
            tag_whitelist: Optional[Iterable[str]] = None,
            bam_fields: Optional[List[str]] = None,
            categorize: Optional[List[str]] = None
    ) -> Generator[pd.DataFrame, None, None]:
        """
        Streaming read BAM in chunks → yield DataFrame chunks (highly recommended for large files)

        Parameters
        ----------
        bam_fpn
        chunk_size
            number of rows per chunk
        tag_whitelist
            only expand these TAGs (None means no TAGs are expanded)
        bam_fields
            only keep these columns (default is deduplicated minimal necessary columns)
        categorize
            convert these columns to category (e.g., ["chrom","cigar"])

        Returns
        -------

        """

        if bam_fields is None:
            bam_fields = list(self._DEFAULT_BAM_FIELDS)
        cat_set = set(categorize or [])

        bam = self.pysam.AlignmentFile(bam_fpn, "rb")
        rows = []
        n = 0
        tag_running_title = 'Chunk'
        for row in self.console._tqdm(
                self._iter_rows(
                    bam=bam,
                    tag_whitelist=tag_whitelist,
                    bam_fields=bam_fields,
                ),
                desc=f"[{tag_running_title} chunk]",
                unit="chunk",
                position=0,
                leave=True,
                dynamic_ncols=False,
            ):
            rows.append(row)
            n += 1
            if n >= chunk_size:
                df = pd.DataFrame.from_records(rows)
                self._tighten_dtypes(df, cat_set)
                yield df
                rows.clear()
                n = 0
        if rows:
            df = pd.DataFrame.from_records(rows)
            self._tighten_dtypes(df, cat_set)
            yield df
        bam.close()

    @staticmethod
    def _tighten_dtypes(
            df: pd.DataFrame,
            cat_set: set,
    ):
        if "pos" in df: df["pos"] = df["pos"].astype("int32", copy=False)
        if "tlen" in df: df["tlen"] = df["tlen"].astype("int32", copy=False)
        if "flag" in df: df["flag"] = df["flag"].astype("int32", copy=False)
        if "mapq" in df: df["mapq"] = df["mapq"].astype("int16", copy=False)
        for b in ("is_read1", "is_read2", "is_reverse"):
            if b in df: df[b] = df[b].astype("bool", copy=False)
        for c in cat_set:
            if c in df: df[c] = df[c].astype("category")

    # /*** -------- One-time full-table processing (reuse for small/medium-scale tasks) -------- ***/
    def todf(
            self,
            chunk_size=1_000_000,
    ) -> pd.DataFrame:
        """Read the full table at once; not recommended for large files — use the chunked API instead."""
        dfs: List[pd.DataFrame] = []
        for df_chunk in self.to_dataframe_chunks(
                bam_fpn=self.bam_fpn,
                chunk_size=chunk_size,
                tag_whitelist=self.tag_whitelist,
                bam_fields=self.bam_fields,
                categorize=self.categorize,
        ):
            dfs.append(df_chunk)
        if not dfs:
            return pd.DataFrame()
        return pd.concat(dfs, ignore_index=True)


if __name__ == "__main__":
    from umiche.path import to

    # @@ /*** -------- based on Reader -------- ***/
    umiche = Reader(
        # bam_fpn=to('data/example.bam'),
        # bam_fpn=to('data/example_bundle1.bam'),
        # bam_fpn=to('data/hgmm_100_STAR_FC_sorted.bam'),
        # to('example/data/assigned_sorted_dedup.bam'),
        # bam_fpn=to('data/simu/monomer/sc/seq_err/permute_0/trimmed/seq_err_0.bam'),
        # bam_fpn=to('data/simu/umi/seq_errs/trimer/permute_0/bam/bipartite/seq_err_0.bam'),
        # bam_fpn=to('data/simu/trimer/pcr8/seq_errs/permute_0/trimmed/seq_err_18.bam'),
        # bam_fpn=to('example/data/deduplicated.bam'),
        # bam_fpn=to('example/data/RM82CLK1_S3_featurecounts_gene_sorted.bam'),

        # bam_fpn=to('data/bone_marrow/example.bam')
        # bam_fpn=to('data/bone_marrow/stripped_XT_sorted.bam')
        # bam_fpn=to('data/bone_marrow/stripped_XT_sorted_little.bam')
        # bam_fpn=to('data/bone_marrow/merge_corrected_sorted_little.bam')
        # bam_fpn=to('data/bone_marrow/merge_corrected_sorted.bam')
        # bam_fpn=to('data/bone_marrow/little.bam')
        # bam_fpn=to('data/bone_marrow/trimmed1.bam')
        # bam_fpn='/mnt/d/Document/Programming/Python/umiche/umiche/data/bone_marrow/trimmed2.bam',

        # ds
        # bam_fpn='/mnt/d/Document/Programming/Python/umiche/umiche/data/r1/trumicount/10xn9k_10c/10xn9k_10c.bam',
        bam_fpn="/mnt/d/Document/Programming/Python/umiche/umiche/data/r1/umitools/umitools.test.RNA-seq.sorted.tagged.bam",
    )

    df = umiche.todf(tags=['MB'])
    # df = umiche.todf(tags=['PO'])
    # df = umiche.todf_trust4(tags=[], filter_by='sequence')
    print(df)
    print(df.columns)
    # print(df.loc[0].values)

    # df = umiche.todf_trust4(tags=['CB', 'UB'], filter_by='sequence')
    # df['r1'] = df.apply(lambda x: x['CB'] + x['UB'], axis=1)
    # print(df[['r1', 'query_name']])
    #
    # from umiche.fastq.Writer import Writer as fastqwriter
    # fastqwriter().togz(
    #     list_2d=df[['r1', 'query_name']].values.tolist(),
    #     sv_fpn=to('data/bone_marrow/trimmed2_r1.fastq.gz'),
    # )

    # from umiche.bam.Writer import Writer as aliwriter
    # aliwriter(df=df).tobam(
    #     tobam_fpn=to('data/bone_marrow/trimmed2.bam'),
    #     tmpl_bam_fpn=to('data/bone_marrow/merge_corrected_sorted_little.bam'),
    #     whitelist=df.index,
    # )
    # print(df.columns)
    # print(df.query_name)
    # df = umiche.todf(tags=['BC', 'XS', 'XT'])
    # df = umiche.todf(tags=['XS', 'XT'])
    # print(df[['query_name', 'BC', 'XS', 'XT']])

    # df = df.loc[df['XS'] == 'Assigned']
    # print(df)


    # @@ /*** -------- based on ReadChunck -------- ***/
    p = ReaderChunk(
        bam_fpn="/mnt/d/Document/Programming/Python/umiche/umiche/data/r1/umitools/umitools.test.RNA-seq.sorted.tagged.bam",
        bam_fields=None,
        tag_whitelist=['MB'],
        categorize=["chrom"],
        verbose=True,
    )
    # df = p.todf()
    # print(df)
    # print(df.columns)

    for i, df_chunk in enumerate(p.to_dataframe_chunks(
            bam_fpn="/mnt/d/Document/Programming/Python/umiche/umiche/data/r1/umitools/umitools.test.RNA-seq.sorted.tagged.bam",
            chunk_size=2_000_000,
            tag_whitelist=None,  # 例：["CB","MB","XF"]
            categorize=["chrom"]  # chrom converts into a category
    )):
        print(f"chunk {i}: shape={df_chunk.shape}")
        # df_chunk

