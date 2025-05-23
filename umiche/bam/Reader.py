__version__ = "v1.0"
__copyright__ = "Copyright 2025"
__license__ = "GPL-3.0"
__developer__ = "Jianfeng Sun"
__maintainer__ = "Jianfeng Sun"
__email__="jianfeng.sunmt@gmail.com"


import time
import pandas as pd
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


if __name__ == "__main__":
    from umiche.path import to

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
        bam_fpn=to('data/bone_marrow/trimmed2.bam')
    )

    # df = umiche.todf(tags=[])
    # df = umiche.todf(tags=['PO'])
    # df = umiche.todf_trust4(tags=[], filter_by='sequence')
    # print(df)
    # print(df.loc[0].values)

    df = umiche.todf_trust4(tags=['CB', 'UB'], filter_by='sequence')
    df['r1'] = df.apply(lambda x: x['CB'] + x['UB'], axis=1)
    print(df[['r1', 'query_name']])

    from umiche.fastq.Writer import Writer as fastqwriter
    fastqwriter().togz(
        list_2d=df[['r1', 'query_name']].values.tolist(),
        sv_fpn=to('data/bone_marrow/trimmed2_r1.fastq.gz'),
    )

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