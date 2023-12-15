__version__ = "v1.0"
__copyright__ = "Copyright 2024"
__license__ = "MIT"
__developer__ = "Jianfeng Sun"
__lab__ = "Cribbslab"

import os
import sys
import time
import numpy as np
import pandas as pd
from umiche.align.Read import read as aliread
from umiche.align.Write import write as aliwrite
from umiche.deduplicate.method.Build import Build as umibuild

from umiche.deduplicate.method.Cluster import cluster as umiclust
from umiche.deduplicate.method.Adjacency import adjacency as umiadj
from umiche.deduplicate.method.Directional import directional as umidirec
from umiche.deduplicate.method.MarkovClustering import markovClustering as umimcl

from umiche.util.Writer import writer as gwriter
from umiche.util.Hamming import hamming
from umiche.util.Number import number as rannum
from umiche.util.Console import Console


class OnePos:

    def __init__(
            self,
            bam_fpn,
            ed_thres,
            method,
            mcl_fold_thres=None,
            inflat_val=2.0,
            exp_val=2,
            iter_num=100,
            umi_='_',
            work_dir='./',
            verbose=False,
    ):
        """

        Parameters
        ----------
        bam_fpn
            str - the full path of a BAM file curated by requirements of different dedup modules
        ed_thres
            int - an edit distance threshold (>1, integer)
        method
            str - a deduplication method (mcl, mcl_val, mcl_ed, cluster, unique, ajacency, directional)
        mode
            str - externally or internally run the module (external by defualt, internal)
        mcl_fold_thres
            float - a mcl fold threshold (1.5 by defualt)
        inflat_val
            float - an inflation value for generating mcl clusters (2.0 by defualt)
        exp_val
            int - an expansion value for generating mcl clusters (2 by defualt)
        iter_num
            int - number of iterations for mcl (100 by defualt)
        sv_fpn
            str - the deduplication file path
        verbose
            bool - print log on the console, (True by default or False)
        """
        self.method = method
        self.bam_fpn = bam_fpn
        self.ed_thres = ed_thres
        self.mcl_fold_thres = mcl_fold_thres
        self.inflat_val = inflat_val
        self.exp_val = exp_val
        self.iter_num = iter_num
        self.work_dir = work_dir
        self.verbose = verbose

        self.umibuild = umibuild
        self.rannum = rannum()
        self.gwriter = gwriter()
        self.console = Console()
        self.console.verbose = self.verbose

        # sys.stdout = open(self.work_dir + self.method + '_log.txt', 'w')


        self.umiclust = umiclust()
        self.umiadj = umiadj()
        self.umidirec = umidirec()
        self.umimcl = umimcl(
            inflat_val=self.inflat_val,
            exp_val=self.exp_val,
            iter_num=self.iter_num,
        )

        self.alireader = aliread(bam_fpn=self.bam_fpn, verbose=self.verbose)
        self.df_bam = self.alireader.todf(tags=[])
        ### @@ self.df_bam
        #           id  ... query_qualities
        # 0          0  ...            None
        # 1          1  ...            None
        # 2          2  ...            None
        # ...      ...  ...             ...
        # 20683  20683  ...            None
        self.console.print('======># of raw reads: {}'.format(self.df_bam.shape[0]))
        self.df_bam = self.df_bam.loc[self.df_bam['reference_id'] != -1]
        self.console.print('======># of reads with qualified chrs: {}'.format(self.df_bam.shape[0]))

        self.df_bam['umi'] = self.df_bam['query_name'].apply(lambda x: x.split(umi_)[1])
        self.console.print('======># of unique umis: {}'.format(self.df_bam['umi'].unique().shape[0]))
        self.console.print('======># of redundant umis: {}'.format(self.df_bam['umi'].shape[0]))
        self.console.print('======>edit distance thres: {}'.format(self.ed_thres))

        self.df_bam['source'] = 1
        ### @@ self.df_bam
        #           id                     query_name  ...        umi  source
        # 0          0   SRR2057595.2985267_ACCGGTTTA  ...  ACCGGTTTA       1
        # 1          1  SRR2057595.13520751_CCAGGTTCT  ...  CCAGGTTCT       1
        # 2          2   SRR2057595.8901432_AGCGGTTAC  ...  AGCGGTTAC       1
        # ...      ...                            ...  ...        ...     ...
        # 20683  20683  SRR2057595.11966225_ACCGGTTGG  ...  ACCGGTTGG       1
        # [20684 rows x 13 columns]
        self.df_bam_gp = self.df_bam.groupby(by=['source'])
        self.pos_gp_keys = self.df_bam_gp.groups.keys()

        self.aliwriter = aliwrite(df=self.df_bam)

        self.console.print('======># of columns in the bam df: {}'.format(len(self.df_bam.columns)))
        self.console.print('======>Columns in the bam df: {}'.format(self.df_bam.columns.tolist()))
        self.console.print('======># of raw reads:\n{}'.format(self.df_bam))

        self.console.print('===>start building umi graphs...')
        umi_graph_build_stime = time.time()
        pos_gps = []
        res_sum = []
        for pos_g in self.pos_gp_keys:
            umi_vignette = self.umibuild(
                df=self.df_bam_gp.get_group(pos_g),
                ed_thres=self.ed_thres,
                verbose=verbose,
            ).data_summary
            # print(umi_vignette)
            # if len(umi_vignette['int_to_umi_dict']) == 1:
            #     print(True)
            #     print(umi_vignette['int_to_umi_dict'])
            #     continue
            # else:
            cc = self.umiclust.cc(umi_vignette['graph_adj'])
            # import json
            # with open('data.json', 'w') as f:
            #     json.dump(cc, f)
            pos_gps.append(pos_g)
            res_sum.append([
                umi_vignette,
                cc,
                umi_vignette['ave_ed'],
                [*umi_vignette['int_to_umi_dict'].keys()],
            ])
            ### @@ self.df['uniq_repr_nodes'] or [*umi_vignette['int_to_umi_dict'].keys()]
            # [0, 1, 2, ..., 1948]
        self.df = pd.DataFrame(
            data=res_sum,
            columns=['vignette', 'cc', 'ave_ed', 'uniq_repr_nodes'],
            index=pos_gps,
        )
        ### @@  self.df
        # vignette  ...   uniq_repr_nodes
        # 1  {'graph_adj': {0: [77, 81, 97, 153, 205, 228, ...  ...  [0, 1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13,...
        # [1 rows x 3 columns]
        self.console.print('===>time for building umi graphs: {:.2f}s'.format(time.time() - umi_graph_build_stime))

        self.df['uniq_umi_len'] = self.df['uniq_repr_nodes'].apply(lambda x: self.length(x))
        ### @@ self.df['uniq_umi_len']
        # 1    1949
        # Name: uniq_umi_len, dtype: int64
        self.console.print('===>start deduplication by the {} method...'.format(self.method))
        if self.method == 'unique':
            dedup_umi_stime = time.time()
            self.df['uniq_sgl_mark'] = self.df['uniq_repr_nodes'].apply(lambda x: 'yes' if len(x) == 1 else 'no')
            self.df = self.df.loc[self.df['uniq_sgl_mark'] == 'no']
            self.console.print('======># of positions with non-single umis: {}'.format(self.df.shape[0]))
            self.console.print('======>finish finding deduplicated UMIs in {:.2f}s'.format(time.time() - dedup_umi_stime))
            # self.console.print('======># of UMIs deduplicated {}'.format(self.df['uniq_umi_len'].loc['yes']))
            self.console.print('======>calculate average edit distances between umis...')
            ### @@ self.df['ave_ed']
            # 1    5.0
            # Name: ave_eds, dtype: float64
            # self.df['ave_ed'] = self.df.apply(lambda x: self.edave(x, by_col='uniq_repr_nodes'), axis=1)
            self.ave_ed_bins = self.df['ave_ed'].value_counts().sort_index().to_frame().reset_index()
            self.console.check("======>bins for average edit distance\n{}".format(self.ave_ed_bins))
            self.df['dedup_uniq_diff_pos'] = self.df.apply(lambda x: self.num_removed_uniq_umis(x, by_col='uniq_repr_nodes'), axis=1)
            self.df['dedup_read_diff_pos'] = self.df.apply(lambda x: self.num_removed_reads(x, by_col='uniq_repr_nodes'), axis=1)
            self.console.print('======># of deduplicated unique umis {}'.format(self.df['dedup_uniq_diff_pos'].sum()))
            self.console.print('======># of deduplicated reads {}'.format(self.df['dedup_read_diff_pos'].sum()))
            self.gwriter.generic(df=self.ave_ed_bins, sv_fpn=self.work_dir + 'uniq_ave_ed_bin.txt', index=True)
            self.df_dedup_sum = self.df[[
                'ave_ed',
                'uniq_umi_len',
                'dedup_uniq_diff_pos',
                'dedup_read_diff_pos',
            ]]
            self.gwriter.generic(
                df=self.df_dedup_sum,
                sv_fpn=self.work_dir + 'uniq_dedup_sum.txt',
                index=True,
                header=True,
            )
            self.console.print('======>start writing deduplicated reads to BAM...')
            dedup_reads_write_stime = time.time()
            self.df['uniq_bam_ids'] = self.df.apply(lambda x: self.bamids(x, by_col='uniq_repr_nodes'), axis=1)
            self.aliwriter.tobam(
                tobam_fpn=self.work_dir + self.method + '_dedup.bam',
                tmpl_bam_fpn=self.bam_fpn,
                whitelist=self.decompose(list_nd=self.df['uniq_bam_ids'].values),
            )
            self.console.print('======>finish writing in {:.2f}s'.format(time.time() - dedup_reads_write_stime))

        elif self.method == 'cluster':
            dedup_umi_stime = time.time()
            self.df['cc_repr_nodes'] = self.df.apply(lambda x: self.umimax(x, by_col='cc'), axis=1)
            ### @@ self.df['cc_repr_nodes']
            # 1    [2]
            # Name: cc_repr_nodes, dtype: object
            self.df['cc_umi_len'] = self.df['cc_repr_nodes'].apply(lambda x: self.length(x))
            ### @@ self.df['cc_umi_len']
            # 1    1
            # Name: cc_umi_len, dtype: int64
            self.console.print('======>finish finding deduplicated umis in {:.2f}s'.format(time.time() - dedup_umi_stime))
            # self.console.print('======># of umis deduplicated to be {}'.format(self.df['cc_umi_len'].loc['yes']))
            self.console.print('======>calculate average edit distances between umis...')
            # self.df['ave_ed'] = self.df.apply(lambda x: self.edave(x, by_col='cc_repr_nodes'), axis=1)
            self.df['dedup_uniq_diff_pos'] = self.df.apply(lambda x: self.num_removed_uniq_umis(x, by_col='cc_repr_nodes'), axis=1)
            self.df['dedup_read_diff_pos'] = self.df.apply(lambda x: self.num_removed_reads(x, by_col='cc_repr_nodes'), axis=1)
            self.console.print('======># of deduplicated unique umis {}'.format(self.df['dedup_uniq_diff_pos'].sum()))
            self.console.print('======># of deduplicated reads {}'.format(self.df['dedup_read_diff_pos'].sum()))
            self.ave_ed_bins = self.df['ave_ed'].value_counts().sort_index().to_frame().reset_index()
            self.console.check("======>bins for average edit distance\n{}".format(self.ave_ed_bins))
            self.gwriter.generic(df=self.ave_ed_bins, sv_fpn=self.work_dir + 'cc_ave_ed_bin.txt', index=True)
            self.df_dedup_sum = self.df[[
                'cc_umi_len',
                'ave_ed',
                'uniq_umi_len',
                'dedup_uniq_diff_pos',
                'dedup_read_diff_pos',
            ]]
            self.gwriter.generic(
                df=self.df_dedup_sum,
                sv_fpn=self.work_dir + 'cc_dedup_sum.txt',
                index=True,
                header=True,
            )
            self.console.print('======>start writing deduplicated reads to BAM...')
            dedup_reads_write_stime = time.time()
            self.df['cc_bam_ids'] = self.df.apply(lambda x: self.bamids(x, by_col='cc_repr_nodes'), axis=1)
            self.aliwriter.tobam(
                tobam_fpn=self.work_dir + self.method + '_dedup.bam',
                tmpl_bam_fpn=self.bam_fpn,
                whitelist=self.decompose(list_nd=self.df['cc_bam_ids'].values),
            )
            self.console.print('======>finish writing in {:.2f}s'.format(time.time() - dedup_reads_write_stime))

        elif self.method == 'adjacency':
            dedup_umi_stime = time.time()
            self.df['adj'] = self.df.apply(
                lambda x: self.umiadj.decompose(
                    cc_sub_dict=self.umiadj.umi_tools(
                        connected_components=x['cc'],
                        df_umi_uniq_val_cnt=x['vignette']['df_umi_uniq_val_cnt'],
                        graph_adj=x['vignette']['graph_adj'],
                    )['clusters'],
                ),
                axis=1,
            )
            self.df['adj_repr_nodes'] = self.df.apply(lambda x: self.umimax(x, by_col='adj'), axis=1)
            self.df['adj_umi_len'] = self.df['adj_repr_nodes'].apply(lambda x: self.length(x))
            self.console.print('======>finish finding deduplicated umis in {:.2f}s'.format(time.time() - dedup_umi_stime))
            # self.console.print('======># of umis deduplicated to be {}'.format(self.df['adj_umi_len'].loc['yes']))
            self.console.print('======>calculate average edit distances between umis...')
#             self.df['ave_ed'] = self.df.apply(lambda x: self.edave(x, by_col='adj_repr_nodes'), axis=1)
            self.df['dedup_uniq_diff_pos'] = self.df.apply(lambda x: self.num_removed_uniq_umis(x, by_col='adj_repr_nodes'), axis=1)
            self.df['dedup_read_diff_pos'] = self.df.apply(lambda x: self.num_removed_reads(x, by_col='adj_repr_nodes'), axis=1)
            self.console.print('======># of deduplicated unique umis {}'.format(self.df['dedup_uniq_diff_pos'].sum()))
            self.console.print('======># of deduplicated reads {}'.format(self.df['dedup_read_diff_pos'].sum()))
            self.ave_ed_bins = self.df['ave_ed'].value_counts().sort_index().to_frame().reset_index()
            self.console.check("======>bins for average edit distance\n{}".format(self.ave_ed_bins))
            self.gwriter.generic(df=self.ave_ed_bins, sv_fpn=self.work_dir + 'adj_ave_ed_bin.txt', index=True)
            self.df_dedup_sum = self.df[[
                'adj_umi_len',
                'ave_ed',
                'uniq_umi_len',
                'dedup_uniq_diff_pos',
                'dedup_read_diff_pos',
            ]]
            self.gwriter.generic(
                df=self.df_dedup_sum,
                sv_fpn=self.work_dir + 'adj_dedup_sum.txt',
                index=True,
                header=True,
            )
            self.console.print('======>start writing deduplicated reads to BAM...')
            dedup_reads_write_stime = time.time()
            self.df['adj_bam_ids'] = self.df.apply(lambda x: self.bamids(x, by_col='adj_repr_nodes'), axis=1)
            self.aliwriter.tobam(
                tobam_fpn=self.work_dir + self.method + '_dedup.bam',
                tmpl_bam_fpn=self.bam_fpn,
                whitelist=self.decompose(list_nd=self.df['adj_bam_ids'].values),
            )
            self.console.print('======>finish writing in {:.2f}s'.format(time.time() - dedup_reads_write_stime))

        elif self.method == 'directional':
            dedup_umi_stime = time.time()
            self.df['direc'] = self.df.apply(
                lambda x: self.umidirec.decompose(
                    cc_sub_dict=self.umidirec.umi_tools(
                        connected_components=x['cc'],
                        df_umi_uniq_val_cnt=x['vignette']['df_umi_uniq_val_cnt'],
                        graph_adj=x['vignette']['graph_adj'],
                    )['clusters'],
                ),
                axis=1,
            )
            self.df['direc_repr_nodes'] = self.df.apply(lambda x: self.umimax(x, by_col='direc'), axis=1)
            self.df['direc_umi_len'] = self.df['direc_repr_nodes'].apply(lambda x: self.length(x))
            self.console.print('======>finish finding deduplicated umis in {:.2f}s'.format(time.time() - dedup_umi_stime))
#             self.console.print('======># of umis deduplicated to be {}'.format(self.df['direc_umi_len'].loc['yes']))
            self.console.print('======>calculate average edit distances between umis...')
#             self.df['ave_ed'] = self.df.apply(lambda x: self.edave(x, by_col='direc_repr_nodes'), axis=1)
            self.df['dedup_uniq_diff_pos'] = self.df.apply(lambda x: self.num_removed_uniq_umis(x, by_col='direc_repr_nodes'), axis=1)
            self.df['dedup_read_diff_pos'] = self.df.apply(lambda x: self.num_removed_reads(x, by_col='direc_repr_nodes'), axis=1)
            self.console.print('======># of deduplicated unique umis {}'.format(self.df['dedup_uniq_diff_pos'].sum()))
            self.console.print('======># of deduplicated reads {}'.format(self.df['dedup_read_diff_pos'].sum()))
            self.ave_ed_bins = self.df['ave_ed'].value_counts().sort_index().to_frame().reset_index()
            self.console.check("======>bins for average edit distance\n{}".format(self.ave_ed_bins))
            self.gwriter.generic(df=self.ave_ed_bins, sv_fpn=self.work_dir + 'direc_ave_ed_bin.txt', index=True)
            self.df_dedup_sum = self.df[[
                'direc_umi_len',
                'ave_ed',
                'uniq_umi_len',
                'dedup_uniq_diff_pos',
                'dedup_read_diff_pos',
            ]]
            self.gwriter.generic(
                df=self.df_dedup_sum,
                sv_fpn=self.work_dir + 'direc_dedup_sum.txt',
                index=True,
                header=True,
            )
            self.console.print('======>start writing deduplicated reads to BAM...')
            dedup_reads_write_stime = time.time()
            self.df['direc_bam_ids'] = self.df.apply(lambda x: self.bamids(x, by_col='direc_repr_nodes'), axis=1)
            self.aliwriter.tobam(
                tobam_fpn=self.work_dir + self.method + '_dedup.bam',
                tmpl_bam_fpn=self.bam_fpn,
                whitelist=self.decompose(list_nd=self.df['direc_bam_ids'].values),
            )
            self.console.print('======>finish writing in {:.2f}s'.format(time.time() - dedup_reads_write_stime))

        elif self.method == 'mcl':
            dedup_umi_stime = time.time()
            self.df['mcl'] = self.df.apply(
                lambda x: self.umimcl.decompose(
                    list_nd=self.umimcl.dfclusters(
                        connected_components=x['cc'],
                        graph_adj=x['vignette']['graph_adj'],
                    )['clusters'].values,
                ),
                axis=1,
            )
            self.df['mcl_repr_nodes'] = self.df.apply(lambda x: self.umimax(x, by_col='mcl'), axis=1)
            self.df['mcl_umi_len'] = self.df['mcl_repr_nodes'].apply(lambda x: self.length(x))
            self.console.print('======>finish finding deduplicated umis in {:.2f}s'.format(time.time() - dedup_umi_stime))
#             self.console.print('======># of umis deduplicated to be {}'.format(self.df['mcl_umi_len'].loc['yes']))
            self.console.print('======>calculate average edit distances between umis...')
#             self.df['ave_ed'] = self.df.apply(lambda x: self.edave(x, by_col='mcl_repr_nodes'), axis=1)
            self.df['dedup_uniq_diff_pos'] = self.df.apply(lambda x: self.num_removed_uniq_umis(x, by_col='mcl_repr_nodes'), axis=1)
            self.df['dedup_read_diff_pos'] = self.df.apply(lambda x: self.num_removed_reads(x, by_col='mcl_repr_nodes'), axis=1)
            self.console.print('======># of deduplicated unique umis {}'.format(self.df['dedup_uniq_diff_pos'].sum()))
            self.console.print('======># of deduplicated reads {}'.format(self.df['dedup_read_diff_pos'].sum()))
            self.ave_ed_bins = self.df['ave_ed'].value_counts().sort_index().to_frame().reset_index()
            self.console.check("======>bins for average edit distance\n{}".format(self.ave_ed_bins))
            self.gwriter.generic(df=self.ave_ed_bins, sv_fpn=self.work_dir + 'mcl_ave_ed_bin.txt', index=True)
            self.df_dedup_sum = self.df[[
                'mcl_umi_len',
                'ave_ed',
                'uniq_umi_len',
                'dedup_uniq_diff_pos',
                'dedup_read_diff_pos',
            ]]
            self.gwriter.generic(
                df=self.df_dedup_sum,
                sv_fpn=self.work_dir + 'mcl_dedup_sum.txt',
                index=True,
                header=True,
            )
            self.console.print('======>start writing deduplicated reads to BAM...')
            dedup_reads_write_stime = time.time()
            self.df['mcl_bam_ids'] = self.df.apply(lambda x: self.bamids(x, by_col='mcl_repr_nodes'), axis=1)
            self.aliwriter.tobam(
                tobam_fpn=self.work_dir + self.method + '_dedup.bam',
                tmpl_bam_fpn=self.bam_fpn,
                whitelist=self.decompose(list_nd=self.df['mcl_bam_ids'].values),
            )
            self.console.print('======>finish writing in {:.2f}s'.format(time.time() - dedup_reads_write_stime))

        elif self.method == 'mcl_val':
            dedup_umi_stime = time.time()
            self.df['mcl_val'] = self.df.apply(
                lambda x: self.umimcl.decompose(
                    list_nd=self.umimcl.maxval_val(
                        df_mcl_ccs=self.umimcl.dfclusters(
                            connected_components=x['cc'],
                            graph_adj=x['vignette']['graph_adj'],
                        ),
                        df_umi_uniq_val_cnt=x['vignette']['df_umi_uniq_val_cnt'],
                        thres_fold=self.mcl_fold_thres,
                    )['clusters'].values,
                ),
                axis=1,
            )
            self.df['mcl_val_repr_nodes'] = self.df.apply(lambda x: self.umimax(x, by_col='mcl_val'), axis=1)
            self.df['mcl_val_umi_len'] = self.df['mcl_val_repr_nodes'].apply(lambda x: self.length(x))
            self.console.print('======>finish finding deduplicated umis in {:.2f}s'.format(time.time() - dedup_umi_stime))
#             self.console.print('======># of umis deduplicated to be {}'.format(self.df['mcl_val_umi_len'].loc['yes']))
            self.console.print('======>calculate average edit distances between umis...')
#             self.df['ave_ed'] = self.df.apply(lambda x: self.edave(x, by_col='mcl_val_repr_nodes'), axis=1)
            self.df['dedup_uniq_diff_pos'] = self.df.apply(lambda x: self.num_removed_uniq_umis(x, by_col='mcl_val_repr_nodes'), axis=1)
            self.df['dedup_read_diff_pos'] = self.df.apply(lambda x: self.num_removed_reads(x, by_col='mcl_val_repr_nodes'), axis=1)
            self.console.print('======># of deduplicated unique umis {}'.format(self.df['dedup_uniq_diff_pos'].sum()))
            self.console.print('======># of deduplicated reads {}'.format(self.df['dedup_read_diff_pos'].sum()))
            self.ave_ed_bins = self.df['ave_ed'].value_counts().sort_index().to_frame().reset_index()
            self.console.check("======>bins for average edit distance\n{}".format(self.ave_ed_bins))
            self.gwriter.generic(df=self.ave_ed_bins, sv_fpn=self.work_dir + 'mcl_val_ave_ed_bin.txt', index=True)
            self.df_dedup_sum = self.df[[
                'mcl_val_umi_len',
                'ave_ed',
                'uniq_umi_len',
                'dedup_uniq_diff_pos',
                'dedup_read_diff_pos',
            ]]
            self.gwriter.generic(
                df=self.df_dedup_sum,
                sv_fpn=self.work_dir + 'mcl_val_dedup_sum.txt',
                index=True,
                header=True,
            )
            self.console.print('======>start writing deduplicated reads to BAM...')
            dedup_reads_write_stime = time.time()
            self.df['mcl_val_bam_ids'] = self.df.apply(lambda x: self.bamids(x, by_col='mcl_val_repr_nodes'), axis=1)
            self.aliwriter.tobam(
                tobam_fpn=self.work_dir + self.method + '_dedup.bam',
                tmpl_bam_fpn=self.bam_fpn,
                whitelist=self.decompose(list_nd=self.df['mcl_val_bam_ids'].values),
            )
            self.console.print('======>finish writing in {:.2f}s'.format(time.time() - dedup_reads_write_stime))

        elif self.method == 'mcl_ed':
            dedup_umi_stime = time.time()
            self.df['mcl_ed'] = self.df.apply(
                lambda x: self.umimcl.decompose(
                    list_nd=self.umimcl.maxval_ed(
                        df_mcl_ccs=self.umimcl.dfclusters(
                            connected_components=x['cc'],
                            graph_adj=x['vignette']['graph_adj'],
                        ),
                        df_umi_uniq_val_cnt=x['vignette']['df_umi_uniq_val_cnt'],
                        umi_uniq_mapped_rev=x['vignette']['int_to_umi_dict'],
                        thres_fold=self.mcl_fold_thres,
                    )['clusters'].values,
                ),
                axis=1,
            )
            self.df['mcl_ed_repr_nodes'] = self.df.apply(lambda x: self.umimax(x, by_col='mcl_ed'), axis=1)
            self.df['mcl_ed_umi_len'] = self.df['mcl_ed_repr_nodes'].apply(lambda x: self.length(x))
            self.console.print('======>finish finding deduplicated umis in {:.2f}s'.format(time.time() - dedup_umi_stime))
#             self.console.print('======># of umis deduplicated to be {}'.format(self.df['mcl_ed_umi_len'].loc['yes']))
            self.console.print('======>calculate average edit distances between umis...')
#             self.df['ave_ed'] = self.df.apply(lambda x: self.edave(x, by_col='mcl_ed_repr_nodes'), axis=1)
            self.df['dedup_uniq_diff_pos'] = self.df.apply(lambda x: self.num_removed_uniq_umis(x, by_col='mcl_ed_repr_nodes'), axis=1)
            self.df['dedup_read_diff_pos'] = self.df.apply(lambda x: self.num_removed_reads(x, by_col='mcl_ed_repr_nodes'), axis=1)
            self.console.print('======># of deduplicated unique umis {}'.format(self.df['dedup_uniq_diff_pos'].sum()))
            self.console.print('======># of deduplicated reads {}'.format(self.df['dedup_read_diff_pos'].sum()))
            self.ave_ed_bins = self.df['ave_ed'].value_counts().sort_index().to_frame().reset_index()
            self.console.check("======>bins for average edit distance\n{}".format(self.ave_ed_bins))
            self.gwriter.generic(df=self.ave_ed_bins, sv_fpn=self.work_dir + 'mcl_ed_ave_ed_bin.txt', index=True)
            self.df_dedup_sum = self.df[[
                'mcl_ed_umi_len',
                'ave_ed',
                'uniq_umi_len',
                'dedup_uniq_diff_pos',
                'dedup_read_diff_pos',
            ]]
            self.gwriter.generic(
                df=self.df_dedup_sum,
                sv_fpn=self.work_dir + 'mcl_ed_dedup_sum.txt',
                index=True,
                header=True,
            )
            self.console.print('======>start writing deduplicated reads to BAM...')
            dedup_reads_write_stime = time.time()
            self.df['mcl_ed_bam_ids'] = self.df.apply(lambda x: self.bamids(x, by_col='mcl_ed_repr_nodes'), axis=1)
            self.aliwriter.tobam(
                tobam_fpn=self.work_dir + self.method + '_dedup.bam',
                tmpl_bam_fpn=self.bam_fpn,
                whitelist=self.decompose(list_nd=self.df['mcl_ed_bam_ids'].values),
            )
            self.console.print('======>finish writing in {:.2f}s'.format(time.time() - dedup_reads_write_stime))

        # sys.stdout.close()

    def num_removed_uniq_umis(self, df_row, by_col):
        """
        It calculates the number of the removed unique UMIs, observed on a given genomic position,
        after CC, Adjacency, Directional, MCL, MCL_val, or MCL_ed.

        Parameters
        ----------
        df_row
            object - a pandas-like df row
        by_col
            str - a column name in question

        Returns
        -------
            int - the sum of deduplicated unique UMI counts per position

        """
        return df_row['uniq_umi_len'] - len(df_row[by_col])

    def num_removed_reads(self, df_row, by_col):
        """
        It calculates the number of the removed reads in total.

        Parameters
        ----------
        df_row
            object - a pandas-like df row
        by_col
            str - a column name in question

        Returns
        -------
            int - the total counts of deduplicated reads per position

        """
        ### @@ set(df_row['uniq_repr_nodes'])
        # {0, 1, 2, ..., 1948}
        ### @@ len(set(df_row['uniq_repr_nodes']))
        # 1949
        ### @@ set(df_row[by_col])
        # {2}
        ### @@ len(set(df_row[by_col]))
        # 1
        diff_nodes = set(df_row['uniq_repr_nodes']) - set(df_row[by_col])
        ### @@ len(diff_nodes)
        # 1948

        if diff_nodes != set():
            # print(diff_nodes)
            umi_val_cnt_dict = df_row['vignette']['df_umi_uniq_val_cnt'].to_dict()
            ### @@ df_row['vignette']['df_umi_uniq_val_cnt'].to_dict()
            # {2: 55, 780: 51, 416: 47, ..., 1948: 1}
            return sum(umi_val_cnt_dict[node] for node in diff_nodes)
        else:
            return 0

    def umimax(
            self,
            df_row,
            by_col,
    ):
        """
        It returns ids of UMIs (i.e., representatives in their groupped UMIs) that has the
        highest count among all reads in their given genomic position.

        Examples
        --------
        {0: [0, 77, 81, ..., 1016], 1: [42, 46, 12], ..., 100: [2, 3, 5]}
        if no. 77 UMI has the highest count among all reads in 0, it will be added
        to umi_maxval_ids.

        Parameters
        ----------
        df_row
            a row of a pandas dataframe
        by_col
            a name of a pandas dataframe column

        Returns
        -------
        a list

        """
        umi_val_cnts = df_row['vignette']['df_umi_uniq_val_cnt']
        ### @@ umi_val_cnts
        # 2       55
        # 780     51
        # 416     47
        #         ..
        # 1948     1
        # Name: count, Length: 1949, dtype: int64
        umi_maxval_ids = []
        for k_c, nodes in df_row[by_col].items():
            self.console.print('key: {} nodes: {}'.format(k_c, nodes))
            ### @@ k_c, nodes
            # 0 [0, 77, 81, ..., 1016]
            ### @@ umi_val_cnts.loc[umi_val_cnts.index.isin(nodes)]
            # 2       55
            # 780     51
            # 416     47
            #         ..
            # 1948     1
            # Name: count, Length: 1949, dtype: int64
            ### @@ umi_val_cnts.loc[umi_val_cnts.index.isin(nodes)].idxmax()
            # 2
            ### @@ umi_val_cnts.loc[umi_val_cnts.index.isin(nodes)].max()
            # 55
            umi_max = umi_val_cnts.loc[umi_val_cnts.index.isin(nodes)].idxmax()
            umi_maxval_ids.append(umi_max)
        return umi_maxval_ids

    def length(self, df_val):
        """

        Parameters
        ----------
        df_val
            a pandas dataframe row

        Returns
        -------
            int - the length of the list

        """
        return len(df_val)

    def decompose(self, list_nd):
        """

        Parameters
        ----------
        x

        Returns
        -------

        """
        print(list_nd)
        print(len(list_nd))
        list_md = []
        for i in list_nd:
            list_md = list_md + i
        self.console.print('======># of the total reads left after deduplication: {}'.format(len(list_md)))
        print(len(list_md))
        return list_md

    def bamids(
            self,
            df_row,
            by_col,
    ):
        """
        It outputs bamids of UMIs that are representative of all nodes in each group.

        Parameters
        ----------
        df_row
            a row of a pandas dataframe
        by_col
            a name of a pandas dataframe column

        Returns
        -------

        """
        uniq_umi_id_to_bam_id_dict = df_row['vignette']['uniq_umi_id_to_bam_id_dict']
        ### @@ uniq_umi_id_to_bam_id_dict
        # {0: 0, 1: 1, 2: 2, ..., 1948: 20609}
        ### @@ len(uniq_umi_id_to_bam_id_dict)
        # 1949
        list_1d = df_row[by_col]
        ### @@ list_1d
        # [2, 780, 416, ..., 1761]
        return [uniq_umi_id_to_bam_id_dict[node] for node in list_1d]

    def edave_deprecated(self, df_row, by_col):
        repr_nodes = df_row[by_col]
        node_len = len(repr_nodes)
        int_to_umi_dict = df_row['vignette']['int_to_umi_dict']
        if node_len != 1:
            ed_list = []
            for i in range(node_len):
                for j in range(i + 1, node_len):
                    ed_list.append(hamming().general(
                        s1=int_to_umi_dict[repr_nodes[i]],
                        s2=int_to_umi_dict[repr_nodes[j]],
                    ))
            return np.ceil(sum(ed_list) / (len(ed_list)))
        else:
            return -1

    def eds_(self, df_row, by_col):
        """"""
        print(df_row.index)
        repr_nodes = df_row[by_col]
        int_to_umi_dict = df_row['vignette']['int_to_umi_dict']
        umi_val_cnts = df_row['vignette']['df_umi_uniq_val_cnt']
        # print(repr_nodes)
        # if len(repr_nodes) == len(np.unique(repr_nodes)):
        #     print(True)
        # else:
        #     print(False)
        node_len = len(repr_nodes)
        if node_len != 1:
            ed_list = []
            for i in range(node_len):
                for j in range(i + 1, node_len):
                    if repr_nodes[i] != repr_nodes[j]:
                        ed_list = ed_list + [hamming().general(
                            int_to_umi_dict[repr_nodes[i]],
                            int_to_umi_dict[repr_nodes[j]])
                        ] * (umi_val_cnts.loc[repr_nodes[i]] * umi_val_cnts.loc[repr_nodes[j]])
            return round(sum(ed_list) / len(ed_list))
        else:
            return -1


if __name__ == "__main__":
    from umiche.path import to

    umiche = OnePos(
        # method='unique',
        # method='cluster',
        method='adjacency',
        # method='directional',
        # method='mcl',
        # method='mcl_val',
        # method='mcl_ed',

        # bam_fpn=to('data/example.bam'),
        bam_fpn=to('data/example_bundle.bam'),
        mcl_fold_thres=1.5,
        inflat_val=1.6,
        exp_val=2,
        iter_num=100,
        ed_thres=1,
        work_dir=to('data/assigned_sorted_dedup.bam'),

        verbose=True,
    )