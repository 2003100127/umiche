__version__ = "v1.0"
__copyright__ = "Copyright 2025"
__license__ = "GPL-3.0"
__developer__ = "Jianfeng Sun"
__maintainer__ = "Jianfeng Sun"
__email__ = "jianfeng.sunmt@gmail.com"


import time

from umiche.bam.Writer import Writer as aliwriter

from umiche.deduplicate.Gadgetry import Gadgetry as umigadgetry

# dedup methods
from umiche.deduplicate.method.Adjacency import Adjacency as umiadj
from umiche.deduplicate.method.Directional import Directional as umidirec
from umiche.deduplicate.method.MarkovClustering import MarkovClustering as umimcl
from umiche.deduplicate.method.Clustering import Clustering as umiclustering
from umiche.deduplicate.method.trimer.SetCover import SetCover as umisc
from umiche.deduplicate.method.trimer.MajorityVote import MajorityVote as umimv

from umiche.util.Writer import Writer as fwriter
from umiche.util.Console import Console


class Tabulate:

    def __init__(
            self,
            df,
            df_bam,
            bam_fpn,
            work_dir,
            heterogeneity,
            verbose=False,
    ):
        self.df = df
        self.df_bam = df_bam
        self.bam_fpn = bam_fpn
        self.work_dir = work_dir
        self.heterogeneity = heterogeneity

        self.aliwriter = aliwriter(df=self.df_bam)

        self.umigadgetry = umigadgetry()

        self.fwriter = fwriter()

        self.console = Console()
        self.console.verbose = verbose

    def set_cover(
            self,
            **kwargs,
    ):
        self.df_umi_uniq = self.df_bam.drop_duplicates(subset=['umi'], keep='first')
        # print(self.df_umi_uniq)
        series_uniq_umi = self.df_umi_uniq.umi
        # print(series_uniq_umi)
        self.umi_to_int_dict = {k: id for id, k in enumerate(series_uniq_umi)}

        dedup_cnt, multimer_umi_solved_by_sc, multimer_umi_not_solved, shortlisted_multimer_umi_list, monomer_umi_lens, multimer_umi_lens = umisc().greedy(
            multimer_list=series_uniq_umi.values,
            recur_len=kwargs['umi_unit_pattern'],
            split_method=kwargs['split_method'],
        )
        # print(monomer_umi_lens)
        # print(multimer_umi_lens)
        # print(dedup_cnt)
        self.df.loc[0, 'dedup_cnt'] = dedup_cnt
        self.df.loc[0, 'num_solved'] = len(multimer_umi_solved_by_sc)
        self.df.loc[0, 'num_not_solved'] = len(multimer_umi_not_solved)
        self.df.loc[0, 'monomer_umi_len'] = ';'.join([str(i) for i in monomer_umi_lens])
        self.df.loc[0, 'multimer_umi_len'] = ';'.join([str(i) for i in multimer_umi_lens])

        sc_bam_ids = []
        for i in shortlisted_multimer_umi_list:
            sc_bam_ids.append(series_uniq_umi.loc[series_uniq_umi.isin([i])].index[0])

        self.console.print('======>start writing deduplicated reads to BAM...')
        dedup_reads_write_stime = time.time()
        # print(self.work_dir)

        import os
        from umiche.util.Folder import Folder as crtfolder
        crtfolder().osmkdir(DIRECTORY=os.path.dirname(kwargs['sv_interm_bam_fpn']))

        self.aliwriter.tobam(
            tobam_fpn=kwargs['sv_interm_bam_fpn'],
            tmpl_bam_fpn=self.bam_fpn,
            whitelist=sc_bam_ids,
        )
        self.console.print('======>finish writing in {:.2f}s'.format(time.time() - dedup_reads_write_stime))
        return self.df

    def majority_vote(
            self,
            **kwargs,
    ):
        self.df_umi_uniq = self.df_bam.drop_duplicates(subset=['umi'], keep='first')
        # print(self.df_umi_uniq)
        series_uniq_umi = self.df_umi_uniq.umi
        # print(series_uniq_umi)
        self.umi_to_int_dict = {k: id for id, k in enumerate(series_uniq_umi)}

        dedup_cnt, uniq_multimer_cnt, shortlisted_multimer_umi_list = umimv().track(
            multimer_list=series_uniq_umi.values,
            recur_len=kwargs['umi_unit_pattern'],
        )
        # print(dedup_cnt)
        self.df.loc[0, 'dedup_cnt'] = dedup_cnt
        sc_bam_ids = []
        for i in shortlisted_multimer_umi_list:
            sc_bam_ids.append(series_uniq_umi.loc[series_uniq_umi.isin([i])].index[0])

        self.console.print('======>start writing deduplicated reads to BAM...')
        dedup_reads_write_stime = time.time()
        # print(self.work_dir)

        import os
        from umiche.util.Folder import Folder as crtfolder
        crtfolder().osmkdir(DIRECTORY=os.path.dirname(kwargs['sv_interm_bam_fpn']))

        self.aliwriter.tobam(
            tobam_fpn=kwargs['sv_interm_bam_fpn'],
            tmpl_bam_fpn=self.bam_fpn,
            whitelist=sc_bam_ids,
        )
        self.console.print('======>finish writing in {:.2f}s'.format(time.time() - dedup_reads_write_stime))
        return self.df

    def unique(self, ):
        dedup_umi_stime = time.time()
        self.df['uniq_sgl_mark'] = self.df['uniq_repr_nodes'].apply(lambda x: 'yes' if len(x) == 1 else 'no')
        self.df = self.df.loc[self.df['uniq_sgl_mark'] == 'no']
        self.console.print('======># of positions with non-single umis: {}'.format(self.df.shape[0]))
        self.console.print('======>finish finding deduplicated UMIs in {:.2f}s'.format(time.time() - dedup_umi_stime))
        # self.console.print('======># of UMIs deduplicated {}'.format(self.df['num_uniq_umis'].loc['yes']))
        self.console.print('======>calculate average edit distances between umis...')
        ### @@ self.df['ave_ed']
        # 1    5.0
        # Name: ave_eds, dtype: float64
        self.df['ave_eds'] = self.df.apply(lambda x: self.umigadgetry.ed_ave(x, by_col='uniq_repr_nodes'), axis=1)
        self.ave_ed_bins = self.df['ave_eds'].value_counts().sort_index().to_frame().reset_index()
        self.console.check("======>bins for average edit distance\n{}".format(self.ave_ed_bins))
        self.df['num_diff_dedup_uniq_umis'] = self.df.apply(
            lambda x: self.umigadgetry.num_removed_uniq_umis(x, by_col='uniq_repr_nodes'), axis=1)
        self.df['num_diff_dedup_reads'] = self.df.apply(
            lambda x: self.umigadgetry.num_removed_reads(x, by_col='uniq_repr_nodes'),
            axis=1)
        self.console.print('======># of deduplicated unique umis {}'.format(self.df['num_diff_dedup_uniq_umis'].sum()))
        self.console.print('======># of deduplicated reads {}'.format(self.df['num_diff_dedup_reads'].sum()))
        if not self.heterogeneity:
            self.fwriter.generic(
                df=self.ave_ed_bins,
                sv_fpn=self.work_dir + 'unique_ave_ed_bin.txt',
                index=True,
                header=True,
            )
            self.fwriter.generic(
                df=self.df[[
                    'ave_ed',
                    'num_uniq_umis',
                    'num_diff_dedup_uniq_umis',
                    'num_diff_dedup_reads',
                ]],
                sv_fpn=self.work_dir + 'unique_dedup_sum.txt',
                index=True,
                header=True,
            )
            self.console.print('======>start writing deduplicated reads to BAM...')
            dedup_reads_write_stime = time.time()
            self.df['uniq_bam_ids'] = self.df.apply(lambda x: self.umigadgetry.bamids(x, by_col='uniq_repr_nodes'), axis=1)
            self.aliwriter.tobam(
                tobam_fpn=self.work_dir + 'unique_dedup.bam',
                tmpl_bam_fpn=self.bam_fpn,
                whitelist=self.umigadgetry.decompose(list_nd=self.df['uniq_bam_ids'].values),
            )
            self.console.print('======>finish writing in {:.2f}s'.format(time.time() - dedup_reads_write_stime))
        return self.df

    def cluster(self, ):
        dedup_umi_stime = time.time()
        self.df['cc_repr_nodes'] = self.df.apply(lambda x: self.umigadgetry.umimax(x, by_col='cc'), axis=1)
        ### @@ self.df['cc_repr_nodes']
        # 1    [2]
        # Name: cc_repr_nodes, dtype: object
        self.df['dedup_cnt'] = self.df['cc_repr_nodes'].apply(lambda x: self.umigadgetry.length(x))
        ### @@ self.df['dedup_cnt']
        # 1    1
        # Name: dedup_cnt, dtype: int64
        self.console.print('======>finish finding deduplicated umis in {:.2f}s'.format(time.time() - dedup_umi_stime))
        # self.console.print('======># of umis deduplicated to be {}'.format(self.df['dedup_cnt'].loc['yes']))
        self.console.print('======>calculate average edit distances between umis...')
        self.df['ave_eds'] = self.df.apply(lambda x: self.umigadgetry.ed_ave(x, by_col='cc_repr_nodes'), axis=1)
        self.df['num_diff_dedup_uniq_umis'] = self.df.apply(
            lambda x: self.umigadgetry.num_removed_uniq_umis(x, by_col='cc_repr_nodes'),
            axis=1)
        self.df['num_diff_dedup_reads'] = self.df.apply(
            lambda x: self.umigadgetry.num_removed_reads(x, by_col='cc_repr_nodes'),
            axis=1)
        self.console.print('======># of deduplicated unique umis {}'.format(self.df['num_diff_dedup_uniq_umis'].sum()))
        self.console.print('======># of deduplicated reads {}'.format(self.df['num_diff_dedup_reads'].sum()))
        self.ave_ed_bins = self.df['ave_eds'].value_counts().sort_index().to_frame().reset_index()
        self.console.check("======>bins for average edit distance\n{}".format(self.ave_ed_bins))
        if not self.heterogeneity:
            self.fwriter.generic(
                df=self.ave_ed_bins,
                sv_fpn=self.work_dir + 'cluster_ave_ed_bin.txt',
                index=True,
                header=True,
            )
            self.fwriter.generic(
                df=self.df[[
                    'dedup_cnt',
                    'ave_ed',
                    'num_uniq_umis',
                    'num_diff_dedup_uniq_umis',
                    'num_diff_dedup_reads',
                ]],
                sv_fpn=self.work_dir + 'cluster_dedup_sum.txt',
                index=True,
                header=True,
            )
            self.console.print('======>start writing deduplicated reads to BAM...')
            dedup_reads_write_stime = time.time()
            self.df['cc_bam_ids'] = self.df.apply(lambda x: self.umigadgetry.bamids(x, by_col='cc_repr_nodes'), axis=1)
            self.aliwriter.tobam(
                tobam_fpn=self.work_dir + 'cluster_dedup.bam',
                tmpl_bam_fpn=self.bam_fpn,
                whitelist=self.umigadgetry.decompose(list_nd=self.df['cc_bam_ids'].values),
            )
            self.console.print('======>finish writing in {:.2f}s'.format(time.time() - dedup_reads_write_stime))
        return self.df

    def adjacency(self, ):
        dedup_umi_stime = time.time()
        umiadj_ob = umiadj()
        self.df['adj'] = self.df.apply(
            lambda x: umiadj_ob.decompose(
                cc_sub_dict=umiadj_ob.umi_tools(
                    connected_components=x['cc'],
                    df_umi_uniq_val_cnt=x['vignette']['df_umi_uniq_val_cnt'],
                    graph_adj=x['vignette']['graph_adj'],
                )['clusters'],
            ),
            axis=1,
        )
        self.df['adj_repr_nodes'] = self.df.apply(lambda x: self.umigadgetry.umimax(x, by_col='adj'), axis=1)
        self.df['dedup_cnt'] = self.df['adj_repr_nodes'].apply(lambda x: self.umigadgetry.length(x))
        self.console.print('======>finish finding deduplicated umis in {:.2f}s'.format(time.time() - dedup_umi_stime))
        # self.console.print('======># of umis deduplicated to be {}'.format(self.df['dedup_cnt'].loc['yes']))
        self.console.print('======>calculate average edit distances between umis...')
        self.df['ave_eds'] = self.df.apply(lambda x: self.umigadgetry.ed_ave(x, by_col='adj_repr_nodes'), axis=1)
        self.df['num_diff_dedup_uniq_umis'] = self.df.apply(
            lambda x: self.umigadgetry.num_removed_uniq_umis(
                x,
                by_col='adj_repr_nodes',
            ),
            axis=1,
        )
        self.df['num_diff_dedup_reads'] = self.df.apply(
            lambda x: self.umigadgetry.num_removed_reads(
                x,
                by_col='adj_repr_nodes',
            ),
            axis=1,
        )
        self.console.print('======># of deduplicated unique umis {}'.format(self.df['num_diff_dedup_uniq_umis'].sum()))
        self.console.print('======># of deduplicated reads {}'.format(self.df['num_diff_dedup_reads'].sum()))
        self.ave_ed_bins = self.df['ave_eds'].value_counts().sort_index().to_frame().reset_index()
        self.console.check("======>bins for average edit distance\n{}".format(self.ave_ed_bins))
        if not self.heterogeneity:
            self.fwriter.generic(
                df=self.ave_ed_bins,
                sv_fpn=self.work_dir + 'adjacency_ave_ed_bin.txt',
                index=True,
                header=True,
            )
            self.fwriter.generic(
                df=self.df[[
                    'dedup_cnt',
                    'ave_ed',
                    'num_uniq_umis',
                    'num_diff_dedup_uniq_umis',
                    'num_diff_dedup_reads',
                ]],
                sv_fpn=self.work_dir + 'adjacency_dedup_sum.txt',
                index=True,
                header=True,
            )
            self.console.print('======>start writing deduplicated reads to BAM...')
            dedup_reads_write_stime = time.time()
            self.df['adj_bam_ids'] = self.df.apply(lambda x: self.umigadgetry.bamids(x, by_col='adj_repr_nodes'), axis=1)
            self.aliwriter.tobam(
                tobam_fpn=self.work_dir + 'adjacency_dedup.bam',
                tmpl_bam_fpn=self.bam_fpn,
                whitelist=self.umigadgetry.decompose(list_nd=self.df['adj_bam_ids'].values),
            )
            self.console.print('======>finish writing in {:.2f}s'.format(time.time() - dedup_reads_write_stime))
        return self.df

    def directional(self, ):
        dedup_umi_stime = time.time()
        umidirec_ob = umidirec(self.heterogeneity)
        # self.df[['count', 'clusters', 'apv', 'disapv']] = self.df.apply(
        #     lambda x: umidirec_ob.umi_tools(
        #         connected_components=x['cc'],
        #         df_umi_uniq_val_cnt=x['vignette']['df_umi_uniq_val_cnt'],
        #         graph_adj=x['vignette']['graph_adj'],
        #     ),
        #     axis=1,
        #     result_type='expand',
        # )
        if self.heterogeneity:
            self.df['count'], self.df['clusters'], self.df['apv'], self.df['disapv'] = zip(
                *self.df.apply(
                    lambda x: umidirec_ob.umi_tools(
                        connected_components=x['cc'],
                        df_umi_uniq_val_cnt=x['vignette']['df_umi_uniq_val_cnt'],
                        graph_adj=x['vignette']['graph_adj'],
                    ),
                    axis=1,
                )
            )
            self.df['direc'] = self.df.apply(
                lambda x: umidirec_ob.decompose(
                    cc_sub_dict=x['clusters'],
                ),
                axis=1,
            )
        else:
            self.df['direc'] = self.df.apply(
                lambda x: umidirec_ob.decompose(
                    cc_sub_dict=umidirec_ob.umi_tools(
                        connected_components=x['cc'],
                        df_umi_uniq_val_cnt=x['vignette']['df_umi_uniq_val_cnt'],
                        graph_adj=x['vignette']['graph_adj'],
                    )['clusters'],
                ),
                axis=1,
            )
        # print(self.df['direc'])
        self.df['direc_repr_nodes'] = self.df.apply(lambda x: self.umigadgetry.umimax(x, by_col='direc'), axis=1)
        # print(self.df['direc_repr_nodes'])
        self.df['dedup_cnt'] = self.df['direc_repr_nodes'].apply(lambda x: self.umigadgetry.length(x))
        self.console.print('======>finish finding deduplicated umis in {:.2f}s'.format(time.time() - dedup_umi_stime))
        # self.console.print('======># of umis deduplicated to be {}'.format(self.df['dedup_cnt'].loc['yes']))
        self.console.print('======>calculate average edit distances between umis...')
        self.df['ave_eds'] = self.df.apply(lambda x: self.umigadgetry.ed_ave(x, by_col='direc_repr_nodes'), axis=1)
        # print(self.df['ave_eds'])
        self.df['num_diff_dedup_uniq_umis'] = self.df.apply(
            lambda x: self.umigadgetry.num_removed_uniq_umis(x, by_col='direc_repr_nodes'),
            axis=1,
        )
        self.df['num_diff_dedup_reads'] = self.df.apply(
            lambda x: self.umigadgetry.num_removed_reads(x, by_col='direc_repr_nodes'),
            axis=1,
        )
        self.console.print('======># of deduplicated unique umis {}'.format(self.df['num_diff_dedup_uniq_umis'].sum()))
        self.console.print('======># of deduplicated reads {}'.format(self.df['num_diff_dedup_reads'].sum()))
        self.ave_ed_bins = self.df['ave_eds'].value_counts().sort_index().to_frame().reset_index()
        # print(self.ave_ed_bins)
        self.console.check("======>bins for average edit distance\n{}".format(self.ave_ed_bins))
        if not self.heterogeneity:
            self.fwriter.generic(
                df=self.ave_ed_bins,
                sv_fpn=self.work_dir + 'directional_ave_ed_bin.txt',
                index=True,
                header=True,
            )
            self.fwriter.generic(
                df=self.df[[
                    'dedup_cnt',
                    'ave_ed',
                    'num_uniq_umis',
                    'num_diff_dedup_uniq_umis',
                    'num_diff_dedup_reads',
                ]],
                sv_fpn=self.work_dir + 'directional_dedup_sum.txt',
                index=True,
                header=True,
            )
            self.console.print('======>start writing deduplicated reads to BAM...')
            dedup_reads_write_stime = time.time()
            self.df['direc_bam_ids'] = self.df.apply(lambda x: self.umigadgetry.bamids(x, by_col='direc_repr_nodes'), axis=1)
            self.aliwriter.tobam(
                tobam_fpn=self.work_dir + 'directional_dedup.bam',
                tmpl_bam_fpn=self.bam_fpn,
                whitelist=self.umigadgetry.decompose(list_nd=self.df['direc_bam_ids'].values),
            )
            self.console.print('======>finish writing in {:.2f}s'.format(time.time() - dedup_reads_write_stime))
        return self.df

    def mcl(
            self,
            **kwargs
    ):
        # print(self.df.columns)
        # print(self.df)
        dedup_umi_stime = time.time()
        umimcl_ob = umimcl(
            inflat_val=kwargs['inflat_val'],
            exp_val=kwargs['exp_val'],
            iter_num=kwargs['iter_num'],
            heterogeneity=self.heterogeneity,
        )
        # print(self.df)
        ### @@ please note that df and df_mcl_res cannot be merged becuase the dimension is not the same.
        if self.heterogeneity:
            df_mcl_res = self.df.apply(
                lambda x: umimcl_ob.dfclusters(
                    connected_components=x['cc'],
                    graph_adj=x['vignette']['graph_adj'],
                ),
                axis=1,
            ).values[0]
            mcl_dict = {'mcl': umimcl_ob.decompose(
                list_nd=df_mcl_res['clusters'].values
            )}
            self.df['mcl'] = self.df.apply(lambda x: mcl_dict['mcl'], axis=1)
            apv_dict = {'apv': [df_mcl_res['apv']]}
            self.df['apv'] = self.df.apply(lambda x: apv_dict['apv'], axis=1)
            # print(self.df['mcl'])
        else:
            self.df['mcl'] = self.df.apply(
                lambda x: umimcl_ob.decompose(
                    list_nd=umimcl_ob.dfclusters(
                        connected_components=x['cc'],
                        graph_adj=x['vignette']['graph_adj'],
                    )['clusters'].values,
                ),
                axis=1,
            )
            # print(self.df['mcl'])
        ### @@ self.df['mcl']
        # 1    {0: [0, 76, 162, 188, 237, 256], 1: [65, 55, 1...
        # Name: mcl, dtype: object
        self.df['mcl_repr_nodes'] = self.df.apply(lambda x: self.umigadgetry.umimax(x, by_col='mcl'), axis=1)
        self.df['dedup_cnt'] = self.df['mcl_repr_nodes'].apply(lambda x: self.umigadgetry.length(x))
        self.console.print('======>finish finding deduplicated umis in {:.2f}s'.format(time.time() - dedup_umi_stime))
        # self.console.print('======># of umis deduplicated to be {}'.format(self.df['dedup_cnt'].loc['yes']))
        self.console.print('======>calculate average edit distances between umis...')
        self.df['ave_eds'] = self.df.apply(lambda x: self.umigadgetry.ed_ave(x, by_col='mcl_repr_nodes'), axis=1)
        self.df['num_diff_dedup_uniq_umis'] = self.df.apply(lambda x: self.umigadgetry.num_removed_uniq_umis(x, by_col='mcl_repr_nodes'),axis=1)
        self.df['num_diff_dedup_reads'] = self.df.apply(lambda x: self.umigadgetry.num_removed_reads(x, by_col='mcl_repr_nodes'),axis=1)
        self.console.print('======># of deduplicated unique umis {}'.format(self.df['num_diff_dedup_uniq_umis'].sum()))
        self.console.print('======># of deduplicated reads {}'.format(self.df['num_diff_dedup_reads'].sum()))
        self.ave_ed_bins = self.df['ave_eds'].value_counts().sort_index().to_frame().reset_index()
        self.console.check("======>bins for average edit distance\n{}".format(self.ave_ed_bins))
        if not self.heterogeneity:
            self.fwriter.generic(
                df=self.ave_ed_bins,
                sv_fpn=self.work_dir + 'mcl_ave_ed_bin.txt',
                index=True,
                header=True,
            )
            self.fwriter.generic(
                df=self.df[[
                    'dedup_cnt',
                    'ave_ed',
                    'num_uniq_umis',
                    'num_diff_dedup_uniq_umis',
                    'num_diff_dedup_reads',
                ]],
                sv_fpn=self.work_dir + 'mcl_dedup_sum.txt',
                index=True,
                header=True,
            )
            self.console.print('======>start writing deduplicated reads to BAM...')
            dedup_reads_write_stime = time.time()
            self.df['mcl_bam_ids'] = self.df.apply(lambda x: self.umigadgetry.bamids(x, by_col='mcl_repr_nodes'), axis=1)
            self.aliwriter.tobam(
                tobam_fpn=self.work_dir + 'mcl_dedup.bam',
                tmpl_bam_fpn=self.bam_fpn,
                whitelist=self.umigadgetry.decompose(list_nd=self.df['mcl_bam_ids'].values),
            )
            self.console.print('======>finish writing in {:.2f}s'.format(time.time() - dedup_reads_write_stime))
        return self.df

    def mcl_cc_all_node_umis(
            self,
            **kwargs,
    ):
        dedup_umi_stime = time.time()
        umimcl_ob = umimcl(
            inflat_val=kwargs['inflat_val'],
            exp_val=kwargs['exp_val'],
            iter_num=kwargs['iter_num'],
            heterogeneity=self.heterogeneity,
        )
        # print(self.df)
        ### @@ please note that df and df_mcl_res cannot be merged becuase the dimension is not the same.
        df_mcl_res = self.df.apply(
            lambda x: umimcl_ob.dfclusters_cc_all_node_umis(
                int_to_umi_dict=x['vignette']['int_to_umi_dict'],
                graph_adj=x['vignette']['graph_adj'],
            ),
            axis=1,
        ).values[0]
        mcl_dict = {'mcl': umimcl_ob.decompose(
            list_nd=df_mcl_res['clusters'].values
        )}
        apv_dict = {'apv': [df_mcl_res['apv']]}
        self.df['mcl'] = self.df.apply(lambda x: mcl_dict['mcl'], axis=1)
        self.df['apv'] = self.df.apply(lambda x: apv_dict['apv'], axis=1)
        # print(self.df['apv'])
        # self.df['mcl'] = self.df.apply(
        #     lambda x: umimcl_ob.decompose(
        #         list_nd=umimcl_ob.dfclusters(
        #             connected_components=x['cc'],
        #             graph_adj=x['vignette']['graph_adj'],
        #         )['clusters'].values,
        #     ),
        #     axis=1,
        # )
        ### @@ self.df['mcl']
        # 1    {0: [0, 76, 162, 188, 237, 256], 1: [65, 55, 1...
        # Name: mcl, dtype: object
        self.df['mcl_repr_nodes'] = self.df.apply(lambda x: self.umigadgetry.umimax(x, by_col='mcl'), axis=1)
        self.df['dedup_cnt'] = self.df['mcl_repr_nodes'].apply(lambda x: self.umigadgetry.length(x))
        self.console.print('======>finish finding deduplicated umis in {:.2f}s'.format(time.time() - dedup_umi_stime))
        # self.console.print('======># of umis deduplicated to be {}'.format(self.df['dedup_cnt'].loc['yes']))
        self.console.print('======>calculate average edit distances between umis...')
        self.df['ave_eds'] = self.df.apply(lambda x: self.umigadgetry.ed_ave(x, by_col='mcl_repr_nodes'), axis=1)
        self.df['num_diff_dedup_uniq_umis'] = self.df.apply(lambda x: self.umigadgetry.num_removed_uniq_umis(x, by_col='mcl_repr_nodes'),axis=1)
        self.df['num_diff_dedup_reads'] = self.df.apply(lambda x: self.umigadgetry.num_removed_reads(x, by_col='mcl_repr_nodes'),axis=1)
        self.console.print('======># of deduplicated unique umis {}'.format(self.df['num_diff_dedup_uniq_umis'].sum()))
        self.console.print('======># of deduplicated reads {}'.format(self.df['num_diff_dedup_reads'].sum()))
        self.ave_ed_bins = self.df['ave_eds'].value_counts().sort_index().to_frame().reset_index()
        self.console.check("======>bins for average edit distance\n{}".format(self.ave_ed_bins))
        if not self.heterogeneity:
            self.fwriter.generic(
                df=self.ave_ed_bins,
                sv_fpn=self.work_dir + 'mcl_ave_ed_bin.txt',
                index=True,
                header=True,
            )
            self.fwriter.generic(
                df=self.df[[
                    'dedup_cnt',
                    'ave_ed',
                    'num_uniq_umis',
                    'num_diff_dedup_uniq_umis',
                    'num_diff_dedup_reads',
                ]],
                sv_fpn=self.work_dir + 'mcl_dedup_sum.txt',
                index=True,
                header=True,
            )
            self.console.print('======>start writing deduplicated reads to BAM...')
            dedup_reads_write_stime = time.time()
            self.df['mcl_bam_ids'] = self.df.apply(lambda x: self.umigadgetry.bamids(x, by_col='mcl_repr_nodes'), axis=1)
            self.aliwriter.tobam(
                tobam_fpn=self.work_dir + 'mcl_dedup.bam',
                tmpl_bam_fpn=self.bam_fpn,
                whitelist=self.umigadgetry.decompose(list_nd=self.df['mcl_bam_ids'].values),
            )
            self.console.print('======>finish writing in {:.2f}s'.format(time.time() - dedup_reads_write_stime))
        return self.df

    def mcl_val(
            self,
            **kwargs,
    ):
        dedup_umi_stime = time.time()
        umimcl_ob = umimcl(
            inflat_val=kwargs['inflat_val'],
            exp_val=kwargs['exp_val'],
            iter_num=kwargs['iter_num'],
            heterogeneity=self.heterogeneity,
        )
        if self.heterogeneity:
            self.df['count'], self.df['clusters'], self.df['apv'], self.df['disapv'] = zip(
                *self.df.apply(
                    lambda x: umimcl_ob.maxval_val(
                        df_mcl_ccs=umimcl_ob.dfclusters(
                            connected_components=x['cc'],
                            graph_adj=x['vignette']['graph_adj'],
                        ),
                        df_umi_uniq_val_cnt=x['vignette']['df_umi_uniq_val_cnt'],
                        thres_fold=kwargs['mcl_fold_thres'],
                    ),
                    axis=1,
                )
            )
            self.df['mcl_val'] = self.df.apply(
                lambda x: umimcl_ob.decompose(
                    list_nd=x['clusters'].values,
                ),
                axis=1,
            )
            # print(self.df['mcl_val'])
        else:
            self.df['mcl_val'] = self.df.apply(
                lambda x: umimcl_ob.decompose(
                    list_nd=umimcl_ob.maxval_val(
                        df_mcl_ccs=umimcl_ob.dfclusters(
                            connected_components=x['cc'],
                            graph_adj=x['vignette']['graph_adj'],
                        ),
                        df_umi_uniq_val_cnt=x['vignette']['df_umi_uniq_val_cnt'],
                        thres_fold=kwargs['mcl_fold_thres'],
                    )['clusters'].values,
                ),
                axis=1,
            )
            # print(self.df['mcl_val'])
        self.df['mcl_val_repr_nodes'] = self.df.apply(lambda x: self.umigadgetry.umimax(x, by_col='mcl_val'), axis=1)
        self.df['dedup_cnt'] = self.df['mcl_val_repr_nodes'].apply(lambda x: self.umigadgetry.length(x))
        self.console.print('======>finish finding deduplicated umis in {:.2f}s'.format(time.time() - dedup_umi_stime))
        # self.console.print('======># of umis deduplicated to be {}'.format(self.df['dedup_cnt'].loc['yes']))
        self.console.print('======>calculate average edit distances between umis...')
        self.df['ave_eds'] = self.df.apply(lambda x: self.umigadgetry.ed_ave(x, by_col='mcl_val_repr_nodes'), axis=1)
        self.df['num_diff_dedup_uniq_umis'] = self.df.apply(
            lambda x: self.umigadgetry.num_removed_uniq_umis(x, by_col='mcl_val_repr_nodes'), axis=1)
        self.df['num_diff_dedup_reads'] = self.df.apply(
            lambda x: self.umigadgetry.num_removed_reads(x, by_col='mcl_val_repr_nodes'), axis=1)
        self.console.print('======># of deduplicated unique umis {}'.format(self.df['num_diff_dedup_uniq_umis'].sum()))
        self.console.print('======># of deduplicated reads {}'.format(self.df['num_diff_dedup_reads'].sum()))
        self.ave_ed_bins = self.df['ave_eds'].value_counts().sort_index().to_frame().reset_index()
        self.console.check("======>bins for average edit distance\n{}".format(self.ave_ed_bins))
        if not self.heterogeneity:
            self.fwriter.generic(
                df=self.ave_ed_bins,
                sv_fpn=self.work_dir + 'mcl_val_ave_ed_bin.txt',
                index=True,
                header=True,
            )
            self.fwriter.generic(
                df=self.df[[
                    'dedup_cnt',
                    'ave_ed',
                    'num_uniq_umis',
                    'num_diff_dedup_uniq_umis',
                    'num_diff_dedup_reads',
                ]],
                sv_fpn=self.work_dir + 'mcl_val_dedup_sum.txt',
                index=True,
                header=True,
            )
            self.console.print('======>start writing deduplicated reads to BAM...')
            dedup_reads_write_stime = time.time()
            self.df['mcl_val_bam_ids'] = self.df.apply(lambda x: self.umigadgetry.bamids(x, by_col='mcl_val_repr_nodes'),
                                                       axis=1)
            self.aliwriter.tobam(
                tobam_fpn=self.work_dir + 'mcl_val_dedup.bam',
                tmpl_bam_fpn=self.bam_fpn,
                whitelist=self.umigadgetry.decompose(list_nd=self.df['mcl_val_bam_ids'].values),
            )
            self.console.print('======>finish writing in {:.2f}s'.format(time.time() - dedup_reads_write_stime))
        return self.df

    def mcl_ed(
            self,
            **kwargs,
    ):
        dedup_umi_stime = time.time()
        umimcl_ob = umimcl(
            inflat_val=kwargs['inflat_val'],
            exp_val=kwargs['exp_val'],
            iter_num=kwargs['iter_num'],
            heterogeneity=self.heterogeneity,
        )
        if self.heterogeneity:
            # self.df[['count', 'clusters', 'apv', 'disapv']]
            self.df['count'], self.df['clusters'], self.df['apv'], self.df['disapv'] = zip(
                *self.df.apply(
                    lambda x: umimcl_ob.maxval_ed(
                        df_mcl_ccs=umimcl_ob.dfclusters(
                            connected_components=x['cc'],
                            graph_adj=x['vignette']['graph_adj'],
                        ),
                        df_umi_uniq_val_cnt=x['vignette']['df_umi_uniq_val_cnt'],
                        int_to_umi_dict=x['vignette']['int_to_umi_dict'],
                        thres_fold=kwargs['mcl_fold_thres'],
                    ),
                    axis=1,
                )
            )
            self.df['mcl_ed'] = self.df.apply(
                lambda x: umimcl_ob.decompose(
                    list_nd=x['clusters'].values,
                    # list_nd=x['clusters'].values[0].tolist(),
                ),
                axis=1,
            )
            # print(self.df['mcl_ed'])
        else:
            self.df['mcl_ed'] = self.df.apply(
                lambda x: umimcl_ob.decompose(
                    list_nd=umimcl_ob.maxval_ed(
                        df_mcl_ccs=umimcl_ob.dfclusters(
                            connected_components=x['cc'],
                            graph_adj=x['vignette']['graph_adj'],
                        ),
                        df_umi_uniq_val_cnt=x['vignette']['df_umi_uniq_val_cnt'],
                        int_to_umi_dict=x['vignette']['int_to_umi_dict'],
                        thres_fold=kwargs['mcl_fold_thres'],
                    )['clusters'].values,
                ),
                axis=1,
            )
        self.df['mcl_ed_repr_nodes'] = self.df.apply(lambda x: self.umigadgetry.umimax(x, by_col='mcl_ed'), axis=1)
        self.df['dedup_cnt'] = self.df['mcl_ed_repr_nodes'].apply(lambda x: self.umigadgetry.length(x))
        # print(self.df['dedup_cnt'])
        self.console.print('======>finish finding deduplicated umis in {:.2f}s'.format(time.time() - dedup_umi_stime))
        # self.console.print('======># of umis deduplicated to be {}'.format(self.df['dedup_cnt'].loc['yes']))
        self.console.print('======>calculate average edit distances between umis...')
        self.df['ave_eds'] = self.df.apply(lambda x: self.umigadgetry.ed_ave(x, by_col='mcl_ed_repr_nodes'), axis=1)
        self.df['num_diff_dedup_uniq_umis'] = self.df.apply(
            lambda x: self.umigadgetry.num_removed_uniq_umis(x, by_col='mcl_ed_repr_nodes'), axis=1)
        self.df['num_diff_dedup_reads'] = self.df.apply(
            lambda x: self.umigadgetry.num_removed_reads(x, by_col='mcl_ed_repr_nodes'),
            axis=1)
        self.console.print('======># of deduplicated unique umis {}'.format(self.df['num_diff_dedup_uniq_umis'].sum()))
        self.console.print('======># of deduplicated reads {}'.format(self.df['num_diff_dedup_reads'].sum()))
        self.ave_ed_bins = self.df['ave_eds'].value_counts().sort_index().to_frame().reset_index()
        self.console.check("======>bins for average edit distance\n{}".format(self.ave_ed_bins))
        if not self.heterogeneity:
            self.fwriter.generic(
                df=self.ave_ed_bins,
                sv_fpn=self.work_dir + 'mcl_ed_ave_ed_bin.txt',
                index=True,
                header=True,
            )
            self.fwriter.generic(
                df=self.df[[
                    'dedup_cnt',
                    'ave_ed',
                    'num_uniq_umis',
                    'num_diff_dedup_uniq_umis',
                    'num_diff_dedup_reads',
                ]],
                sv_fpn=self.work_dir + 'mcl_ed_dedup_sum.txt',
                index=True,
                header=True,
            )
            self.console.print('======>start writing deduplicated reads to BAM...')
            dedup_reads_write_stime = time.time()
            self.df['mcl_ed_bam_ids'] = self.df.apply(lambda x: self.umigadgetry.bamids(x, by_col='mcl_ed_repr_nodes'), axis=1)
            self.aliwriter.tobam(
                tobam_fpn=self.work_dir + 'mcl_ed_dedup.bam',
                tmpl_bam_fpn=self.bam_fpn,
                whitelist=self.umigadgetry.decompose(list_nd=self.df['mcl_ed_bam_ids'].values),
            )
            self.console.print('======>finish writing in {:.2f}s'.format(time.time() - dedup_reads_write_stime))
        return self.df

    def clustering_umi_seq_onehot(
            self,
            clustering_method,
            **kwargs
    ):
        dedup_umi_stime = time.time()
        self.umiclustering = umiclustering(
            clustering_method=clustering_method,
            **kwargs,
        )
        df_clustering_res = self.df.apply(
            lambda x: self.umiclustering.dfclusters(
                connected_components=x['cc'],
                graph_adj=x['vignette']['graph_adj'],
                int_to_umi_dict=x['vignette']['int_to_umi_dict'],
            ),
            axis=1,
        ).values[0]
        mcl_dict = {clustering_method: self.umiclustering.decompose(
            list_nd=df_clustering_res['clusters'].values
        )}
        apv_dict = {'apv': [df_clustering_res['apv']]}
        self.df[clustering_method] = self.df.apply(lambda x: mcl_dict[clustering_method], axis=1)
        self.df['apv'] = self.df.apply(lambda x: apv_dict['apv'], axis=1)
        # print(self.df['apv'])
        self.df[clustering_method + '_repr_nodes'] = self.df.apply(lambda x: self.umigadgetry.umimax(x, by_col=clustering_method), axis=1)
        self.df['dedup_cnt'] = self.df[clustering_method + '_repr_nodes'].apply(lambda x: self.umigadgetry.length(x))
        self.console.print('======>finish finding deduplicated umis in {:.2f}s'.format(time.time() - dedup_umi_stime))
        # self.console.print('======># of umis deduplicated to be {}'.format(self.df['dedup_cnt'].loc['yes']))
        self.console.print('======>calculate average edit distances between umis...')
        self.df['ave_eds'] = self.df.apply(lambda x: self.umigadgetry.ed_ave(x, by_col=clustering_method + '_repr_nodes'), axis=1)
        self.df['num_diff_dedup_uniq_umis'] = self.df.apply(
            lambda x: self.umigadgetry.num_removed_uniq_umis(x, by_col=clustering_method + '_repr_nodes'),
            axis=1)
        self.df['num_diff_dedup_reads'] = self.df.apply(
            lambda x: self.umigadgetry.num_removed_reads(x, by_col=clustering_method + '_repr_nodes'),
            axis=1)
        self.console.print('======># of deduplicated unique umis {}'.format(self.df['num_diff_dedup_uniq_umis'].sum()))
        self.console.print('======># of deduplicated reads {}'.format(self.df['num_diff_dedup_reads'].sum()))
        self.ave_ed_bins = self.df['ave_eds'].value_counts().sort_index().to_frame().reset_index()
        self.console.check("======>bins for average edit distance\n{}".format(self.ave_ed_bins))
        if not self.heterogeneity:
            self.fwriter.generic(
                df=self.ave_ed_bins,
                sv_fpn=self.work_dir + clustering_method + '_ave_ed_bin.txt',
                index=True,
                header=True,
            )
            self.fwriter.generic(
                df=self.df[[
                    'dedup_cnt',
                    'ave_ed',
                    'num_uniq_umis',
                    'num_diff_dedup_uniq_umis',
                    'num_diff_dedup_reads',
                ]],
                sv_fpn=self.work_dir + clustering_method + '_dedup_sum.txt',
                index=True,
                header=True,
            )
            self.console.print('======>start writing deduplicated reads to BAM...')
            dedup_reads_write_stime = time.time()
            self.df[clustering_method + '_bam_ids'] = self.df.apply(lambda x: self.umigadgetry.bamids(x, by_col=clustering_method + '_repr_nodes'), axis=1)
            self.aliwriter.tobam(
                tobam_fpn=self.work_dir + clustering_method + '_dedup.bam',
                tmpl_bam_fpn=self.bam_fpn,
                whitelist=self.umigadgetry.decompose(list_nd=self.df[clustering_method + '_bam_ids'].values),
            )
            self.console.print('======>finish writing in {:.2f}s'.format(time.time() - dedup_reads_write_stime))
        return self.df