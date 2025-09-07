__version__ = "v1.0"
__copyright__ = "Copyright 2025"
__license__ = "GPL-3.0"
__developer__ = "Jianfeng Sun"
__maintainer__ = "Jianfeng Sun"
__email__ = "jianfeng.sunmt@gmail.com"

import os
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

from umiche.deduplicate.method.PCRArtefact import PCRArtefact as pcrartefact

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

        self.bam_fn = os.path.splitext(os.path.basename(self.bam_fpn))[0]

        self.aliwriter = aliwriter(df=self.df_bam)

        self.umigadgetry = umigadgetry()
        self.pcrartefact = pcrartefact(verbose=verbose)

        self.fwriter = fwriter()

        self.console = Console()
        self.verbose = verbose
        self.console.verbose = self.verbose

    def set_cover_depracated(
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

    @Console.vignette()
    def set_cover(
            self,
            **kwargs,
    ):
        self.df['clusters'] = self.df.apply(
            lambda x: umisc(verbose=self.verbose).greedy(
                multimer_list=x['vignette']['umi'],
                recur_len=kwargs['umi_unit_pattern'],
                # split_method='split_to_all',
                split_method='split_to_',
            )['clusters'],
            axis=1,
        )
        # print(self.df['clusters'])
        self.df['repr_node_map'] = self.df.apply(lambda x: self.umigadgetry.umimax(x, by_col='clusters', met='map'), axis=1)
        # print(self.df['repr_node_map'])
        self.df['repr_nodes'] = self.df.apply(lambda x: self.umigadgetry.umimax(x, by_col='clusters'), axis=1)
        # print(self.df['repr_nodes'])
        self.df['dedup_cnt'] = self.df['repr_nodes'].apply(lambda x: self.umigadgetry.length(x))
        # print(self.df['dedup_cnt'])

        self.df_bam = self.pcrartefact.denote(
            df_bam=self.df_bam,
            df_condition=self.df,
            condition_set=kwargs['granul_lvl_list'],
            umi_col=kwargs['umi_col'],
            new_col='UMI_mapped',
            pd_col='PD',
            inplace=False,
        )
        # self.df_bam.to_csv(self.work_dir + 'self.df_bam.txt', sep='\t', index=False, header=True)
        # print(self.df_bam)
        # gp_df = self.df_bam.groupby(by=['spikeUMI'])
        # keys = gp_df.groups.keys()
        # for key in keys:
        #     df_sub = gp_df.get_group(key)
        #     print(key, df_sub['UMI_mapped'].unique().shape[0])

        self.console.print('======>calculate average edit distances between umis...')
        self.df['ave_ed'] = self.df.apply(lambda x: self.umigadgetry.ed_ave(x, by_col='repr_nodes'), axis=1)
        self.df['num_diff_dedup_uniq_umis'] = self.df.apply(
            lambda x: self.umigadgetry.num_removed_uniq_umis(
                x,
                by_col='repr_nodes',
            ),
            axis=1,
        )
        self.df['num_diff_dedup_reads'] = self.df.apply(
            lambda x: self.umigadgetry.num_removed_reads(
                x,
                by_col='repr_nodes',
            ),
            axis=1,
        )
        self.console.print('======># of deduplicated unique umis {}'.format(self.df['num_diff_dedup_uniq_umis'].sum()))
        self.console.print('======># of deduplicated reads {}'.format(self.df['num_diff_dedup_reads'].sum()))
        self.ave_ed_bins = self.df['ave_ed'].value_counts().sort_index().to_frame().reset_index()
        self.console.check("======>bins for average edit distance\n{}".format(self.ave_ed_bins))
        if not self.heterogeneity:

            from umiche.util.Folder import Folder as crtfolder
            crtfolder().osmkdir(DIRECTORY=self.work_dir + '/set_cover/')

            self.fwriter.generic(
                df=self.ave_ed_bins,
                sv_fpn=self.work_dir + '/set_cover/' + self.bam_fn + '_' + str(kwargs['token']) + '_ave_ed_bin.txt',
                index=True,
                header=True,
            )
            # print(self.df)
            # print(self.df.columns)
            self.fwriter.generic(
                df=self.df[[
                    'dedup_cnt',
                    'ave_ed',
                    'num_uniq_umis',
                    'num_diff_dedup_uniq_umis',
                    'num_diff_dedup_reads',
                ]],
                sv_fpn=self.work_dir + '/set_cover/' + self.bam_fn + '_' + str(kwargs['token']) + '_dedup_sum.txt',
                index=True,
                header=True,
            )
            self.console.print('======>start writing deduplicated reads to BAM...')
            dedup_reads_write_stime = time.time()
            self.df['sv_bam_ids'] = self.df.apply(lambda x: self.umigadgetry.bamids(x, by_col='repr_nodes'), axis=1)
            self.aliwriter.tobam(
                tobam_fpn=self.work_dir + '/set_cover/' + self.bam_fn + '_' + str(kwargs['token']) + '_dedup.bam',
                tmpl_bam_fpn=self.bam_fpn,
                whitelist=self.umigadgetry.decompose(list_nd=self.df['sv_bam_ids'].values),
            )
            self.console.print('======>finish writing in {:.2f}s'.format(time.time() - dedup_reads_write_stime))
        return {
            'df': self.df,
            'df_bam': self.df_bam,
        }

    def majority_vote_deprecated(
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

        from umiche.util.Folder import Folder as crtfolder
        crtfolder().osmkdir(DIRECTORY=os.path.dirname(kwargs['sv_interm_bam_fpn']))

        self.aliwriter.tobam(
            tobam_fpn=kwargs['sv_interm_bam_fpn'],
            tmpl_bam_fpn=self.bam_fpn,
            whitelist=sc_bam_ids,
        )
        self.console.print('======>finish writing in {:.2f}s'.format(time.time() - dedup_reads_write_stime))
        return self.df

    @Console.vignette()
    def majority_vote(
            self,
            **kwargs,
    ):
        self.df['clusters'] = self.df.apply(
            lambda x: umimv(verbose=self.verbose).track(
                multimer_list=x['vignette']['umi'],
                recur_len=kwargs['umi_unit_pattern'],
            )['clusters'],
            axis=1,
        )
        print(self.df['clusters'])
        self.df['repr_node_map'] = self.df.apply(lambda x: self.umigadgetry.umimax(x, by_col='clusters', met='map'), axis=1)
        print(self.df['repr_node_map'])
        self.df['repr_nodes'] = self.df.apply(lambda x: self.umigadgetry.umimax(x, by_col='clusters'), axis=1)
        print(self.df['repr_nodes'])
        self.df['dedup_cnt'] = self.df['repr_nodes'].apply(lambda x: self.umigadgetry.length(x))
        print(self.df['dedup_cnt'])

        self.df_bam = self.pcrartefact.denote(
            df_bam=self.df_bam,
            df_condition=self.df,
            condition_set=kwargs['granul_lvl_list'],
            umi_col=kwargs['umi_col'],
            new_col='UMI_mapped',
            pd_col='PD',
            inplace=False
        )
        # self.df_bam.to_csv(self.work_dir + 'self.df_bam.txt', sep='\t', index=False, header=True)
        # print(self.df_bam)
        # gp_df = self.df_bam.groupby(by=['spikeUMI'])
        # keys = gp_df.groups.keys()
        # for key in keys:
        #     df_sub = gp_df.get_group(key)
        #     print(key, df_sub['UMI_mapped'].unique().shape[0])

        self.console.print('======>calculate average edit distances between umis...')
        self.df['ave_ed'] = self.df.apply(lambda x: self.umigadgetry.ed_ave(x, by_col='repr_nodes'), axis=1)
        self.df['num_diff_dedup_uniq_umis'] = self.df.apply(
            lambda x: self.umigadgetry.num_removed_uniq_umis(
                x,
                by_col='repr_nodes',
            ),
            axis=1,
        )
        self.df['num_diff_dedup_reads'] = self.df.apply(
            lambda x: self.umigadgetry.num_removed_reads(
                x,
                by_col='repr_nodes',
            ),
            axis=1,
        )
        self.console.print('======># of deduplicated unique umis {}'.format(self.df['num_diff_dedup_uniq_umis'].sum()))
        self.console.print('======># of deduplicated reads {}'.format(self.df['num_diff_dedup_reads'].sum()))
        self.ave_ed_bins = self.df['ave_ed'].value_counts().sort_index().to_frame().reset_index()
        self.console.check("======>bins for average edit distance\n{}".format(self.ave_ed_bins))
        if not self.heterogeneity:
            from umiche.util.Folder import Folder as crtfolder
            crtfolder().osmkdir(DIRECTORY=self.work_dir + '/majority_vote/')

            self.fwriter.generic(
                df=self.ave_ed_bins,
                sv_fpn=self.work_dir + '/majority_vote/' + self.bam_fn + '_' + str(kwargs['token']) + '_ave_ed_bin.txt',
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
                sv_fpn=self.work_dir + '/majority_vote/' + self.bam_fn + '_' + str(kwargs['token']) + '_dedup_sum.txt',
                index=True,
                header=True,
            )
            self.console.print('======>start writing deduplicated reads to BAM...')
            dedup_reads_write_stime = time.time()
            self.df['sv_bam_ids'] = self.df.apply(lambda x: self.umigadgetry.bamids(x, by_col='repr_nodes'), axis=1)
            self.aliwriter.tobam(
                tobam_fpn=self.work_dir + '/majority_vote/' + self.bam_fn + '_' + str(kwargs['token']) + '_dedup.bam',
                tmpl_bam_fpn=self.bam_fpn,
                whitelist=self.umigadgetry.decompose(list_nd=self.df['sv_bam_ids'].values),
            )
            self.console.print('======>finish writing in {:.2f}s'.format(time.time() - dedup_reads_write_stime))
        return {
            'df': self.df,
            'df_bam': self.df_bam,
        }

    @Console.vignette()
    def unique(
            self,
            **kwargs,
    ):
        dedup_umi_stime = time.time()
        self.df['repr_nodes'] = self.df['uniq_repr_nodes']
        self.df['clusters'] = self.df['uniq_repr_nodes'].apply(lambda x: {i: [e] for i, e in enumerate(x)})
        print(self.df['clusters'])
        self.df['repr_node_map'] = self.df.apply(lambda x: self.umigadgetry.umimax(x, by_col='clusters', met='map'), axis=1)

        self.df['uniq_sgl_mark'] = self.df['repr_nodes'].apply(lambda x: 'yes' if len(x) == 1 else 'no')
        self.df = self.df.loc[self.df['uniq_sgl_mark'] == 'no']
        self.console.print('======># of positions with non-single umis: {}'.format(self.df.shape[0]))
        self.console.print('======>finish finding deduplicated UMIs in {:.2f}s'.format(time.time() - dedup_umi_stime))
        # self.console.print('======># of UMIs deduplicated {}'.format(self.df['num_uniq_umis'].loc['yes']))


        self.df_bam = self.pcrartefact.denote(
            df_bam=self.df_bam,
            df_condition=self.df,
            condition_set=kwargs['granul_lvl_list'],
            umi_col=kwargs['umi_col'],
            new_col='UMI_mapped',
            pd_col='PD',
            inplace=False,
        )
        # self.df_bam.to_csv(self.work_dir + 'self.df_bam.txt', sep='\t', index=False, header=True)
        # print(self.df_bam)
        # gp_df = self.df_bam.groupby(by=['spikeUMI'])
        # keys = gp_df.groups.keys()
        # for key in keys:
        #     df_sub = gp_df.get_group(key)
        # print(key, df_sub['UMI_mapped'].unique().shape[0])


        self.console.print('======>calculate average edit distances between umis...')
        ### @@ self.df['ave_ed']
        # 1    5.0
        # Name: ave_eds, dtype: float64
        self.df['ave_ed'] = self.df.apply(lambda x: self.umigadgetry.ed_ave(x, by_col='repr_nodes'), axis=1)
        self.ave_ed_bins = self.df['ave_ed'].value_counts().sort_index().to_frame().reset_index()
        self.console.check("======>bins for average edit distance\n{}".format(self.ave_ed_bins))
        self.df['num_diff_dedup_uniq_umis'] = self.df.apply(
            lambda x: self.umigadgetry.num_removed_uniq_umis(x, by_col='repr_nodes'), axis=1)
        self.df['num_diff_dedup_reads'] = self.df.apply(
            lambda x: self.umigadgetry.num_removed_reads(x, by_col='repr_nodes'),
            axis=1)
        self.console.print('======># of deduplicated unique umis {}'.format(self.df['num_diff_dedup_uniq_umis'].sum()))
        self.console.print('======># of deduplicated reads {}'.format(self.df['num_diff_dedup_reads'].sum()))
        if not self.heterogeneity:
            from umiche.util.Folder import Folder as crtfolder
            crtfolder().osmkdir(DIRECTORY=self.work_dir + '/unique/')

            self.fwriter.generic(
                df=self.ave_ed_bins,
                sv_fpn=self.work_dir + '/unique/' + self.bam_fn + '_' + str(kwargs['token']) + '_ave_ed_bin.txt',
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
                sv_fpn=self.work_dir + '/unique/' + self.bam_fn + '_' + str(kwargs['token']) + '_dedup_sum.txt',
                index=True,
                header=True,
            )
            self.console.print('======>start writing deduplicated reads to BAM...')
            dedup_reads_write_stime = time.time()
            self.df['uniq_bam_ids'] = self.df.apply(lambda x: self.umigadgetry.bamids(x, by_col='repr_nodes'), axis=1)
            self.aliwriter.tobam(
                tobam_fpn=self.work_dir + '/unique/' + self.bam_fn + '_' + str(kwargs['token']) + '_dedup.bam',
                tmpl_bam_fpn=self.bam_fpn,
                whitelist=self.umigadgetry.decompose(list_nd=self.df['uniq_bam_ids'].values),
            )
            self.console.print('======>finish writing in {:.2f}s'.format(time.time() - dedup_reads_write_stime))
        return {
            'df': self.df,
            'df_bam': self.df_bam,
        }

    @Console.vignette()
    def cluster(
            self,
            **kwargs,
    ):
        dedup_umi_stime = time.time()
        self.df['clusters'] = self.df['cc']
        self.df['repr_node_map'] = self.df.apply(lambda x: self.umigadgetry.umimax(x, by_col='clusters', met='map'), axis=1)
        self.df['repr_nodes'] = self.df.apply(lambda x: self.umigadgetry.umimax(x, by_col='cc'), axis=1)
        ### @@ self.df['repr_nodes']
        # 1    [2]
        # Name: repr_nodes, dtype: object
        self.df['dedup_cnt'] = self.df['repr_nodes'].apply(lambda x: self.umigadgetry.length(x))
        ### @@ self.df['dedup_cnt']
        # 1    1
        # Name: dedup_cnt, dtype: int64
        self.console.print('======>finish finding deduplicated umis in {:.2f}s'.format(time.time() - dedup_umi_stime))
        # self.console.print('======># of umis deduplicated to be {}'.format(self.df['dedup_cnt'].loc['yes']))


        self.df_bam = self.pcrartefact.denote(
            df_bam=self.df_bam,
            df_condition=self.df,
            condition_set=kwargs['granul_lvl_list'],
            umi_col=kwargs['umi_col'],
            new_col='UMI_mapped',
            pd_col='PD',
            inplace=False
        )
        # self.df_bam.to_csv(self.work_dir + 'self.df_bam.txt', sep='\t', index=False, header=True)
        # print(self.df_bam)
        # gp_df = self.df_bam.groupby(by=['spikeUMI'])
        # keys = gp_df.groups.keys()
        # for key in keys:
        #     df_sub = gp_df.get_group(key)
        # print(key, df_sub['UMI_mapped'].unique().shape[0])


        self.console.print('======>calculate average edit distances between umis...')
        self.df['ave_ed'] = self.df.apply(lambda x: self.umigadgetry.ed_ave(x, by_col='repr_nodes'), axis=1)
        self.df['num_diff_dedup_uniq_umis'] = self.df.apply(
            lambda x: self.umigadgetry.num_removed_uniq_umis(x, by_col='repr_nodes'),
            axis=1)
        self.df['num_diff_dedup_reads'] = self.df.apply(
            lambda x: self.umigadgetry.num_removed_reads(x, by_col='repr_nodes'),
            axis=1)
        self.console.print('======># of deduplicated unique umis {}'.format(self.df['num_diff_dedup_uniq_umis'].sum()))
        self.console.print('======># of deduplicated reads {}'.format(self.df['num_diff_dedup_reads'].sum()))
        self.ave_ed_bins = self.df['ave_ed'].value_counts().sort_index().to_frame().reset_index()
        self.console.check("======>bins for average edit distance\n{}".format(self.ave_ed_bins))
        if not self.heterogeneity:
            from umiche.util.Folder import Folder as crtfolder
            crtfolder().osmkdir(DIRECTORY=self.work_dir + '/cluster/')

            self.fwriter.generic(
                df=self.ave_ed_bins,
                sv_fpn=self.work_dir + '/cluster/' + self.bam_fn + '_' + str(kwargs['token']) + '_ave_ed_bin.txt',
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
                sv_fpn=self.work_dir + '/cluster/' + self.bam_fn + '_' + str(kwargs['token']) + '_dedup_sum.txt',
                index=True,
                header=True,
            )
            self.console.print('======>start writing deduplicated reads to BAM...')
            dedup_reads_write_stime = time.time()
            self.df['cc_bam_ids'] = self.df.apply(lambda x: self.umigadgetry.bamids(x, by_col='repr_nodes'), axis=1)
            self.aliwriter.tobam(
                tobam_fpn=self.work_dir + '/cluster/' + self.bam_fn + '_' + str(kwargs['token']) + '_dedup.bam',
                tmpl_bam_fpn=self.bam_fpn,
                whitelist=self.umigadgetry.decompose(list_nd=self.df['cc_bam_ids'].values),
            )
            self.console.print('======>finish writing in {:.2f}s'.format(time.time() - dedup_reads_write_stime))
        return {
            'df': self.df,
            'df_bam': self.df_bam,
        }

    @Console.vignette()
    def adjacency(
            self,
            **kwargs,
    ):
        # print('adj kwargs:', kwargs)
        dedup_umi_stime = time.time()
        umiadj_ob = umiadj(verbose=False)
        self.df['clusters'] = self.df.apply(
            lambda x: umiadj_ob.decompose(
                cc_sub_dict=umiadj_ob.umi_tools(
                    connected_components=x['cc'],
                    df_umi_uniq_val_cnt=x['vignette']['df_umi_uniq_val_cnt'],
                    graph_adj=x['vignette']['graph_adj'],
                )['clusters'],
            ),
            axis=1,
        )
        # print(self.df['clusters'])
        self.df['repr_node_map'] = self.df.apply(lambda x: self.umigadgetry.umimax(x, by_col='clusters', met='map'), axis=1)
        self.df['repr_nodes'] = self.df.apply(lambda x: self.umigadgetry.umimax(x, by_col='clusters'), axis=1)
        self.df['dedup_cnt'] = self.df['repr_nodes'].apply(lambda x: self.umigadgetry.length(x))
        # print(self.df['dedup_cnt'])
        # print(self.df)
        # print(self.df.columns)

        self.df_bam = self.pcrartefact.denote(
            df_bam=self.df_bam,
            df_condition=self.df,
            condition_set=kwargs['granul_lvl_list'],
            umi_col=kwargs['umi_col'],
            new_col='UMI_mapped',
            pd_col='PD',
            inplace=False
        )
        # self.df_bam.to_csv(self.work_dir + 'self.df_bam.txt', sep='\t', index=False, header=True)
        # print(self.df_bam)
#         gp_df = self.df_bam.groupby(by=['spikeUMI'])
#         keys = gp_df.groups.keys()
#         for key in keys:
#             df_sub = gp_df.get_group(key)
#             print(key, df_sub['UMI_mapped'].unique().shape[0])

        self.console.print('======>finish finding deduplicated umis in {:.2f}s'.format(time.time() - dedup_umi_stime))
        # self.console.print('======># of umis deduplicated to be {}'.format(self.df['dedup_cnt'].loc['yes']))
        self.console.print('======>calculate average edit distances between umis...')
        self.df['ave_ed'] = self.df.apply(lambda x: self.umigadgetry.ed_ave(x, by_col='repr_nodes'), axis=1)
        self.df['num_diff_dedup_uniq_umis'] = self.df.apply(
            lambda x: self.umigadgetry.num_removed_uniq_umis(
                x,
                by_col='repr_nodes',
            ),
            axis=1,
        )
        self.df['num_diff_dedup_reads'] = self.df.apply(
            lambda x: self.umigadgetry.num_removed_reads(
                x,
                by_col='repr_nodes',
            ),
            axis=1,
        )
        self.console.print('======># of deduplicated unique umis {}'.format(self.df['num_diff_dedup_uniq_umis'].sum()))
        self.console.print('======># of deduplicated reads {}'.format(self.df['num_diff_dedup_reads'].sum()))
        self.ave_ed_bins = self.df['ave_ed'].value_counts().sort_index().to_frame().reset_index()
        self.console.check("======>bins for average edit distance\n{}".format(self.ave_ed_bins))
        if not self.heterogeneity:
            from umiche.util.Folder import Folder as crtfolder
            crtfolder().osmkdir(DIRECTORY=self.work_dir + '/adjacency/')

            self.fwriter.generic(
                df=self.ave_ed_bins,
                sv_fpn=self.work_dir + '/adjacency/' + self.bam_fn + '_' + str(kwargs['token']) + '_ave_ed_bin.txt',
                index=True,
                header=True,
            )
            # print(self.df)
            # print(self.df.columns)
            self.fwriter.generic(
                df=self.df[[
                    'dedup_cnt',
                    'ave_ed',
                    'num_uniq_umis',
                    'num_diff_dedup_uniq_umis',
                    'num_diff_dedup_reads',
                ]],
                sv_fpn=self.work_dir + '/adjacency/' + self.bam_fn + '_' + str(kwargs['token']) + '_dedup_sum.txt',
                index=True,
                header=True,
            )
            self.console.print('======>start writing deduplicated reads to BAM...')
            dedup_reads_write_stime = time.time()
            self.df['adj_bam_ids'] = self.df.apply(lambda x: self.umigadgetry.bamids(x, by_col='repr_nodes'), axis=1)
            self.aliwriter.tobam(
                tobam_fpn=self.work_dir + '/adjacency/' + self.bam_fn + '_' + str(kwargs['token']) + '_dedup.bam',
                tmpl_bam_fpn=self.bam_fpn,
                whitelist=self.umigadgetry.decompose(list_nd=self.df['adj_bam_ids'].values),
            )
            self.console.print('======>finish writing in {:.2f}s'.format(time.time() - dedup_reads_write_stime))
        return {
            'df': self.df,
            'df_bam': self.df_bam,
        }

    @Console.vignette()
    def directional(
            self,
            **kwargs,
    ):
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
            self.df['clusters'] = self.df.apply(
                lambda x: umidirec_ob.decompose(
                    cc_sub_dict=umidirec_ob.umi_tools(
                        connected_components=x['cc'],
                        df_umi_uniq_val_cnt=x['vignette']['df_umi_uniq_val_cnt'],
                        graph_adj=x['vignette']['graph_adj'],
                    )['clusters'],
                ),
                axis=1,
            )
        print(self.df['clusters'])
        self.df['repr_node_map'] = self.df.apply(lambda x: self.umigadgetry.umimax(x, by_col='clusters', met='map'), axis=1)
        self.df['repr_nodes'] = self.df.apply(lambda x: self.umigadgetry.umimax(x, by_col='clusters'), axis=1)
        # print(self.df['repr_nodes'])
        self.df['dedup_cnt'] = self.df['repr_nodes'].apply(lambda x: self.umigadgetry.length(x))
        print(self.df['dedup_cnt'])
        self.console.print('======>finish finding deduplicated umis in {:.2f}s'.format(time.time() - dedup_umi_stime))
        # self.console.print('======># of umis deduplicated to be {}'.format(self.df['dedup_cnt'].loc['yes']))

        self.df_bam = self.pcrartefact.denote(
            df_bam=self.df_bam,
            df_condition=self.df,
            condition_set=kwargs['granul_lvl_list'],
            umi_col=kwargs['umi_col'],
            new_col='UMI_mapped',
            pd_col='PD',
            inplace=False
        )
        # self.df_bam.to_csv(self.work_dir + 'self.df_bam.txt', sep='\t', index=False, header=True)
        # print(self.df_bam)
        # gp_df = self.df_bam.groupby(by=['spikeUMI'])
        # keys = gp_df.groups.keys()
        # for key in keys:
        #     df_sub = gp_df.get_group(key)
        # print(key, df_sub['UMI_mapped'].unique().shape[0])

        self.console.print('======>calculate average edit distances between umis...')
        self.df['ave_ed'] = self.df.apply(lambda x: self.umigadgetry.ed_ave(x, by_col='repr_nodes'), axis=1)
        # print(self.df['ave_ed'])
        self.df['num_diff_dedup_uniq_umis'] = self.df.apply(
            lambda x: self.umigadgetry.num_removed_uniq_umis(x, by_col='repr_nodes'),
            axis=1,
        )
        self.df['num_diff_dedup_reads'] = self.df.apply(
            lambda x: self.umigadgetry.num_removed_reads(x, by_col='repr_nodes'),
            axis=1,
        )
        self.console.print('======># of deduplicated unique umis {}'.format(self.df['num_diff_dedup_uniq_umis'].sum()))
        self.console.print('======># of deduplicated reads {}'.format(self.df['num_diff_dedup_reads'].sum()))
        self.ave_ed_bins = self.df['ave_ed'].value_counts().sort_index().to_frame().reset_index()
        # print(self.ave_ed_bins)
        self.console.check("======>bins for average edit distance\n{}".format(self.ave_ed_bins))
        if not self.heterogeneity:
            from umiche.util.Folder import Folder as crtfolder
            crtfolder().osmkdir(DIRECTORY=self.work_dir + '/directional/')

            self.fwriter.generic(
                df=self.ave_ed_bins,
                sv_fpn=self.work_dir + '/directional/' + self.bam_fn + '_' + str(kwargs['token']) + '_ave_ed_bin.txt',
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
                sv_fpn=self.work_dir + '/directional/' + self.bam_fn + '_' + str(kwargs['token']) + '_dedup_sum.txt',
                index=True,
                header=True,
            )
            self.console.print('======>start writing deduplicated reads to BAM...')
            dedup_reads_write_stime = time.time()
            self.df['direc_bam_ids'] = self.df.apply(lambda x: self.umigadgetry.bamids(x, by_col='repr_nodes'), axis=1)
            self.aliwriter.tobam(
                tobam_fpn=self.work_dir + '/directional/' + self.bam_fn + '_' + str(kwargs['token']) + '_dedup.bam',
                tmpl_bam_fpn=self.bam_fpn,
                whitelist=self.umigadgetry.decompose(list_nd=self.df['direc_bam_ids'].values),
            )
            self.console.print('======>finish writing in {:.2f}s'.format(time.time() - dedup_reads_write_stime))
        return {
            'df': self.df,
            'df_bam': self.df_bam,
        }

    @Console.vignette()
    def umicountr(
            self,
            clustering_method,
            **kwargs,
    ):
        print('clustering method: {}'.format(clustering_method))
        dedup_umi_stime = time.time()
        umiadj_ob = umiadj()
        if clustering_method == 'adj':
            dedup_func = umiadj_ob.umicountr
        elif clustering_method == 'adj_direc':
            dedup_func = umiadj_ob.umicountr_directional
        elif clustering_method == 'adj_singleton':
            dedup_func = umiadj_ob.umicountr_singleton
        else:
            dedup_func = umiadj_ob.umicountr

        self.df['clusters'] = self.df.apply(
            lambda x: umiadj_ob.decompose(
                cc_sub_dict=dedup_func(
                    connected_components=x['cc'],
                    df_umi_uniq_val_cnt=x['vignette']['df_umi_uniq_val_cnt'],
                    graph_adj=x['vignette']['graph_adj'],
                )['clusters'],
            ),
            axis=1,
        )
        print(self.df['clusters'])
        # self.df['assigned'] = self.df.apply(
        #     lambda x: umiadj_ob.umicountr(
        #             connected_components=x['cc'],
        #             df_umi_uniq_val_cnt=x['vignette']['df_umi_uniq_val_cnt'],
        #             graph_adj=x['vignette']['graph_adj'],
        #         )['assigned'],
        #     axis=1,
        # )
        # print(self.df['assigned'])
        self.df['repr_node_map'] = self.df.apply(lambda x: self.umigadgetry.umimax(x, by_col='clusters', met='map'), axis=1)
        self.df['repr_nodes'] = self.df.apply(lambda x: self.umigadgetry.umimax(x, by_col='clusters'), axis=1)
        self.df['dedup_cnt'] = self.df['repr_nodes'].apply(lambda x: self.umigadgetry.length(x))
        print(self.df['dedup_cnt'])
        self.console.print('======>finish finding deduplicated umis in {:.2f}s'.format(time.time() - dedup_umi_stime))
        # self.console.print('======># of umis deduplicated to be {}'.format(self.df['dedup_cnt'].loc['yes']))

        self.df_bam = self.pcrartefact.denote(
            df_bam=self.df_bam,
            df_condition=self.df,
            condition_set=kwargs['granul_lvl_list'],
            umi_col=kwargs['umi_col'],
            new_col='UMI_mapped',
            pd_col='PD',
            inplace=False
        )
        # self.df_bam.to_csv(self.work_dir + 'self.df_bam.txt', sep='\t', index=False, header=True)
        # print(self.df_bam)
        # gp_df = self.df_bam.groupby(by=['spikeUMI'])
        # keys = gp_df.groups.keys()
        # for key in keys:
        #     df_sub = gp_df.get_group(key)
        # print(key, df_sub['UMI_mapped'].unique().shape[0])


        self.console.print('======>calculate average edit distances between umis...')
        self.df['ave_ed'] = self.df.apply(lambda x: self.umigadgetry.ed_ave(x, by_col='repr_nodes'), axis=1)
        self.df['num_diff_dedup_uniq_umis'] = self.df.apply(
            lambda x: self.umigadgetry.num_removed_uniq_umis(
                x,
                by_col='repr_nodes',
            ),
            axis=1,
        )
        self.df['num_diff_dedup_reads'] = self.df.apply(
            lambda x: self.umigadgetry.num_removed_reads(
                x,
                by_col='repr_nodes',
            ),
            axis=1,
        )
        self.console.print('======># of deduplicated unique umis {}'.format(self.df['num_diff_dedup_uniq_umis'].sum()))
        self.console.print('======># of deduplicated reads {}'.format(self.df['num_diff_dedup_reads'].sum()))
        self.ave_ed_bins = self.df['ave_ed'].value_counts().sort_index().to_frame().reset_index()
        self.console.check("======>bins for average edit distance\n{}".format(self.ave_ed_bins))
        if not self.heterogeneity:
            from umiche.util.Folder import Folder as crtfolder
            crtfolder().osmkdir(DIRECTORY=self.work_dir + '/' + str(clustering_method) + '/')

            self.fwriter.generic(
                df=self.ave_ed_bins,
                sv_fpn=self.work_dir + '/' + str(clustering_method) + '/' + self.bam_fn + '_' + str(kwargs['token']) + '_ave_ed_bin.txt',
                index=True,
                header=True,
            )
            # print(self.df)
            # print(self.df.columns)
            self.fwriter.generic(
                df=self.df[[
                    'dedup_cnt',
                    'ave_ed',
                    'num_uniq_umis',
                    'num_diff_dedup_uniq_umis',
                    'num_diff_dedup_reads',
                ]],
                sv_fpn=self.work_dir + '/' + str(clustering_method) + '/' + self.bam_fn + '_' + str(kwargs['token']) + '_dedup_sum.txt',
                index=True,
                header=True,
            )
            self.console.print('======>start writing deduplicated reads to BAM...')
            dedup_reads_write_stime = time.time()
            self.df['adj_bam_ids'] = self.df.apply(lambda x: self.umigadgetry.bamids(x, by_col='repr_nodes'), axis=1)
            self.aliwriter.tobam(
                tobam_fpn=self.work_dir + '/' + str(clustering_method) + '/' + self.bam_fn + '_' + str(kwargs['token']) + '_dedup.bam',
                tmpl_bam_fpn=self.bam_fpn,
                whitelist=self.umigadgetry.decompose(list_nd=self.df['adj_bam_ids'].values),
            )
            self.console.print('======>finish writing in {:.2f}s'.format(time.time() - dedup_reads_write_stime))
        return {
            'df': self.df,
            'df_bam': self.df_bam,
        }

    @Console.vignette()
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
            self.df['clusters'] = self.df.apply(
                lambda x: umimcl_ob.decompose(
                    list_nd=umimcl_ob.dfclusters(
                        connected_components=x['cc'],
                        graph_adj=x['vignette']['graph_adj'],
                    )['clusters'].values,
                ),
                axis=1,
            )
            # print(self.df['clusters'])
        ### @@ self.df['clusters']
        # 1    {0: [0, 76, 162, 188, 237, 256], 1: [65, 55, 1...
        # Name: mcl, dtype: object
        self.df['repr_node_map'] = self.df.apply(lambda x: self.umigadgetry.umimax(x, by_col='clusters', met='map'), axis=1)
        self.df['repr_nodes'] = self.df.apply(lambda x: self.umigadgetry.umimax(x, by_col='clusters'), axis=1)
        self.df['dedup_cnt'] = self.df['repr_nodes'].apply(lambda x: self.umigadgetry.length(x))
        self.console.print('======>finish finding deduplicated umis in {:.2f}s'.format(time.time() - dedup_umi_stime))
        # self.console.print('======># of umis deduplicated to be {}'.format(self.df['dedup_cnt'].loc['yes']))

        self.df_bam = self.pcrartefact.denote(
            df_bam=self.df_bam,
            df_condition=self.df,
            condition_set=kwargs['granul_lvl_list'],
            umi_col=kwargs['umi_col'],
            new_col='UMI_mapped',
            pd_col='PD',
            inplace=False
        )
        # self.df_bam.to_csv(self.work_dir + 'self.df_bam.txt', sep='\t', index=False, header=True)
        # print(self.df_bam)
        # gp_df = self.df_bam.groupby(by=['spikeUMI'])
        # keys = gp_df.groups.keys()
        # for key in keys:
        #     df_sub = gp_df.get_group(key)
        # print(key, df_sub['UMI_mapped'].unique().shape[0])


        self.console.print('======>calculate average edit distances between umis...')
        self.df['ave_ed'] = self.df.apply(lambda x: self.umigadgetry.ed_ave(x, by_col='repr_nodes'), axis=1)
        self.df['num_diff_dedup_uniq_umis'] = self.df.apply(lambda x: self.umigadgetry.num_removed_uniq_umis(x, by_col='repr_nodes'),axis=1)
        self.df['num_diff_dedup_reads'] = self.df.apply(lambda x: self.umigadgetry.num_removed_reads(x, by_col='repr_nodes'),axis=1)
        self.console.print('======># of deduplicated unique umis {}'.format(self.df['num_diff_dedup_uniq_umis'].sum()))
        self.console.print('======># of deduplicated reads {}'.format(self.df['num_diff_dedup_reads'].sum()))
        self.ave_ed_bins = self.df['ave_ed'].value_counts().sort_index().to_frame().reset_index()
        self.console.check("======>bins for average edit distance\n{}".format(self.ave_ed_bins))
        if not self.heterogeneity:
            from umiche.util.Folder import Folder as crtfolder
            crtfolder().osmkdir(DIRECTORY=self.work_dir + '/mcl/')

            self.fwriter.generic(
                df=self.ave_ed_bins,
                sv_fpn=self.work_dir + '/mcl/' + self.bam_fn + '_' + str(kwargs['token']) + '_ave_ed_bin.txt',
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
                sv_fpn=self.work_dir + '/mcl/' + self.bam_fn + '_' + str(kwargs['token']) + '_dedup_sum.txt',
                index=True,
                header=True,
            )
            self.console.print('======>start writing deduplicated reads to BAM...')
            dedup_reads_write_stime = time.time()
            self.df['mcl_bam_ids'] = self.df.apply(lambda x: self.umigadgetry.bamids(x, by_col='repr_nodes'), axis=1)
            self.aliwriter.tobam(
                tobam_fpn=self.work_dir + '/mcl/' + self.bam_fn + '_' + str(kwargs['token']) + '_dedup.bam',
                tmpl_bam_fpn=self.bam_fpn,
                whitelist=self.umigadgetry.decompose(list_nd=self.df['mcl_bam_ids'].values),
            )
            self.console.print('======>finish writing in {:.2f}s'.format(time.time() - dedup_reads_write_stime))
        return {
            'df': self.df,
            'df_bam': self.df_bam,
        }

    @Console.vignette()
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
        self.df['ave_ed'] = self.df.apply(lambda x: self.umigadgetry.ed_ave(x, by_col='mcl_repr_nodes'), axis=1)
        self.df['num_diff_dedup_uniq_umis'] = self.df.apply(lambda x: self.umigadgetry.num_removed_uniq_umis(x, by_col='mcl_repr_nodes'),axis=1)
        self.df['num_diff_dedup_reads'] = self.df.apply(lambda x: self.umigadgetry.num_removed_reads(x, by_col='mcl_repr_nodes'),axis=1)
        self.console.print('======># of deduplicated unique umis {}'.format(self.df['num_diff_dedup_uniq_umis'].sum()))
        self.console.print('======># of deduplicated reads {}'.format(self.df['num_diff_dedup_reads'].sum()))
        self.ave_ed_bins = self.df['ave_ed'].value_counts().sort_index().to_frame().reset_index()
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
        return {
            'df': self.df,
            'df_bam': self.df_bam,
        }

    @Console.vignette()
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
            self.df['clusters'] = self.df.apply(
                lambda x: umimcl_ob.decompose(
                    list_nd=x['clusters'].values,
                ),
                axis=1,
            )
            # print(self.df['clusters'])
        else:
            self.df['clusters'] = self.df.apply(
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
            # print(self.df['clusters'])

        self.df['repr_node_map'] = self.df.apply(lambda x: self.umigadgetry.umimax(x, by_col='clusters', met='map'), axis=1)
        self.df['repr_nodes'] = self.df.apply(lambda x: self.umigadgetry.umimax(x, by_col='clusters'), axis=1)
        self.df['dedup_cnt'] = self.df['repr_nodes'].apply(lambda x: self.umigadgetry.length(x))
        self.console.print('======>finish finding deduplicated umis in {:.2f}s'.format(time.time() - dedup_umi_stime))
        # self.console.print('======># of umis deduplicated to be {}'.format(self.df['dedup_cnt'].loc['yes']))


        self.df_bam = self.pcrartefact.denote(
            df_bam=self.df_bam,
            df_condition=self.df,
            condition_set=kwargs['granul_lvl_list'],
            umi_col=kwargs['umi_col'],
            new_col='UMI_mapped',
            pd_col='PD',
            inplace=False
        )
        # self.df_bam.to_csv(self.work_dir + 'self.df_bam.txt', sep='\t', index=False, header=True)
        # print(self.df_bam)
        # gp_df = self.df_bam.groupby(by=['spikeUMI'])
        # keys = gp_df.groups.keys()
        # for key in keys:
        #     df_sub = gp_df.get_group(key)
        # print(key, df_sub['UMI_mapped'].unique().shape[0])


        self.console.print('======>calculate average edit distances between umis...')
        self.df['ave_ed'] = self.df.apply(lambda x: self.umigadgetry.ed_ave(x, by_col='repr_nodes'), axis=1)
        self.df['num_diff_dedup_uniq_umis'] = self.df.apply(
            lambda x: self.umigadgetry.num_removed_uniq_umis(x, by_col='repr_nodes'), axis=1)
        self.df['num_diff_dedup_reads'] = self.df.apply(
            lambda x: self.umigadgetry.num_removed_reads(x, by_col='repr_nodes'), axis=1)
        self.console.print('======># of deduplicated unique umis {}'.format(self.df['num_diff_dedup_uniq_umis'].sum()))
        self.console.print('======># of deduplicated reads {}'.format(self.df['num_diff_dedup_reads'].sum()))
        self.ave_ed_bins = self.df['ave_ed'].value_counts().sort_index().to_frame().reset_index()
        self.console.check("======>bins for average edit distance\n{}".format(self.ave_ed_bins))
        if not self.heterogeneity:
            from umiche.util.Folder import Folder as crtfolder
            crtfolder().osmkdir(DIRECTORY=self.work_dir + '/mcl_val/')

            self.fwriter.generic(
                df=self.ave_ed_bins,
                sv_fpn=self.work_dir + '/mcl_val/' + self.bam_fn + '_' + str(kwargs['token']) + '_ave_ed_bin.txt',
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
                sv_fpn=self.work_dir + '/mcl_val/' + self.bam_fn + '_' + str(kwargs['token']) + '_dedup_sum.txt',
                index=True,
                header=True,
            )
            self.console.print('======>start writing deduplicated reads to BAM...')
            dedup_reads_write_stime = time.time()
            self.df['mcl_val_bam_ids'] = self.df.apply(lambda x: self.umigadgetry.bamids(x, by_col='repr_nodes'),
                                                       axis=1)
            self.aliwriter.tobam(
                tobam_fpn=self.work_dir + '/mcl_val/' + self.bam_fn + '_' + str(kwargs['token']) + '_dedup.bam',
                tmpl_bam_fpn=self.bam_fpn,
                whitelist=self.umigadgetry.decompose(list_nd=self.df['mcl_val_bam_ids'].values),
            )
            self.console.print('======>finish writing in {:.2f}s'.format(time.time() - dedup_reads_write_stime))
        return {
            'df': self.df,
            'df_bam': self.df_bam,
        }

    @Console.vignette()
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
            self.df['clusters'] = self.df.apply(
                lambda x: umimcl_ob.decompose(
                    list_nd=x['clusters'].values,
                    # list_nd=x['clusters'].values[0].tolist(),
                ),
                axis=1,
            )
            # print(self.df['clusters'])
        else:
            self.df['clusters'] = self.df.apply(
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

        self.df['repr_node_map'] = self.df.apply(lambda x: self.umigadgetry.umimax(x, by_col='clusters', met='map'), axis=1)
        self.df['repr_nodes'] = self.df.apply(lambda x: self.umigadgetry.umimax(x, by_col='clusters'), axis=1)
        self.df['dedup_cnt'] = self.df['repr_nodes'].apply(lambda x: self.umigadgetry.length(x))
        # print(self.df['dedup_cnt'])
        self.console.print('======>finish finding deduplicated umis in {:.2f}s'.format(time.time() - dedup_umi_stime))
        # self.console.print('======># of umis deduplicated to be {}'.format(self.df['dedup_cnt'].loc['yes']))


        self.df_bam = self.pcrartefact.denote(
            df_bam=self.df_bam,
            df_condition=self.df,
            condition_set=kwargs['granul_lvl_list'],
            umi_col=kwargs['umi_col'],
            new_col='UMI_mapped',
            pd_col='PD',
            inplace=False
        )
        # self.df_bam.to_csv(self.work_dir + 'self.df_bam.txt', sep='\t', index=False, header=True)
        # print(self.df_bam)
        # gp_df = self.df_bam.groupby(by=['spikeUMI'])
        # keys = gp_df.groups.keys()
        # for key in keys:
        #     df_sub = gp_df.get_group(key)
        # print(key, df_sub['UMI_mapped'].unique().shape[0])


        self.console.print('======>calculate average edit distances between umis...')
        self.df['ave_ed'] = self.df.apply(lambda x: self.umigadgetry.ed_ave(x, by_col='repr_nodes'), axis=1)
        self.df['num_diff_dedup_uniq_umis'] = self.df.apply(
            lambda x: self.umigadgetry.num_removed_uniq_umis(x, by_col='repr_nodes'), axis=1)
        self.df['num_diff_dedup_reads'] = self.df.apply(
            lambda x: self.umigadgetry.num_removed_reads(x, by_col='repr_nodes'),
            axis=1)
        self.console.print('======># of deduplicated unique umis {}'.format(self.df['num_diff_dedup_uniq_umis'].sum()))
        self.console.print('======># of deduplicated reads {}'.format(self.df['num_diff_dedup_reads'].sum()))
        self.ave_ed_bins = self.df['ave_ed'].value_counts().sort_index().to_frame().reset_index()
        self.console.check("======>bins for average edit distance\n{}".format(self.ave_ed_bins))
        if not self.heterogeneity:
            from umiche.util.Folder import Folder as crtfolder
            crtfolder().osmkdir(DIRECTORY=self.work_dir + '/mcl_ed/')

            self.fwriter.generic(
                df=self.ave_ed_bins,
                sv_fpn=self.work_dir + '/mcl_ed/' + self.bam_fn + '_' + str(kwargs['token']) + '_ave_ed_bin.txt',
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
                sv_fpn=self.work_dir + '/mcl_ed/' + self.bam_fn + '_' + str(kwargs['token']) + '_dedup_sum.txt',
                index=True,
                header=True,
            )
            self.console.print('======>start writing deduplicated reads to BAM...')
            dedup_reads_write_stime = time.time()
            self.df['mcl_ed_bam_ids'] = self.df.apply(lambda x: self.umigadgetry.bamids(x, by_col='repr_nodes'), axis=1)
            self.aliwriter.tobam(
                tobam_fpn=self.work_dir + '/mcl_ed/' + self.bam_fn + '_' + str(kwargs['token']) + '_dedup.bam',
                tmpl_bam_fpn=self.bam_fpn,
                whitelist=self.umigadgetry.decompose(list_nd=self.df['mcl_ed_bam_ids'].values),
            )
            self.console.print('======>finish writing in {:.2f}s'.format(time.time() - dedup_reads_write_stime))
        return {
            'df': self.df,
            'df_bam': self.df_bam,
        }

    @Console.vignette()
    def clustering_umi_seq_onehot(
            self,
            clustering_method,
            **kwargs,
    ):
        print('clustering method: {}'.format(clustering_method))
        dedup_umi_stime = time.time()
        self.umiclustering = umiclustering(
            clustering_method=clustering_method,
            **kwargs,
        )
        self.df['clusters'] = self.df.apply(
            lambda x: self.umiclustering.decompose(
                list_nd=self.umiclustering.dfclusters(
                    connected_components=x['cc'],
                    graph_adj=x['vignette']['graph_adj'],
                    int_to_umi_dict=x['vignette']['int_to_umi_dict'],
                )['clusters'].values,
            ),
            axis=1,
        )
        print(self.df['clusters'])
        self.df['repr_node_map'] = self.df.apply(lambda x: self.umigadgetry.umimax(x, by_col='clusters', met='map'), axis=1)
        self.df['repr_nodes'] = self.df.apply(lambda x: self.umigadgetry.umimax(x, by_col='clusters'), axis=1)
        self.df['dedup_cnt'] = self.df['repr_nodes'].apply(lambda x: self.umigadgetry.length(x))
        self.console.print('======>finish finding deduplicated umis in {:.2f}s'.format(time.time() - dedup_umi_stime))
        # self.console.print('======># of umis deduplicated to be {}'.format(self.df['dedup_cnt'].loc['yes']))

        # print(self.df['repr_nodes'])

        self.df_bam = self.pcrartefact.denote(
            df_bam=self.df_bam,
            df_condition=self.df,
            condition_set=kwargs['granul_lvl_list'],
            umi_col=kwargs['umi_col'],
            new_col='UMI_mapped',
            pd_col='PD',
            inplace=False
        )
        # self.df_bam.to_csv(self.work_dir + 'self.df_bam.txt', sep='\t', index=False, header=True)
        # print(self.df_bam)
        # gp_df = self.df_bam.groupby(by=['spikeUMI'])
        # keys = gp_df.groups.keys()
        # for key in keys:
        #     df_sub = gp_df.get_group(key)
        #     print(key, df_sub['UMI_mapped'].unique().shape[0])

        self.console.print('======>calculate average edit distances between umis...')
        self.df['ave_ed'] = self.df.apply(lambda x: self.umigadgetry.ed_ave(x, by_col='repr_nodes'), axis=1)
        self.df['num_diff_dedup_uniq_umis'] = self.df.apply(
            lambda x: self.umigadgetry.num_removed_uniq_umis(x, by_col='repr_nodes'),
            axis=1)
        self.df['num_diff_dedup_reads'] = self.df.apply(
            lambda x: self.umigadgetry.num_removed_reads(x, by_col='repr_nodes'),
            axis=1)
        self.console.print('======># of deduplicated unique umis {}'.format(self.df['num_diff_dedup_uniq_umis'].sum()))
        self.console.print('======># of deduplicated reads {}'.format(self.df['num_diff_dedup_reads'].sum()))
        self.ave_ed_bins = self.df['ave_ed'].value_counts().sort_index().to_frame().reset_index()
        self.console.check("======>bins for average edit distance\n{}".format(self.ave_ed_bins))
        if not self.heterogeneity:
            from umiche.util.Folder import Folder as crtfolder
            crtfolder().osmkdir(DIRECTORY=self.work_dir + '/' + str(clustering_method) + '/')
            self.fwriter.generic(
                df=self.ave_ed_bins,
                sv_fpn=self.work_dir + '/' + str(clustering_method) + '/' + self.bam_fn + '_' + str(kwargs['token']) + '_ave_ed_bin.txt',
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
                sv_fpn=self.work_dir + '/' + str(clustering_method) + '/' + self.bam_fn + '_' + str(kwargs['token']) + '_dedup_sum.txt',
                index=True,
                header=True,
            )
            self.console.print('======>start writing deduplicated reads to BAM...')
            dedup_reads_write_stime = time.time()
            self.df[clustering_method + '_bam_ids'] = self.df.apply(lambda x: self.umigadgetry.bamids(x, by_col='repr_nodes'), axis=1)
            self.aliwriter.tobam(
                tobam_fpn=self.work_dir + '/' + str(clustering_method) + '/' + self.bam_fn + '_' + str(kwargs['token']) + '_dedup.bam',
                tmpl_bam_fpn=self.bam_fpn,
                whitelist=self.umigadgetry.decompose(list_nd=self.df[clustering_method + '_bam_ids'].values),
            )
            self.console.print('======>finish writing in {:.2f}s'.format(time.time() - dedup_reads_write_stime))
        return {
            'df': self.df,
            'df_bam': self.df_bam,
        }

    @Console.vignette()
    def starsolo(
            self,
            **kwargs,
    ):
        # print('adj kwargs:', kwargs)
        from umiche.deduplicate.method.starsolo.STARsolo import STARsolo

        dedup_umi_stime = time.time()
        starsolo_ob = STARsolo(verbose=True)
        self.df['clusters'] = self.df.apply(
            lambda x: starsolo_ob.dedup_from_graph(
                connected_components=x['cc'],
                counts=x['vignette']['df_umi_uniq_val_cnt'].to_dict(),
                graph_adj=x['vignette']['graph_adj'],
            )['clusters'],
            axis=1,
        )
        print(self.df['clusters'])
        self.df['repr_node_map'] = self.df.apply(lambda x: self.umigadgetry.umimax(x, by_col='clusters', met='map'), axis=1)
        self.df['repr_nodes'] = self.df.apply(lambda x: self.umigadgetry.umimax(x, by_col='clusters'), axis=1)
        self.df['dedup_cnt'] = self.df['repr_nodes'].apply(lambda x: self.umigadgetry.length(x))
        # print(self.df['dedup_cnt'])
        # print(self.df)
        # print(self.df.columns)

        self.df_bam = self.pcrartefact.denote(
            df_bam=self.df_bam,
            df_condition=self.df,
            condition_set=kwargs['granul_lvl_list'],
            umi_col=kwargs['umi_col'],
            new_col='UMI_mapped',
            pd_col='PD',
            inplace=False
        )
        # self.df_bam.to_csv(self.work_dir + 'self.df_bam.txt', sep='\t', index=False, header=True)
        # print(self.df_bam)
        # gp_df = self.df_bam.groupby(by=['spikeUMI'])
        # keys = gp_df.groups.keys()
        # for key in keys:
        #     df_sub = gp_df.get_group(key)
        #     print(key, df_sub['UMI_mapped'].unique().shape[0])
        self.console.print('======>finish finding deduplicated umis in {:.2f}s'.format(time.time() - dedup_umi_stime))
        # self.console.print('======># of umis deduplicated to be {}'.format(self.df['dedup_cnt'].loc['yes']))
        self.console.print('======>calculate average edit distances between umis...')
        self.df['ave_ed'] = self.df.apply(lambda x: self.umigadgetry.ed_ave(x, by_col='repr_nodes'), axis=1)
        self.df['num_diff_dedup_uniq_umis'] = self.df.apply(
            lambda x: self.umigadgetry.num_removed_uniq_umis(
                x,
                by_col='repr_nodes',
            ),
            axis=1,
        )
        self.df['num_diff_dedup_reads'] = self.df.apply(
            lambda x: self.umigadgetry.num_removed_reads(
                x,
                by_col='repr_nodes',
            ),
            axis=1,
        )
        self.console.print('======># of deduplicated unique umis {}'.format(self.df['num_diff_dedup_uniq_umis'].sum()))
        self.console.print('======># of deduplicated reads {}'.format(self.df['num_diff_dedup_reads'].sum()))
        self.ave_ed_bins = self.df['ave_ed'].value_counts().sort_index().to_frame().reset_index()
        self.console.check("======>bins for average edit distance\n{}".format(self.ave_ed_bins))
        if not self.heterogeneity:
            from umiche.util.Folder import Folder as crtfolder
            crtfolder().osmkdir(DIRECTORY=self.work_dir + '/starsolo/')
            self.fwriter.generic(
                df=self.ave_ed_bins,
                sv_fpn=self.work_dir + '/starsolo/' + self.bam_fn + '_' + str(kwargs['token']) + '_ave_ed_bin.txt',
                index=True,
                header=True,
            )
            # print(self.df)
            # print(self.df.columns)
            self.fwriter.generic(
                df=self.df[[
                    'dedup_cnt',
                    'ave_ed',
                    'num_uniq_umis',
                    'num_diff_dedup_uniq_umis',
                    'num_diff_dedup_reads',
                ]],
                sv_fpn=self.work_dir + '/starsolo/' + self.bam_fn + '_' + str(kwargs['token']) + '_dedup_sum.txt',
                index=True,
                header=True,
            )
            self.console.print('======>start writing deduplicated reads to BAM...')
            dedup_reads_write_stime = time.time()
            self.df['adj_bam_ids'] = self.df.apply(lambda x: self.umigadgetry.bamids(x, by_col='repr_nodes'), axis=1)
            self.aliwriter.tobam(
                tobam_fpn=self.work_dir + '/starsolo/' + self.bam_fn + '_' + str(kwargs['token']) + '_dedup.bam',
                tmpl_bam_fpn=self.bam_fpn,
                whitelist=self.umigadgetry.decompose(list_nd=self.df['adj_bam_ids'].values),
            )
            self.console.print('======>finish writing in {:.2f}s'.format(time.time() - dedup_reads_write_stime))
        return {
            'df': self.df,
            'df_bam': self.df_bam,
        }

    @Console.vignette()
    def gencore(
            self,
            **kwargs,
    ):
        # print('adj kwargs:', kwargs)
        from umiche.deduplicate.method.gencore.Gencore import Gencore

        dedup_umi_stime = time.time()
        gencore_ob = Gencore(verbose=True)
        self.df['clusters'] = self.df.apply(
            lambda x: gencore_ob.fit_graph(
                connected_components=x['cc'],
                counts=x['vignette']['df_umi_uniq_val_cnt'].to_dict(),
                graph_adj=x['vignette']['graph_adj'],
                int_to_umi_dict=x['vignette']['int_to_umi_dict'],
            )['clusters'],
            axis=1,
        )
        # print(self.df['clusters'])
        self.df['repr_node_map'] = self.df.apply(lambda x: self.umigadgetry.umimax(x, by_col='clusters', met='map'), axis=1)
        self.df['repr_nodes'] = self.df.apply(lambda x: self.umigadgetry.umimax(x, by_col='clusters'), axis=1)
        self.df['dedup_cnt'] = self.df['repr_nodes'].apply(lambda x: self.umigadgetry.length(x))
        # print(self.df['dedup_cnt'])
        # print(self.df)
        # print(self.df.columns)

        self.df_bam = self.pcrartefact.denote(
            df_bam=self.df_bam,
            df_condition=self.df,
            condition_set=kwargs['granul_lvl_list'],
            umi_col=kwargs['umi_col'],
            new_col='UMI_mapped',
            pd_col='PD',
            inplace=False
        )
        # self.df_bam.to_csv(self.work_dir + 'self.df_bam.txt', sep='\t', index=False, header=True)
        # print(self.df_bam)
        # gp_df = self.df_bam.groupby(by=['spikeUMI'])
        # keys = gp_df.groups.keys()
        # for key in keys:
        #     df_sub = gp_df.get_group(key)
        #     print(key, df_sub['UMI_mapped'].unique().shape[0])

        self.console.print('======>finish finding deduplicated umis in {:.2f}s'.format(time.time() - dedup_umi_stime))
        # self.console.print('======># of umis deduplicated to be {}'.format(self.df['dedup_cnt'].loc['yes']))
        self.console.print('======>calculate average edit distances between umis...')
        self.df['ave_ed'] = self.df.apply(lambda x: self.umigadgetry.ed_ave(x, by_col='repr_nodes'), axis=1)
        self.df['num_diff_dedup_uniq_umis'] = self.df.apply(
            lambda x: self.umigadgetry.num_removed_uniq_umis(
                x,
                by_col='repr_nodes',
            ),
            axis=1,
        )
        self.df['num_diff_dedup_reads'] = self.df.apply(
            lambda x: self.umigadgetry.num_removed_reads(
                x,
                by_col='repr_nodes',
            ),
            axis=1,
        )
        self.console.print('======># of deduplicated unique umis {}'.format(self.df['num_diff_dedup_uniq_umis'].sum()))
        self.console.print('======># of deduplicated reads {}'.format(self.df['num_diff_dedup_reads'].sum()))
        self.ave_ed_bins = self.df['ave_ed'].value_counts().sort_index().to_frame().reset_index()
        self.console.check("======>bins for average edit distance\n{}".format(self.ave_ed_bins))
        if not self.heterogeneity:
            from umiche.util.Folder import Folder as crtfolder
            crtfolder().osmkdir(DIRECTORY=self.work_dir + '/gencore/')

            self.fwriter.generic(
                df=self.ave_ed_bins,
                sv_fpn=self.work_dir + '/gencore/' + self.bam_fn + '_' + str(kwargs['token']) + '_ave_ed_bin.txt',
                index=True,
                header=True,
            )
            # print(self.df)
            # print(self.df.columns)
            self.fwriter.generic(
                df=self.df[[
                    'dedup_cnt',
                    'ave_ed',
                    'num_uniq_umis',
                    'num_diff_dedup_uniq_umis',
                    'num_diff_dedup_reads',
                ]],
                sv_fpn=self.work_dir + '/gencore/' + self.bam_fn + '_' + str(kwargs['token']) + '_dedup_sum.txt',
                index=True,
                header=True,
            )
            self.console.print('======>start writing deduplicated reads to BAM...')
            dedup_reads_write_stime = time.time()
            self.df['adj_bam_ids'] = self.df.apply(lambda x: self.umigadgetry.bamids(x, by_col='repr_nodes'), axis=1)
            self.aliwriter.tobam(
                tobam_fpn=self.work_dir + '/gencore/' + self.bam_fn + '_' + str(kwargs['token']) + '_dedup.bam',
                tmpl_bam_fpn=self.bam_fpn,
                whitelist=self.umigadgetry.decompose(list_nd=self.df['adj_bam_ids'].values),
            )
            self.console.print('======>finish writing in {:.2f}s'.format(time.time() - dedup_reads_write_stime))
        return {
            'df': self.df,
            'df_bam': self.df_bam,
        }

    def dropest(
            self,
            **kwargs,
    ):
        from umiche.deduplicate.method.dropest.DropEst import dropest_caller, write_dedup_bam
        df_sum = dropest_caller(
            self.df_bam,
            method='bayesian',
            cell_col=kwargs['granul_lvl_list'][0],
            gene_col=kwargs['granul_lvl_list'][1],
            # gene_col=kwargs['granul_lvl_list'][0],
            umi_col=kwargs['umi_col'],
            qual_col='umi_qual',
            return_type='summary',
        )
        # print(df_sum.head())
        if not self.heterogeneity:
            from umiche.util.Folder import Folder as crtfolder
            crtfolder().osmkdir(DIRECTORY=self.work_dir + '/dropest/')
            self.fwriter.generic(
                df=df_sum,
                sv_fpn=self.work_dir + '/dropest/' + self.bam_fn + '_' + str(kwargs['token']) + '_dedup_sum.txt',
                index=True,
                header=True,
            )
            write_dedup_bam(
                df=self.df_bam,
                src_bam_fpn=self.bam_fpn,
                out_bam_fpn=self.work_dir + '/dropest/' + self.bam_fn + '_' + str(kwargs['token']) + '_dedup.bam',
                method='bayesian',
                cell_col=kwargs['granul_lvl_list'][0],
                gene_col=kwargs['granul_lvl_list'][1],
                umi_col=kwargs['umi_col'],
                qual_col='umi_qual',
                keep_tag=None,
                threads=2,
            )
        return {
            'df': self.df,
            'df_bam': self.df_bam,
        }

    def irescue(
            self,
            **kwargs,
    ):
        from umiche.deduplicate.method.irescue.Irescue import Irescue

        irescuer = Irescue(
            cell_col=kwargs['granul_lvl_list'][0],
            feature_col=kwargs['granul_lvl_list'][1],
            # feature_col=kwargs['granul_lvl_list'][0],
            umi_col=kwargs['umi_col'],
            read_id_col=None,
            max_hd=kwargs['max_hd'],
            em_cycles=100,
            em_tol=1e-5,
            dump_ec=True,
            no_umi=False,
        )
        counts_long, df_sum = irescuer.fit_transform(self.df_bam)

        print(df_sum.head())
        import pandas as pd, ast
        def second_index_level(idx: pd.Index) -> pd.Index:
            if isinstance(idx, pd.MultiIndex) and idx.nlevels > 1:
                return idx.get_level_values(1)
            # for (a,b) or "(a, b)"
            out = []
            for x in idx:
                if isinstance(x, tuple) and len(x) >= 2:
                    out.append(x[1])
                elif isinstance(x, str) and x.startswith("("):
                    try:
                        t = ast.literal_eval(x)
                        out.append(t[1] if isinstance(t, tuple) and len(t) >= 2 else x)
                    except Exception:
                        out.append(x)
                else:
                    out.append(x)
            return pd.Index(out)
        df_sum.index = second_index_level(df_sum.index)
        print(df_sum.head())

        if not self.heterogeneity:
            from umiche.util.Folder import Folder as crtfolder
            crtfolder().osmkdir(DIRECTORY=self.work_dir + '/irescue/')
            self.fwriter.generic(
                df=df_sum,
                sv_fpn=self.work_dir + '/irescue/' + self.bam_fn + '_' + str(kwargs['token']) + '_dedup_sum.txt',
                index=True,
                header=True,
            )
            from umiche.deduplicate.method.irescue.Irescue import write_dedup_bam
            write_dedup_bam(
                df=self.df_bam,
                src_bam_fpn=self.bam_fpn,
                out_bam_fpn=self.work_dir + '/irescue/' + self.bam_fn + '_' + str(kwargs['token']) + '_dedup.bam',
                cell_col=kwargs['granul_lvl_list'][0],
                feature_col=kwargs['granul_lvl_list'][1],
                # feature_col=kwargs['granul_lvl_list'][0],
                umi_col=kwargs['umi_col'],
                read_col="read",
                read_id_col=None,
                max_hd=kwargs['max_hd'],
                no_umi=False,
                mark_only=False,
                duplicate_tag="PD",
                threads=4,
            )
            # print("wrote dedup BAM:", kept, "kept,", dups, "duplicates")
        return {
            'df': self.df,
            'df_bam': self.df_bam,
        }

    def umis(
            self,
            **kwargs,
    ):
        from umiche.deduplicate.method.umis.UMIS import UMIS, write_dedup_bam

        tag_map = {
            'cell': kwargs['granul_lvl_list'][0],
            'umi': kwargs['umi_col'],
            'gene': kwargs['granul_lvl_list'][0],
            'ref_name': 'chrom',
            'nh': 'NH',
            'pos': 'pos',
            'is_unmapped': 'is_unmapped',
        }

        counter = UMIS(tag_map=tag_map, min_evidence=1.0, weighted=True, positional=False)
        df_sum = counter.count_df(self.df_bam, drop_zeros=True)
        print(df_sum.head())
        import pandas as pd, ast
        def second_index_level(idx: pd.Index) -> pd.Index:
            if isinstance(idx, pd.MultiIndex) and idx.nlevels > 1:
                return idx.get_level_values(1)
            # for (a,b) or "(a, b)"
            out = []
            for x in idx:
                if isinstance(x, tuple) and len(x) >= 2:
                    out.append(x[1])
                elif isinstance(x, str) and x.startswith("("):
                    try:
                        t = ast.literal_eval(x)
                        out.append(t[1] if isinstance(t, tuple) and len(t) >= 2 else x)
                    except Exception:
                        out.append(x)
                else:
                    out.append(x)
            return pd.Index(out)

        df_sum.index = second_index_level(df_sum.index)
        print(df_sum.head())

        if not self.heterogeneity:
            from umiche.util.Folder import Folder as crtfolder
            crtfolder().osmkdir(DIRECTORY=self.work_dir + '/umis/')

            self.fwriter.generic(
                df=df_sum,
                sv_fpn=self.work_dir + '/umis/' + self.bam_fn + '_' + str(kwargs['token']) + '_dedup_sum.txt',
                index=True,
                header=True,
            )

            write_dedup_bam(
                df=self.df_bam,
                out_bam_fpn=self.work_dir + '/umis/' + self.bam_fn + '_' + str(kwargs['token']) + '_dedup.bam',
                template_bam_path=self.bam_fpn,
                tag_map=tag_map,
                min_evidence=1.0,
                weighted=True,
                positional=False,
                mode="filter",
                duplicate_tag="PD",
            )

        return {
            'df': self.df,
            'df_bam': self.df_bam,
        }