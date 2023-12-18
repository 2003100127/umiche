__version__ = "v1.0"
__copyright__ = "Copyright 2024"
__license__ = "MIT"
__developer__ = "Jianfeng Sun"
__maintainer__ = "Jianfeng Sun"
__email__="jianfeng.sunmt@gmail.com"
__lab__ = "Cribbslab"

import time

from umiche.align.Write import write as aliwrite

from umiche.deduplicate.Gadgetry import Gadgetry as umigadgetry

# dedup methods
from umiche.deduplicate.method.Adjacency import Adjacency as umiadj
from umiche.deduplicate.method.Directional import Directional as umidirec
from umiche.deduplicate.method.MarkovClustering import MarkovClustering as umimcl
from umiche.deduplicate.trimer.SetCoverOptimization import setCoverOptimization as umiscp

from umiche.util.Writer import writer as gwriter
from umiche.util.Console import Console


class Tabulate:

    def __init__(
            self,
            df,
            df_bam,
            bam_fpn,
            work_dir,
            verbose=False,
    ):
        self.df = df
        self.df_bam = df_bam
        self.bam_fpn = bam_fpn
        self.work_dir = work_dir

        self.aliwriter = aliwrite(df=self.df_bam)

        self.umigadgetry = umigadgetry()

        self.umiadj = umiadj()
        self.umidirec = umidirec()
        self.umiscp = umiscp()

        self.gwriter = gwriter()

        self.console = Console()
        self.console.verbose = verbose

    def set_cover(self, ):
        dedup_umi_stime = time.time()
        self.dedup_num = self.umiscp.count(
            inbam=self.bam_fpn,
            tag=self.pos_tag,
            sep='_',
        )
        print(self.dedup_num)
        # print(self.df_bam.columns)
        # print(self.df_bam)
        return

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
        self.df['num_diff_dedup_uniq_umis'] = self.df.apply(lambda x: self.umigadgetry.num_removed_uniq_umis(x, by_col='uniq_repr_nodes'), axis=1)
        self.df['num_diff_dedup_reads'] = self.df.apply(lambda x: self.umigadgetry.num_removed_reads(x, by_col='uniq_repr_nodes'),
                                                       axis=1)
        self.console.print('======># of deduplicated unique umis {}'.format(self.df['num_diff_dedup_uniq_umis'].sum()))
        self.console.print('======># of deduplicated reads {}'.format(self.df['num_diff_dedup_reads'].sum()))
        self.gwriter.generic(
            df=self.ave_ed_bins,
            sv_fpn=self.work_dir + 'unique_ave_ed_bin.txt',
            index=True,
        )
        self.gwriter.generic(
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
        self.df['cc_umi_len'] = self.df['cc_repr_nodes'].apply(lambda x: self.umigadgetry.length(x))
        ### @@ self.df['cc_umi_len']
        # 1    1
        # Name: cc_umi_len, dtype: int64
        self.console.print('======>finish finding deduplicated umis in {:.2f}s'.format(time.time() - dedup_umi_stime))
        # self.console.print('======># of umis deduplicated to be {}'.format(self.df['cc_umi_len'].loc['yes']))
        self.console.print('======>calculate average edit distances between umis...')
        self.df['ave_eds'] = self.df.apply(lambda x: self.umigadgetry.ed_ave(x, by_col='cc_repr_nodes'), axis=1)
        self.df['num_diff_dedup_uniq_umis'] = self.df.apply(lambda x: self.umigadgetry.num_removed_uniq_umis(x, by_col='cc_repr_nodes'),
                                                       axis=1)
        self.df['num_diff_dedup_reads'] = self.df.apply(lambda x: self.umigadgetry.num_removed_reads(x, by_col='cc_repr_nodes'),
                                                       axis=1)
        self.console.print('======># of deduplicated unique umis {}'.format(self.df['num_diff_dedup_uniq_umis'].sum()))
        self.console.print('======># of deduplicated reads {}'.format(self.df['num_diff_dedup_reads'].sum()))
        self.ave_ed_bins = self.df['ave_eds'].value_counts().sort_index().to_frame().reset_index()
        self.console.check("======>bins for average edit distance\n{}".format(self.ave_ed_bins))
        self.gwriter.generic(
            df=self.ave_ed_bins,
            sv_fpn=self.work_dir + 'cluster_ave_ed_bin.txt',
            index=True,
        )
        self.gwriter.generic(
            df=self.df[[
                'cc_umi_len',
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
        self.df['adj_repr_nodes'] = self.df.apply(lambda x: self.umigadgetry.umimax(x, by_col='adj'), axis=1)
        self.df['adj_umi_len'] = self.df['adj_repr_nodes'].apply(lambda x: self.umigadgetry.length(x))
        self.console.print('======>finish finding deduplicated umis in {:.2f}s'.format(time.time() - dedup_umi_stime))
        # self.console.print('======># of umis deduplicated to be {}'.format(self.df['adj_umi_len'].loc['yes']))
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
        self.gwriter.generic(
            df=self.ave_ed_bins,
            sv_fpn=self.work_dir + 'adjacency_ave_ed_bin.txt',
            index=True,
        )
        self.gwriter.generic(
            df=self.df[[
                'adj_umi_len',
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
        self.df['direc_repr_nodes'] = self.df.apply(lambda x: self.umigadgetry.umimax(x, by_col='direc'), axis=1)
        print(self.df['direc_repr_nodes'])
        self.df['direc_umi_len'] = self.df['direc_repr_nodes'].apply(lambda x: self.umigadgetry.length(x))
        self.console.print('======>finish finding deduplicated umis in {:.2f}s'.format(time.time() - dedup_umi_stime))
        #             self.console.print('======># of umis deduplicated to be {}'.format(self.df['direc_umi_len'].loc['yes']))
        self.console.print('======>calculate average edit distances between umis...')
        self.df['ave_eds'] = self.df.apply(lambda x: self.umigadgetry.ed_ave(x, by_col='direc_repr_nodes'), axis=1)
        print(self.df['ave_eds'])
        self.df['num_diff_dedup_uniq_umis'] = self.df.apply(
            lambda x: self.umigadgetry.num_removed_uniq_umis(x, by_col='direc_repr_nodes'), axis=1)
        self.df['num_diff_dedup_reads'] = self.df.apply(lambda x: self.umigadgetry.num_removed_reads(x, by_col='direc_repr_nodes'),
                                                       axis=1)
        self.console.print('======># of deduplicated unique umis {}'.format(self.df['num_diff_dedup_uniq_umis'].sum()))
        self.console.print('======># of deduplicated reads {}'.format(self.df['num_diff_dedup_reads'].sum()))
        self.ave_ed_bins = self.df['ave_eds'].value_counts().sort_index().to_frame().reset_index()
        print(self.ave_ed_bins)
        self.console.check("======>bins for average edit distance\n{}".format(self.ave_ed_bins))
        self.gwriter.generic(
            df=self.ave_ed_bins,
            sv_fpn=self.work_dir + 'directional_ave_ed_bin.txt',
            index=True,
        )
        self.gwriter.generic(
            df=self.df[[
                'direc_umi_len',
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
            inflat_val=2.0,
            exp_val=2,
            iter_num=100,

    ):
        dedup_umi_stime = time.time()
        self.umimcl = umimcl(
            inflat_val=inflat_val,
            exp_val=exp_val,
            iter_num=iter_num,
        )
        self.df['mcl'] = self.df.apply(
            lambda x: self.umimcl.decompose(
                list_nd=self.umimcl.dfclusters(
                    connected_components=x['cc'],
                    graph_adj=x['vignette']['graph_adj'],
                )['clusters'].values,
            ),
            axis=1,
        )
        self.df['mcl_repr_nodes'] = self.df.apply(lambda x: self.umigadgetry.umimax(x, by_col='mcl'), axis=1)
        self.df['mcl_umi_len'] = self.df['mcl_repr_nodes'].apply(lambda x: self.umigadgetry.length(x))
        self.console.print('======>finish finding deduplicated umis in {:.2f}s'.format(time.time() - dedup_umi_stime))
        #             self.console.print('======># of umis deduplicated to be {}'.format(self.df['mcl_umi_len'].loc['yes']))
        self.console.print('======>calculate average edit distances between umis...')
        self.df['ave_eds'] = self.df.apply(lambda x: self.umigadgetry.ed_ave(x, by_col='mcl_repr_nodes'), axis=1)
        self.df['num_diff_dedup_uniq_umis'] = self.df.apply(lambda x: self.umigadgetry.num_removed_uniq_umis(x, by_col='mcl_repr_nodes'),
                                                       axis=1)
        self.df['num_diff_dedup_reads'] = self.df.apply(lambda x: self.umigadgetry.num_removed_reads(x, by_col='mcl_repr_nodes'),
                                                       axis=1)
        self.console.print('======># of deduplicated unique umis {}'.format(self.df['num_diff_dedup_uniq_umis'].sum()))
        self.console.print('======># of deduplicated reads {}'.format(self.df['num_diff_dedup_reads'].sum()))
        self.ave_ed_bins = self.df['ave_eds'].value_counts().sort_index().to_frame().reset_index()
        self.console.check("======>bins for average edit distance\n{}".format(self.ave_ed_bins))
        self.gwriter.generic(
            df=self.ave_ed_bins,
            sv_fpn=self.work_dir + 'mcl_ave_ed_bin.txt',
            index=True,
        )
        self.gwriter.generic(
            df=self.df[[
                'mcl_umi_len',
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
            inflat_val=2.0,
            exp_val=2,
            iter_num=100,
            mcl_fold_thres=2,
    ):
        self.umimcl = umimcl(
            inflat_val=inflat_val,
            exp_val=exp_val,
            iter_num=iter_num,
        )
        dedup_umi_stime = time.time()
        self.df['mcl_val'] = self.df.apply(
            lambda x: self.umimcl.decompose(
                list_nd=self.umimcl.maxval_val(
                    df_mcl_ccs=self.umimcl.dfclusters(
                        connected_components=x['cc'],
                        graph_adj=x['vignette']['graph_adj'],
                    ),
                    df_umi_uniq_val_cnt=x['vignette']['df_umi_uniq_val_cnt'],
                    thres_fold=mcl_fold_thres,
                )['clusters'].values,
            ),
            axis=1,
        )
        self.df['mcl_val_repr_nodes'] = self.df.apply(lambda x: self.umigadgetry.umimax(x, by_col='mcl_val'), axis=1)
        self.df['mcl_val_umi_len'] = self.df['mcl_val_repr_nodes'].apply(lambda x: self.umigadgetry.length(x))
        self.console.print('======>finish finding deduplicated umis in {:.2f}s'.format(time.time() - dedup_umi_stime))
        #             self.console.print('======># of umis deduplicated to be {}'.format(self.df['mcl_val_umi_len'].loc['yes']))
        self.console.print('======>calculate average edit distances between umis...')
        self.df['ave_eds'] = self.df.apply(lambda x: self.umigadgetry.ed_ave(x, by_col='mcl_val_repr_nodes'), axis=1)
        self.df['num_diff_dedup_uniq_umis'] = self.df.apply(
            lambda x: self.umigadgetry.num_removed_uniq_umis(x, by_col='mcl_val_repr_nodes'), axis=1)
        self.df['num_diff_dedup_reads'] = self.df.apply(lambda x: self.umigadgetry.num_removed_reads(x, by_col='mcl_val_repr_nodes'),
                                                       axis=1)
        self.console.print('======># of deduplicated unique umis {}'.format(self.df['num_diff_dedup_uniq_umis'].sum()))
        self.console.print('======># of deduplicated reads {}'.format(self.df['num_diff_dedup_reads'].sum()))
        self.ave_ed_bins = self.df['ave_eds'].value_counts().sort_index().to_frame().reset_index()
        self.console.check("======>bins for average edit distance\n{}".format(self.ave_ed_bins))
        self.gwriter.generic(
            df=self.ave_ed_bins,
            sv_fpn=self.work_dir + 'mcl_val_ave_ed_bin.txt',
            index=True,
        )
        self.gwriter.generic(
            df=self.df[[
                'mcl_val_umi_len',
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
        self.df['mcl_val_bam_ids'] = self.df.apply(lambda x: self.umigadgetry.bamids(x, by_col='mcl_val_repr_nodes'), axis=1)
        self.aliwriter.tobam(
            tobam_fpn=self.work_dir + 'mcl_val_dedup.bam',
            tmpl_bam_fpn=self.bam_fpn,
            whitelist=self.umigadgetry.decompose(list_nd=self.df['mcl_val_bam_ids'].values),
        )
        self.console.print('======>finish writing in {:.2f}s'.format(time.time() - dedup_reads_write_stime))
        return self.df

    def mcl_ed(
            self,
            inflat_val=2.0,
            exp_val=2,
            iter_num=100,
            mcl_fold_thres=2,
    ):
        self.umimcl = umimcl(
            inflat_val=inflat_val,
            exp_val=exp_val,
            iter_num=iter_num,
        )
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
                    thres_fold=mcl_fold_thres,
                )['clusters'].values,
            ),
            axis=1,
        )
        self.df['mcl_ed_repr_nodes'] = self.df.apply(lambda x: self.umigadgetry.umimax(x, by_col='mcl_ed'), axis=1)
        self.df['mcl_ed_umi_len'] = self.df['mcl_ed_repr_nodes'].apply(lambda x: self.umigadgetry.length(x))
        self.console.print('======>finish finding deduplicated umis in {:.2f}s'.format(time.time() - dedup_umi_stime))
        #             self.console.print('======># of umis deduplicated to be {}'.format(self.df['mcl_ed_umi_len'].loc['yes']))
        self.console.print('======>calculate average edit distances between umis...')
        self.df['ave_eds'] = self.df.apply(lambda x: self.umigadgetry.ed_ave(x, by_col='mcl_ed_repr_nodes'), axis=1)
        self.df['num_diff_dedup_uniq_umis'] = self.df.apply(
            lambda x: self.umigadgetry.num_removed_uniq_umis(x, by_col='mcl_ed_repr_nodes'), axis=1)
        self.df['num_diff_dedup_reads'] = self.df.apply(lambda x: self.umigadgetry.num_removed_reads(x, by_col='mcl_ed_repr_nodes'),
                                                       axis=1)
        self.console.print('======># of deduplicated unique umis {}'.format(self.df['num_diff_dedup_uniq_umis'].sum()))
        self.console.print('======># of deduplicated reads {}'.format(self.df['num_diff_dedup_reads'].sum()))
        self.ave_ed_bins = self.df['ave_eds'].value_counts().sort_index().to_frame().reset_index()
        self.console.check("======>bins for average edit distance\n{}".format(self.ave_ed_bins))
        self.gwriter.generic(
            df=self.ave_ed_bins,
            sv_fpn=self.work_dir + 'mcl_ed_ave_ed_bin.txt',
            index=True,
        )
        self.gwriter.generic(
            df=self.df[[
                'mcl_ed_umi_len',
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
