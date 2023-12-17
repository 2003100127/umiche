__version__ = "v1.0"
__copyright__ = "Copyright 2024"
__license__ = "MIT"
__developer__ = "Jianfeng Sun"
__maintainer__ = "Jianfeng Sun"
__email__="jianfeng.sunmt@gmail.com"
__lab__ = "Cribbslab"

import time
import pandas as pd

from umiche.align.Read import read as aliread
from umiche.align.Write import write as aliwrite

from umiche.deduplicate.Build import Build as umibuild
from umiche.deduplicate.Gadgetry import Gadgetry as umigadgetry

# dedup methods
from umiche.deduplicate.method.Cluster import Cluster as umiclust
from umiche.deduplicate.method.Adjacency import Adjacency as umiadj
from umiche.deduplicate.method.Directional import Directional as umidirec
from umiche.deduplicate.method.MarkovClustering import MarkovClustering as umimcl

from umiche.util.Number import number as rannum
from umiche.util.Writer import writer as gwriter
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
        self.umigadgetry = umigadgetry()
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

        self.df['num_uniq_umis'] = self.df['uniq_repr_nodes'].apply(lambda x: self.umigadgetry.length(x))
        ### @@ self.df['num_uniq_umis']
        # 1    1949
        # Name: num_uniq_umis, dtype: int64
        self.console.print('===>start deduplication by the {} method...'.format(self.method))
        # sys.stdout.close()

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
        # self.df['ave_ed'] = self.df.apply(lambda x: self.edave(x, by_col='uniq_repr_nodes'), axis=1)
        self.ave_ed_bins = self.df['ave_ed'].value_counts().sort_index().to_frame().reset_index()
        self.console.check("======>bins for average edit distance\n{}".format(self.ave_ed_bins))
        self.df['num_diff_dedup_uniq_umis'] = self.df.apply(lambda x: self.umigadgetry.num_removed_uniq_umis(x, by_col='uniq_repr_nodes'), axis=1)
        self.df['num_diff_dedup_reads'] = self.df.apply(lambda x: self.umigadgetry.num_removed_reads(x, by_col='uniq_repr_nodes'),
                                                       axis=1)
        self.console.print('======># of deduplicated unique umis {}'.format(self.df['num_diff_dedup_uniq_umis'].sum()))
        self.console.print('======># of deduplicated reads {}'.format(self.df['num_diff_dedup_reads'].sum()))
        self.gwriter.generic(df=self.ave_ed_bins, sv_fpn=self.work_dir + 'uniq_ave_ed_bin.txt', index=True)
        self.df_dedup_sum = self.df[[
            'ave_ed',
            'num_uniq_umis',
            'num_diff_dedup_uniq_umis',
            'num_diff_dedup_reads',
        ]]
        self.gwriter.generic(
            df=self.df_dedup_sum,
            sv_fpn=self.work_dir + 'uniq_dedup_sum.txt',
            index=True,
            header=True,
        )
        self.console.print('======>start writing deduplicated reads to BAM...')
        dedup_reads_write_stime = time.time()
        self.df['uniq_bam_ids'] = self.df.apply(lambda x: self.umigadgetry.bamids(x, by_col='uniq_repr_nodes'), axis=1)
        self.aliwriter.tobam(
            tobam_fpn=self.work_dir + self.method + '_dedup.bam',
            tmpl_bam_fpn=self.bam_fpn,
            whitelist=self.umigadgetry.decompose(list_nd=self.df['uniq_bam_ids'].values),
        )
        self.console.print('======>finish writing in {:.2f}s'.format(time.time() - dedup_reads_write_stime))
        return

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
        # self.df['ave_ed'] = self.df.apply(lambda x: self.edave(x, by_col='cc_repr_nodes'), axis=1)
        self.df['num_diff_dedup_uniq_umis'] = self.df.apply(lambda x: self.umigadgetry.num_removed_uniq_umis(x, by_col='cc_repr_nodes'),
                                                       axis=1)
        self.df['num_diff_dedup_reads'] = self.df.apply(lambda x: self.umigadgetry.num_removed_reads(x, by_col='cc_repr_nodes'),
                                                       axis=1)
        self.console.print('======># of deduplicated unique umis {}'.format(self.df['num_diff_dedup_uniq_umis'].sum()))
        self.console.print('======># of deduplicated reads {}'.format(self.df['num_diff_dedup_reads'].sum()))
        self.ave_ed_bins = self.df['ave_ed'].value_counts().sort_index().to_frame().reset_index()
        self.console.check("======>bins for average edit distance\n{}".format(self.ave_ed_bins))
        self.gwriter.generic(df=self.ave_ed_bins, sv_fpn=self.work_dir + 'cc_ave_ed_bin.txt', index=True)
        self.df_dedup_sum = self.df[[
            'cc_umi_len',
            'ave_ed',
            'num_uniq_umis',
            'num_diff_dedup_uniq_umis',
            'num_diff_dedup_reads',
        ]]
        self.gwriter.generic(
            df=self.df_dedup_sum,
            sv_fpn=self.work_dir + 'cc_dedup_sum.txt',
            index=True,
            header=True,
        )
        self.console.print('======>start writing deduplicated reads to BAM...')
        dedup_reads_write_stime = time.time()
        self.df['cc_bam_ids'] = self.df.apply(lambda x: self.umigadgetry.bamids(x, by_col='cc_repr_nodes'), axis=1)
        self.aliwriter.tobam(
            tobam_fpn=self.work_dir + self.method + '_dedup.bam',
            tmpl_bam_fpn=self.bam_fpn,
            whitelist=self.umigadgetry.decompose(list_nd=self.df['cc_bam_ids'].values),
        )
        self.console.print('======>finish writing in {:.2f}s'.format(time.time() - dedup_reads_write_stime))
        return

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
        #             self.df['ave_ed'] = self.df.apply(lambda x: self.edave(x, by_col='adj_repr_nodes'), axis=1)
        self.df['num_diff_dedup_uniq_umis'] = self.df.apply(lambda x: self.umigadgetry.num_removed_uniq_umis(x, by_col='adj_repr_nodes'),
                                                       axis=1)
        self.df['num_diff_dedup_reads'] = self.df.apply(lambda x: self.umigadgetry.num_removed_reads(x, by_col='adj_repr_nodes'),
                                                       axis=1)
        self.console.print('======># of deduplicated unique umis {}'.format(self.df['num_diff_dedup_uniq_umis'].sum()))
        self.console.print('======># of deduplicated reads {}'.format(self.df['num_diff_dedup_reads'].sum()))
        self.ave_ed_bins = self.df['ave_ed'].value_counts().sort_index().to_frame().reset_index()
        self.console.check("======>bins for average edit distance\n{}".format(self.ave_ed_bins))
        self.gwriter.generic(df=self.ave_ed_bins, sv_fpn=self.work_dir + 'adj_ave_ed_bin.txt', index=True)
        self.df_dedup_sum = self.df[[
            'adj_umi_len',
            'ave_ed',
            'num_uniq_umis',
            'num_diff_dedup_uniq_umis',
            'num_diff_dedup_reads',
        ]]
        self.gwriter.generic(
            df=self.df_dedup_sum,
            sv_fpn=self.work_dir + 'adj_dedup_sum.txt',
            index=True,
            header=True,
        )
        self.console.print('======>start writing deduplicated reads to BAM...')
        dedup_reads_write_stime = time.time()
        self.df['adj_bam_ids'] = self.df.apply(lambda x: self.umigadgetry.bamids(x, by_col='adj_repr_nodes'), axis=1)
        self.aliwriter.tobam(
            tobam_fpn=self.work_dir + self.method + '_dedup.bam',
            tmpl_bam_fpn=self.bam_fpn,
            whitelist=self.umigadgetry.decompose(list_nd=self.df['adj_bam_ids'].values),
        )
        self.console.print('======>finish writing in {:.2f}s'.format(time.time() - dedup_reads_write_stime))
        return

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
        self.df['direc_umi_len'] = self.df['direc_repr_nodes'].apply(lambda x: self.umigadgetry.length(x))
        self.console.print('======>finish finding deduplicated umis in {:.2f}s'.format(time.time() - dedup_umi_stime))
        #             self.console.print('======># of umis deduplicated to be {}'.format(self.df['direc_umi_len'].loc['yes']))
        self.console.print('======>calculate average edit distances between umis...')
        #             self.df['ave_ed'] = self.df.apply(lambda x: self.edave(x, by_col='direc_repr_nodes'), axis=1)
        self.df['num_diff_dedup_uniq_umis'] = self.df.apply(
            lambda x: self.umigadgetry.num_removed_uniq_umis(x, by_col='direc_repr_nodes'), axis=1)
        self.df['num_diff_dedup_reads'] = self.df.apply(lambda x: self.umigadgetry.num_removed_reads(x, by_col='direc_repr_nodes'),
                                                       axis=1)
        self.console.print('======># of deduplicated unique umis {}'.format(self.df['num_diff_dedup_uniq_umis'].sum()))
        self.console.print('======># of deduplicated reads {}'.format(self.df['num_diff_dedup_reads'].sum()))
        self.ave_ed_bins = self.df['ave_ed'].value_counts().sort_index().to_frame().reset_index()
        self.console.check("======>bins for average edit distance\n{}".format(self.ave_ed_bins))
        self.gwriter.generic(df=self.ave_ed_bins, sv_fpn=self.work_dir + 'direc_ave_ed_bin.txt', index=True)
        self.df_dedup_sum = self.df[[
            'direc_umi_len',
            'ave_ed',
            'num_uniq_umis',
            'num_diff_dedup_uniq_umis',
            'num_diff_dedup_reads',
        ]]
        self.gwriter.generic(
            df=self.df_dedup_sum,
            sv_fpn=self.work_dir + 'direc_dedup_sum.txt',
            index=True,
            header=True,
        )
        self.console.print('======>start writing deduplicated reads to BAM...')
        dedup_reads_write_stime = time.time()
        self.df['direc_bam_ids'] = self.df.apply(lambda x: self.umigadgetry.bamids(x, by_col='direc_repr_nodes'), axis=1)
        self.aliwriter.tobam(
            tobam_fpn=self.work_dir + self.method + '_dedup.bam',
            tmpl_bam_fpn=self.bam_fpn,
            whitelist=self.umigadgetry.decompose(list_nd=self.df['direc_bam_ids'].values),
        )
        self.console.print('======>finish writing in {:.2f}s'.format(time.time() - dedup_reads_write_stime))

    def mcl(self, ):
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
        self.df['mcl_repr_nodes'] = self.df.apply(lambda x: self.umigadgetry.umimax(x, by_col='mcl'), axis=1)
        self.df['mcl_umi_len'] = self.df['mcl_repr_nodes'].apply(lambda x: self.umigadgetry.length(x))
        self.console.print('======>finish finding deduplicated umis in {:.2f}s'.format(time.time() - dedup_umi_stime))
        #             self.console.print('======># of umis deduplicated to be {}'.format(self.df['mcl_umi_len'].loc['yes']))
        self.console.print('======>calculate average edit distances between umis...')
        #             self.df['ave_ed'] = self.df.apply(lambda x: self.edave(x, by_col='mcl_repr_nodes'), axis=1)
        self.df['num_diff_dedup_uniq_umis'] = self.df.apply(lambda x: self.umigadgetry.num_removed_uniq_umis(x, by_col='mcl_repr_nodes'),
                                                       axis=1)
        self.df['num_diff_dedup_reads'] = self.df.apply(lambda x: self.umigadgetry.num_removed_reads(x, by_col='mcl_repr_nodes'),
                                                       axis=1)
        self.console.print('======># of deduplicated unique umis {}'.format(self.df['num_diff_dedup_uniq_umis'].sum()))
        self.console.print('======># of deduplicated reads {}'.format(self.df['num_diff_dedup_reads'].sum()))
        self.ave_ed_bins = self.df['ave_ed'].value_counts().sort_index().to_frame().reset_index()
        self.console.check("======>bins for average edit distance\n{}".format(self.ave_ed_bins))
        self.gwriter.generic(df=self.ave_ed_bins, sv_fpn=self.work_dir + 'mcl_ave_ed_bin.txt', index=True)
        self.df_dedup_sum = self.df[[
            'mcl_umi_len',
            'ave_ed',
            'num_uniq_umis',
            'num_diff_dedup_uniq_umis',
            'num_diff_dedup_reads',
        ]]
        self.gwriter.generic(
            df=self.df_dedup_sum,
            sv_fpn=self.work_dir + 'mcl_dedup_sum.txt',
            index=True,
            header=True,
        )
        self.console.print('======>start writing deduplicated reads to BAM...')
        dedup_reads_write_stime = time.time()
        self.df['mcl_bam_ids'] = self.df.apply(lambda x: self.umigadgetry.bamids(x, by_col='mcl_repr_nodes'), axis=1)
        self.aliwriter.tobam(
            tobam_fpn=self.work_dir + self.method + '_dedup.bam',
            tmpl_bam_fpn=self.bam_fpn,
            whitelist=self.umigadgetry.decompose(list_nd=self.df['mcl_bam_ids'].values),
        )
        self.console.print('======>finish writing in {:.2f}s'.format(time.time() - dedup_reads_write_stime))

    def mcl_val(self, ):
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
        self.df['mcl_val_repr_nodes'] = self.df.apply(lambda x: self.umigadgetry.umimax(x, by_col='mcl_val'), axis=1)
        self.df['mcl_val_umi_len'] = self.df['mcl_val_repr_nodes'].apply(lambda x: self.umigadgetry.length(x))
        self.console.print('======>finish finding deduplicated umis in {:.2f}s'.format(time.time() - dedup_umi_stime))
        #             self.console.print('======># of umis deduplicated to be {}'.format(self.df['mcl_val_umi_len'].loc['yes']))
        self.console.print('======>calculate average edit distances between umis...')
        #             self.df['ave_ed'] = self.df.apply(lambda x: self.edave(x, by_col='mcl_val_repr_nodes'), axis=1)
        self.df['num_diff_dedup_uniq_umis'] = self.df.apply(
            lambda x: self.umigadgetry.num_removed_uniq_umis(x, by_col='mcl_val_repr_nodes'), axis=1)
        self.df['num_diff_dedup_reads'] = self.df.apply(lambda x: self.umigadgetry.num_removed_reads(x, by_col='mcl_val_repr_nodes'),
                                                       axis=1)
        self.console.print('======># of deduplicated unique umis {}'.format(self.df['num_diff_dedup_uniq_umis'].sum()))
        self.console.print('======># of deduplicated reads {}'.format(self.df['num_diff_dedup_reads'].sum()))
        self.ave_ed_bins = self.df['ave_ed'].value_counts().sort_index().to_frame().reset_index()
        self.console.check("======>bins for average edit distance\n{}".format(self.ave_ed_bins))
        self.gwriter.generic(df=self.ave_ed_bins, sv_fpn=self.work_dir + 'mcl_val_ave_ed_bin.txt', index=True)
        self.df_dedup_sum = self.df[[
            'mcl_val_umi_len',
            'ave_ed',
            'num_uniq_umis',
            'num_diff_dedup_uniq_umis',
            'num_diff_dedup_reads',
        ]]
        self.gwriter.generic(
            df=self.df_dedup_sum,
            sv_fpn=self.work_dir + 'mcl_val_dedup_sum.txt',
            index=True,
            header=True,
        )
        self.console.print('======>start writing deduplicated reads to BAM...')
        dedup_reads_write_stime = time.time()
        self.df['mcl_val_bam_ids'] = self.df.apply(lambda x: self.umigadgetry.bamids(x, by_col='mcl_val_repr_nodes'), axis=1)
        self.aliwriter.tobam(
            tobam_fpn=self.work_dir + self.method + '_dedup.bam',
            tmpl_bam_fpn=self.bam_fpn,
            whitelist=self.umigadgetry.decompose(list_nd=self.df['mcl_val_bam_ids'].values),
        )
        self.console.print('======>finish writing in {:.2f}s'.format(time.time() - dedup_reads_write_stime))
        return

    def mcl_ed(self, ):
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
        self.df['mcl_ed_repr_nodes'] = self.df.apply(lambda x: self.umigadgetry.umimax(x, by_col='mcl_ed'), axis=1)
        self.df['mcl_ed_umi_len'] = self.df['mcl_ed_repr_nodes'].apply(lambda x: self.umigadgetry.length(x))
        self.console.print('======>finish finding deduplicated umis in {:.2f}s'.format(time.time() - dedup_umi_stime))
        #             self.console.print('======># of umis deduplicated to be {}'.format(self.df['mcl_ed_umi_len'].loc['yes']))
        self.console.print('======>calculate average edit distances between umis...')
        #             self.df['ave_ed'] = self.df.apply(lambda x: self.edave(x, by_col='mcl_ed_repr_nodes'), axis=1)
        self.df['num_diff_dedup_uniq_umis'] = self.df.apply(
            lambda x: self.umigadgetry.num_removed_uniq_umis(x, by_col='mcl_ed_repr_nodes'), axis=1)
        self.df['num_diff_dedup_reads'] = self.df.apply(lambda x: self.umigadgetry.num_removed_reads(x, by_col='mcl_ed_repr_nodes'),
                                                       axis=1)
        self.console.print('======># of deduplicated unique umis {}'.format(self.df['num_diff_dedup_uniq_umis'].sum()))
        self.console.print('======># of deduplicated reads {}'.format(self.df['num_diff_dedup_reads'].sum()))
        self.ave_ed_bins = self.df['ave_ed'].value_counts().sort_index().to_frame().reset_index()
        self.console.check("======>bins for average edit distance\n{}".format(self.ave_ed_bins))
        self.gwriter.generic(df=self.ave_ed_bins, sv_fpn=self.work_dir + 'mcl_ed_ave_ed_bin.txt', index=True)
        self.df_dedup_sum = self.df[[
            'mcl_ed_umi_len',
            'ave_ed',
            'num_uniq_umis',
            'num_diff_dedup_uniq_umis',
            'num_diff_dedup_reads',
        ]]
        self.gwriter.generic(
            df=self.df_dedup_sum,
            sv_fpn=self.work_dir + 'mcl_ed_dedup_sum.txt',
            index=True,
            header=True,
        )
        self.console.print('======>start writing deduplicated reads to BAM...')
        dedup_reads_write_stime = time.time()
        self.df['mcl_ed_bam_ids'] = self.df.apply(lambda x: self.umigadgetry.bamids(x, by_col='mcl_ed_repr_nodes'), axis=1)
        self.aliwriter.tobam(
            tobam_fpn=self.work_dir + self.method + '_dedup.bam',
            tmpl_bam_fpn=self.bam_fpn,
            whitelist=self.umigadgetry.decompose(list_nd=self.df['mcl_ed_bam_ids'].values),
        )
        self.console.print('======>finish writing in {:.2f}s'.format(time.time() - dedup_reads_write_stime))


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
        work_dir=to('data/'),

        verbose=True,
    )