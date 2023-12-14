__version__ = "v1.0"
__copyright__ = "Copyright 2023"
__license__ = "MIT"
__developer__ = "Jianfeng Sun"
__lab__ = "cribbslab"

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
            sv_fpn='./dedup.bam',
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
        self.sv_fpn = sv_fpn
        self.verbose = verbose

        self.umibuild = umibuild
        self.rannum = rannum()
        self.gwriter = gwriter()
        self.console = Console()
        self.console.verbose = self.verbose

        self.dirname = os.path.dirname(self.sv_fpn) + '/'
        # sys.stdout = open(self.dirname + self.method + '_log.txt', 'w')


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
        self.console.print('======># of raw reads: {}'.format(self.df_bam.shape[0]))
        self.df_bam = self.df_bam.loc[self.df_bam['reference_id'] != -1]
        self.console.print('======># of reads with qualified chrs: {}'.format(self.df_bam.shape[0]))

        self.df_bam['umi'] = self.df_bam['query_name'].apply(lambda x: x.split(umi_)[1])
        self.console.print('======># of unique umis: {}'.format(self.df_bam['umi'].unique().shape[0]))
        self.console.print('======># of redundant umis: {}'.format(self.df_bam['umi'].shape[0]))
        self.console.print('======>edit distance thres: {}'.format(self.ed_thres))

        self.df_bam['source'] = 1
        self.df_bam_gp = self.df_bam.groupby(by=['source'])
        self.pos_gp_keys = self.df_bam_gp.groups.keys()

        self.aliwriter = aliwrite(df=self.df_bam)

        self.console.print('======># of columns in the bam df: {}'.format(len(self.df_bam.columns)))
        self.console.print('======>Columns in the bam df: {}'.format(self.df_bam.columns.tolist()))
        self.console.print('======># of raw reads:\n{}'.format(self.df_bam))

        self.console.print('===>start building umi graphs...')
        umi_graph_build_stime = time.time()
        gps = []
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
            import json
            with open('data.json', 'w') as f:
                json.dump(cc, f)
            gps.append(pos_g)
            res_sum.append([
                umi_vignette,
                cc,
                [*umi_vignette['int_to_umi_dict'].keys()],
            ])
        self.df = pd.DataFrame(
            data=res_sum,
            columns=['vignette', 'cc', 'uniq_repr_nodes'],
            index=gps,
        )
        print(self.df.columns)
        print(self.df)

        self.console.print('===>time for building umi graphs: {:.2f}s'.format(time.time() - umi_graph_build_stime))

        self.df['uniq_umi_len'] = self.df['uniq_repr_nodes'].apply(lambda x: self.length(x))

        self.console.print('===>start deduplication by the {} method...'.format(self.method))
        if self.method == 'unique':
            dedup_umi_stime = time.time()
            # self.df['uniq_sgl_mark'] = self.df['uniq_repr_nodes'].apply(lambda x: self.markSingleUMI(x))
            # self.df = self.df.loc[self.df['uniq_sgl_mark'] == 'no']
            self.console.print('======># of positions with non-single umis: {}'.format(self.df.shape[0]))
            self.console.print('======>finish finding deduplicated umis in {:.2f}s'.format(time.time() - dedup_umi_stime))
            self.console.print('======># of umis deduplicated to be {}'.format(self.df['uniq_umi_len'].loc['yes']))
            self.console.print('======>calculate average edit distances between umis...')
            dedup_umi_edave_stime = time.time()
            self.df['ave_eds'] = self.df.apply(lambda x: self.edave(x, by_col='uniq_repr_nodes'), axis=1)
            self.console.print('======>finish calculating ave eds in {:.2f}s'.format(time.time() - dedup_umi_edave_stime))
            ave_ed_bins = self.df['ave_eds'].value_counts().sort_index()
            self.console.check(ave_ed_bins)
            self.df['dedup_uniq_diff_pos'] = self.df.apply(lambda x: self.diffDedupUniqCountPos(x, by_col='uniq_repr_nodes'), axis=1)
            self.df['dedup_read_diff_pos'] = self.df.apply(lambda x: self.diffDedupReadCountPos(x, by_col='uniq_repr_nodes'), axis=1)
            self.console.print('======># of deduplicated unique umis {} on the basis of the unique method'.format(self.df['dedup_uniq_diff_pos'].sum()))
            self.console.print('======># of deduplicated reads {} on the basis of the unique method'.format(self.df['dedup_read_diff_pos'].sum()))
            self.gwriter.generic(df=ave_ed_bins, sv_fpn=self.dirname + 'uniq_ave_ed_pos_bin.txt', index=True)
            self.df_dedup_sum = self.df[[
                'ave_eds',
                'uniq_umi_len',
                'dedup_uniq_diff_pos',
                'dedup_read_diff_pos',
            ]]
            self.gwriter.generic(
                df=self.df_dedup_sum,
                sv_fpn=self.dirname + 'uniq_dedup_sum.txt',
                index=True,
                header=True,
            )
            self.console.print('======>start writing deduplicated reads to BAM...')
            dedup_reads_write_stime = time.time()
            self.df['uniq_bam_ids'] = self.df.apply(lambda x: self.bamids(x, by_col='uniq_repr_nodes'), axis=1)
            self.aliwriter.tobam(
                tobam_fpn=self.dirname + self.method + '_dedup.bam',
                tmpl_bam_fpn=self.bam_fpn,
                whitelist=self.decompose(list_nd=self.df['uniq_bam_ids'].values),
            )
            self.console.print('======>finish writing in {:.2f}s'.format(time.time() - dedup_reads_write_stime))

        elif self.method == 'cluster':
            dedup_umi_stime = time.time()
            self.df['cc_repr_nodes'] = self.df.apply(lambda x: self.umimax(x, by_col='cc'), axis=1)
            self.df['cc_umi_len'] = self.df['cc_repr_nodes'].apply(lambda x: self.length(x))
            self.console.print('======>finish finding deduplicated umis in {:.2f}s'.format(time.time() - dedup_umi_stime))
            self.console.print('======># of umis deduplicated to be {}'.format(self.df['cc_umi_len'].loc['yes']))
            self.console.print('======>calculate average edit distances between umis...')
            dedup_umi_edave_stime = time.time()
            self.df['ave_eds'] = self.df.apply(lambda x: self.edave(x, by_col='cc_repr_nodes'), axis=1)
            self.console.print('======>finish calculating ave eds in {:.2f}s'.format(time.time() - dedup_umi_edave_stime))
            self.df['dedup_uniq_diff_pos'] = self.df.apply(lambda x: self.diffDedupUniqCountPos(x, by_col='cc_repr_nodes'), axis=1)
            self.df['dedup_read_diff_pos'] = self.df.apply(lambda x: self.diffDedupReadCountPos(x, by_col='cc_repr_nodes'), axis=1)
            self.console.print('======># of deduplicated unique umis {} on the basis of the unique method'.format(self.df['dedup_uniq_diff_pos'].sum()))
            self.console.print('======># of deduplicated reads {} on the basis of the unique method'.format(self.df['dedup_read_diff_pos'].sum()))
            ave_ed_bins = self.df['ave_eds'].value_counts().sort_index()
            self.console.check(ave_ed_bins)
            self.gwriter.generic(df=ave_ed_bins, sv_fpn=self.dirname + 'cc_ave_ed_pos_bin.txt', index=True)
            self.df_dedup_sum = self.df[[
                'cc_umi_len',
                'ave_eds',
                'uniq_umi_len',
                'dedup_uniq_diff_pos',
                'dedup_read_diff_pos',
            ]]
            self.gwriter.generic(
                df=self.df_dedup_sum,
                sv_fpn=self.dirname + 'cc_dedup_sum.txt',
                index=True,
                header=True,
            )
            self.console.print('======>start writing deduplicated reads to BAM...')
            dedup_reads_write_stime = time.time()
            self.df['cc_bam_ids'] = self.df.apply(lambda x: self.bamids(x, by_col='cc_repr_nodes'), axis=1)
            self.aliwriter.tobam(
                tobam_fpn=self.dirname + self.method + '_dedup.bam',
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
            self.console.print('======># of umis deduplicated to be {}'.format(self.df['adj_umi_len'].loc['yes']))
            self.console.print('======>calculate average edit distances between umis...')
            dedup_umi_edave_stime = time.time()
            self.df['ave_eds'] = self.df.apply(lambda x: self.edave(x, by_col='adj_repr_nodes'), axis=1)
            self.console.print('======>finish calculating ave eds in {:.2f}s'.format(time.time() - dedup_umi_edave_stime))
            self.df['dedup_uniq_diff_pos'] = self.df.apply(lambda x: self.diffDedupUniqCountPos(x, by_col='adj_repr_nodes'), axis=1)
            self.df['dedup_read_diff_pos'] = self.df.apply(lambda x: self.diffDedupReadCountPos(x, by_col='adj_repr_nodes'), axis=1)
            self.console.print('======># of deduplicated unique umis {} on the basis of the unique method'.format(self.df['dedup_uniq_diff_pos'].sum()))
            self.console.print('======># of deduplicated reads {} on the basis of the unique method'.format(self.df['dedup_read_diff_pos'].sum()))
            ave_ed_bins = self.df['ave_eds'].value_counts().sort_index()
            self.console.check(ave_ed_bins)
            self.gwriter.generic(df=ave_ed_bins, sv_fpn=self.dirname + 'adj_ave_ed_pos_bin.txt', index=True)
            self.df_dedup_sum = self.df[[
                'adj_umi_len',
                'ave_eds',
                'uniq_umi_len',
                'dedup_uniq_diff_pos',
                'dedup_read_diff_pos',
            ]]
            self.gwriter.generic(
                df=self.df_dedup_sum,
                sv_fpn=self.dirname + 'adj_dedup_sum.txt',
                index=True,
                header=True,
            )
            self.console.print('======>start writing deduplicated reads to BAM...')
            dedup_reads_write_stime = time.time()
            self.df['adj_bam_ids'] = self.df.apply(lambda x: self.bamids(x, by_col='adj_repr_nodes'), axis=1)
            self.aliwriter.tobam(
                tobam_fpn=self.dirname + self.method + '_dedup.bam',
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
            self.console.print('======># of umis deduplicated to be {}'.format(self.df['direc_umi_len'].loc['yes']))
            self.console.print('======>calculate average edit distances between umis...')
            dedup_umi_edave_stime = time.time()
            self.df['ave_eds'] = self.df.apply(lambda x: self.edave(x, by_col='direc_repr_nodes'), axis=1)
            self.console.print('======>finish calculating ave eds in {:.2f}s'.format(time.time() - dedup_umi_edave_stime))
            self.df['dedup_uniq_diff_pos'] = self.df.apply(lambda x: self.diffDedupUniqCountPos(x, by_col='direc_repr_nodes'), axis=1)
            self.df['dedup_read_diff_pos'] = self.df.apply(lambda x: self.diffDedupReadCountPos(x, by_col='direc_repr_nodes'), axis=1)
            self.console.print('======># of deduplicated unique umis {} on the basis of the unique method'.format(self.df['dedup_uniq_diff_pos'].sum()))
            self.console.print('======># of deduplicated reads {} on the basis of the unique method'.format(self.df['dedup_read_diff_pos'].sum()))
            ave_ed_bins = self.df['ave_eds'].value_counts().sort_index()
            self.console.check(ave_ed_bins)
            self.gwriter.generic(df=ave_ed_bins, sv_fpn=self.dirname + 'direc_ave_ed_pos_bin.txt', index=True)
            self.df_dedup_sum = self.df[[
                'direc_umi_len',
                'ave_eds',
                'uniq_umi_len',
                'dedup_uniq_diff_pos',
                'dedup_read_diff_pos',
            ]]
            self.gwriter.generic(
                df=self.df_dedup_sum,
                sv_fpn=self.dirname + 'direc_dedup_sum.txt',
                index=True,
                header=True,
            )
            self.console.print('======>start writing deduplicated reads to BAM...')
            dedup_reads_write_stime = time.time()
            self.df['direc_bam_ids'] = self.df.apply(lambda x: self.bamids(x, by_col='direc_repr_nodes'), axis=1)
            self.aliwriter.tobam(
                tobam_fpn=self.dirname + self.method + '_dedup.bam',
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
            self.console.print('======># of umis deduplicated to be {}'.format(self.df['mcl_umi_len'].loc['yes']))
            self.console.print('======>calculate average edit distances between umis...')
            dedup_umi_edave_stime = time.time()
            self.df['ave_eds'] = self.df.apply(lambda x: self.edave(x, by_col='mcl_repr_nodes'), axis=1)
            self.console.print('======>finish calculating ave eds in {:.2f}s'.format(time.time() - dedup_umi_edave_stime))
            self.df['dedup_uniq_diff_pos'] = self.df.apply(lambda x: self.diffDedupUniqCountPos(x, by_col='mcl_repr_nodes'), axis=1)
            self.df['dedup_read_diff_pos'] = self.df.apply(lambda x: self.diffDedupReadCountPos(x, by_col='mcl_repr_nodes'), axis=1)
            self.console.print('======># of deduplicated unique umis {} on the basis of the unique method'.format(self.df['dedup_uniq_diff_pos'].sum()))
            self.console.print('======># of deduplicated reads {} on the basis of the unique method'.format(self.df['dedup_read_diff_pos'].sum()))
            ave_ed_bins = self.df['ave_eds'].value_counts().sort_index()
            self.console.check(ave_ed_bins)
            self.gwriter.generic(df=ave_ed_bins, sv_fpn=self.dirname + 'mcl_ave_ed_pos_bin.txt', index=True)
            self.df_dedup_sum = self.df[[
                'mcl_umi_len',
                'ave_eds',
                'uniq_umi_len',
                'dedup_uniq_diff_pos',
                'dedup_read_diff_pos',
            ]]
            self.gwriter.generic(
                df=self.df_dedup_sum,
                sv_fpn=self.dirname + 'mcl_dedup_sum.txt',
                index=True,
                header=True,
            )
            self.console.print('======>start writing deduplicated reads to BAM...')
            dedup_reads_write_stime = time.time()
            self.df['mcl_bam_ids'] = self.df.apply(lambda x: self.bamids(x, by_col='mcl_repr_nodes'), axis=1)
            self.aliwriter.tobam(
                tobam_fpn=self.dirname + self.method + '_dedup.bam',
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
            self.console.print('======># of umis deduplicated to be {}'.format(self.df['mcl_val_umi_len'].loc['yes']))
            self.console.print('======>calculate average edit distances between umis...')
            dedup_umi_edave_stime = time.time()
            self.df['ave_eds'] = self.df.apply(lambda x: self.edave(x, by_col='mcl_val_repr_nodes'), axis=1)
            self.console.print('======>finish calculating ave eds in {:.2f}s'.format(time.time() - dedup_umi_edave_stime))
            self.df['dedup_uniq_diff_pos'] = self.df.apply(lambda x: self.diffDedupUniqCountPos(x, by_col='mcl_val_repr_nodes'), axis=1)
            self.df['dedup_read_diff_pos'] = self.df.apply(lambda x: self.diffDedupReadCountPos(x, by_col='mcl_val_repr_nodes'), axis=1)
            self.console.print('======># of deduplicated unique umis {} on the basis of the unique method'.format(self.df['dedup_uniq_diff_pos'].sum()))
            self.console.print('======># of deduplicated reads {} on the basis of the unique method'.format(self.df['dedup_read_diff_pos'].sum()))
            ave_ed_bins = self.df['ave_eds'].value_counts().sort_index()
            self.console.check(ave_ed_bins)
            self.gwriter.generic(df=ave_ed_bins, sv_fpn=self.dirname + 'mcl_val_ave_ed_pos_bin.txt', index=True)
            self.df_dedup_sum = self.df[[
                'mcl_val_umi_len',
                'ave_eds',
                'uniq_umi_len',
                'dedup_uniq_diff_pos',
                'dedup_read_diff_pos',
            ]]
            self.gwriter.generic(
                df=self.df_dedup_sum,
                sv_fpn=self.dirname + 'mcl_val_dedup_sum.txt',
                index=True,
                header=True,
            )
            self.console.print('======>start writing deduplicated reads to BAM...')
            dedup_reads_write_stime = time.time()
            self.df['mcl_val_bam_ids'] = self.df.apply(lambda x: self.bamids(x, by_col='mcl_val_repr_nodes'), axis=1)
            self.aliwriter.tobam(
                tobam_fpn=self.dirname + self.method + '_dedup.bam',
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
            self.console.print('======># of umis deduplicated to be {}'.format(self.df['mcl_ed_umi_len'].loc['yes']))
            self.console.print('======>calculate average edit distances between umis...')
            dedup_umi_edave_stime = time.time()
            self.df['ave_eds'] = self.df.apply(lambda x: self.edave(x, by_col='mcl_ed_repr_nodes'), axis=1)
            self.console.print('======>finish calculating ave eds in {:.2f}s'.format(time.time() - dedup_umi_edave_stime))
            self.df['dedup_uniq_diff_pos'] = self.df.apply(lambda x: self.diffDedupUniqCountPos(x, by_col='mcl_ed_repr_nodes'), axis=1)
            self.df['dedup_read_diff_pos'] = self.df.apply(lambda x: self.diffDedupReadCountPos(x, by_col='mcl_ed_repr_nodes'), axis=1)
            self.console.print('======># of deduplicated unique umis {} on the basis of the unique method'.format(self.df['dedup_uniq_diff_pos'].sum()))
            self.console.print('======># of deduplicated reads {} on the basis of the unique method'.format(self.df['dedup_read_diff_pos'].sum()))
            ave_ed_bins = self.df['ave_eds'].value_counts().sort_index()
            self.console.check(ave_ed_bins)
            self.gwriter.generic(df=ave_ed_bins, sv_fpn=self.dirname + 'mcl_ed_ave_ed_pos_bin.txt', index=True)
            self.df_dedup_sum = self.df[[
                'mcl_ed_umi_len',
                'ave_eds',
                'uniq_umi_len',
                'dedup_uniq_diff_pos',
                'dedup_read_diff_pos',
            ]]
            self.gwriter.generic(
                df=self.df_dedup_sum,
                sv_fpn=self.dirname + 'mcl_ed_dedup_sum.txt',
                index=True,
                header=True,
            )
            self.console.print('======>start writing deduplicated reads to BAM...')
            dedup_reads_write_stime = time.time()
            self.df['mcl_ed_bam_ids'] = self.df.apply(lambda x: self.bamids(x, by_col='mcl_ed_repr_nodes'), axis=1)
            self.aliwriter.tobam(
                tobam_fpn=self.dirname + self.method + '_dedup.bam',
                tmpl_bam_fpn=self.bam_fpn,
                whitelist=self.decompose(list_nd=self.df['mcl_ed_bam_ids'].values),
            )
            self.console.print('======>finish writing in {:.2f}s'.format(time.time() - dedup_reads_write_stime))

        # sys.stdout.close()

    def diffDedupUniqCountPos(self, df_row, by_col):
        """

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

    def diffDedupReadCountPos(self, df_row, by_col):
        """

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
        diff_nodes = set(df_row['uniq_repr_nodes']) - set(df_row[by_col])
        if diff_nodes != set():
            # print(diff_nodes)
            umi_val_cnt_dict = df_row['vignette']['df_umi_uniq_val_cnt'].to_dict()
            # print(umi_val_cnt_dict)
            return sum(umi_val_cnt_dict[node] for node in diff_nodes)
        else:
            return 0

    def length(self, df_val):
        """

        Parameters
        ----------
        df_val
            list - a python list

        Returns
        -------
            int - the length of the list

        """
        return len(df_val)

    def markSingleUMI(self, df_val):
        if len(df_val) == 1:
            return 'yes'
        else:
            return 'no'

    def correct(self, umi):
        vernier = [i for i in range(36) if i % 3 == 0]
        umi_trimers = [umi[v: v+3] for v in vernier]
        # umi_trimers = textwrap.wrap(umi, 3)
        t = []
        for umi_trimer in umi_trimers:
            s = set(umi_trimer)
            if len(s) == 3:
                rand_index = self.rannum.uniform(low=0, high=3, num=1, use_seed=False)[0]
                t.append(umi_trimer[rand_index])
            elif len(s) == 2:
                sdict = {umi_trimer.count(i): i for i in s}
                t.append(sdict[2])
            else:
                t.append(umi_trimer[0])
        return ''.join(t)

    def decompose(self, list_nd):
        """

        Parameters
        ----------
        x

        Returns
        -------

        """
        list_md = []
        for i in list_nd:
            list_md = list_md + i
        self.console.print('======># of the total reads left after deduplication: {}'.format(len(list_md)))
        return list_md

    def bamids(self, df_row, by_col):
        """"""
        bam_id_maps = df_row['vignette']['umi_bam_ids']
        list_1d = df_row[by_col]
        return [bam_id_maps[node] for node in list_1d]

    def umimax(self, df_row, by_col):
        umi_val_cnts = df_row['vignette']['df_umi_uniq_val_cnt']
        umi_cc = []
        for k_c, nodes in df_row[by_col].items():
            # self.console.print('cc: ', x['cc'])
            # self.console.print('vc: ', umi_val_cnts)
            # self.console.print('nodes: ',nodes)
            # self.console.print('val_cnts: ', umi_val_cnts.loc[umi_val_cnts.index.isin(nodes)].max())
            umi_max = umi_val_cnts.loc[umi_val_cnts.index.isin(nodes)].idxmax()
            umi_cc.append(umi_max)
            # self.console.print('val_cnts1: ',)
        return umi_cc

    def edave(self, df_row, by_col):
        repr_nodes = df_row[by_col]
        umi_maps = df_row['vignette']['int_to_umi_dict']
        node_len = len(repr_nodes)
        if node_len != 1:
            ed_list = []
            for i in range(node_len):
                for j in range(i + 1, node_len):
                    ed_list.append(hamming().general(
                        s1=umi_maps[repr_nodes[i]],
                        s2=umi_maps[repr_nodes[j]],
                    ))
            return np.ceil(sum(ed_list) / (len(ed_list)))
        else:
            return -1

    def eds_(self, df_row, by_col):
        """"""
        print(df_row.index)
        repr_nodes = df_row[by_col]
        umi_maps = df_row['vignette']['int_to_umi_dict']
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
                            umi_maps[repr_nodes[i]],
                            umi_maps[repr_nodes[j]])
                        ] * (umi_val_cnts.loc[repr_nodes[i]] * umi_val_cnts.loc[repr_nodes[j]])
            return round(sum(ed_list) / len(ed_list))
        else:
            return -1

    def evaluate(self, ):
        return


if __name__ == "__main__":
    from umiche.path import to

    umiche = OnePos(
        # method='unique',
        # method='cluster',
        # method='adjacency',
        method='directional',
        # method='mcl',
        # method='mcl_val',
        # method='mcl_ed',

        bam_fpn=to('data/example.bam'),
        # bam_fpn=to('data/example_bundle.bam'),
        mcl_fold_thres=1.5,
        inflat_val=1.6,
        exp_val=2,
        iter_num=100,
        ed_thres=1,
        sv_fpn=to('data/assigned_sorted_dedup.bam'),

        verbose=True,
    )