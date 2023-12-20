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

from umiche.deduplicate.Build import Build as umibuild
from umiche.deduplicate.Gadgetry import Gadgetry as umigadgetry
from umiche.deduplicate.Tabulate import Tabulate as umitab

# dedup methods
from umiche.deduplicate.method.Cluster import Cluster as umiclust

from umiche.util.Number import number as rannum
from umiche.util.Writer import writer as gwriter
from umiche.util.Console import Console


class OnePos:

    def __init__(
            self,
            bam_fpn,
            ed_thres,
            mcl_fold_thres=None,
            inflat_val=2.0,
            exp_val=2,
            iter_num=100,
            umi_='_',
            work_dir='./',

            heterogeneity=False,

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
        self.bam_fpn = bam_fpn
        self.ed_thres = ed_thres
        self.mcl_fold_thres = mcl_fold_thres
        self.inflat_val = inflat_val
        self.exp_val = exp_val
        self.iter_num = iter_num
        self.work_dir = work_dir
        self.heterogeneity = heterogeneity
        self.verbose = verbose

        self.umiclust = umiclust()

        self.umibuild = umibuild
        self.umigadgetry = umigadgetry()
        self.rannum = rannum()
        self.gwriter = gwriter()
        self.console = Console()
        self.console.verbose = self.verbose

        # sys.stdout = open(self.work_dir + self.method + '_log.txt', 'w')

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
        # self.console.print('===>start deduplication by the {} method...'.format(self.method))
        # sys.stdout.close()

    def unique(self, ) -> pd.DataFrame:
        return umitab(
            df=self.df,
            df_bam=self.df_bam,
            bam_fpn=self.bam_fpn,
            work_dir=self.work_dir,
            heterogeneity=self.heterogeneity,
            verbose=False,
        ).unique()

    def cluster(self, ) -> pd.DataFrame:
        return umitab(
            df=self.df,
            df_bam=self.df_bam,
            bam_fpn=self.bam_fpn,
            work_dir=self.work_dir,
            heterogeneity=self.heterogeneity,
            verbose=False,
        ).cluster()


    def adjacency(self, ) -> pd.DataFrame:
        return umitab(
            df=self.df,
            df_bam=self.df_bam,
            bam_fpn=self.bam_fpn,
            work_dir=self.work_dir,
            heterogeneity=self.heterogeneity,
            verbose=False,
        ).adjacency()

    def directional(self, ) -> pd.DataFrame:
        return umitab(
            df=self.df,
            df_bam=self.df_bam,
            bam_fpn=self.bam_fpn,
            work_dir=self.work_dir,
            heterogeneity=self.heterogeneity,
            verbose=False,
        ).directional()

    def mcl(self, ) -> pd.DataFrame:
        return umitab(
            df=self.df,
            df_bam=self.df_bam,
            bam_fpn=self.bam_fpn,
            work_dir=self.work_dir,
            heterogeneity=self.heterogeneity,
            verbose=False,
        ).mcl(
            inflat_val=self.inflat_val,
            exp_val=self.exp_val,
            iter_num=self.iter_num,
        )

    def mcl_val(self, ) -> pd.DataFrame:
        return umitab(
            df=self.df,
            df_bam=self.df_bam,
            bam_fpn=self.bam_fpn,
            work_dir=self.work_dir,
            heterogeneity=self.heterogeneity,
            verbose=False,
        ).mcl_val(
            inflat_val=self.inflat_val,
            exp_val=self.exp_val,
            iter_num=self.iter_num,
            mcl_fold_thres=self.mcl_fold_thres,
        )

    def mcl_ed(self, ) -> pd.DataFrame:
        return umitab(
            df=self.df,
            df_bam=self.df_bam,
            bam_fpn=self.bam_fpn,
            work_dir=self.work_dir,
            heterogeneity=self.heterogeneity,
            verbose=False,
        ).mcl_ed(
            inflat_val=self.inflat_val,
            exp_val=self.exp_val,
            iter_num=self.iter_num,
            mcl_fold_thres=self.mcl_fold_thres,
        )


if __name__ == "__main__":
    from umiche.path import to

    umiche = OnePos(
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

    print(umiche.unique())
    print(umiche.cluster())
    print(umiche.adjacency())
    print(umiche.directional())
    print(umiche.mcl())
    print(umiche.mcl_val())
    print(umiche.mcl_ed())
