__version__ = "v1.0"
__copyright__ = "Copyright 2025"
__license__ = "GPL-3.0"
__developer__ = "Jianfeng Sun"
__maintainer__ = "Jianfeng Sun"
__email__="jianfeng.sunmt@gmail.com"


import sys
import pandas as pd

from umiche.deduplicate.Gadgetry import Gadgetry as umigadgetry
from umiche.deduplicate.Tabulate import Tabulate as umitab

from umiche.util.Console import Console
sys.setrecursionlimit(15000000)


class SeqtechSimu:

    def __init__(
            self,
            bam_fpn,
            df_bam,
            ed_thres,
            umi_tag,
            token,
            granul_lvl_list,
            build_method='graph',
            work_dir='./',
            verbose=False,

            **kwargs,
    ):
        self.bam_fpn = bam_fpn
        self.df_bam = df_bam
        self.ed_thres = ed_thres
        self.umi_tag = umi_tag
        self.token = token
        self.granul_lvl_list = granul_lvl_list
        self.work_dir = work_dir
        self.verbose = verbose
        self.kwargs = kwargs

        self.umigadgetry = umigadgetry()

        self.console = Console()
        self.console.verbose = self.verbose

        # sys.stdout = open(self.work_dir + self.method + '_log.txt', 'w')

        self.console.df_column_summary(df=self.df_bam, title='df_bam')

        self.console.print('======># of raw reads: {}'.format(self.df_bam.shape[0]))
        self.console.print('======># of overall unique UMIs: {}'.format(self.df_bam[umi_tag].unique().shape[0]))

        # self.df_bam = self.df_bam.loc[self.df_bam['reference_id'] != -1]
        # self.console.print('======># of reads with qualified chrs: {}'.format(self.df_bam.shape[0]))

        self.df_bam_gp = self.df_bam.groupby(by=self.granul_lvl_list)
        self.pos_gp_keys = self.df_bam_gp.groups.keys()
        self.console.print('======># of conditions given by the granularity level list: {}'.format(len(self.pos_gp_keys)))
        pos_gps = []
        res_sum = []
        if build_method == 'graph':
            for pos_g in self.console._tqdm(
                    self.pos_gp_keys,
                    # total=bam_in.count(until_eof=True),
                    desc="[UMI graph build]",
                    unit="reads",
                    position=0,
                    leave=True,
                    dynamic_ncols=False,
                ):
                # print(self.df_bam_gp.get_group(pos_g))
                from umiche.bam.Build import Graph as umibuildg
                self.umibuildg = umibuildg
                umi_vignette = self.umibuildg(
                    df=self.df_bam_gp.get_group(pos_g),
                    ed_thres=self.ed_thres,
                    umi_tag=self.umi_tag,
                    verbose=False,
                ).data_summary
                # print(umi_vignette)
                from umiche.deduplicate.method.Cluster import Cluster as umiclust
                self.umiclust = umiclust()
                cc = self.umiclust.cc(umi_vignette['graph_adj'])
                pos_gps.append(pos_g)
                res_sum.append([
                    umi_vignette,
                    cc,
                    umi_vignette['ave_ed'],
                    [*umi_vignette['int_to_umi_dict'].keys()],
                ])
                # ### @@ self.df['uniq_repr_nodes'] or [*umi_vignette['int_to_umi_dict'].keys()]
                # # [0, 1, 2, ..., 1948]
            self.df = pd.DataFrame(
                data=res_sum,
                columns=['vignette', 'cc', 'ave_ed', 'uniq_repr_nodes'],
                index=pos_gps,
            )
            # print(self.df)
            self.df['num_uniq_umis'] = self.df['uniq_repr_nodes'].apply(lambda x: self.umigadgetry.length(x))

        if build_method == 'basic':
            # print(1)
            for pos_g in self.console._tqdm(
                    self.pos_gp_keys,
                    # total=bam_in.count(until_eof=True),
                    desc="[UMI graph build]",
                    unit="reads",
                    position=0,
                    leave=True,
                    dynamic_ncols=False,
                ):
                # print(self.df_bam_gp.get_group(pos_g))
                from umiche.bam.Build import Basic as umibuildb
                self.umibuildb = umibuildb
                umi_vignette = self.umibuildb(
                    df=self.df_bam_gp.get_group(pos_g),
                    umi_tag=self.umi_tag,
                    verbose=False,
                ).data_summary
                # print(umi_vignette)
                pos_gps.append(pos_g)
                res_sum.append([
                    umi_vignette,
                    [*umi_vignette['int_to_umi_dict'].keys()],
                ])
                # print(umi_vignette['int_to_umi_dict'])
                # @@ umi_vignette['int_to_umi_dict']
                # {0: 'ACTGAGTG', 1: 'AGTGGACA', 2: 'AAAGGCCC', ..., 2627: 'AACTTGGT', 2628: 'GGTATATG'}
                # print([*umi_vignette['int_to_umi_dict'].keys()])
                # @@ [*umi_vignette['int_to_umi_dict'].keys()]
                # [0, 1, 2, ..., 2628]
                # print(len([*umi_vignette['int_to_umi_dict'].keys()]))
                # @@ len([*umi_vignette['int_to_umi_dict'].keys()])
                # 2629
            self.df = pd.DataFrame(
                data=res_sum,
                columns=['vignette', 'uniq_repr_nodes'],
                index=pos_gps,
            )
            self.console.df_column_summary(df=self.df, title='df')
            self.df['num_uniq_umis'] = self.df['uniq_repr_nodes'].apply(lambda x: self.umigadgetry.length(x))

    def unique(self, ) -> pd.DataFrame:
        return umitab(
            df=self.df,
            df_bam=self.df_bam,
            bam_fpn=self.bam_fpn,
            work_dir=self.work_dir,
            heterogeneity=False,
            verbose=self.verbose,
        ).unique(
            umi_col=self.umi_tag,
            token=self.token,
            granul_lvl_list=self.granul_lvl_list,
        )

    def cluster(self, ) -> pd.DataFrame:
        return umitab(
            df=self.df,
            df_bam=self.df_bam,
            bam_fpn=self.bam_fpn,
            work_dir=self.work_dir,
            heterogeneity=False,
            verbose=self.verbose,
        ).cluster(
            umi_col=self.umi_tag,
            token=self.token,
            granul_lvl_list=self.granul_lvl_list,
        )

    def adjacency(self, ) -> pd.DataFrame:
        return umitab(
            df=self.df,
            df_bam=self.df_bam,
            bam_fpn=self.bam_fpn,
            work_dir=self.work_dir,
            heterogeneity=False,
            verbose=self.verbose,
        ).adjacency(
            umi_col=self.umi_tag,
            token=self.token,
            granul_lvl_list=self.granul_lvl_list,
        )

    def umicountr_adj(self, ) -> pd.DataFrame:
        return umitab(
            df=self.df,
            df_bam=self.df_bam,
            bam_fpn=self.bam_fpn,
            work_dir=self.work_dir,
            heterogeneity=False,
            verbose=self.verbose,
        ).umicountr(
            clustering_method='adj',
            umi_col=self.umi_tag,
            token=self.token,
            granul_lvl_list=self.granul_lvl_list,
        )

    def umicountr_adj_direc(self, ) -> pd.DataFrame:
        return umitab(
            df=self.df,
            df_bam=self.df_bam,
            bam_fpn=self.bam_fpn,
            work_dir=self.work_dir,
            heterogeneity=False,
            verbose=self.verbose,
        ).umicountr(
            clustering_method='adj_direc',
            umi_col=self.umi_tag,
            token=self.token,
            granul_lvl_list=self.granul_lvl_list,
        )

    def umicountr_adj_singleton(self, ) -> pd.DataFrame:
        return umitab(
            df=self.df,
            df_bam=self.df_bam,
            bam_fpn=self.bam_fpn,
            work_dir=self.work_dir,
            heterogeneity=False,
            verbose=self.verbose,
        ).umicountr(
            clustering_method='adj_singleton',
            umi_col=self.umi_tag,
            token=self.token,
            granul_lvl_list=self.granul_lvl_list,
        )

    def directional(self, ) -> pd.DataFrame:
        return umitab(
            df=self.df,
            df_bam=self.df_bam,
            bam_fpn=self.bam_fpn,
            work_dir=self.work_dir,
            heterogeneity=False,
            verbose=self.verbose,
        ).directional(
            umi_col=self.umi_tag,
            token=self.token,
            granul_lvl_list=self.granul_lvl_list,
        )

    def mcl(self, ) -> pd.DataFrame:
        return umitab(
            df=self.df,
            df_bam=self.df_bam,
            bam_fpn=self.bam_fpn,
            work_dir=self.work_dir,
            heterogeneity=False,
            verbose=self.verbose,
        ).mcl(
            umi_col=self.umi_tag,
            token=self.token,
            granul_lvl_list=self.granul_lvl_list,
            **self.kwargs
        )

    def mcl_ed(self, ) -> pd.DataFrame:
        return umitab(
            df=self.df,
            df_bam=self.df_bam,
            bam_fpn=self.bam_fpn,
            work_dir=self.work_dir,
            heterogeneity=False,
            verbose=self.verbose,
        ).mcl_ed(
            umi_col=self.umi_tag,
            token=self.token,
            granul_lvl_list=self.granul_lvl_list,
            **self.kwargs
        )

    def mcl_val(self, ) -> pd.DataFrame:
        return umitab(
            df=self.df,
            df_bam=self.df_bam,
            bam_fpn=self.bam_fpn,
            work_dir=self.work_dir,
            heterogeneity=False,
            verbose=self.verbose,
        ).mcl_val(
            umi_col=self.umi_tag,
            token=self.token,
            granul_lvl_list=self.granul_lvl_list,
            **self.kwargs
        )

    def dbscan_seq_onehot(self, ) -> pd.DataFrame:
        return umitab(
            df=self.df,
            df_bam=self.df_bam,
            bam_fpn=self.bam_fpn,
            work_dir=self.work_dir,
            heterogeneity=False,
            verbose=self.verbose,
        ).clustering_umi_seq_onehot(
            clustering_method='dbscan',
            umi_col=self.umi_tag,
            token=self.token,
            granul_lvl_list=self.granul_lvl_list,
            **self.kwargs
        )

    def birch_seq_onehot(self, ) -> pd.DataFrame:
        return umitab(
            df=self.df,
            df_bam=self.df_bam,
            bam_fpn=self.bam_fpn,
            work_dir=self.work_dir,
            heterogeneity=False,
            verbose=self.verbose,
        ).clustering_umi_seq_onehot(
            clustering_method='birch',
            umi_col=self.umi_tag,
            token=self.token,
            granul_lvl_list=self.granul_lvl_list,
            **self.kwargs
        )

    def hdbscan_seq_onehot(self, ) -> pd.DataFrame:
        return umitab(
            df=self.df,
            df_bam=self.df_bam,
            bam_fpn=self.bam_fpn,
            work_dir=self.work_dir,
            heterogeneity=False,
            verbose=self.verbose,
        ).clustering_umi_seq_onehot(
            clustering_method='hdbscan',
            umi_col=self.umi_tag,
            token=self.token,
            granul_lvl_list=self.granul_lvl_list,
            **self.kwargs
        )

    def aprop_seq_onehot(self, ) -> pd.DataFrame:
        return umitab(
            df=self.df,
            df_bam=self.df_bam,
            bam_fpn=self.bam_fpn,
            work_dir=self.work_dir,
            heterogeneity=False,
            verbose=self.verbose,
        ).clustering_umi_seq_onehot(
            clustering_method='aprop',
            umi_col=self.umi_tag,
            token=self.token,
            granul_lvl_list=self.granul_lvl_list,
            **self.kwargs
        )

    def set_cover(self, ) -> pd.DataFrame:
        return umitab(
            df=self.df,
            df_bam=self.df_bam,
            bam_fpn=self.bam_fpn,
            work_dir=self.work_dir,
            heterogeneity=False,
            verbose=self.verbose,
        ).set_cover(
            umi_col=self.umi_tag,
            token=self.token,
            granul_lvl_list=self.granul_lvl_list,
            **self.kwargs
        )

    def majority_vote(self, ) -> pd.DataFrame:
        return umitab(
            df=self.df,
            df_bam=self.df_bam,
            bam_fpn=self.bam_fpn,
            work_dir=self.work_dir,
            heterogeneity=False,
            verbose=self.verbose,
        ).majority_vote(
            umi_col=self.umi_tag,
            token=self.token,
            granul_lvl_list=self.granul_lvl_list,
            **self.kwargs
        )

    def starsolo(self, ) -> pd.DataFrame:
        return umitab(
            df=self.df,
            df_bam=self.df_bam,
            bam_fpn=self.bam_fpn,
            work_dir=self.work_dir,
            heterogeneity=False,
            verbose=self.verbose,
        ).starsolo(
            umi_col=self.umi_tag,
            token=self.token,
            granul_lvl_list=self.granul_lvl_list,
        )

    def gencore(self, ) -> pd.DataFrame:
        return umitab(
            df=self.df,
            df_bam=self.df_bam,
            bam_fpn=self.bam_fpn,
            work_dir=self.work_dir,
            heterogeneity=False,
            verbose=self.verbose,
        ).gencore(
            umi_col=self.umi_tag,
            token=self.token,
            granul_lvl_list=self.granul_lvl_list,
        )

    def dropest(self, ) -> pd.DataFrame:
        return umitab(
            df=self.df,
            df_bam=self.df_bam,
            bam_fpn=self.bam_fpn,
            work_dir=self.work_dir,
            heterogeneity=False,
            verbose=self.verbose,
        ).dropest(
            umi_col=self.umi_tag,
            token=self.token,
            granul_lvl_list=self.granul_lvl_list,
        )

    def irescue(self, ) -> pd.DataFrame:
        return umitab(
            df=self.df,
            df_bam=self.df_bam,
            bam_fpn=self.bam_fpn,
            work_dir=self.work_dir,
            heterogeneity=False,
            verbose=self.verbose,
        ).irescue(
            umi_col=self.umi_tag,
            token=self.token,
            granul_lvl_list=self.granul_lvl_list,
            max_hd=self.ed_thres,
        )

    def umis(self, ) -> pd.DataFrame:
        return umitab(
            df=self.df,
            df_bam=self.df_bam,
            bam_fpn=self.bam_fpn,
            work_dir=self.work_dir,
            heterogeneity=False,
            verbose=self.verbose,
        ).umis(
            umi_col=self.umi_tag,
            token=self.token,
            granul_lvl_list=self.granul_lvl_list,
        )


if __name__ == "__main__":
    # @@ /*** -------------------- ***/
    # calc homotrimer UMIs
    # @@ /*** -------------------- ***/

    seqtechs = [
        'vasaseq',
        'pipseq',
        'scifiseq',
        'scrbseq',
        'shareseq',
        'snareseq',
        'splitseq',
        'strtseqc1',
        'strtseq2i',
        'quartzseq2',
        'petipseq',
        'marsseq2',
        'pairedseq',
        'issaacseq',
        'indrop',
        'dropseq',
        '10xv2',
        '10xv3',
        'celseq2',
        'flashseq',
    ]

    for seqtech in seqtechs:
        print("seqtech: {}".format(seqtech))
        bam_fpn = "/mnt/d/Document/Programming/Python/umiche/umiche/data/r1/seqtech/simu/" + seqtech + "/" + seqtech + ".bam"
        from umiche.bam.Reader import ReaderChunk
        df_bam = ReaderChunk(
            bam_fpn=bam_fpn,
            bam_fields=['contig', 'pos', 'CIGAR', 'read'],
            tag_whitelist=['CB', 'MB', 'XF'],
            verbose=True,
        ).todf(chunk_size=1_000_000)
        Console(True).df_column_summary(df=df_bam)
        df_bam = df_bam.reset_index(drop=True)
        # df_bam = df_bam.dropna(subset=["XT"]).reset_index(drop=True)
        # Console(True).df_column_summary(df=df_bam)
        # df_bam = df_bam[df_bam["XS"] == "Assigned"]
        # Console(True).df_column_summary(df=df_bam)
        print(df_bam)
        print(df_bam['XF'].value_counts())

        from umiche.deduplicate.method.trimer.Expand import Expand
        df_bam["htUMI"] = df_bam["MB"].apply(Expand().homotrimer_umi)
        print(df_bam["htUMI"])

        # @@ downsampling pct
        df_bam = df_bam.groupby("XF", group_keys=False).sample(frac=0.1, random_state=42)

        # @@ downsampling num
        # df_bam = (
        #     df_bam.groupby("XF", group_keys=False)
        #     .apply(lambda g: g.sample(n=min(len(g), 1000), random_state=42))
        #     .reset_index(drop=True)
        # )

        mut_rates = [
            # 0,
        #     0.0005,
        #     0.001, 0.01,
            0.05,
            0.1,
        #     0.2, 0.3,
        #     # 0.4, 0.5
        ]
        for mut_rate in mut_rates:

            print("===>mut_rate: {}".format(mut_rate))

            from umiche.deduplicate.method.trimer.Error import Error
            df_bam["htUMI_" + str(mut_rate)] = df_bam["htUMI"].apply(lambda umi: Error().mutated(umi, mut_rate=mut_rate, mode="normal"))

            # p = SeqtechSimu(
            #     bam_fpn=bam_fpn,
            #     df_bam=df_bam,
            #     ed_thres=1,
            #     # umi_tag="MB",
            #     umi_tag="htUMI_" + str(mut_rate),
            #     token=str(mut_rate),
            #     # granul_lvl_list=['CB', 'XF'],
            #     granul_lvl_list=['XF'],
            #
            #     build_method='basic',  # basic graph
            #
            #     work_dir='/mnt/d/Document/Programming/Python/umiche/umiche/data/r1/seqtech/simu/' + seqtech + '/ed1_ht_spl_pct10/',
            #
            #     mcl_fold_thres=1.5,
            #     inflat_val=1.6,
            #     exp_val=2,
            #     iter_num=100,
            #
            #     umi_unit_pattern=3,
            #
            #     verbose=True,  # False True
            # )
            # res_sc = p.set_cover()
            # res_mv = p.majority_vote()

            p = SeqtechSimu(
                bam_fpn=bam_fpn,
                df_bam=df_bam,
                ed_thres=1,
                # umi_tag="MB",
                umi_tag="htUMI_" + str(mut_rate),
                token=str(mut_rate),
                granul_lvl_list=['XF'],

                build_method='graph',  # basic graph

                work_dir='/mnt/d/Document/Programming/Python/umiche/umiche/data/r1/seqtech/simu/' + seqtech + '/ed1_ht_spl_pct10/',

                mcl_fold_thres=1.5,
                inflat_val=1.6,
                exp_val=2,
                iter_num=100,

                umi_unit_pattern=1,

                verbose=True,  # False True
            )
            # res_unique = p.unique()
            # res_cluster = p.cluster()
            # res_adjacency = p.adjacency()
            # res_directional = p.directional()
            # res_umicountr_adj = p.umicountr_adj()
            # res_umicountr_adj_direc = p.umicountr_adj_direc()
            # res_umicountr_adj_singleton = p.umicountr_adj_singleton()
            res_mcl = p.mcl()
            res_mcl_ed = p.mcl_ed()
            res_mcl_val = p.mcl_val()
            # res_dbscan_seq_onehot = p.dbscan_seq_onehot()
            # res_aprop_seq_onehot = p.aprop_seq_onehot()
            # res_birch_seq_onehot = p.birch_seq_onehot()
            # res_starsolo = p.starsolo()
            # res_gencore = p.gencore()
            # # res_dropest = p.dropest()
            # res_irescue = p.irescue()
            # res_umis = p.umis()