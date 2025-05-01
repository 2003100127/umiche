__version__ = "v1.0"
__copyright__ = "Copyright 2025"
__license__ = "GPL-3.0"
__developer__ = "Jianfeng Sun"
__maintainer__ = "Jianfeng Sun"
__email__="jianfeng.sunmt@gmail.com"


import time
import numpy as np
import pandas as pd
from umiche.network.Edge import Edge as guuedge
from umiche.util.Hamming import Hamming
from umiche.util.Console import Console


class Build:

    def __init__(
            self,
            df,
            ed_thres,
            verbose=False,
    ):
        """

        Parameters
        ----------
        df
        ed_thres
        """
        self.df = df
        # print(df)
        self.hamming = Hamming()
        self.guuedge = guuedge(verbose=False)
        self.console = Console()
        self.console.verbose = verbose

        self.char_to_int = {'A': 0, 'C': 1, 'G': 2, 'T': 3}

        umi_keymap_stime = time.time()
        self.df_umi_uniq = df.drop_duplicates(subset=['umi'], keep='first')
        ### @@ self.df_umi_uniq
        #           id                     query_name  ...        umi  source
        # 0          0   SRR2057595.2985267_ACCGGTTTA  ...  ACCGGTTTA       1
        # 1          1  SRR2057595.13520751_CCAGGTTCT  ...  CCAGGTTCT       1
        # 2          2   SRR2057595.8901432_AGCGGTTAC  ...  AGCGGTTAC       1
        # ...      ...                            ...  ...        ...     ...
        # 20609  20609   SRR2057595.5197759_TAGGTTCCT  ...  TAGGTTCCT       1
        # [1949 rows x 13 columns]
        self.umi_to_bamid_dict = pd.Series(self.df_umi_uniq.index, index=self.df_umi_uniq.umi.values).to_dict()
        ### @@ self.umi_to_bamid_dict
        # {'ACCGGTTTA': 0, 'CCAGGTTCT': 1, 'AGCGGTTAC': 2, ..., 'TAGGTTCCT': 20609}
        self.uniq_umis = self.df_umi_uniq['umi'].values
        # print(self.uniq_umis)
        ### @@ self.uniq_umis
        # ['ACCGGTTTA' 'CCAGGTTCT' 'AGCGGTTAC', ..., 'TAGGTTCCT']
        # self.uniq_umis1 = self.df['umi'].unique()
        self.uniq_umi_num = self.uniq_umis.shape[0]
        self.console.print('=========># of unique UMIs: {}'.format(self.uniq_umi_num))
        self.umi_to_int_dict = {k: id for id, k in enumerate(self.uniq_umis)}
        ### @@ self.umi_to_int_dict
        # {'ACCGGTTTA': 0, 'CCAGGTTCT': 1, 'AGCGGTTAC': 2, ..., 'TAGGTTCCT': 1948}
        self.int_to_umi_dict = {id: k for k, id in self.umi_to_int_dict.items()}
        ### @@ self.int_to_umi_dict
        # {0: 'ACCGGTTTA', 1: 'CCAGGTTCT', 2: 'AGCGGTTAC', ..., 1948: 'TAGGTTCCT'}
        self.df_umi_uniq_val_cnt = self.df['umi'].value_counts(ascending=False)
        ### @@ self.df_umi_uniq_val_cnt
        # umi
        # AGCGGTTAC    55
        # TACGGTTAT    51
        # AATGGTTTC    47
        #              ..
        # TAGGTTCCT     1
        df_umi_uniq_val_cnt_ids = self.df_umi_uniq_val_cnt.index
        ### @@ df_umi_uniq_val_cnt_ids
        # Index(['AGCGGTTAC', 'TACGGTTAT', 'AATGGTTTC',  ... 'TAGGTTCCT'], dtype='object', name='umi', length=1949)
        self.df_umi_uniq_val_cnt.index = [self.umi_to_int_dict[i] for i in df_umi_uniq_val_cnt_ids]
        ### @@ self.df_umi_uniq_val_cnt
        # 2       55
        # 780     51
        # 416     47
        #         ..
        # 1948     1
        # Name: count, Length: 1949, dtype: int64
        self.console.print('=========>umi keymap time: {:.3f}s'.format(time.time() - umi_keymap_stime))

        self.uniq_umi_id_to_bam_id_dict = {}
        # self.umi_bam_ids1 = {}
        for k, v in self.int_to_umi_dict.items():
            self.uniq_umi_id_to_bam_id_dict[k] = self.umi_to_bamid_dict[v]
            ### @@ self.umi_to_bamid_dict[v]
            # 0
            # 1
            # 2
            # ...
            # 20609
            ### @@ slow scheme
            # self.umi_bam_ids1[k] = df.loc[df['umi'].isin([v])]['id'].values[0]
        ed_list_stime = time.time()
        self.ed_list = self.ed_list_()
        ### @@ self.ed_list
        # [[0, 1, 4],
        # [0, 2, 3],
        # [0, 3, 2],
        # ...
        # [1947, 1948, 4]]
        self.console.print('=========>edit distance list construction time: {:.3f}s'.format(time.time() - ed_list_stime))
        # self.ed_list = self.calc_ed()
        # print(self.ed_list)
        self.df_eds = pd.DataFrame(self.ed_list, columns=['node_1', 'node_2', 'ed'])
        if len(self.df_eds['ed']):
            self.ave_ed = np.ceil(sum(self.df_eds['ed']) / (len(self.df_eds['ed'])))
        else:
            self.ave_ed = -1
            ### @@ self.ave_ed
        # 5.0
        ### @@ self.df_eds
        #          node_1  node_2  ed
        # 0             0       1   4
        # 1             0       2   3
        # 2             0       3   2
        # ...         ...     ...  ..
        # 1898325    1947    1948   4
        # [1898326 rows x 3 columns]
        self.df_ed_sel = self.df_eds.loc[self.df_eds['ed'] == ed_thres]
        # print(self.df_ed_sel['ed'].value_counts())

        ### @@ self.df_ed_sel
        #          node_1  node_2  ed
        # 76            0      77   1
        # 80            0      81   1
        # 96            0      97   1
        # ...         ...     ...  ..
        # 1898293    1940    1944   1
        # [14200 rows x 3 columns]
        self.console.print('=========>edit distance list construction time: {:.3f}s'.format(time.time() - ed_list_stime))

        edge_list_stime = time.time()
        self.edge_list = self.guuedge.fromdf(self.df_ed_sel, to_tuple=False)
        ### @@  self.edge_list
        # [[0, 77], [0, 81], [0, 97], ..., [1940, 1944]]
        self.console.print('=========>edge list construction time: {:.3f}s'.format(time.time() - edge_list_stime))

        graph_adj_stime = time.time()
        self.graph_adj = {i: [] for i in [*self.umi_to_int_dict.values()]}
        ### @@ self.graph_adj
        # {0: [], 1: [], 2: [], ..., 1948: []}
        self.guuedge.graph = self.edge_list
        self.graph_adj_edges = self.guuedge.to_adj_dict()
        ### @@ self.graph_adj_edges
        # {0: [77, 81, 97, 153, 205, 228, 239, 433, 518, 562, 602, 728, 791, 839, 930],
        # 1: [111, 328, 348, 472, 494, 779, 789, 927, 945, 946, 1029, 1055, 1210, 1345, 1398],
        # 2: [156, 163, 168, 215, 282, 422, 480, 565, 765, 831, 843, 1022, 1166, 1229, 1474],
        # ...,
        # 1947: [143, 1066, 1080, 1183, 1199, 1308, 1318, 1606, 1635, 1640, 1682, 1751, 1906, 1935],
        # 1948: [55, 151, 241, 609, 1025, 1032, 1060, 1191, 1230, 1307, 1489, 1490, 1713, 1851, 1922]}
        ### @@ len(self.graph_adj_edges)
        # 1949
        for k, v in self.graph_adj_edges.items():
            self.graph_adj[k] += v
        ### @@ self.graph_adj
        # {0: [77, 81, 97, 153, 205, 228, 239, 433, 518, 562, 602, 728, 791, 839, 930],
        # 1: [111, 328, 348, 472, 494, 779, 789, 927, 945, 946, 1029, 1055, 1210, 1345, 1398],
        # 2: [156, 163, 168, 215, 282, 422, 480, 565, 765, 831, 843, 1022, 1166, 1229, 1474],
        # ...,
        # 1948: [55, 151, 241, 609, 1025, 1032, 1060, 1191, 1230, 1307, 1489, 1490, 1713, 1851, 1922]}

        # from umiche.deduplicate.method.Cluster import cluster as umimonoclust
        # cc = umimonoclust().cc(self.graph_adj)
        # print(cc)
        # print(len(cc))
        # fff = False
        # if len(self.int_to_umi_dict) > 200:
        #     print('!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!')
        #     print(self.edge_list)
        #     print(self.graph_adj)
        #     print(self.int_to_umi_dict)
        #     print(self.df_umi_uniq_val_cnt)
        #     fff = True
        # else:
        #     fff = False
            # break
        self.console.print('===>graph adjacency list construction time: {:.3f}s'.format(time.time() - graph_adj_stime))
        self.data_summary = {
            'graph_adj': self.graph_adj,
            'int_to_umi_dict': self.int_to_umi_dict,
            'df_umi_uniq_val_cnt': self.df_umi_uniq_val_cnt,
            'uniq_umi_id_to_bam_id_dict': self.uniq_umi_id_to_bam_id_dict,
            'ave_ed': self.ave_ed if len(self.int_to_umi_dict) != 1 else -1,
            # 'fff': fff,
        }

    def calc_ed(self, ):
        nda_umi_db = np.array([list(umi) for umi in self.uniq_umis])
        nda_umi_db = np.vectorize(self.char_to_int.get)(nda_umi_db)
        df = pd.DataFrame()
        for i in range(self.uniq_umi_num):
            # r_ids = np.arange(i, self.uniq_umi_num)
            # l_map_ids = [self.umi_to_int_dict[self.uniq_umis[i]]] * (self.uniq_umi_num - i)
            # r_map_ids = [self.umi_to_int_dict[self.uniq_umis[r]] for r in r_ids]
            dd = self.np_vectorize(
                # l_map_ids=l_map_ids,
                # r_map_ids=r_map_ids,
                db_ref=nda_umi_db[i:],
                umi=self.uniq_umis[i],
            )
            df = pd.concat([df, dd], axis=0)
        return df

    def np_vectorize(
            self,
            # l_map_ids,
            # r_map_ids,
            db_ref,
            umi,
    ):
        # stime = time.time()
        df = pd.DataFrame()
        si_arr = np.sum(db_ref != [self.char_to_int[s] for s in umi], axis=1)
        # print(si_arr)
        # print("==================>time: {:.5f}s".format(time.time() - stime))
        # df[0] = l_map_ids
        # df[1] = r_map_ids
        df['ed'] = si_arr
        # print(df)
        # print(cnt_dict)
        # with open(self.sv_fp + seq1 + '.json', 'w') as fp:
        #     json.dump(cnt_dict, fp)
        return df

    def ed_list_(self, ):
        eds = []
        for i in range(self.uniq_umi_num):
            for j in range(i + 1, self.uniq_umi_num):
                l = self.uniq_umis[i]
                r = self.uniq_umis[j]
                # if self.umi_to_int_dict[l] == 31:
                #     print(l)
                # if self.umi_to_int_dict[r] == 50:
                #     print(r)
                eds.append([
                    self.umi_to_int_dict[l],
                    self.umi_to_int_dict[r],
                    self.hamming.general(l, r),
                ])
        # print(len(eds))
        return eds