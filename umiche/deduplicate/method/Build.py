__version__ = "v1.0"
__copyright__ = "Copyright 2023"
__license__ = "MIT"
__developer__ = "Jianfeng Sun"
__lab__ = "cribbslab"

import time
import numpy as np
import pandas as pd
from umiche.path import to
from umiche.util.Number import number as rannum
from umiche.util.Hamming import hamming
from umiche.network.Edge import edge as guuedge
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
        self.hamming = hamming()
        self.guuedge = guuedge()
        self.console = Console()
        self.console.verbose = verbose

        self.char_to_int = {'A': 0, 'C': 1, 'G': 2, 'T': 3}

        umi_keymap_stime = time.time()
        self.uniq_umis = self.df['umi'].unique()
        self.uniq_umi_num = self.uniq_umis.shape[0]
        self.console.print('===># of unique UMIs: {}'.format(self.uniq_umi_num))

        self.umi_to_int_dict = {k: id for id, k in enumerate(self.uniq_umis)}
        # print(self.umi_to_int_dict)
        self.int_to_umi_dict = {id: k for k, id in self.umi_to_int_dict.items()}
        # print(self.int_to_umi_dict)
        self.df_umi_uniq_val_cnt = self.df['umi'].value_counts(ascending=False)
        # print(self.df_umi_uniq_val_cnt)
        df_umi_uniq_val_cnt_indexes = self.df_umi_uniq_val_cnt.index
        self.df_umi_uniq_val_cnt.index = [self.umi_to_int_dict[i] for i in df_umi_uniq_val_cnt_indexes]
        # print(self.df_umi_uniq_val_cnt)
        self.console.print('===>umi keymap time: {:.3f}s'.format(time.time() - umi_keymap_stime))

        self.umi_bam_ids = {}
        for k, v in self.int_to_umi_dict.items():
            self.umi_bam_ids[k] = df.loc[df['umi'].isin([v])]['id'].values[0]

        ed_list_stime = time.time()
        self.ed_list = self.ed_list_()
        # print(self.ed_list)
        self.console.print('===>edit distance list construction time: {:.3f}s'.format(time.time() - ed_list_stime))
        ed_list_stime = time.time()
        self.ed_list = self.calc_ed()
        print(self.ed_list)
        self.df_eds = pd.DataFrame(self.ed_list, columns=['node_1', 'node_2', 'ed'])
        self.df_ed_sel = self.df_eds.loc[self.df_eds['ed'] == ed_thres]
        self.console.print('===>edit distance list construction time: {:.3f}s'.format(time.time() - ed_list_stime))

        edge_list_stime = time.time()
        self.edge_list = self.guuedge.fromdf(self.df_ed_sel, to_tuple=False)
        self.console.print('===>edge list construction time: {:.3f}s'.format(time.time() - edge_list_stime))

        graph_adj_stime = time.time()
        self.graph_adj = {i: [] for i in [*self.umi_to_int_dict.values()]}
        self.guuedge.graph = self.edge_list
        self.graph_adj_edges = self.guuedge.toAdjacencyDict()
        for k, v in self.graph_adj_edges.items():
            self.graph_adj[k] += v
        # from umiche.deduplicate.monomer.Cluster import cluster as umimonoclust
        # cc = umimonoclust().cc(self.graph_adj)
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
            'umi_bam_ids': self.umi_bam_ids,
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

    def pcrnum(self, x):
        """

        Parameters
        ----------
        x

        Returns
        -------

        """
        c = x.split('_')[0].split('-')
        if c[1] == 'init':
            return -1
        else:
            return c[2]