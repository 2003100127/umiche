__version__ = "v1.0"
__copyright__ = "Copyright 2024"
__license__ = "MIT"
__developer__ = "Jianfeng Sun"
__maintainer__ = "Jianfeng Sun"
__email__="jianfeng.sunmt@gmail.com"
__lab__ = "Cribbslab"

import numpy as np
import pandas as pd
from sklearn.cluster import DBSCAN as skdbscan
from sklearn.cluster import Birch
from umiche.util.Single import Single as dnasgl


class DBSCAN:

    def __init__(self, ):
        self.dnasgl = dnasgl()
        self.nt_to_int_dict = self.dnasgl.todict(nucleotides=self.dnasgl.get(universal=True), reverse=False)
        self.int_to_nt_dict = self.dnasgl.todict(nucleotides=self.dnasgl.get(universal=True), reverse=True)

    def dfclusters(
            self,
            connected_components,
            df_umi_uniq_val_cnt,
            int_to_umi_dict,
    ):
        df_ccs = pd.DataFrame({'cc_vertices': [*connected_components.values()]})
        print(df_ccs)
        ### @@ df_ccs
        #              cc_vertices
        # 0  [A, B, C, D, E, F, G]
        df_ccs['cc_max_id'] = df_ccs['cc_vertices'].apply(lambda cc: self.maxid(df_umi_uniq_val_cnt, cc))
        ### @@ df_ccs['cc_max_id']
        # 0    A
        # Name: cc_max_id, dtype: object
        df_ccs['cc_max_id2seq'] = df_ccs['cc_max_id'].apply(lambda x: int_to_umi_dict[x])
        print(df_ccs['cc_max_id2seq'])
        # 0    AGATCTCGCA
        # Name: cc_max_id2seq, dtype: object
        vertex_onehot = df_ccs['cc_max_id2seq'].apply(lambda x: self.onehot(x)).values.tolist()
        df_vertex_onehot = pd.DataFrame(vertex_onehot)
        print(df_vertex_onehot)
        # d = skdbscan(eps=2.5, min_samples=1).fit(df_vertex_onehot)
        d =  Birch(threshold=1.8, n_clusters=None).fit(df_vertex_onehot)
        asd = np.unique(d.labels_)
        asdas = np.array(d.labels_)
        labels = d.labels_
        print(d.labels_)
        print(asd)
        print(asdas)
        n_clusters_ = len(set(labels)) - (1 if -1 in labels else 0)
        n_noise_ = list(labels).count(-1)
        print(n_clusters_)
        print(n_noise_)
        # return len(asd)
        return len(asd) + len(asdas[asdas == -1]), len(asdas[asdas == -1])

    def maxid(
            self,
            df_umi_uniq_val_cnt,
            cc,
    ):
        """

        Parameters
        ----------
        df_umi_uniq_val_cnt
        cc

        Returns
        -------

        """
        xx = [*df_umi_uniq_val_cnt.loc[df_umi_uniq_val_cnt.index.isin(cc)].sort_values(ascending=False).to_dict().keys()]
        return xx[0]

    def dusters(self, df_umis):
        """
        # df_umis = pd.DataFrame.from_dict(df_umis, orient='index', columns=['umi'])

        Parameters
        ----------
        df_umis

        Returns
        -------

        """
        arr_umis = df_umis['umi'].apply(lambda x: self.onehot(x)).values.tolist()
        d = skdbscan(eps=1.6, min_samples=1).fit(arr_umis)
        # d =  Birch(threshold=1.8, n_clusters=None).fit(cccc)
        labels = d.labels_
        n_clusters_ = len(set(labels)) - (1 if -1 in labels else 0)
        return n_clusters_

    def onehot(self, x):
        umi_ltr2digi = [self.nt_to_int_dict[i] for i in list(x)]
        ids_first_pos = [i*4 for i in range(len(x))]
        ids_to_be_one = [i+j for i,j in zip(umi_ltr2digi, ids_first_pos)]
        # print(ids_to_be_one)
        one_hot = np.zeros(len(x)*4)
        one_hot[ids_to_be_one] = 1
        return one_hot.astype(int)


if __name__ == "__main__":
    from umiche.deduplicate.method.Cluster import Cluster as umimonoclust

    p = DBSCAN()

    ### @@ data from UMI-tools
    # graph_adj = {
    #     'A': ['B', 'C', 'D'],
    #     'B': ['A', 'C'],
    #     'C': ['A', 'B'],
    #     'D': ['A', 'E', 'F'],
    #     'E': ['D'],
    #     'F': ['D'],
    # }
    # print("An adjacency list of a graph:\n{}".format(graph_adj))
    #
    # node_val_sorted = pd.Series({
    #     'A': 456,
    #     'E': 90,
    #     'D': 72,
    #     'B': 2,
    #     'C': 2,
    #     'F': 1,
    # })
    # print("Counts sorted:\n{}".format(node_val_sorted))

    ### @@ data from mine
    graph_adj = {
        'A': ['B', 'C', 'D'],
        'B': ['A', 'C'],
        'C': ['A', 'B'],
        'D': ['A',],
        'E': [ 'G'],
        'F': [ 'G'],
        'G': ['E', 'F'],
    }
    print("An adjacency list of a graph:\n{}".format(graph_adj))

    node_val_sorted = pd.Series({
        'A': 120,
        'D': 90,
        'E': 50,
        'G': 5,
        'B': 2,
        'C': 2,
        'F': 1,
    })
    print("Counts sorted:\n{}".format(node_val_sorted))

    int_to_umi_dict = {
        'A': 'AGATCTCGCA',
        'B': 'AGATCCCGCA',
        'C': 'AGATCACGCA',
        'D': 'AGATCGCGCA',
        'E': 'AGATCGCGGA',
        'F': 'AGATCGCGTA',
        'G': 'TGATCGCGAA',
    }

    ccs = umimonoclust().cc(graph_adj=graph_adj)
    print("Connected components:\n{}".format(ccs))

    res = p.dfclusters(
        connected_components=ccs,
        df_umi_uniq_val_cnt=node_val_sorted,
        int_to_umi_dict=int_to_umi_dict,
    )
    print(res)

    # df_uniq_umi = pd.DataFrame.from_dict(int_to_umi_dict, orient='index', columns=['umi'])
    # print(df_uniq_umi)
    # res = p.dusters(df_uniq_umi)
    # print(res)