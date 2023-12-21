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
from sklearn.cluster import Birch as skbirch

from umiche.deduplicate.method.ReformKit import ReformKit as refkit



class Clustering:

    def __init__(
            self,
            clustering_method='dbscan',
    ):
        self.refkit = refkit()
        self.clustering_method = clustering_method

    def dfclusters(
            self,
            connected_components,
            graph_adj,
            int_to_umi_dict,
    ):
        """

        Parameters
        ----------
        connected_components
            connected components in dict format:
            {
                'cc0': [...] # nodes,
                'cc1': [...],
                'cc2': [...],
                ...
                'ccn': [...],
            }
            e.g.
            {
                0: ['A', 'B', 'C', 'D', 'E', 'F'],
            }
        graph_adj
            the adjacency list of a graph

        Returns
        -------
            a pandas dataframe
            each connected component is decomposed into more connected subcomponents.

        """
        ### @@ graph_cc_adj
        # {'A': ['B', 'C', 'D'], 'B': ['A', 'C'], 'C': ['A', 'B'], 'D': ['A'], 'E': ['G'], 'F': ['G'], 'G': ['E', 'F']}

        # When an adjacency list of a graph is shown as above, we have the output in the following.
        df_ccs = pd.DataFrame({'cc_vertices': [*connected_components.values()]})
        ### @@ df_ccs
        #     cc_vertices
        # 0  [A, B, C, D]
        # 1     [E, G, F]
        df_ccs['umi'] = df_ccs['cc_vertices'].apply(lambda x: [int_to_umi_dict[node] for node in x])
        ### @@ df_ccs['umi']
        # 0    [AGATCTCGCA, AGATCCCGCA, AGATCACGCA, AGATCGCGCA]
        # 1                [AGATCGCGGA, TGATCGCGAA, AGATCGCGTA]
        # Name: umi, dtype: object
        df_ccs['onehot'] = df_ccs['umi'].apply(lambda umi_arr: [self.refkit.onehot(umi=umi) for umi in umi_arr])
        ### @@ df_ccs['onehot']
        # 0    [[1, 0, 0, 0, 0, 0, 0, 1, 1, 0, 0, 0, 0, 1, 0,...
        # 1    [[1, 0, 0, 0, 0, 0, 0, 1, 1, 0, 0, 0, 0, 1, 0,...
        # Name: onehot, dtype: object

        clustering_ins = self.tool[self.clustering_method]
        df_ccs['clustering_clusters'] = df_ccs['onehot'].apply(lambda onehot_2d_arrs: [
            clustering_ins.fit(onehot_2d_arrs).labels_
        ])
        print(df_ccs['mcl_clusters'])

        # df_ccs['graph_cc_adj'] = df_ccs['cc_vertices'].apply(lambda x: self.refkit.graph_cc_adj(x, graph_adj))
        # ### @@ df_ccs['graph_cc_adj']
        # # 0    {'A': ['B', 'C', 'D'], 'B': ['A', 'C'], 'C': [...
        # # 1            {'E': ['G'], 'G': ['E', 'F'], 'F': ['G']}
        # # Name: graph_cc_adj, dtype: object
        # df_ccs['nt_to_int_map'] = df_ccs['graph_cc_adj'].apply(lambda x: self.refkit.keymap(graph_adj=x, reverse=False))
        # df_ccs['int_to_nt_map'] = df_ccs['graph_cc_adj'].apply(lambda x: self.refkit.keymap(graph_adj=x, reverse=True))

        ### @@ nt_to_int_map
        # {'A': 0, 'B': 1, 'C': 2, 'D': 3, 'E': 4, 'F': 5}
        ### @@ int_to_nt_map
        # {0: 'A', 1: 'B', 2: 'C', 3: 'D', 4: 'E', 5: 'F'}
        # df_ccs['mcl_clusters'] = df_ccs['cc_adj_mat'].apply(lambda x: self.refkit.cluster(x))
        # df_ccs['clusters'] = df_ccs.apply(lambda x: self.refkit.key2node(list_2d=x['mcl_clusters'], keymap=x['int_to_nt_map']), axis=1)
        # ### @@ mcl_clusters
        # # [(0, 1, 2), (3, 4, 5)]
        # ### @@ clusters
        # # [['A', 'B', 'C'], ['D', 'E', 'F']]
        # df_ccs['clust_num'] = df_ccs['clusters'].apply(lambda x: len(x))
        ### @@ clust_num
        # 2
        return df_ccs

    @property
    def tool(self, ):
        return {
            'dbscan': skdbscan(eps=1.5, min_samples=1),
            'birch': skbirch(threshold=1.8, n_clusters=None),
        }

    def dfclusters1(
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
        df_ccs['cc_max_id'] = df_ccs['cc_vertices'].apply(lambda cc: self.refkit.maxid(df_umi_uniq_val_cnt, cc))
        ### @@ df_ccs['cc_max_id']
        # 0    A
        # Name: cc_max_id, dtype: object
        df_ccs['cc_max_id2seq'] = df_ccs['cc_max_id'].apply(lambda x: int_to_umi_dict[x])
        print(df_ccs['cc_max_id2seq'])
        # 0    AGATCTCGCA
        # Name: cc_max_id2seq, dtype: object
        vertex_onehot = df_ccs['cc_max_id2seq'].apply(lambda umi: self.refkit.onehot(umi=umi)).values.tolist()
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



    def dusters(self, df_umis):
        """
        # df_umis = pd.DataFrame.from_dict(df_umis, orient='index', columns=['umi'])

        Parameters
        ----------
        df_umis

        Returns
        -------

        """
        arr_umis = df_umis['umi'].apply(lambda umi: self.refkit.onehot(umi=umi)).values.tolist()
        d = skdbscan(eps=1.6, min_samples=1).fit(arr_umis)
        # d =  Birch(threshold=1.8, n_clusters=None).fit(cccc)
        labels = d.labels_
        n_clusters_ = len(set(labels)) - (1 if -1 in labels else 0)
        return n_clusters_


if __name__ == "__main__":
    from umiche.deduplicate.method.Cluster import Cluster as umimonoclust

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
    # graph_adj = {
    #     'A': ['B', 'C', 'D'],
    #     'B': ['A', 'C'],
    #     'C': ['A', 'B'],
    #     'D': ['A', 'E', 'F'],
    #     'E': ['D', 'G'],
    #     'F': ['D', 'G'],
    #     'G': ['E', 'F'],
    # }
    graph_adj = {
        'A': ['B', 'C', 'D'],
        'B': ['A', 'C'],
        'C': ['A', 'B'],
        'D': ['A',],
        'E': ['G'],
        'F': ['G'],
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

    p = Clustering(clustering_method='dbscan')

    res = p.dfclusters(
        connected_components=ccs,
        graph_adj=graph_adj,
        # df_umi_uniq_val_cnt=node_val_sorted,
        int_to_umi_dict=int_to_umi_dict,
    )
    # print(res)

    # df_uniq_umi = pd.DataFrame.from_dict(int_to_umi_dict, orient='index', columns=['umi'])
    # print(df_uniq_umi)
    # res = p.dusters(df_uniq_umi)
    # print(res)