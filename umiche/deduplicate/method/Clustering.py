__version__ = "v1.0"
__copyright__ = "Copyright 2025"
__license__ = "GPL-3.0"
__developer__ = "Jianfeng Sun"
__maintainer__ = "Jianfeng Sun"
__email__="jianfeng.sunmt@gmail.com"


import numpy as np
import pandas as pd
from sklearn.cluster import HDBSCAN as skhdbscan
from sklearn.cluster import DBSCAN as skdbscan
from sklearn.cluster import Birch as skbirch
from sklearn.cluster import AffinityPropagation as skaprop

from umiche.deduplicate.method.ReformKit import ReformKit as refkit
from umiche.network.Adjacency import Adjacency as netadj
from umiche.util.Console import Console


class Clustering:

    def __init__(
            self,
            clustering_method='dbscan',
            heterogeneity=None,
            verbose=True,
            **kwargs
    ):
        self.refkit = refkit()
        self.netadj = netadj()

        self.clustering_method = clustering_method
        self.heterogeneity = heterogeneity

        self.kwargs = kwargs
        if 'dbscan_eps' in self.kwargs.keys():
            self.dbscan_eps = self.kwargs['dbscan_eps']
        else:
            self.dbscan_eps = None
        if 'dbscan_min_spl' in self.kwargs.keys():
            self.dbscan_min_spl = self.kwargs['dbscan_min_spl']
        else:
            self.dbscan_min_spl = None
        if 'birch_thres' in self.kwargs.keys():
            self.birch_thres = self.kwargs['birch_thres']
        else:
            self.birch_thres = None
        if 'birch_n_clusters' in self.kwargs.keys():
            self.birch_n_clusters = self.kwargs['birch_n_clusters']
        else:
            self.birch_n_clusters = None
        if 'hdbscan_min_spl' in self.kwargs.keys():
            self.hdbscan_min_spl = self.kwargs['hdbscan_min_spl']
        else:
            self.hdbscan_min_spl = None
        if 'aprop_preference' in self.kwargs.keys():
            self.aprop_preference = self.kwargs['aprop_preference']
        else:
            self.aprop_preference = None
        if 'aprop_random_state' in self.kwargs.keys():
            self.aprop_random_state = self.kwargs['aprop_random_state']
        else:
            self.aprop_random_state = None

        self.console = Console()
        self.console.verbose = verbose

    @property
    def tool(self, ):
        return {
            'dbscan': skdbscan(eps=self.dbscan_eps, min_samples=self.dbscan_min_spl),
            'birch': skbirch(threshold=self.birch_thres, n_clusters=None),
            'hdbscan': skhdbscan(min_samples=self.hdbscan_min_spl),
            'aprop': skaprop(preference=self.aprop_preference, random_state=self.aprop_random_state),
        }

    def adj_to_edge_list(self, graph):
        self.netadj.graph = graph
        return self.netadj.to_edge_list()

    def tovertex(self, x):
        clustering_cluster_arr = x['clustering_clusters'][0]
        cc_vertex_arr = x['cc_vertices']
        uniq_cls_arr = np.unique(clustering_cluster_arr)
        clusters = [[] for _ in range(len(uniq_cls_arr))]
        for i, cls in enumerate(clustering_cluster_arr):
            clusters[cls].append(cc_vertex_arr[i])
        return clusters

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
        # print(df_ccs['clustering_clusters'])
        df_ccs['clusters'] = df_ccs.apply(lambda x: self.tovertex(x), axis=1)
        # print(df_ccs['clusters'])
        df_ccs['clust_num'] = df_ccs['clusters'].apply(lambda x: len(x))
        # print(df_ccs['clust_num'])
        df_ccs['graph_cc_adj'] = df_ccs['cc_vertices'].apply(lambda x: self.refkit.graph_cc_adj(x, graph_adj))
        # print(df_ccs['graph_cc_adj'])
        df_ccs['edge_list'] = df_ccs['graph_cc_adj'].apply(lambda graph: self.adj_to_edge_list(graph=graph))
        # print(df_ccs['edge_list'])
        df_ccs['apv'] = df_ccs['edge_list'].apply(lambda edge_list: [list(el) for el in edge_list])
        # print(df_ccs['apv'])
        return df_ccs

    def dfclusters_adj_mat(
            self,
            connected_components,
            graph_adj,
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
        df_ccs['graph_cc_adj'] = df_ccs['cc_vertices'].apply(lambda x: self.refkit.graph_cc_adj(x, graph_adj))
        # print(df_ccs['graph_cc_adj'])
        df_ccs['cc_adj_mat'] = df_ccs.apply(
            # lambda x: self.matrix(
            #     graph_adj=x['graph_cc_adj'],
            #     key_map=x['nt_to_int_map'],
            # ),
            lambda x: netadj(graph=x['graph_cc_adj']).to_matrix(),
            axis=1,
        )
        # print(df_ccs['cc_adj_mat'])
        clustering_ins = self.tool[self.clustering_method]
        df_ccs['clustering_clusters'] = df_ccs['cc_adj_mat'].apply(lambda adj_mat: [
            clustering_ins.fit(adj_mat).labels_
        ])
        # print(df_ccs['clustering_clusters'])
        df_ccs['clusters'] = df_ccs.apply(lambda x: self.tovertex(x), axis=1)
        # print(df_ccs['clusters'])
        df_ccs['clust_num'] = df_ccs['clusters'].apply(lambda x: len(x))
        # print(df_ccs['clust_num'])
        df_ccs['edge_list'] = df_ccs['graph_cc_adj'].apply(lambda graph: self.adj_to_edge_list(graph=graph))
        # print(df_ccs['edge_list'])
        df_ccs['apv'] = df_ccs['edge_list'].apply(lambda edge_list: [list(el) for el in edge_list])
        # print(df_ccs['apv'])
        return df_ccs

    def dfclusters_cc_fuse(
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
        d =  skbirch(threshold=1.8, n_clusters=None).fit(df_vertex_onehot)
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

    def dfclusters_cc_all_node_umis(
            self,
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
        df_ccs = pd.DataFrame()
        df_ccs.loc[0, 'method'] = 'dfclusters_cc_all_node_umis'
        onehot_2d_arrs = [list(self.refkit.onehot(umi=umi)) for k, umi in int_to_umi_dict.items()]
        d1 = {'dfclusters_cc_all_node_umis': onehot_2d_arrs}
        d2 = {'dfclusters_cc_all_node_umis': [*int_to_umi_dict.keys()]}
        df_ccs['onehot'] = df_ccs['method'].apply(lambda x: d1[x])
        df_ccs['cc_vertices'] = df_ccs['method'].apply(lambda x: d2[x])
        # print(df_ccs)
        clustering_ins = self.tool[self.clustering_method]
        df_ccs['clustering_clusters'] = df_ccs['onehot'].apply(lambda onehot_2d_arrs: [
            clustering_ins.fit(onehot_2d_arrs).labels_
        ])
        # print(df_ccs['clustering_clusters'])
        df_ccs['clusters'] = df_ccs.apply(lambda x: self.tovertex(x), axis=1)
        # print(df_ccs['clusters'])
        df_ccs['clust_num'] = df_ccs['clusters'].apply(lambda x: len(x))
        # print(df_ccs['clust_num'])
        df_ccs['graph_cc_adj'] = df_ccs['cc_vertices'].apply(lambda x: self.refkit.graph_cc_adj(x, graph_adj))
        # print(df_ccs['graph_cc_adj'])
        df_ccs['edge_list'] = df_ccs['graph_cc_adj'].apply(lambda graph: self.adj_to_edge_list(graph=graph))
        # print(df_ccs['edge_list'])
        df_ccs['apv'] = df_ccs['edge_list'].apply(lambda edge_list: [list(el) for el in edge_list])
        # print(df_ccs['apv'])
        return df_ccs

    def decompose(
            self,
            list_nd,
    ):
        """

        Parameters
        ----------
        df

        Returns
        -------
        {

        }

        """
        # print(list_nd)
        list_md = []
        for i in list_nd:
            list_md = list_md + i
        res = {}
        for i, cc_sub_each_mcl in enumerate(list_md):
            res[i] = cc_sub_each_mcl
        # print(res)
        return res


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

    p = Clustering(
        clustering_method='dbscan',
        dbscan_eps=1.5,
        dbscan_min_spl=1,
        birch_thres=1.8,
        birch_n_clusters=None,
    )

    df = p.dfclusters(
        connected_components=ccs,
        graph_adj=graph_adj,
        # df_umi_uniq_val_cnt=node_val_sorted,
        int_to_umi_dict=int_to_umi_dict,
    )
    print(df)
    df_decomposed = p.decompose(list_nd=df['clusters'].values)
    print("deduplicated clusters decomposed:\n{}".format(df_decomposed))

    # df = p.dfclusters_adj_mat(
    #     connected_components=ccs,
    #     graph_adj=graph_adj,
    # )
    # print(df)
    # df_decomposed = p.decompose(list_nd=df['clusters'].values)
    # print("deduplicated clusters decomposed:\n{}".format(df_decomposed))

    # df = p.dfclusters_cc_all_node_umis(
    #     graph_adj=graph_adj,
    #     int_to_umi_dict=int_to_umi_dict,
    # )
    # print(df)
    # df_decomposed = p.decompose(list_nd=df['clusters'].values)
    # print("deduplicated clusters decomposed:\n{}".format(df_decomposed))