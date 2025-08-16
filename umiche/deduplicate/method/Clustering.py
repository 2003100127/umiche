__version__ = "v1.0"
__copyright__ = "Copyright 2025"
__license__ = "GPL-3.0"
__developer__ = "Jianfeng Sun"
__maintainer__ = "Jianfeng Sun"
__email__="jianfeng.sunmt@gmail.com"

from typing import Dict, List, Tuple

import numpy as np
import pandas as pd
from collections import defaultdict
from tqdm.auto import tqdm
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

        self.clustering_method = clustering_method.lower()
        self.heterogeneity = heterogeneity

        self.kwargs = kwargs

        if self.clustering_method == 'dbscan':
            self.dbscan_eps = self.kwargs.get("dbscan_eps", 1.5)
            self.dbscan_min_spl = self.kwargs.get("dbscan_min_spl", 3)
        elif self.clustering_method == 'birch':
            self.birch_thres = self.kwargs.get("birch_thres", 0.5)
            self.birch_n_clusters = self.kwargs.get("birch_n_clusters", None)
        elif self.clustering_method == 'aprop':
            self.aprop_preference = self.kwargs.get("aprop_preference", None)
            self.aprop_random_state = self.kwargs.get("aprop_random_state", 0)
        elif self.clustering_method == 'hdbscan':
            self.hdbscan_min_spl = self.kwargs.get("hdbscan_min_spl", 3)
        else:
            raise ValueError(f"Unsupported clustering method: {self.clustering_method}.")

        self.console = Console()
        self.console.verbose = verbose

    @property
    def tool(self, ):
        return {
            'dbscan': skdbscan(
                # eps=self.dbscan_eps,
                # min_samples=self.dbscan_min_spl,
                eps=getattr(self, 'dbscan_eps', 1.5),
                min_samples=getattr(self, 'dbscan_min_spl', 3),
            ),
            'birch': skbirch(
                # threshold=self.birch_thres,
                # n_clusters=self.birch_n_clusters,
                threshold=getattr(self, 'birch_thres', 0.5),
                n_clusters=getattr(self, 'birch_n_clusters', None),
            ),
            'aprop': skaprop(
                # preference=self.aprop_preference,
                # random_state=self.aprop_random_state,
                preference=getattr(self, 'aprop_preference', None),
                random_state=getattr(self, 'aprop_random_state', 0),
            ),
            'hdbscan': skhdbscan(
                # min_samples=self.hdbscan_min_spl,
                min_samples=getattr(self, 'hdbscan_min_spl', 3),
            ),
        }

    @Console.vignette()
    def clustering(
            self,
            connected_components: Dict[int, List[str]],
            graph_adj: Dict[str, List[str]],
            int_to_umi_dict: Dict[str, str],
            df_umi_uniq_val_cnt=None,
    ):
        """
        参照 dfclusters() 的逐 CC 流程（不用 pandas）：
          - 对每个 CC：取节点 -> UMI -> onehot，调用 self.tool[self.clustering_method] 聚类
          - 整理为子簇列表，选代表，生成 assigned 与 clusters
          - apv 为该 CC 子图的边列表（与 dfclusters 一致）
        返回：
          assigned_all, {
            'count': 去重簇数,
            'clusters': { 'cc_0': {'node_<rep>': [rep, ...]}, ... },
            'apv':      { 'cc_0': [[u,v], ...], ... }
          }
        """
        clustering_ins = self.tool[self.clustering_method]

        clusters_all: Dict[str, Dict[str, List[str]]] = {}
        apv_all: Dict[str, List[List[str]]] = {}
        assigned_all: Dict[str, str] = {}

        # 可能是 pandas.Series；统一成 dict 的 getter
        def _count_of(node: str) -> int:
            if df_umi_uniq_val_cnt is None:
                return 0
            try:
                # 支持 Series/Dict 两类
                return int(df_umi_uniq_val_cnt[node])
            except Exception:
                return 0

        cc_items = list(connected_components.items())
        with tqdm(
                total=len(cc_items),
                desc=f"[{self.clustering_method}] CCs",
                unit="cc",
                position=0,
                leave=True,
                dynamic_ncols=True,
        ) as p_cc:
            for i, cc_nodes in cc_items:
                p_cc.set_postfix_str(f"cc={i} nodes={len(cc_nodes)}")
                cc_key = f"cc_{i}"

                cc_sub, apv_sub, assigned_sub = self.clustering_(
                    cc=cc_nodes,
                    graph_adj=graph_adj,
                    int_to_umi_dict=int_to_umi_dict,
                    clustering_ins=clustering_ins,
                    count_getter=_count_of,
                    _tqdm_position=1,
                )

                clusters_all[cc_key] = cc_sub
                apv_all[cc_key] = apv_sub
                assigned_all.update(assigned_sub)

                p_cc.update(1)

        dedup_count = sum(len(v) for v in clusters_all.values())

        return {
            "count": dedup_count,
            "clusters": clusters_all,
            "apv": apv_all,
            "assigned": assigned_all,
        }

    @Console.vignette()
    def clustering_(
            self,
            cc: List[str],
            graph_adj: Dict[str, List[str]],
            int_to_umi_dict: Dict[str, str],
            clustering_ins,
            count_getter=lambda _: 0,
            _tqdm_position=1,
    ) -> Tuple[Dict[str, List[str]], List[List[str]], Dict[str, str]]:
        """
        单个 CC 的核心实现（参照 dfclusters，但去掉 pandas）：
          1) 子图邻接 -> 节点序
          2) 节点 UMI -> onehot 作为特征
          3) 用 self.tool[...] 的聚类器 .fit(onehot) 获取 labels_
          4) labels -> 子簇（label==-1 作为 singleton）
          5) 每簇选代表（count 最大，tie 按字典序），生成 clusters 与 assigned
          6) apv = CC 子图边列表
        """
        # 1) /*** graph_cc_adj ***/
        graph_cc_adj = self.refkit.graph_cc_adj(cc, graph_adj)
        nodes = list(graph_cc_adj.keys())

        # 2) /*** onehot ***/
        onehot_2d = []
        with tqdm(
                total=len(nodes),
                desc="  ↳ UMI → onehot encodings",
                unit="node",
                position=_tqdm_position,
                leave=False,
                dynamic_ncols=True,
        ) as p_enc:
            for n in nodes:
                onehot_2d.append(list(self.refkit.onehot(umi=int_to_umi_dict[n])))
                p_enc.update(1)

        # 3) /*** call clustering method ***/
        tqdm.write(f"  ↳ [{getattr(self, 'clustering_method', '?')}] fitting on {len(nodes)} nodes.")
        labels = clustering_ins.fit(onehot_2d).labels_.tolist()
        tqdm.write(f"  ↳ [{getattr(self, 'clustering_method', '?')}] fit done. final labels={labels}")

        # 4) /*** treat -1 (means noise) as separate cluster, singleton  ***/
        lab2members: Dict[int, List[str]] = defaultdict(list)
        for n, lab in zip(nodes, labels):
            lab2members[lab].append(n)
        subclusters: List[List[str]] = []
        for lab, mems in lab2members.items():
            if lab == -1:
                for n in mems:
                    subclusters.append([n])
            else:
                subclusters.append(mems)

        # 5) repr node + clusters + assigned
        cc_sub: Dict[str, List[str]] = {}
        assigned_sub: Dict[str, str] = {}

        with tqdm(
                total=len(subclusters),
                desc="  ↳ assemble clusters",
                unit="sub",
                position=_tqdm_position,
                leave=False,
                dynamic_ncols=True,
        ) as p_asm:
            for clust in subclusters:
                # repr node：highest count in a cluster，
                rep = sorted(clust, key=lambda n: (-count_getter(n), n))[0]
                key = f"node_{rep}"
                # repr node in the first pos，others in descending order
                members_sorted = sorted(clust, key=lambda x: (x != rep, x))
                cc_sub[key] = members_sorted
                for n in clust:
                    if n != rep:
                        assigned_sub[n] = rep
                p_asm.update(1)

        # 6) apv：子图边列表（与 dfclusters 相同）
        edge_list = self.adj_to_edge_list(graph=graph_cc_adj)
        apv_edges = [list(e) for e in edge_list]

        return cc_sub, apv_edges, assigned_sub

    def adj_to_edge_list(self, graph):
        self.netadj.graph = graph
        return self.netadj.to_edge_list()

    def tovertex(self, x):
        # print(x)
        clustering_cluster_arr = x['clustering_clusters'][0]
        # print("['clustering_clusters'][0]",clustering_cluster_arr)
        cc_vertex_arr = x['cc_vertices']
        # print(cc_vertex_arr)
        uniq_cls_arr = np.unique(clustering_cluster_arr)
        # print(uniq_cls_arr)
        clusters = [[] for _ in range(len(uniq_cls_arr))]
        # print(clusters)
        for i, cls in enumerate(clustering_cluster_arr):
            # print(cls)
            if cls != -1:
                # print('ss', cc_vertex_arr[i])
                clusters[cls].append(cc_vertex_arr[i])
            else:
                clusters.append([cc_vertex_arr[i]])
        clusters = [sublist for sublist in clusters if sublist]
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
        # print(int_to_umi_dict)
        # print(connected_components)
        # print(dcnt)
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
        # print('clustering_clusters:',df_ccs['clustering_clusters'])
        df_ccs['clusters'] = df_ccs.apply(lambda x: self.tovertex(x), axis=1)
        # print(df_ccs['clusters'])
        # print(self.decompose(df_ccs['clusters'].values))
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
    graph_adj = {
        'A': ['B', 'C', 'D'],
        'B': ['A', 'C'],
        'C': ['A', 'B'],
        'D': ['A', 'E', 'F'],
        'E': ['D', 'G'],
        'F': ['D', 'G'],
        'G': ['E', 'F'],
    }
    # graph_adj = {
    #     'A': ['B', 'C', 'D'],
    #     'B': ['A', 'C'],
    #     'C': ['A', 'B'],
    #     'D': ['A'],
    #     'E': ['G'],
    #     'F': ['G'],
    #     'G': ['E', 'F'],
    # }
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

    ccs = umimonoclust().cc(graph_adj=graph_adj)
    print("Connected components:\n{}".format(ccs))

    int_to_umi_dict = {
        'A': 'AGATCTCGCA',
        'B': 'AGATCCCGCA',
        'C': 'AGATCACGCA',
        'D': 'AGATCTGGCA',
        'E': 'AGATCTGGGA',
        'F': 'AGATCTGGCT',
        'G': 'AGATCTGGGT',
    }
    print("int_to_umi_dict:\n{}".format(ccs))

    p = Clustering(
        clustering_method='dbscan',
        dbscan_eps=1.4,
        dbscan_min_spl=1,
        # birch_thres=1.8,
        # birch_n_clusters=None,
        verbose=False, # True False
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

    res = p.clustering(
        connected_components=ccs,
        graph_adj=graph_adj,
        int_to_umi_dict=int_to_umi_dict,
        df_umi_uniq_val_cnt=node_val_sorted,
    )

    print("assigned:", res["assigned"])
    print("count:", res["count"])
    print("clusters:", res['clusters'])
    print("apv:", res["apv"])


    df = p.dfclusters_adj_mat(
        connected_components=ccs,
        graph_adj=graph_adj,
    )
    print(df)
    df_decomposed = p.decompose(list_nd=df['clusters'].values)
    print("deduplicated clusters decomposed:\n{}".format(df_decomposed))

    df = p.dfclusters_cc_all_node_umis(
        graph_adj=graph_adj,
        int_to_umi_dict=int_to_umi_dict,
    )
    print(df)
    df_decomposed = p.decompose(list_nd=df['clusters'].values)
    print("deduplicated clusters decomposed:\n{}".format(df_decomposed))