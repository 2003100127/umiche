__version__ = "v1.0"
__copyright__ = "Copyright 2025"
__license__ = "GPL-3.0"
__developer__ = "Jianfeng Sun"
__maintainer__ = "Jianfeng Sun"
__email__="jianfeng.sunmt@gmail.com"


from typing import Dict

# import numpy as np
import pandas as pd
import markov_clustering as mc
from umiche.network.CC import CC as gbfscc
from umiche.deduplicate.method.ReformKit import ReformKit as refkit
from umiche.util.Hamming import Hamming
from umiche.network.Adjacency import Adjacency as netadj
from umiche.util.Console import Console


class MarkovClustering:

    def __init__(
            self,
            inflat_val,
            exp_val,
            iter_num,
            verbose=True,
            **kwargs,
    ):
        self.inflat_val = inflat_val
        self.exp_val = exp_val
        self.iter_num = iter_num
        self.kwargs = kwargs
        if 'heterogeneity' in self.kwargs.keys():
            self.heterogeneity = self.kwargs['heterogeneity']
        else:
            self.heterogeneity = False

        self.netadj = netadj()
        self.gbfscc = gbfscc()
        self.refkit = refkit()

        self.console = Console()
        self.console.verbose = verbose

    def adj_to_edge_list(self, graph):
        self.netadj.graph = graph
        return self.netadj.to_edge_list()

    def clustering(
            self,
            connected_components,
            df_umi_uniq_val_cnt,
            graph_adj,
    ):
        """"""
        clusters_all = {}
        apv_all = {}
        assigned_all = {}

        cc_index = 0
        for cc_id, cc_nodes in connected_components.items():
            cc_key = f"cc_{cc_index}"
            cc_index += 1

            cc_sub, apv_sub, assign_sub = self.clustering_(
                df_umi_uniq_val_cnt=df_umi_uniq_val_cnt,
                cc=cc_nodes,
                graph_adj=graph_adj,
            )
            clusters_all[cc_key] = cc_sub
            apv_all[cc_key] = apv_sub
            assigned_all.update(assign_sub)

        dedup_count = sum(len(v) for v in clusters_all.values())

        return {
            "count": dedup_count,
            "clusters": clusters_all,
            "apv": apv_all,
            "assigned": assigned_all,
        }

    def clustering_(
            self,
            df_umi_uniq_val_cnt,
            cc,
            graph_adj,
    ):
        """
        for single one connected component

        Parameters
        ----------
        df_umi_uniq_val_cnt
        cc
        graph_adj

        Returns
        -------
            cc_clusters : {'node_<rep>': [rep, members…]}
            apv_edges   : [[u,v], ...]   (Edge list = approval)
            assigned_sub: raw → representative

        """
        # /*** 1. graph_cc_adj ***/
        graph_cc_adj = self.refkit.graph_cc_adj(cc, graph_adj)
        # key2int = self.refkit.keymap(graph_adj=graph_cc_adj, reverse=False)
        int2key = self.refkit.keymap(graph_adj=graph_cc_adj, reverse=True)

        # /*** 2. cc_adj_mat ***/
        cc_adj_mat = netadj(graph=graph_cc_adj).to_matrix()
        # cc_adj_mat = _graph_to_matrix(graph_cc_adj, key2int)
        mcl_clusters = self.cluster(cc_adj_mat=cc_adj_mat)

        # /*** 3. turn nodes ***/
        clusters_nt = [[int2key[i] for i in tup] for tup in mcl_clusters]

        # /*** 4. clusters dict & assigned ***/
        assigned_sub = {}
        cc_sub = {}
        for clust in clusters_nt:
            # rep = highest read count in this cluster
            rep = sorted(clust, key=lambda n: (-df_umi_uniq_val_cnt[n], n))[0]
            key = f"node_{rep}"
            cc_sub[key] = sorted(clust, key=lambda x: (x != rep, x))
            for n in clust:
                if n != rep:
                    assigned_sub[n] = rep

        # /*** 5. approval edge list ***/
        edge_list = self.adj_to_edge_list(graph=graph_cc_adj)
        apv_edges = [list(e) for e in edge_list]

        return cc_sub, apv_edges, assigned_sub

    # ================================================================
    #  Max-value-based collapsing
    # ================================================================
    def clustering_val(
            self,
            connected_components,
            df_umi_uniq_val_cnt,
            graph_adj,
            thres_fold: float = 2.0,
    ):
        clusters_all, apv_all, dis_all, assigned_all = {}, {}, {}, {}
        cc_idx = 0

        for cc_nodes in connected_components.values():
            # /*** mcl clusters ***/
            graph_cc_adj = self.refkit.graph_cc_adj(cc_nodes, graph_adj)
            int2key = self.refkit.keymap(graph_adj=graph_cc_adj, reverse=True)
            mat = netadj(graph=graph_cc_adj).to_matrix()
            mcl_int = self.cluster(cc_adj_mat=mat)  # list(tuple(int))
            mcl_clusts = [[int2key[i] for i in tup] for tup in mcl_int]

            # /*** turn mcl_clusts to dict format ***/
            mcl_dict = {}
            for sub in mcl_clusts:
                # key of the repr node (highest UMI count)
                rep = max(sub, key=lambda n: (df_umi_uniq_val_cnt[n], -ord(n[0])))
                mcl_dict[f"node_{rep}"] = sub
            # single CC dict
            mcl_cc_dict = {f"cc_{cc_idx}": mcl_dict}

            # /*** max-value merge ***/
            cc_sub, apv, disapv, assign_sub = self.clustering_val_(
                mcl_cc_dict,
                df_umi_uniq_val_cnt,
                thres_fold
            )
            tag = f"cc_{cc_idx}"
            clusters_all[tag] = cc_sub
            apv_all[tag] = apv
            dis_all[tag] = disapv
            assigned_all.update(assign_sub)
            cc_idx += 1

        dedup_cnt = sum(len(v) for v in clusters_all.values())
        return {
            "count": dedup_cnt,
            "clusters": clusters_all,
            "apv": apv_all,
            "disapv": dis_all,
            "assigned": assigned_all,
        }

    # ================================================================
    #
    # ================================================================
    def clustering_val_(
            self,
            mcl_clusters,
            df_umi_uniq_val_cnt,
            thres_fold: float,
    ):
        # /*** parse input ***/
        if isinstance(mcl_clusters, dict):
            sub_clusters = []
            for cc_dict in mcl_clusters.values():
                sub_clusters.extend(cc_dict.values())
        else: # if list
            sub_clusters = list(mcl_clusters)

        # /*** 1. get repr nodes from each sub-cluster ***/
        rep_to_members, reps, weights = {}, [], {}
        for sub in sub_clusters:
            rep = max(sub, key=lambda n: (df_umi_uniq_val_cnt[n], -ord(n[0])))
            reps.append(rep)
            weights[rep] = int(df_umi_uniq_val_cnt[rep])
            rep_to_members[rep] = sub[:]

        # /*** 2. add edge between repr nodes according to fold ***/
        graph = {r: set() for r in reps}
        apv_edges, dis_edges = [], []
        for i, u in enumerate(reps):
            for v in reps[i + 1:]:
                wu, wv = weights[u], weights[v]
                if (wu >= thres_fold * wv - 1) or (wv >= thres_fold * wu - 1):
                    graph[u].add(v);
                    graph[v].add(u)
                    apv_edges.append([u, v])
                else:
                    dis_edges.append([u, v])

        # /*** 3. cc in repr graph ***/
        groups, visited = [], set()
        for r in reps:
            if r in visited: continue
            q, comp = [r], []
            visited.add(r)
            while q:
                x = q.pop()
                comp.append(x)
                for y in graph[x]:
                    if y not in visited:
                        visited.add(y);
                        q.append(y)
            groups.append(comp)

        # /*** 4. clusters & assigned /***
        clusters_cc, assigned_sub = {}, {}
        for g in groups:
            # repr = highest UMI count
            rep_final = max(g, key=lambda n: (weights[n], -ord(n[0])))
            members = []
            for rep in g:
                members.extend(rep_to_members[rep])
            # cluster list
            uniq_members = list(dict.fromkeys(members))
            uniq_members.sort(key=lambda x: (x != rep_final, x))
            clusters_cc[f"node_{rep_final}"] = uniq_members
            for n in uniq_members:
                if n != rep_final:
                    assigned_sub[n] = rep_final

        return clusters_cc, apv_edges, dis_edges, assigned_sub

    # ================================================================
    #  merging mcl clusters based on edit/Hamming dist
    # ================================================================
    def clustering_ed(
            self,
            connected_components,
            df_umi_uniq_val_cnt,
            graph_adj,
            int_to_umi_dict,
            thres_fold: int = 1,
    ):
        """

        Parameters
        ----------
        connected_components
        df_umi_uniq_val_cnt
        graph_adj
        int_to_umi_dict
        thres_fold

        Returns
        -------
            dict:
              {
                "count":    dedup cnt,
                "clusters": { "cc_0": {"node_A":[...], ...}, ... },
                "apv":      { "cc_0": [[A,D],...], ... },
                "disapv":   { "cc_0": [[A,E],...], ... },
                "assigned": { raw → representative, ... },
              }

        """
        clusters_all, apv_all, dis_all, assigned_all = {}, {}, {}, {}
        cc_idx = 0

        for cc_nodes in connected_components.values():
            # ------- 1. Markov clustering -------
            g_cc = self.refkit.graph_cc_adj(cc_nodes, graph_adj)
            int2k = self.refkit.keymap(graph_adj=g_cc, reverse=True)
            mat_cc = netadj(graph=g_cc).to_matrix()
            mcl_int = self.cluster(cc_adj_mat=mat_cc)
            mcl_lis = [[int2k[i] for i in tup] for tup in mcl_int]

            # ------- 2. turn to dict -------
            mcl_dict = {}
            for sub in mcl_lis:
                rep = max(sub, key=lambda n: (df_umi_uniq_val_cnt[n], -ord(n[0])))
                mcl_dict[f"node_{rep}"] = sub
            mcl_cc_dict = {f"cc_{cc_idx}": mcl_dict}

            # ------- 3. Hamming-based merge -------
            cc_sub, apv, dis, assign_sub = self.clustering_ed_(
                mcl_cc_dict,
                df_umi_uniq_val_cnt,
                int_to_umi_dict,
                thres_fold,
            )
            tag = f"cc_{cc_idx}"
            clusters_all[tag] = cc_sub
            apv_all[tag] = apv
            dis_all[tag] = dis
            assigned_all.update(assign_sub)
            cc_idx += 1

        dedup_cnt = sum(len(v) for v in clusters_all.values())
        return {
            "count": dedup_cnt,
            "clusters": clusters_all,
            "apv": apv_all,
            "disapv": dis_all,
            "assigned": assigned_all,
        }

    # ================================================================
    #  Hamming-based representatives merging based on a mcl CC
    # ================================================================
    def clustering_ed_(
            self,
            mcl_clusters,
            df_umi_uniq_val_cnt,
            int_to_umi_dict,
            thres_fold: int,
    ):
        """
        Return
        ------
        cc_clusters : {'node_<rep>': [rep, ...]}
        apv_edges   : [[u,v]...]  (Hamming ≤ thres_fold)
        dis_edges   : [[u,v]...]  (Hamming  > thres_fold)
        assigned    : raw → representative
        """
        # /*** 0. parse ***/
        if isinstance(mcl_clusters, dict):
            sub_clusters = []
            for cc_dict in mcl_clusters.values():
                sub_clusters.extend(cc_dict.values())
        else:
            sub_clusters = list(mcl_clusters)

        # /*** 1. get repr nodes ***/
        rep_to_members, reps, wgt = {}, [], {}
        for sub in sub_clusters:
            rep = max(sub, key=lambda n: (df_umi_uniq_val_cnt[n], -ord(n[0])))
            reps.append(rep)
            wgt[rep] = int(df_umi_uniq_val_cnt[rep])
            rep_to_members[rep] = sub[:]

        # /*** 2. repr nodes based on Hamming ***/
        graph = {r: set() for r in reps}
        apv_edges, dis_edges = [], []
        for i, u in enumerate(reps):
            seq_u = int_to_umi_dict[u]
            for v in reps[i + 1:]:
                seq_v = int_to_umi_dict[v]
                ed = Hamming().general(seq_u, seq_v)
                if ed <= thres_fold:
                    graph[u].add(v);
                    graph[v].add(u)
                    apv_edges.append([u, v])
                else:
                    dis_edges.append([u, v])

        # /*** 3. CCs ***/
        groups, seen = [], set()
        for r in reps:
            if r in seen: continue
            q, comp = [r], []
            seen.add(r)
            while q:
                x = q.pop();
                comp.append(x)
                for y in graph[x]:
                    if y not in seen:
                        seen.add(y);
                        q.append(y)
            groups.append(comp)

        # /*** 4. clusters & assigned ***/
        clusters_cc, assigned = {}, {}
        for g in groups:
            rep_final = max(g, key=lambda n: (wgt[n], -ord(n[0])))
            members = []
            for rep in g:
                members.extend(rep_to_members[rep])
            uniq = list(dict.fromkeys(members))
            uniq.sort(key=lambda x: (x != rep_final, x))
            clusters_cc[f"node_{rep_final}"] = uniq
            for n in uniq:
                if n != rep_final:
                    assigned[n] = rep_final

        return clusters_cc, apv_edges, dis_edges, assigned

    def dfclusters(
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
        # print(connected_components)
        # print([*connected_components.values()])
        df_ccs = pd.DataFrame({'cc_vertices': [*connected_components.values()]})
        # print(df_ccs['cc_vertices'])
        df_ccs['graph_cc_adj'] = df_ccs['cc_vertices'].apply(lambda x: self.refkit.graph_cc_adj(x, graph_adj))
        ### @@ graph_cc_adj
        # {'A': ['B', 'C', 'D'], 'B': ['A', 'C'], 'C': ['A', 'B'], 'D': ['A', 'E', 'F'], 'E': ['D'], 'F': ['D']}
        df_ccs['nt_to_int_map'] = df_ccs['graph_cc_adj'].apply(lambda x: self.refkit.keymap(graph_adj=x, reverse=False))
        df_ccs['int_to_nt_map'] = df_ccs['graph_cc_adj'].apply(lambda x: self.refkit.keymap(graph_adj=x, reverse=True))
        ### @@ nt_to_int_map
        # {'A': 0, 'B': 1, 'C': 2, 'D': 3, 'E': 4, 'F': 5}
        ### @@ int_to_nt_map
        # {0: 'A', 1: 'B', 2: 'C', 3: 'D', 4: 'E', 5: 'F'}
        df_ccs['cc_adj_mat'] = df_ccs.apply(
            # lambda x: self.matrix(
            #     graph_adj=x['graph_cc_adj'],
            #     key_map=x['nt_to_int_map'],
            # ),
            lambda x: netadj(graph=x['graph_cc_adj']).to_matrix(),
            axis=1,
        )
        ### @@ cc_adj_mat
        # [[0. 1. 1. 1. 0. 0.]
        #  [1. 0. 1. 0. 0. 0.]
        #  [1. 1. 0. 0. 0. 0.]
        #  [1. 0. 0. 0. 1. 1.]
        #  [0. 0. 0. 1. 0. 0.]
        #  [0. 0. 0. 1. 0. 0.]]
        df_ccs['mcl_clusters'] = df_ccs['cc_adj_mat'].apply(lambda x: self.cluster(x))
        df_ccs['clusters'] = df_ccs.apply(lambda x: self.refkit.key2node(list_2d=x['mcl_clusters'], keymap=x['int_to_nt_map']), axis=1)
        # print(df_ccs['clusters'].values.tolist())
        ### @@ mcl_clusters
        # [(0, 1, 2), (3, 4, 5)]
        ### @@ clusters
        # [['A', 'B', 'C'], ['D', 'E', 'F']]
        df_ccs['clust_num'] = df_ccs['clusters'].apply(lambda x: len(x))
        # df_ccs['clust_num'] = df_ccs['clusters'].apply(lambda x: print(x) if len(x) > 1 else 1)
        ### @@ clust_num
        # 2
        df_ccs['edge_list'] = df_ccs['graph_cc_adj'].apply(lambda graph: self.adj_to_edge_list(graph=graph))
        # print(df_ccs['edge_list'])
        df_ccs['apv'] = df_ccs['edge_list'].apply(lambda edge_list: [list(el) for el in edge_list])
        # print(df_ccs['apv'])
        # df_ccs.apply(lambda x: self.refkit.breakpoint(x, connected_components) if x['clust_num'] > 1 else 1, axis=1)
        return df_ccs

    def cluster(
            self,
            cc_adj_mat,
    ):
        """

        Parameters
        ----------
        cc_adj_mat
            [[0. 1. 1. 1. 0. 0.]
             [1. 0. 1. 0. 0. 0.]
             [1. 1. 0. 0. 0. 0.]
             [1. 0. 0. 0. 1. 1.]
             [0. 0. 0. 1. 0. 0.]
             [0. 0. 0. 1. 0. 0.]]

             # for {'A': ['B', 'C', 'D'], 'B': ['A', 'C'], 'C': ['A', 'B'], 'D': ['A', 'E', 'F'], 'E': ['D'], 'F': ['D']}


        Returns
        -------

        """
        result = mc.run_mcl(
            cc_adj_mat,
            inflation=self.inflat_val,
            expansion=self.exp_val,
            iterations=int(self.iter_num),
        )
        clusters = mc.get_clusters(result)
        # print(clusters)
        return clusters

    def maxval_val(
            self,
            df_mcl_ccs,
            df_umi_uniq_val_cnt,
            thres_fold,
    ):
        """

        Parameters
        ----------
        df_mcl_ccs
        df_umi_uniq_val_cnt
        thres_fold

        Returns
        -------

        """
        # print(df_umi_uniq_val_cnt)
        df_mcl_ccs['mscmv_val'] = df_mcl_ccs['clusters'].apply(
            lambda x: self.maxval_val_(
                mcl_clusters_per_cc=x,
                df_umi_uniq_val_cnt=df_umi_uniq_val_cnt,
                thres_fold=thres_fold,
            )
        )
        df_mcl_ccs['mscmv_val_len'] = df_mcl_ccs['mscmv_val'].apply(lambda x: x[0])
        df_mcl_ccs['mscmv_val_clusters'] = df_mcl_ccs['mscmv_val'].apply(lambda x: x[1])
        df_mcl_ccs['mscmv_val_apv'] = df_mcl_ccs['mscmv_val'].apply(lambda x: x[2])
        df_mcl_ccs['mscmv_val_disapv'] = df_mcl_ccs['mscmv_val'].apply(lambda x: x[3])
        # print(df_mcl_ccs['mscmv_val_len'].sum())
        # print(df_mcl_ccs[['mscmv_val_clusters', 'mscmv_val_disapv', ]])
        # print(df_mcl_ccs['mscmv_val_disapv'],)
        # return {
        #     'count': df_mcl_ccs['mscmv_val_len'],
        #     'clusters': df_mcl_ccs['mscmv_val_clusters'],
        #     'apv': df_mcl_ccs['mscmv_val_apv'],
        #     'disapv': df_mcl_ccs['mscmv_val_disapv'],
        # }
        if self.heterogeneity:
            return (
                df_mcl_ccs['mscmv_val_len'],
                df_mcl_ccs['mscmv_val_clusters'],
                df_mcl_ccs['mscmv_val_apv'],
                df_mcl_ccs['mscmv_val_disapv'],
            )
        else:
            return {
                "count":  df_mcl_ccs['mscmv_val_len'],
                "clusters": df_mcl_ccs['mscmv_val_clusters'],
                "apv": df_mcl_ccs['mscmv_val_apv'],
                "disapv": df_mcl_ccs['mscmv_val_disapv'],
            }

    def maxval_val_(
            self,
            mcl_clusters_per_cc,
            df_umi_uniq_val_cnt,
            thres_fold,
    ):
        """

        Parameters
        ----------
        mcl_clusters_per_cc
        df_umi_uniq_val_cnt
        thres_fold

        Returns
        -------

        """
        mcl_sub_clust_max_val_graph = {}
        mcl_sub_clust_max_val_weights = {}
        for clust in mcl_clusters_per_cc:
            cc_clust_sorted = self.refkit.sort_vals(df_umi_uniq_val_cnt, cc=clust)
            nodes = [*cc_clust_sorted.keys()]
            weights = [*cc_clust_sorted.values()]
            mcl_sub_clust_max_val_graph[nodes[0]] = set()
            mcl_sub_clust_max_val_weights[nodes[0]] = weights[0]
        # print(mcl_sub_clust_max_val_graph)
        approval = []
        disapproval = []
        mscmv_nodes = [*mcl_sub_clust_max_val_weights.keys()]
        mscmv_weights = [*mcl_sub_clust_max_val_weights.values()]
        mscmv_len = len(mscmv_nodes)
        for i in range(mscmv_len):
            for j in range(i + 1, mscmv_len):
                node_i = mscmv_nodes[i]
                node_j = mscmv_nodes[j]
                node_weight_i = mscmv_weights[i]
                node_weight_j = mscmv_weights[j]
                if node_weight_i > thres_fold * node_weight_j - 1:
                    mcl_sub_clust_max_val_graph[node_i].add(node_j)
                    mcl_sub_clust_max_val_graph[node_j].add(node_i)
                    approval.append([node_i, node_j])
                elif node_weight_j > thres_fold * node_weight_i - 1:
                    mcl_sub_clust_max_val_graph[node_i].add(node_j)
                    mcl_sub_clust_max_val_graph[node_j].add(node_i)
                    approval.append([node_i, node_j])
                else:
                    disapproval.append([node_i, node_j])
        # print(mcl_sub_clust_max_val_graph)
        clusters = list(self.gbfscc.deque(mcl_sub_clust_max_val_graph))
        return len(clusters), clusters, approval, disapproval

    def maxval_ed(
            self,
            df_mcl_ccs,
            df_umi_uniq_val_cnt,
            int_to_umi_dict,
            thres_fold,
    ):
        """

        Parameters
        ----------
        df_mcl_ccs
        df_umi_uniq_val_cnt
        int_to_umi_dict
        thres_fold

        Returns
        -------

        """
        df_mcl_ccs['mscmv_ed'] = df_mcl_ccs['clusters'].apply(
            lambda x: self.maxval_ed_(
                mcl_clusters_per_cc=x,
                df_umi_uniq_val_cnt=df_umi_uniq_val_cnt,
                int_to_umi_dict=int_to_umi_dict,
                thres_fold=thres_fold,
            )
        )
        df_mcl_ccs['mscmv_ed_len'] = df_mcl_ccs['mscmv_ed'].apply(lambda x: x[0])
        df_mcl_ccs['mscmv_ed_clusters'] = df_mcl_ccs['mscmv_ed'].apply(lambda x: x[1])
        df_mcl_ccs['mscmv_ed_apv'] = df_mcl_ccs['mscmv_ed'].apply(lambda x: x[2])
        df_mcl_ccs['mscmv_ed_disapv'] = df_mcl_ccs['mscmv_ed'].apply(lambda x: x[3])
        # print(df_mcl_ccs['mscmv_ed_len'].sum())
        # print(df_mcl_ccs[['mscmv_ed_clusters', 'mscmv_ed_disapv', ]])
        # print(df_mcl_ccs['mscmv_ed_clusters'])
        # print('mscmv_ed_apv', df_mcl_ccs['mscmv_ed_apv'])
        # print('mscmv_ed_disapv', df_mcl_ccs['mscmv_ed_disapv'])
        if self.heterogeneity:
            return (
                df_mcl_ccs['mscmv_ed_len'],
                df_mcl_ccs['mscmv_ed_clusters'],
                df_mcl_ccs['mscmv_ed_apv'],
                df_mcl_ccs['mscmv_ed_disapv']
            )
        else:
            return {
                'count': df_mcl_ccs['mscmv_ed_len'],
                'clusters': df_mcl_ccs['mscmv_ed_clusters'],
                'apv': df_mcl_ccs['mscmv_ed_apv'],
                'disapv': df_mcl_ccs['mscmv_ed_disapv'],
            }

    def maxval_ed_(
            self,
            mcl_clusters_per_cc,
            df_umi_uniq_val_cnt,
            int_to_umi_dict,
            thres_fold,
    ):
        """
        # for k1, v1 in mcl_sub_clust_max_val_weights.items():
        #     for k2, v2 in mcl_sub_clust_max_val_weights.items():
        #         if k1 != k2:
        #             edh = Hamming().general(
        #                 int_to_umi_dict[k1],
        #                 int_to_umi_dict[k2],
        #             )
        #             if edh <= thres_fold:
        #                 mcl_sub_clust_max_val_graph[k1].add(k2)
        #                 mcl_sub_clust_max_val_graph[k2].add(k1)
        #                 approval.append([k1, k2])
        #             else:
        #                 disapproval.append([k1, k2])

        Parameters
        ----------
        mcl_clusters_per_cc
            [['A', 'B', 'C'], ['D', 'E', 'F', 'G']]
        df_umi_uniq_val_cnt
            A    120
            D     90
            E     50
            G      5
            B      2
            C      2
            F      1
            dtype: int64
        int_to_umi_dict
            {'A': 'AGATCTCGCA', 'B': 'AGATCCCGCA', 'C': 'AGATCACGCA', 'D': 'AGATCGCGCA', 'E': 'AGATCGCGGA', 'F': 'AGATCGCGTA', 'G': 'TGATCGCGAA'}
        thres_fold
            1

        Returns
        -------

        """
        ### @@ mcl_clusters_per_cc
        # [['A', 'B', 'C'], ['D', 'E', 'F', 'G']]
        ### @@ df_umi_uniq_val_cnt
        # A    120
        # D     90
        # E     50
        # G      5
        # B      2
        # C      2
        # F      1
        # dtype: int64
        ### @@ int_to_umi_dict
        # {'A': 'AGATCTCGCA', 'B': 'AGATCCCGCA', 'C': 'AGATCACGCA', 'D': 'AGATCGCGCA', 'E': 'AGATCGCGGA', 'F': 'AGATCGCGTA', 'G': 'TGATCGCGAA'}
        ### @@ thres_fold
        # 1
        mcl_sub_clust_max_val_graph = {}
        mcl_sub_clust_max_val_weights = {}
        for clust in mcl_clusters_per_cc:
            cc_clust_sorted = self.refkit.sort_vals(df_umi_uniq_val_cnt, cc=clust)
            nodes = [*cc_clust_sorted.keys()]
            weights = [*cc_clust_sorted.values()]
            mcl_sub_clust_max_val_graph[nodes[0]] = set()
            mcl_sub_clust_max_val_weights[nodes[0]] = weights[0]
        # print(mcl_sub_clust_max_val_graph)
        approval = []
        disapproval = []
        mscmv_nodes = [*mcl_sub_clust_max_val_graph.keys()]
        ### @@ mscmv_nodes
        # ['A', 'D']
        mscmv_len = len(mscmv_nodes)
        for i in range(mscmv_len):
            for j in range(i + 1, mscmv_len):
                node_i = mscmv_nodes[i]
                node_j = mscmv_nodes[j]
                edh = Hamming().general(
                    int_to_umi_dict[node_i],
                    int_to_umi_dict[node_j],
                )
                if edh <= thres_fold:
                    mcl_sub_clust_max_val_graph[node_i].add(node_j)
                    mcl_sub_clust_max_val_graph[node_j].add(node_i)
                    approval.append([node_i, node_j])
                else:
                    disapproval.append([node_i, node_j])
        ### @@ mcl_sub_clust_max_val_graph
        # {'A': {'D'}, 'D': {'A'}}
        clusters = list(self.gbfscc.deque(mcl_sub_clust_max_val_graph))
        ### @@ clusters
        # [['A', 'D']]
        return len(clusters), clusters, approval, disapproval

    def decompose(
            self,
            list_nd,
    ) -> Dict:
        """

        Parameters
        ----------
        df

        Returns
        -------

        """
        ### @@ list_nd
        # [list([['A', 'B', 'C'], ['D', 'E', 'F', 'G']])]
        list_md = []
        for i in list_nd:
            ### @@ i
            # [['A', 'B', 'C'], ['D', 'E', 'F', 'G']]
            list_md = list_md + i
        ### @@ list_md
        # [['A', 'B', 'C'], ['D', 'E', 'F', 'G']]
        res = {}
        for i, cc_sub_each_mcl in enumerate(list_md):
            res[i] = cc_sub_each_mcl
        ### @@ res
        # {0: ['A', 'B', 'C'], 1: ['D', 'E', 'F', 'G']}
        return res

    def get_full_subcc(
            self,
            ccs_dict : Dict,
            mcl_ccs_dict : Dict,
    ) -> Dict:
        d = {}
        for ccid2, cc in ccs_dict.items():
            d[ccid2] = []
            for i in cc:
                for ccid1, mcl_cc in mcl_ccs_dict.items():
                    if i in mcl_cc:
                        d[ccid2] = d[ccid2] + mcl_cc
        return d

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
        df_ccs['graph_cc_adj'] = df_ccs['cc_vertices'].apply(lambda x: self.refkit.graph_cc_adj(x, graph_adj))
        ### @@ graph_cc_adj
        # {'A': ['B', 'C', 'D'], 'B': ['A', 'C'], 'C': ['A', 'B'], 'D': ['A', 'E', 'F'], 'E': ['D'], 'F': ['D']}
        df_ccs['nt_to_int_map'] = df_ccs['graph_cc_adj'].apply(lambda x: self.refkit.keymap(graph_adj=x, reverse=False))
        df_ccs['int_to_nt_map'] = df_ccs['graph_cc_adj'].apply(lambda x: self.refkit.keymap(graph_adj=x, reverse=True))
        ### @@ nt_to_int_map
        # {'A': 0, 'B': 1, 'C': 2, 'D': 3, 'E': 4, 'F': 5}
        ### @@ int_to_nt_map
        # {0: 'A', 1: 'B', 2: 'C', 3: 'D', 4: 'E', 5: 'F'}
        df_ccs['cc_adj_mat'] = df_ccs.apply(
            # lambda x: self.matrix(
            #     graph_adj=x['graph_cc_adj'],
            #     key_map=x['nt_to_int_map'],
            # ),
            lambda x: netadj(graph=x['graph_cc_adj']).to_matrix(),
            axis=1,
        )
        ### @@ cc_adj_mat
        # [[0. 1. 1. 1. 0. 0.]
        #  [1. 0. 1. 0. 0. 0.]
        #  [1. 1. 0. 0. 0. 0.]
        #  [1. 0. 0. 0. 1. 1.]
        #  [0. 0. 0. 1. 0. 0.]
        #  [0. 0. 0. 1. 0. 0.]]
        df_ccs['mcl_clusters'] = df_ccs['cc_adj_mat'].apply(lambda x: self.cluster(x))
        df_ccs['clusters'] = df_ccs.apply(
            lambda x: self.refkit.key2node(list_2d=x['mcl_clusters'], keymap=x['int_to_nt_map']), axis=1)
        ### @@ mcl_clusters
        # [(0, 1, 2), (3, 4, 5)]
        ### @@ clusters
        # [['A', 'B', 'C'], ['D', 'E', 'F']]
        df_ccs['clust_num'] = df_ccs['clusters'].apply(lambda x: len(x))
        ### @@ clust_num
        # 2

        df_ccs['edge_list'] = df_ccs['graph_cc_adj'].apply(lambda graph: self.adj_to_edge_list(graph=graph))
        # print(df_ccs['edge_list'])
        df_ccs['apv'] = df_ccs['edge_list'].apply(lambda edge_list: [list(el) for el in edge_list])
        # print(df_ccs['apv'])
        return df_ccs


if __name__ == "__main__":
    from umiche.deduplicate.method.Cluster import Cluster as umimonoclust

    p = MarkovClustering(
        inflat_val=1.6,
        exp_val=2,
        iter_num=100,
    )

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

    ## @@@ /*** mcl clustering ***/
    res_mcl = p.clustering(
        connected_components=ccs,
        df_umi_uniq_val_cnt=node_val_sorted,
        graph_adj=graph_adj,
    )
    print("assigned mapping: ", res_mcl["assigned"])
    print("dedup count: ", res_mcl["count"])
    print("clusters: ", res_mcl["clusters"])
    print("approval edges: ", res_mcl["apv"])

    ## @@@ /*** mcl dfclusters ***/
    df = p.dfclusters(
        connected_components=ccs,
        graph_adj=graph_adj,
    )
    print(df)
    # print(df.columns)
    # print("vertices: {}".format(df.loc[0, 'cc_vertices']))
    # print("graph_cc_adj: {}".format(df.loc[0, 'graph_cc_adj']))
    # print("nt_to_int_map: {}".format(df.loc[0, 'nt_to_int_map']))
    # print("int_to_nt_map: {}".format(df.loc[0, 'int_to_nt_map']))
    # print("cc_adj_mat: {}".format(df.loc[0, 'cc_adj_mat']))
    # print("mcl_clusters: {}".format(df.loc[0, 'mcl_clusters']))
    # print("clusters: {}".format(df.loc[0, 'clusters']))
    # print("clust_num: {}".format(df.loc[0, 'clust_num']))
    import numpy as np
    print(df['clusters'].tolist())
    # df_mcl_decomposed = p.decompose(list_nd=df['clusters'].values)
    # print("deduplicated clusters decomposed (mcl):\n{}".format(df_mcl_decomposed))

    ## @@@ /*** mcl_val mcl_val ***/
    # df_mcl_val = p.maxval_val(
    #     df_mcl_ccs=df,
    #     df_umi_uniq_val_cnt=node_val_sorted,
    #     thres_fold=2,
    # )
    # print(df_mcl_val)
    # dedup_count = df_mcl_val['count']
    # dedup_clusters = df_mcl_val['clusters']
    # print("deduplicated count (mcl_val):\n{}".format(dedup_count))
    # print("deduplicated clusters (mcl_val):\n{}".format(dedup_clusters))
    #
    # df_mcl_val = p.decompose(list_nd=df_mcl_val['clusters'].values)
    # print("deduplicated clusters decomposed (mcl_val):\n{}".format(df_mcl_val))

    ## @@@ /*** mcl_val clustering_val ***/
    # res_val = p.clustering_val(
    #     connected_components=ccs,
    #     df_umi_uniq_val_cnt=node_val_sorted,
    #     graph_adj=graph_adj,
    #     thres_fold=2,
    # )
    # print("assigned (val): ", res_val["assigned"])
    # print("dedup count: ", res_val["count"])
    # print("clusters: ", res_val["clusters"])
    # print("approval edges: ", res_val["apv"])

    ## @@@ /*** mcl_ed mcl_ed ***/
    df_mcl_ed = p.maxval_ed(
        df_mcl_ccs=df,
        df_umi_uniq_val_cnt=node_val_sorted,
        thres_fold=1,
        int_to_umi_dict=int_to_umi_dict,
    )
    dedup_count = df_mcl_ed['count']
    dedup_clusters = df_mcl_ed['clusters']
    print('approval: {}'.format(df_mcl_ed['apv']))

    print("deduplicated count (mcl_ed):\n{}".format(dedup_count))
    print("deduplicated clusters (mcl_ed):\n{}".format(dedup_clusters))

    df_mcl_ed = p.decompose(list_nd=df_mcl_ed['clusters'].values)
    print("deduplicated clusters decomposed (mcl_ed):\n{}".format(df_mcl_ed))


    ## @@@ /*** mcl_ed mcl_ed ***/
    res_ed = p.clustering_ed(
        connected_components=ccs,
        df_umi_uniq_val_cnt=node_val_sorted,
        graph_adj=graph_adj,
        int_to_umi_dict=int_to_umi_dict,
        thres_fold=1,
    )

    print("assigned:", res_ed["assigned"])
    print("dedup count:", res_ed["count"])
    print("clusters:", res_ed["clusters"])
    print("apv edges:", res_ed["apv"])
