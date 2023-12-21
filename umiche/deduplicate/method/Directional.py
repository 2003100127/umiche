__version__ = "v1.0"
__copyright__ = "Copyright 2024"
__license__ = "MIT"
__developer__ = "Jianfeng Sun"
__maintainer__ = "Jianfeng Sun"
__email__="jianfeng.sunmt@gmail.com"
__lab__ = "Cribbslab"

import sys
import pandas as pd
sys.setrecursionlimit(15000000)


class Directional:

    def __init__(
            self,
            heterogeneity=False,
    ):
        self.heterogeneity = heterogeneity

    def umi_tools(
            self,
            connected_components,
            df_umi_uniq_val_cnt,
            graph_adj,
    ):
        """

        Parameters
        ----------
        connected_components
            {0: ['A', 'B', 'C', 'D', 'E', 'F']}
        df_umi_uniq_val_cnt
            A    456
            E     90
            D     72
            B      2
            C      2
            F      1
            dtype: int64
        graph_adj
            {'A': ['B', 'C', 'D'], 'B': ['A', 'C'], 'C': ['A', 'B'], 'D': ['A', 'E', 'F'], 'E': ['D'], 'F': ['D']}
        Returns
        -------

        """
        cc_sub_cnt = []
        cc_subs = {}
        cc_apvs = {}
        cc_disapvs = {}
        for i, cc in connected_components.items():
            cc_sub, apv_node_nbr, disapv_node_nbr = self.umi_tools_(
                df_umi_uniq_val_cnt=df_umi_uniq_val_cnt,
                cc=cc,
                graph_adj=graph_adj,
            )
            # print(cc_sub)
            cc_sub_cnt.append(len(cc_sub))
            cc_subs['cc_' + str(i)] = cc_sub
            cc_apvs['cc_' + str(i)] = apv_node_nbr
            cc_disapvs['cc_' + str(i)] = disapv_node_nbr
        # print(sum(cc_sub_cnt))
        # print(cc_subs)
        # print(cc_apvs)
        # print(cc_disapvs)
        if self.heterogeneity:
            return (
                sum(cc_sub_cnt),
                cc_subs,
                cc_apvs,
                cc_disapvs,
            )
        else:
            return {
                "count": sum(cc_sub_cnt),
                "clusters": cc_subs,
                "apv": cc_apvs,
                "disapv": cc_disapvs,
            }

    def umi_tools_(
            self,
            df_umi_uniq_val_cnt,
            cc,
            graph_adj,
    ):
        """

        Parameters
        ----------
        df_umi_uniq_val_cnt
            A    456
            E     90
            D     72
            B      2
            C      2
            F      1
            dtype: int64
        cc
            {0: ['A', 'B', 'C', 'D', 'E', 'F']}
        graph_adj
            {'A': ['B', 'C', 'D'], 'B': ['A', 'C'], 'C': ['A', 'B'], 'D': ['A', 'E', 'F'], 'E': ['D'], 'F': ['D']}

        Returns
        -------

        """
        cc_node_sorted = df_umi_uniq_val_cnt.loc[df_umi_uniq_val_cnt.index.isin(cc)].sort_values(ascending=False).to_dict()
        ### @@ cc_umi_sorted
        # {'A': 456, 'E': 90, 'D': 72, 'B': 2, 'C': 2, 'F': 1}
        nodes = [*cc_node_sorted.keys()]
        # print(nodes)
        ### @@ cc_sorted
        # ['A', 'E', 'D', 'B', 'C', 'F']
        node_cp = nodes.copy()
        node_set_remaining = set(node_cp)
        ### @@ node_set_remaining
        # {'C', 'F', 'E', 'B', 'D', 'A'}
        cc_sub = {}
        apv_node_nbr = {}
        disapv_node_nbr = {}
        while nodes:
            e = nodes.pop(0)
            if e in node_set_remaining:
                seen, apv, disapv = self.dfs(
                    node=e,
                    node_val_sorted=cc_node_sorted,
                    node_set_remaining=node_set_remaining,
                    graph_adj=graph_adj,
                )
                ### @@ e, seen
                # A {'C', 'D', 'F', 'A', 'B'}
                # E {'E'}
                cc_sub['node_' + str(e)] = list(seen)
                apv_node_nbr['node_' + str(e)] = apv
                disapv_node_nbr['node_' + str(e)] = disapv
                node_set_remaining = node_set_remaining - seen
                # print('remaining: {}'.format(node_set_remaining))
                # print('disapproval {}'.format(disapv))
                ### @@ print('disapproval {}'.format(disapv))
                # disapproval []
                # disapproval [[183, 103]]
                # disapproval [[131, 4], [131, 147]]
                # ...
                # disapproval [[133, 194]]
                # disapproval []
            else:
                continue
        ### @@ disapv_node_nbr
        # {'node_0': []}
        # {'node_36': [[183, 103]]}
        # {'node_29': [[131, 4], [131, 147]], 'node_4': []}
        # {'node_7': []}
        # {'node_28': [[8, 57]]}
        # ...
        # {'node_59': [[133, 194]]}
        # {'node_63': []}
        return cc_sub, apv_node_nbr, disapv_node_nbr

    def dfs(
            self,
            node,
            node_val_sorted,
            node_set_remaining,
            graph_adj,
    ):
        """

        Parameters
        ----------
        node
            'A'
        node_val_sorted
            A    456
            E     90
            D     72
            B      2
            C      2
            F      1
            dtype: int64
        node_set_remaining
            {'F', 'A', 'D', 'C', 'E', 'B'}
            {'E'}
        graph_adj
            {'A': ['B', 'C', 'D'], 'B': ['A', 'C'], 'C': ['A', 'B'], 'D': ['A', 'E', 'F'], 'E': ['D'], 'F': ['D']}
        Returns
        -------

        """
        visited = set()
        approval = []
        disapproval = []
        g = graph_adj
        def search(node):
            visited.add(node)
            # print(visited)
            for neighbor in g[node]:
                # print(neighbor)
                if neighbor not in visited:
                    if neighbor in node_set_remaining:
                        if node_val_sorted[node] >= 2 * node_val_sorted[neighbor] - 1:
                            approval.append([node, neighbor])
                            search(neighbor)
                        else:
                            disapproval.append([node, neighbor])
        search(node)
        ### @@ approval
        # {'cc_0': {'node_A': [['A', 'B'], ['A', 'C'], ['A', 'D'], ['D', 'F']], 'node_E': []}}
        ### @@ disapproval
        # {'cc_0': {'node_A': [['B', 'C'], ['D', 'E']], 'node_E': []}}
        return visited, approval, disapproval

    def decompose(
            self,
            cc_sub_dict,
    ):
        """

        Parameters
        ----------
        cc_sub_dict
            input: {'cc_0': {'node_A': ['C', 'A', 'D', 'B', 'F'], 'node_E': ['E']}}
        Returns
        -------
            output: {0: ['C', 'A', 'D', 'B', 'F'], 1: ['E']}

        """
        cc_cnt = 0
        ccs = {}
        for k1, v1 in cc_sub_dict.items():
            for k2, v2 in v1.items():
                ccs[cc_cnt] = v2
                cc_cnt += 1
        return ccs


if __name__ == "__main__":
    import pandas as pd
    from umiche.deduplicate.method.Cluster import Cluster as umiclust

    p = Directional()

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

    ccs = umiclust().cc(graph_adj=graph_adj)
    print("Connected components:\n{}".format(ccs))

    dedup_res = p.umi_tools(
        connected_components=ccs,
        df_umi_uniq_val_cnt=node_val_sorted,
        graph_adj=graph_adj
    )
    dedup_count = dedup_res['count']
    dedup_clusters = dedup_res['clusters']
    print("deduplicated count:\n{}".format(dedup_count))
    print("deduplicated clusters:\n{}".format(dedup_clusters))

    dedup_clusters_dc = p.decompose(dedup_clusters)
    print("deduplicated clusters decomposed:\n{}".format(dedup_clusters_dc))

    print(dedup_res['apv'])
    print(dedup_res['disapv'])
    # {'cc_0': {'node_A': [['A', 'B'], ['A', 'C']], 'node_D': [['D', 'F']], 'node_E': [['E', 'G']]}}
    # {'cc_0': {'node_A': [['B', 'C'], ['A', 'D']], 'node_D': [['D', 'E'], ['F', 'G']], 'node_E': []}}