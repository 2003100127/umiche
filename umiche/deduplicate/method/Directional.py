__version__ = "v1.0"
__copyright__ = "Copyright 2024"
__license__ = "MIT"
__developer__ = "Jianfeng Sun"
__lab__ = "Cribbslab"

import sys
sys.setrecursionlimit(15000000)


class Directional:

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
        df_umi_uniq_val_cnt
        graph_adj

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
                    node_val_sorted=node_val_sorted,
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
            else:
                continue
        # print(disapv_node_nbr)
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
        node_val_sorted
        node_set_remaining
        graph_adj

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

    def formatApvsDisapv(
            self,
            cc_dict,
    ):
        """
        the format of input to the directional method in umi-tools

        Parameters
        ----------
        cc_dict

        Returns
        -------

        """
        cc_lvl_df = pd.Series(cc_dict)
        return cc_lvl_df.apply(lambda x: self.dictTo2d(x))

    def dictTo2d(
            self,
            x,
    ):
        """

        Parameters
        ----------
        x

        Returns
        -------

        """
        node_list_3d = [*x.values()]
        res_2d = []
        for i in node_list_3d:
            res_2d = res_2d + i
        return res_2d

    def formatCCS(
            self,
            cc_dict,
    ):
        """

        Parameters
        ----------
        cc_dict

        Returns
        -------

        """
        cc_lvl_df = pd.Series(cc_dict)
        return cc_lvl_df.apply(lambda x: [*x.values()])

    def decompose(
            self,
            cc_sub_dict,
    ):
        """

        Parameters
        ----------
        cc_sub_dict

        Returns
        -------

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
    from umiche.deduplicate.method.Cluster import cluster as umimonoclust

    p = Directional()

    graph_adj = {
        'A': ['B', 'C', 'D'],
        'B': ['A', 'C'],
        'C': ['A', 'B'],
        'D': ['A', 'E', 'F'],
        'E': ['D'],
        'F': ['D'],
    }
    print("An adjacency list of a graph:\n{}".format(graph_adj))

    node_val_sorted = pd.Series({
        'A': 456,
        'E': 90,
        'D': 72,
        'B': 2,
        'C': 2,
        'F': 1,
    })
    print("Counts sorted:\n{}".format(node_val_sorted))

    ccs = umimonoclust().cc(graph_adj=graph_adj)
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
