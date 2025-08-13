__version__ = "v1.0"
__copyright__ = "Copyright 2025"
__license__ = "GPL-3.0"
__developer__ = "Jianfeng Sun"
__maintainer__ = "Jianfeng Sun"
__email__="jianfeng.sunmt@gmail.com"

from typing import Dict

from collections import defaultdict
from umiche.util.Console import Console


class Adjacency:

    def __init__(
            self,
            verbose=False,
    ):
        self.console = Console()
        self.console.verbose = verbose

    def umi_tools(
            self,
            connected_components,
            df_umi_uniq_val_cnt,
            graph_adj,
    ) -> Dict:
        """
        Examples
        --------
        umi_tools adjacency wrap

        Parameters
        ----------
        connected_components
        df_umi_uniq_val_cnt
        graph_adj

        Returns
        -------

        """
        tcl = []
        cc_subs = {}
        assigned = {}
        for i, cc in connected_components.items():
            self.console.print('======>connected_components: {}'.format(cc))
            step, cc_sub, sub_assign = self.umi_tools_(
                df_umi_uniq_val_cnt=df_umi_uniq_val_cnt,
                cc=cc,
                graph_adj=graph_adj,
            )
            tcl.append(step)
            cc_subs['cc_' + str(i)] = cc_sub
            assigned.update(sub_assign)
        return {
            'count': sum(tcl),
            'clusters': cc_subs,
            'assigned': assigned,
        }

    def umi_tools_(
            self,
            df_umi_uniq_val_cnt,
            cc,
            graph_adj,
    ) -> Dict:
        """
        umi_tools adjacency

        Parameters
        ----------
        df_umi_uniq_val_cnt
            unique umi counts
        cc
            connected_components
        graph_adj
            the adjacency list of a graph

        Returns
        -------

        """
        cc_umi_sorted = df_umi_uniq_val_cnt.loc[df_umi_uniq_val_cnt.index.isin(cc)].sort_values(ascending=False).to_dict()
        ### @@ cc_umi_sorted
        # {'A': 456, 'E': 90, 'D': 72, 'B': 2, 'C': 2, 'F': 1}
        cc_sorted = [*cc_umi_sorted.keys()]
        ### @@ cc_sorted
        # ['A', 'E', 'D', 'B', 'C', 'F']
        visited = set()
        step = 1
        subcomponents = {}
        cc_set = set(cc_sorted)
        while cc_sorted:
            e = cc_sorted.pop(0)
            subcomponents[e] = []
            for node in graph_adj[e]:
                if node not in visited:
                    subcomponents[e].append(node)
                    self.console.print('=========>subcomponents: {}'.format(subcomponents))
            ### @@ e, subcomponents[e]
            # A ['B', 'C', 'D']
            # E []
            # D ['F']
            visited.add(e)
            ### @@ e, visited
            # A {'A'}
            # E {'B', 'A', 'E', 'C', 'D'}
            # D {'B', 'A', 'E', 'C', 'D'}
            visited.update(graph_adj[e])
            ### @@ e, visited
            # A {'B', 'D', 'A', 'C'}
            # E {'B', 'A', 'E', 'C', 'D'}
            # D {'B', 'F', 'A', 'E', 'C', 'D'}
            subcomponents[e] = graph_adj[e]
            ### @@ e, subcomponents[e]
            # A ['B', 'C', 'D']
            # E ['D']
            # D ['A', 'E', 'F']
            self.console.print('======>The ccurent ele popping out: {} {}'.format(e, visited))
            if visited == cc_set:
                # print(step)
                break
            else:
                step += 1
        ### @@ subcomponents
        # {'A': ['B', 'C', 'D'], 'E': ['D'], 'D': ['A', 'E', 'F']}
        vertex = [*subcomponents.keys()]
        ### @@ vertex
        # ['A', 'E', 'D']
        cc_sub = {}
        for k, v in subcomponents.items():
            cc_sub['node_' + str(k)] = [k]
            for i in v:
                if i not in vertex:
                    cc_sub['node_' + str(k)].append(i)

        # ## /*** assigned ***/
        vertex = list(subcomponents.keys())
        cc_sub = {}
        assigned_sub = {}

        for rep, neigh in subcomponents.items():
            key = f'node_{rep}'
            cc_sub[key] = [rep]
            for n in neigh:
                if n not in vertex:
                    cc_sub[key].append(n)
                    assigned_sub[n] = rep

        return step, cc_sub, assigned_sub

    def decompose(self, cc_sub_dict):
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

    def umicountr(
            self,
            graph_adj,
            connected_components,
            df_umi_uniq_val_cnt,
    ) -> Dict:
        """"""
        self.console.print("===>UMIcountR adjacency method is being used.")
        assigned, taken = {}, set()
        umis_uniq_ordered = df_umi_uniq_val_cnt.index.tolist()
        for u in umis_uniq_ordered:
            if u in taken:
                continue
            taken.add(u)
            for v in graph_adj[u]:
                if v in taken:
                    continue
                assigned[v] = u
                taken.add(v)

        clusters_out = {}
        for cc_id, cc_nodes in connected_components.items():
            cc_sub = defaultdict(list)
            for n in cc_nodes:
                rep = assigned.get(n, n)
                cc_sub[f'node_{rep}'].append(n)
            for rep_key, members in cc_sub.items():
                rep_node = rep_key.replace("node_", "")
                members.sort(key=lambda x: (x != rep_node, x))  # rep first
            clusters_out[f'cc_{cc_id}'] = dict(cc_sub)

        dedup_count = sum(len(sub) for sub in clusters_out.values())

        return {
            'count': dedup_count,
            'clusters': clusters_out,
            'assigned': assigned,
        }

    def umicountr_directional(
            self,
            graph_adj,
            connected_components,
            df_umi_uniq_val_cnt,
    ) -> Dict:
        """
        merged when count ≤ 0.5×rep

        Parameters
        ----------
        graph_adj
        connected_components
        df_umi_uniq_val_cnt

        Returns
        -------

        """
        self.console.print("===>UMIcountR adjacency_directional method is being used.")
        umis_cnt_dict = df_umi_uniq_val_cnt.to_dict()  # UMI → read_count
        assigned, taken = {}, set()
        umis_uniq_ordered = df_umi_uniq_val_cnt.index.tolist()
        for u in umis_uniq_ordered:
            if u in taken:
                continue
            taken.add(u)
            cup = umis_cnt_dict[u]
            for v in graph_adj[u]:
                if v in taken:
                    continue
                if umis_cnt_dict[v] <= 0.5 * cup:
                    assigned[v] = u
                    taken.add(v)

        clusters_out = {}
        for cc_id, cc_nodes in connected_components.items():
            cc_sub = defaultdict(list)
            for n in cc_nodes:
                rep = assigned.get(n, n)
                cc_sub[f'node_{rep}'].append(n)
            for key, members in cc_sub.items():
                rep_node = key.replace("node_", "")
                members.sort(key=lambda x: (x != rep_node, x))
            clusters_out[f'cc_{cc_id}'] = dict(cc_sub)

        dedup_count = sum(len(sub) for sub in clusters_out.values())
        return {
            'count': dedup_count,
            'clusters': clusters_out,
            'assigned': assigned,
        }

    def umicountr_singleton(
            self,
            graph_adj,
            connected_components,
            df_umi_uniq_val_cnt,
    ) -> Dict:
        """仅合并 count==1 的直接邻居到代表。"""
        self.console.print("===>UMIcountR adjacency_singleton method is being used.")
        umis_uniq_ordered = df_umi_uniq_val_cnt.index.tolist()
        umis_cnt_dict = df_umi_uniq_val_cnt.to_dict()
        assigned, taken = {}, set()
        for u in umis_uniq_ordered:
            if u in taken:
                continue
            taken.add(u)
            for v in graph_adj[u]:
                if v in taken:
                    continue
                if umis_cnt_dict[v] == 1:
                    assigned[v] = u
                    taken.add(v)

        clusters_out = {}
        for cc_id, cc_nodes in connected_components.items():
            cc_sub = defaultdict(list)
            for n in cc_nodes:
                rep = assigned.get(n, n)
                cc_sub[f'node_{rep}'].append(n)
            for key, members in cc_sub.items():
                rep_node = key.replace("node_", "")
                members.sort(key=lambda x: (x != rep_node, x))
            clusters_out[f'cc_{cc_id}'] = dict(cc_sub)

        dedup_count = sum(len(sub) for sub in clusters_out.values())

        return {
            'count': dedup_count,
            'clusters': clusters_out,
            'assigned': assigned,
        }


if __name__ == "__main__":
    import pandas as pd
    from umiche.deduplicate.method.Cluster import Cluster as umiclust

    p = Adjacency(
        verbose=True,
    )

    # ### @@ data from UMI-tools
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
    dedup_assigned = dedup_res['assigned']
    print("deduplicated count:\n{}".format(dedup_count))
    print("deduplicated clusters:\n{}".format(dedup_clusters))
    print("deduplicated assigned:\n{}".format(dedup_assigned))

    dedup_clusters_dc = p.decompose(dedup_clusters)
    print("deduplicated clusters decomposed:\n{}".format(dedup_clusters_dc))

    # ## /*** umicountr ***/
    int_to_umi_dict = {
        'A': 'AGATCTCGCA',
        'B': 'AGATCCCGCA',
        'C': 'AGATCACGCA',
        'D': 'AGATCTGGCA',
        'E': 'AGATCTGGGA',
        'F': 'AGATCTGGCT',
        'G': 'AGATCTGGGT',
    }

    tt = p.umicountr(
        graph_adj=graph_adj,
        connected_components=ccs,
        df_umi_uniq_val_cnt=node_val_sorted,
    )
    print(tt)

    tt = p.umicountr_directional(
        graph_adj=graph_adj,
        connected_components=ccs,
        df_umi_uniq_val_cnt=node_val_sorted,
    )
    print(tt)

    tt = p.umicountr_singleton(
        graph_adj=graph_adj,
        connected_components=ccs,
        df_umi_uniq_val_cnt=node_val_sorted,
    )
    print(tt)