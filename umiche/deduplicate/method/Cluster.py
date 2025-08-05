__version__ = "v1.0"
__copyright__ = "Copyright 2025"
__license__ = "GPL-3.0"
__developer__ = "Jianfeng Sun"
__maintainer__ = "Jianfeng Sun"
__email__="jianfeng.sunmt@gmail.com"

from typing import List, Dict

from umiche.network.CC import CC as gbfscc
from collections import defaultdict
from umiche.util.Console import Console
from collections import deque


class Cluster:

    def __init__(
            self,
            verbose=True,
    ):
        self.console = Console()
        self.console.verbose = verbose

    def cc(
            self,
            graph_adj,
    ):
        """

        Parameters
        ----------
        graph_adj

        Returns
        -------

        """
        connected_components = list(gbfscc().deque(graph_adj))
        return {i: cc for i, cc in enumerate(connected_components)}

    def ccnx(
            self,
            edge_list,
    ):
        """

        Parameters
        ----------
        edge_list

        Returns
        -------

        """
        import networkx as nx
        G = nx.Graph()
        for edge in edge_list:
            G.add_edge(edge[0], edge[1])
        return {i: G.subgraph(cc).nodes() for i, cc in enumerate(nx.connected_components(G))}

    def umicountr(
            self,
            graph_adj,
            connected_components,
            df_umi_uniq_val_cnt,
    ) -> Dict:
        """

        Parameters
        ----------
        graph_adj
        connected_components
        df_umi_uniq_val_cnt

        Returns
        -------
           {
          'count': ,
          'clusters': {
              'cc_0': {
                  'node_A': ['A', 'B', 'C', ...],
                  ...
              },
              'cc_1': {...},
              ...
          }
        }
        """
        umis_cnt_dict = df_umi_uniq_val_cnt.to_dict()
        umis_uniq_ordered = df_umi_uniq_val_cnt.index.tolist()

        assigned, seen = {}, set()
        for u in umis_uniq_ordered:
            if u in seen:
                continue
            comp = []
            q = deque([u])
            seen.add(u)
            while q:
                x = q.popleft()
                comp.append(x)
                for y in graph_adj[x]:
                    if y not in seen:
                        seen.add(y)
                        q.append(y)
            rep = sorted(comp, key=lambda s: (-umis_cnt_dict[s], s))[0]
            for w in comp:
                if w != rep:
                    assigned[w] = rep

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

        dedup_count = sum(len(cc_sub) for cc_sub in clusters_out.values())

        return {
            'count': dedup_count,
            'clusters': clusters_out,
            'assigned': assigned,
        }


if __name__ == "__main__":
    p = Cluster()
    # print(p.cc({0: [1,], 1: [0, 2], 2: [1]}))
    # print(p.cc({0: []}))
    graph_adj_mclumi = {
     'A': ['B', 'C', 'D'],
     'B': ['A', 'C'],
     'C': ['A', 'B'],
     'D': ['A', 'E', 'F'],
     'E': ['D', 'G'],
     'F': ['D', 'G'],
     'G': ['E', 'F'],
    }
    edge_list = [('B', 'A'), ('D', 'A'), ('C', 'B'), ('F', 'D'), ('C', 'A'), ('G', 'F'), ('E', 'D'), ('G', 'E')]

    ccs = p.cc(graph_adj=graph_adj_mclumi)
    print(ccs)
    print(p.ccnx(edge_list=edge_list))

    import pandas as pd
    node_val_sorted = pd.Series({
        'A': 120,
        'D': 90,
        'E': 50,
        'G': 5,
        'B': 2,
        'C': 2,
        'F': 1,
    })
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

    print(p.umicountr(
        graph_adj=graph_adj_mclumi,
        connected_components=ccs,
        df_umi_uniq_val_cnt=node_val_sorted,
    ))