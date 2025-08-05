__version__ = "v1.0"
__copyright__ = "Copyright 2025"
__license__ = "GPL-3.0"
__developer__ = "Jianfeng Sun"
__maintainer__ = "Jianfeng Sun"
__email__="jianfeng.sunmt@gmail.com"

from typing import List, Dict

from umiche.network.CC import CC as gbfscc
from umiche.deduplicate.method.ReformKit import ReformKit
from umiche.util.Console import Console


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
            dist_func,
            umis: List[str],
            ed_thres: int,
    ) -> Dict[str, str]:
        """
        每个分量内，选最丰度（ties 用字典序）的 UMI 为代表，分量里的其它全部并到代表。
        这实现了图示中的“传递合并”（即到代表的 HD 可 > editham）。
        """
        from collections import Counter, deque
        counts = Counter(umis)
        umis_uniq_ordered = list(counts.keys())
        nbrs = ReformKit().neighbor_graph(
            dist_func=dist_func,
            umis=umis_uniq_ordered,
            ed_thres=ed_thres,
        )
        assigned: Dict[str, str] = {}
        seen = set()
        for u in umis_uniq_ordered:
            if u in seen:
                continue
            comp = []
            q = deque([u])
            seen.add(u)
            while q:
                x = q.popleft()
                comp.append(x)
                for y in nbrs[x]:
                    if y not in seen:
                        seen.add(y)
                        q.append(y)
            # 分量代表：最高丰度，tie 按字典序
            rep = sorted(comp, key=lambda s: (-counts[s], s))[0]
            for w in comp:
                if w != rep:
                    assigned[w] = rep
        return assigned


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
    print(p.cc(graph_adj=graph_adj_mclumi))
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
        'D': 'AGATCTGGCA',  # old 'AGATCGCGCA'
        'E': 'AGATCTGGGA',
        'F': 'AGATCTGGCT',
        'G': 'AGATCTGGGT',
    }
    umis = []
    for i in node_val_sorted.index:
        print(i)
        print(node_val_sorted.loc[i])
        umis += [int_to_umi_dict[i]] * node_val_sorted.loc[i]

    from umiche.util.Hamming import Hamming

    print(p.umicountr(
        dist_func=Hamming().umicountr,
        umis=umis,
        # d=int_to_umi_dict,
        ed_thres=1,
    ))