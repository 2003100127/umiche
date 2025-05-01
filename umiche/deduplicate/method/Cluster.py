__version__ = "v1.0"
__copyright__ = "Copyright 2025"
__license__ = "GPL-3.0"
__developer__ = "Jianfeng Sun"
__maintainer__ = "Jianfeng Sun"
__email__="jianfeng.sunmt@gmail.com"


from umiche.network.CC import CC as gbfscc
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
