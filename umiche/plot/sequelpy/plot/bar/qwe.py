__version__ = "v1.0"
__copyright__ = "Copyright 2024"
__license__ = "MIT"
__developer__ = "Jianfeng Sun"
__maintainer__ = "Jianfeng Sun"
__email__="jianfeng.sunmt@gmail.com"
__lab__ = "Cribbslab"

import networkx as nx
import matplotlib.pyplot as plt
from umiche.plot.Element import Element as pele
from umiche.network.Adjacency import Adjacency as netadj


class Graph:

    def __init__(
            self,
            graph,
    ):
        self.pele = pele()
        self.netadj = netadj()
        self.graph = graph
        self.netadj.graph = self.graph
        self.color_list = self.pele.color(which='ccs4', is_random=True)

    def draw(
            self,
    ):

        el = self.netadj.to_edge_list()
        ccs = [['A:120', 'B:2', 'C:2', 'D:90', 'E:50', 'F:1']]  # cc
        ccs = [['A:120', 'B:2', 'C:2'], ['D:90', 'E:50', 'F:1'], ['G:1']]  # adj
        ccs = [['A:120', 'B:2', 'C:2'], ['D:90', 'F:1', ], ['E:50', 'G:1']]  # direc
        # ccs = [['A:120', 'B:2', 'C:2'], ['D:90', 'E:50', 'F:1', 'G:1']] # mcl
        # ccs = [['A:120', 'B:2', 'C:2', 'D:90', 'E:50', 'F:1', 'G:1']] # mcl_ed
        # ccs = [['A:120', 'B:2', 'C:2'], ['D:90', 'E:50', 'F:1', 'G:1']] # mcl_val
        ncs = []
        for i, cc in enumerate(ccs):
            for j in cc:
                ncs.append(self.color_list[i])
        print(ncs)

        # relationships = pd.DataFrame({'from': ['10', '10', '10', '5000', '10000'],
        #                               'to':   ['100', '500', '1000', '500', '500']})

        G = nx.Graph()
        G.add_nodes_from(['A:120', 'B:2', 'C:2', 'D:90', 'E:50', 'F:1', 'G:1'])
        G.add_edges_from(el)
        # explicitly set positions
        pos = {
            'A:120': (0, 0),
            'B:2': (-1, 1),
            'C:2': (1, 1),
            'D:90': (0, -2),
            'E:50': (-1, -4.5),
            'F:1': (1, -4.5),
            'G:1': (0, -6),
        }

        options = {
            # "font_size": 36,
            "font_size": 16,
            "node_size": 3000,
            "node_color": "bisque",
            "edgecolors": ncs,
            # "edgecolors": ['firebrick', 'firebrick', 'firebrick', 'royalblue', 'springgreen', 'royalblue', 'springgreen'],
            "linewidths": 5,
            "width": 5,

        }
        print(G.nodes())
        nx.draw_networkx(G, pos, **options)

        # Set margins for the axes so that nodes aren't clipped
        ac = plt.gca()
        ac.margins(0.20)
        plt.axis("off")
        plt.show()

        # # Create graph object
        # G = nx.Graph(el)

        # # Manually set each nodes size
        # node_sizes = [4560, 20, 20, 720, 900, 10,]

        # # Draw graph
        # nx.draw(G, with_labels=True, node_size=node_sizes)

        # plt.show()


if __name__ == "__main__":
    import pandas as pd

    graph_adj = {
        'A:120': ['B:2', 'C:2', 'D:90'],
        'B:2': ['A:120', 'C:2'],
        'C:2': ['A:120', 'B:2'],
        'D:90': ['A:120', 'E:50', 'F:1'],
        'E:50': ['D:90', 'G:1'],
        'F:1': ['D:90', 'G:1'],
        'G:1': ['E:50', 'F:1'],
    }

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

    p = Graph(graph_adj)

    fig, ax = plt.subplots(figsize=(7, 7))

    print(p.draw())