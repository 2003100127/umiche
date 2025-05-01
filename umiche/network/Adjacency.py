__version__ = "v1.0"
__copyright__ = "Copyright 2025"
__license__ = "GPL-3.0"
__developer__ = "Jianfeng Sun"
__maintainer__ = "Jianfeng Sun"
__email__="jianfeng.sunmt@gmail.com"


import numpy as np
from umiche.util.Console import Console


class Adjacency:

    def __init__(
            self,
            graph=None,
            verbose=True,
    ):
        self._graph = graph
        self._graph_mapped = None

        self.console = Console()
        self.console.verbose = verbose

    @property
    def graph(self, ):
        self.console.print('===>The current graph: {}'.format(self._graph))
        return self._graph

    @graph.setter
    def graph(
            self,
            value,
    ):
        """

        Parameters
        ----------
        value

        Returns
        -------

        """
        self._graph = value

    @property
    def graph_mapped(self, ):
        if self._graph == None:
            raise 'Please set and input your graph'
        else:
            return self.map(self._graph)

    @property
    def glen(self, ):
        # print('the number of egdes in the current graph is {}'.format(len(self._graph)))
        return len(self._graph)

    @property
    def key_mapped(self, ):
        return {[*self._graph.keys()][k]: k for k in range(len([*self._graph.keys()]))}

    def map(
            self,
            graph,
    ):
        """

        Parameters
        ----------
        graph

        Returns
        -------

        """
        self.console.print('===>the graph is being mapped')
        self.console.print('======>key map: {}'.format(self.key_mapped))
        if isinstance(graph, dict):
            self.console.print('======>the graph is a dict')
            g_mapped = {}
            for k, vals in self._graph.items():
                g_mapped[self.key_mapped[k]] = []
                for val in vals:
                    g_mapped[self.key_mapped[k]].append(self.key_mapped[val])
            self.console.print('======>the mapped graph: {}'.format(g_mapped))
            return g_mapped

    def dict(self, ):
        return self._graph

    def set(self, ):
        adj_set = {}
        for k, vals in self._graph.items():
            adj_set[k] = set(vals)
        return adj_set

    def list(self, ):
        return [*self._graph.values()]

    def to_matrix(self, ):
        adj_mat = np.zeros(shape=[self.glen, self.glen])
        for k, vals in self._graph.items():
            for val in vals:
                adj_mat[self.key_mapped[k], self.key_mapped[val]] = 1
        return adj_mat

    def hash(self, ):
        return

    def to_edge_list(
            self,
            rr=True,
    ):
        """

        Parameters
        ----------
        rr
            if it is not symmetry.
            False
            [(0, 1), (0, 2), (0, 3), (1, 0), (1, 2), (2, 0), (2, 1), (3, 0), (3, 4), (3, 5), (4, 3), (4, 6), (5, 3), (5, 6), (6, 4), (6, 5)]
            True
            [(2, 1), (6, 5), (4, 3), (2, 0), (6, 4), (3, 0), (1, 0), (5, 3)]

        Returns
        -------

        """
        edges = []
        for k, vals in self._graph.items():
            for val in vals:
                edges.append((k, val))
        if rr:
            repeat = []
            edge_set = set(edges)
            while edges:
                edge = edges.pop(0)
                if tuple(reversed(edge)) in edges:
                    repeat.append(edge)
            edges = list(edge_set.difference(set(repeat)))
        return edges


if __name__ == "__main__":
    # ### @@ data from UMI-tools
    # graph_adj_dict = {
    #     'A': ['B', 'C'],
    #     'B': ['A', 'C'],
    #     'C': ['A', 'B'],
    #     'D': ['E', 'F'],
    #     'E': ['D'],
    #     'F': ['D'],
    # }

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

    p = Adjacency(graph_adj)

    p.graph = p.graph_mapped

    print(p.graph)

    print(p.dict())

    print(p.set())

    print(p.list())

    print(p.to_matrix())

    print(p.to_edge_list(rr=False))
    # print(p.to_edge_list(rr=True))