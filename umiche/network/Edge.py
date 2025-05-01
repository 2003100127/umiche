__version__ = "v1.0"
__copyright__ = "Copyright 2025"
__license__ = "GPL-3.0"
__developer__ = "Jianfeng Sun"
__maintainer__ = "Jianfeng Sun"
__email__="jianfeng.sunmt@gmail.com"


import numpy as np
from umiche.util.Console import Console


class Edge:

    def __init__(
            self,
            graph=None,
            verbose=False,
    ):
        self._graph = graph
        self._graph_mapped = None

        self.console = Console()
        self.console.verbose = verbose

    @property
    def graph(self, ):
        self.console.print('======>the current graph is {}'.format(self._graph))
        return self._graph

    @property
    def glen(self, ):
        self.console.print('======>the number of egdes in the current graph is {}'.format(len(self._graph)))
        return len(self._graph)

    @property
    def nodes(self, ):
        g_np = np.array(self._graph)
        # print(self._graph)
        if g_np.size != 0:
            l = np.array(g_np)[:, 0]
            r = np.array(g_np)[:, 1]
            return list(np.unique(np.concatenate((l, r))))
        else:
            return []

    @property
    def key_mapped(self, ):
        return {k: id for id, k in enumerate(self.nodes)}

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
            raise 'please set your graph using an edge list'
        else:
            return self.map(self._graph)

    @property
    def rvredanduncy(self, ):
        repeat = []
        edges = self._graph.copy()
        edge_set = set(self._graph)
        while edges:
            edge = edges.pop(0)
            if tuple(reversed(edge)) in edges:
                repeat.append(edge)
        edges = list(edge_set.difference(set(repeat)))
        # print(self._graph)
        return edges

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
        # print('===>the graph is being mapped')
        # print('======>key map: {}'.format(self.key_mapped))
        g_mapped = []
        for i, edge in enumerate(graph):
            g_mapped.append((self.key_mapped[edge[0]], self.key_mapped[edge[1]]))
        # print('======>the mapped graph: {}'.format(g_mapped))
        return g_mapped

    def to_adj_dict(self, ):
        adj_list = {}
        # scan the arrays edge_u and edge_v
        # print('nodes are {}'.format(self.nodes))
        for i in self.nodes:
            adj_list[i] = []
        for i in range(self.glen):
            l = self._graph[i][0]
            r = self._graph[i][1]
            adj_list[l].append(r)
            adj_list[r].append(l)
        # print(adj_list)
        return adj_list

    def fromlist(
            self,
            list_2d,
    ):
        """

        Parameters
        ----------
        list_2d

        Returns
        -------

        """
        tuples = []
        for e in list_2d:
            tuples.append((e[0], e[1]))
        return tuples

    def fromdf(
            self,
            df,
            col_name1='node_1',
            col_name2='node_2',
            to_tuple=False,
    ):
        """

        Parameters
        ----------
        df
        col_name1
        col_name2
        to_tuple

        Returns
        -------

        """
        edge_list = df[[col_name1, col_name2]].values.tolist()
        if to_tuple:
            return self.fromlist(edge_list)
        else:
            return edge_list


if __name__ == "__main__":
    graph_edge_list = [
        ('A', 'B'),
        ('A', 'C'),
        ('B', 'C'),
        ('A', 'D'),
        ('D', 'E'),
        ('D', 'F'),
        ('D', 'F'),
    ]

    p = Edge(graph_edge_list)

    # p.graph = p.graph_mapped

    print(p.graph)
    print(p.key_mapped)

    p.graph = p.rvredanduncy

    print(p.graph)
    print(p.key_mapped)

    # print(p.rvredanduncy)

    print(p.to_adj_dict())
    print(p.graph_mapped)

    # print(p.nodes)