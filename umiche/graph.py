__version__ = "v1.0"
__copyright__ = "Copyright 2025"
__license__ = "GPL-3.0"
__developer__ = "Jianfeng Sun"
__maintainer__ = "Jianfeng Sun"
__email__="jianfeng.sunmt@gmail.com"


from typing import Dict

from umiche.network.Adjacency import Adjacency
from umiche.network.Edge import Edge
from umiche.network.CC import CC


def adjacency(
        graph_adj : Dict,
):
    return Adjacency(
        graph=graph_adj,
    )


def edge(
        graph_edge_list,
):
    return Edge(
        graph=graph_edge_list
    )


def cc(
        graph_adj : Dict,
        method : str = 'deque',
        verbose : bool = True,
):
    if method == "deque":
        return CC(
            verbose=verbose,
        ).deque(
            graph=graph_adj,
        )
    if method == "set":
        return CC(
            verbose=verbose,
        ).set(
            graph=graph_adj,
        )


if __name__ == "__main__":
    graph_adj_umitools = {
        'A': ['B', 'C', 'D'],
        'B': ['A', 'C'],
        'C': ['A', 'B'],
        'D': ['A', 'E', 'F'],
        'E': ['D'],
        'F': ['D'],
    }
    graph_adj_mclumi = {
        'A': ['B', 'C', 'D'],
        'B': ['A', 'C'],
        'C': ['A', 'B'],
        'D': ['A', 'E', 'F'],
        'E': ['D', 'G'],
        'F': ['D', 'G'],
        'G': ['E', 'F'],
    }

    # adj = adjacency(
    #     # graph_adj=graph_adj_umitools,
    #     graph_adj=graph_adj_mclumi,
    # )
    # print(adj.graph_mapped)
    # print(adj.graph)
    # print(adj.dict())
    # print(adj.set())
    # print(adj.list())
    # print(adj.to_matrix())
    # print(adj.to_edge_list(rr=False))
    # print(adj.to_edge_list(rr=True))

    # eg = edge(
    #     graph_edge_list=adj.to_edge_list(rr=False),
    # )
    # print(eg.graph)
    # print(eg.nodes)
    # print(eg.key_mapped)
    # print(eg.rvredanduncy)
    # eg.graph = eg.rvredanduncy
    # print(eg.graph)
    # print(eg.key_mapped)
    # print(eg.rvredanduncy)
    # print(eg.to_adj_dict())
    # print(eg.graph_mapped)

    connected_components = cc(
        graph_adj=graph_adj_mclumi,
        method='deque',
        verbose=True
    )
    print(list(connected_components))

    connected_components = cc(
        graph_adj=graph_adj_mclumi,
        method='set',
        verbose=True
    )
    print(connected_components)