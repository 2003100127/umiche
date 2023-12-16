from igraph import *


def toEdgeList(graph_adj, rr=True):
    edges = []
    for k, vals in graph_adj.items():
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

# g = Graph.Erdos_Renyi(30,0.3)
# print(g)
g = Graph()
g.add_vertices(['A','B','C','D','E','F',])
graph_adj = {
        'A': ['B', 'C', 'D'],
        'B': ['A', 'C'],
        'C': ['A', 'B'],
        'D': ['A', 'E', 'F'],
        'E': ['D'],
        'F': ['D'],
    }
print()
g.add_edges(toEdgeList(graph_adj))
print(g)
comms = g.community_multilevel()

print(comms)
print(comms[0])
print(comms[1])
plot(comms, mark_groups = True)