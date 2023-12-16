import matplotlib.pyplot as plt
import pandas as pd
import numpy as np
import networkx as nx
import matplotlib.colors as mcolors
color = mcolors.CSS4_COLORS
# colors = np.random.permutation([*color.keys()]).tolist()
colors = ['firebrick', 'royalblue', 'springgreen', 'burlywood', 'oldlace', 'lightgoldenrodyellow', 'grey', 'cyan', 'crimson', 'mediumvioletred', 'maroon', 'mediumturquoise', 'teal', 'azure', 'palevioletred', 'mediumslateblue', 'olivedrab', 'darkred', 'dimgrey', 'lightgreen', 'blueviolet', 'chartreuse', 'whitesmoke', 'mediumspringgreen', 'yellowgreen', 'green', 'darkorchid', 'lawngreen', 'black', 'darkgoldenrod', 'moccasin', 'seashell', 'orange', 'slateblue', 'hotpink', 'tan', 'darkolivegreen', 'mediumorchid', 'darkgreen', 'navajowhite', 'khaki', 'paleturquoise', 'darkgray', 'darkorange', 'salmon', 'floralwhite', 'lightgray', 'darkblue', 'lightslategray', 'darkslateblue', 'goldenrod', 'lightcyan', 'papayawhip', 'wheat', 'lightslategrey', 'darksalmon', 'skyblue', 'violet', 'coral', 'mediumseagreen', 'darkmagenta', 'palegreen', 'magenta', 'lavenderblush', 'darkkhaki', 'darkturquoise', 'lightsalmon', 'blanchedalmond', 'deepskyblue', 'mediumaquamarine', 'mintcream', 'lightyellow', 'linen', 'peachpuff', 'aliceblue', 'gold', 'ghostwhite', 'silver', 'bisque', 'aquamarine', 'lightblue', 'lavender', 'antiquewhite', 'dimgray', 'lightskyblue', 'pink', 'beige', 'indianred', 'rosybrown', 'rebeccapurple', 'lightgrey', 'powderblue', 'chocolate', 'darkgrey', 'purple', 'darkslategrey', 'cornsilk', 'turquoise', 'blue', 'greenyellow', 'red', 'lime', 'ivory', 'royalblue', 'palegoldenrod', 'lightsteelblue', 'slategray', 'gray', 'brown', 'mistyrose', 'cornflowerblue', 'tomato', 'indigo', 'snow', 'darkslategray', 'lightpink', 'darkviolet', 'limegreen', 'honeydew', 'yellow', 'sienna', 'mediumblue', 'fuchsia', 'forestgreen', 'orchid', 'plum', 'lemonchiffon', 'gainsboro', 'sandybrown', 'slategrey', 'mediumpurple', 'saddlebrown', 'lightcoral', 'midnightblue', 'navy', 'thistle', 'dodgerblue', 'lightseagreen', 'darkseagreen', 'cadetblue', 'seagreen', 'darkcyan', 'steelblue', 'white', 'peru', 'deeppink', 'olive', 'firebrick', 'aqua', 'orangered']
print(colors)
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


graph_adj = {
        'A:120': ['B:2', 'C:2', 'D:90'],
        'B:2': ['A:120', 'C:2'],
        'C:2': ['A:120', 'B:2'],
        'D:90': ['A:120', 'E:50', 'F:1'],
        'E:50': ['D:90', 'G:1'],
        'F:1': ['D:90', 'G:1'],
        'G:1': ['E:50', 'F:1'],
    }

el = toEdgeList(graph_adj)
ccs = [['A:120', 'B:2', 'C:2', 'D:90', 'E:50', 'F:1']] # cc
ccs = [['A:120', 'B:2', 'C:2'], ['D:90', 'E:50', 'F:1'], ['G:1']] # adj
ccs = [['A:120', 'B:2', 'C:2'], ['D:90', 'F:1', ], ['E:50', 'G:1']] # direc
# ccs = [['A:120', 'B:2', 'C:2'], ['D:90', 'E:50', 'F:1', 'G:1']] # mcl
# ccs = [['A:120', 'B:2', 'C:2', 'D:90', 'E:50', 'F:1', 'G:1']] # mcl_ed
# ccs = [['A:120', 'B:2', 'C:2'], ['D:90', 'E:50', 'F:1', 'G:1']] # mcl_val
ncs = []
for i, cc in enumerate(ccs):
    for j in cc:
        ncs.append(colors[i])
print(ncs)
fig, ax = plt.subplots(figsize=(7, 7))

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
#
# # Manually set each nodes size
# node_sizes = [4560, 20, 20, 720, 900, 10,]
#
# # Draw graph
# nx.draw(G, with_labels=True, node_size=node_sizes)

# plt.show()
