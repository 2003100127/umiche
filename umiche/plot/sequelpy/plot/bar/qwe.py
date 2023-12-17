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
            which_color="ccs4",
    ):
        self.pele = pele()
        self.netadj = netadj()
        self.graph = graph
        self.netadj.graph = self.graph
        self.color_list = self.pele.color(which=which_color, is_random=True)

    def draw(
            self,
            ccs,
            ax,
    ):
        el = self.netadj.to_edge_list()
        color_per_subcc = [self.color_list[i] for i, cc in enumerate(ccs)]
        # print(color_per_subcc)

        # relationships = pd.DataFrame({'from': ['10', '10', '10', '5000', '10000'],
        #                               'to':   ['100', '500', '1000', '500', '500']})

        G = nx.Graph()
        G.add_nodes_from(self.graph.keys())
        G.add_edges_from(el)
        print(G.nodes())

        # explicitly set positions
        pos = {
            'A:120': (0, 0),
            'B:2': (-1, 1),
            'C:2': (1, 1),
            'D:90': (0, -2),
            'E:50': (-1, -4.5),
            'F:1': (1, -4.5),
            'G:5': (0, -6),
        }

        options = {
            # "font_size": 36,
            "font_size": 16,
            "node_size": 3000,
            "node_color": "bisque",
            "edgecolors": color_per_subcc,
            # "edgecolors": ['firebrick', 'firebrick', 'firebrick', 'royalblue', 'springgreen', 'royalblue', 'springgreen'],
            "linewidths": 5,
            "width": 5,

        }
        nx.draw_networkx(G, pos, ax=ax, **options)

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
        'E:50': ['D:90', 'G:5'],
        'F:1': ['D:90', 'G:5'],
        'G:5': ['E:50', 'F:1'],
    }

    # graph_adj = {
    #     'A': ['B', 'C', 'D'],
    #     'B': ['A', 'C'],
    #     'C': ['A', 'B'],
    #     'D': ['A', 'E', 'F'],
    #     'E': ['D', 'G'],
    #     'F': ['D', 'G'],
    #     'G': ['E', 'F'],
    # }
    # print("An adjacency list of a graph:\n{}".format(graph_adj))

    node_val_sorted = pd.Series({
        'A:120': 120,
        'D:90': 90,
        'E:50': 50,
        'G:5': 5,
        'B:2': 2,
        'C:2': 2,
        'F:1': 1,
    })
    print("Counts sorted:\n{}".format(node_val_sorted))

    # ccs = [['A:120', 'B:2', 'C:2', 'D:90', 'E:50', 'F:1']]  # cc
    # ccs = [['A:120', 'B:2', 'C:2'], ['D:90', 'E:50', 'F:1'], ['G:5']]  # adj
    # ccs = [['A:120', 'B:2', 'C:2'], ['D:90', 'F:1', ], ['E:50', 'G:5']]  # direc
    # ccs = [['A:120', 'B:2', 'C:2'], ['D:90', 'E:50', 'F:1', 'G:5']] # mcl
    # ccs = [['A:120', 'B:2', 'C:2', 'D:90', 'E:50', 'F:1', 'G:5']] # mcl_ed
    # ccs = [['A:120', 'B:2', 'C:2'], ['D:90', 'E:50', 'F:1', 'G:5']] # mcl_val

    ### @@@ ******Connected components******
    from umiche.deduplicate.method.Cluster import Cluster as umiclust
    ccs = umiclust().cc(graph_adj=graph_adj)
    print("Connected components:\n{}".format(ccs))

    ### @@@ ******UMI-tools Adjacency******
    from umiche.deduplicate.method.Adjacency import Adjacency as umiadj
    from umiche.deduplicate.method.Directional import Directional as umidirec
    from umiche.deduplicate.method.MarkovClustering import MarkovClustering as umimcl
    dedup_res_adj = umiadj().umi_tools(
        connected_components=ccs,
        df_umi_uniq_val_cnt=node_val_sorted,
        graph_adj=graph_adj,
    )
    dedup_res_adj_dc = umiadj().decompose(dedup_res_adj['clusters'])
    print("deduplicated clusters (UMI-tools Adjacency):\n{}".format(dedup_res_adj_dc))

    ### @@@ ******UMI-tools Directional******
    from umiche.deduplicate.method.Directional import Directional as umidirec
    from umiche.deduplicate.method.MarkovClustering import MarkovClustering as umimcl
    dedup_res_direc = umidirec().umi_tools(
        connected_components=ccs,
        df_umi_uniq_val_cnt=node_val_sorted,
        graph_adj=graph_adj,
    )
    dedup_res_direc_dc = umidirec().decompose(dedup_res_direc['clusters'])
    print("deduplicated clusters (UMI-tools Directional):\n{}".format(dedup_res_direc_dc))

    ### @@@ ******MCL******
    from umiche.deduplicate.method.MarkovClustering import MarkovClustering as umimcl
    mcl = umimcl(
        inflat_val=1.6,
        exp_val=2,
        iter_num=100,
    )
    df_mcl = mcl.dfclusters(
        connected_components=ccs,
        graph_adj=graph_adj,
    )
    dedup_res_mcl_dc = mcl.decompose(list_nd=df_mcl['clusters'].values)
    print("deduplicated clusters (MCL):\n{}".format(dedup_res_mcl_dc))

    ### @@@ ******MCL mcl_val******
    df_mcl_val = mcl.maxval_val(
        df_mcl_ccs=df_mcl,
        df_umi_uniq_val_cnt=node_val_sorted,
        thres_fold=2,
    )
    dedup_res_mcl_val_dc = mcl.decompose(list_nd=df_mcl_val['clusters'].values)
    print("deduplicated clusters decomposed (mcl_val):\n{}".format(dedup_res_mcl_val_dc))
    dedup_res_mcl_val_dc_full = mcl.get_full_subcc(ccs_dict=dedup_res_mcl_val_dc, mcl_ccs_dict=dedup_res_mcl_dc)
    print("deduplicated clusters decomposed full list(mcl_val):\n{}".format(dedup_res_mcl_val_dc_full))

    ### @@@ ******MCL mcl_ed******
    int_to_umi_dict = {
        'A:120': 'AGATCTCGCA',
        'B:2': 'AGATCCCGCA',
        'C:2': 'AGATCACGCA',
        'D:90': 'AGATCGCGCA',
        'E:50': 'AGATCGCGGA',
        'F:1': 'AGATCGCGTA',
        'G:5': 'TGATCGCGAA',
    }
    df_mcl_ed = mcl.maxval_ed(
        df_mcl_ccs=df_mcl,
        df_umi_uniq_val_cnt=node_val_sorted,
        thres_fold=1,
        umi_uniq_mapped_rev=int_to_umi_dict,
    )
    dedup_res_mcl_ed_dc = mcl.decompose(list_nd=df_mcl_ed['clusters'].values)
    print("deduplicated clusters decomposed (mcl_ed):\n{}".format(dedup_res_mcl_ed_dc))
    dedup_res_mcl_ed_dc_full = mcl.get_full_subcc(ccs_dict=dedup_res_mcl_ed_dc, mcl_ccs_dict=dedup_res_mcl_dc)
    print("deduplicated clusters decomposed full list(mcl_ed):\n{}".format(dedup_res_mcl_ed_dc_full))

    fig, ax = plt.subplots(nrows=2, ncols=3, figsize=(12, 10))
    p = Graph(graph=graph_adj)
    p.draw(ccs, ax=ax[0, 1])
    plt.show()
    # print()