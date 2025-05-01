__version__ = "v1.0"
__copyright__ = "Copyright 2025"
__license__ = "GPL-3.0"
__developer__ = "Jianfeng Sun"
__maintainer__ = "Jianfeng Sun"
__email__="jianfeng.sunmt@gmail.com"


import networkx as nx
import seaborn as sns
import matplotlib.pyplot as plt
from umiche.plot.gadget.Element import Element as pele
from umiche.network.Adjacency import Adjacency as netadj


class Graph:

    def __init__(
            self,
            graph,
            which_color="tableau",
    ):
        self.pele = pele()
        self.netadj = netadj()
        self.graph = graph
        self.netadj.graph = self.graph
        self.color_list = self.pele.color(which=which_color, is_random=True)

        sns.set(font="Helvetica")
        sns.set_style("ticks")

    def draw(
            self,
            cc_dict,
            title,
            ax,
    ):
        el = self.netadj.to_edge_list()
        color_per_subcc = {}
        for cc_i, (cc_id, cc) in enumerate(cc_dict.items()):
            for node_j, node in enumerate(cc):
                color_per_subcc[node] = self.color_list[cc_i]
        print(color_per_subcc)
        color_per_subcc = [color_per_subcc[k] for k in self.graph.keys()]

        G = nx.Graph()
        G.add_nodes_from(self.graph.keys())
        G.add_edges_from(el)
        print(G.nodes())

        # pos = {
        #     'A:120': (0, 0),
        #     'B:2': (-1, 1),
        #     'C:2': (1, 1),
        #     'D:90': (0, -2.5),
        #     'E:50': (-1, -4.5),
        #     'F:1': (1, -4.5),
        #     'G:5': (0, -6),
        # }

        pos = {
            'A': (0, 0),
            'B': (-1, 1),
            'C': (1, 1),
            'D': (0, -2.5),
            'E': (-1, -4.5),
            'F': (1, -4.5),
        }

        options = {
            "font_size": 14,
            "node_size": 2000,
            "node_color": "white", # bisque
            "edgecolors": color_per_subcc, # node edge color
            "edge_color": 'dimgray', # simply edge color
            "linewidths": 3,
            "width": 4,
        }
        nx.draw_networkx(G, pos, ax=ax, **options)
        ax.set_title(title, fontsize=18)
        ax.set_axis_off()


if __name__ == "__main__":
    import pandas as pd

    # graph_adj = {
    #     'A:120': ['B:2', 'C:2', 'D:90'],
    #     'B:2': ['A:120', 'C:2'],
    #     'C:2': ['A:120', 'B:2'],
    #     'D:90': ['A:120', 'E:50', 'F:1'],
    #     'E:50': ['D:90', 'G:5'],
    #     'F:1': ['D:90', 'G:5'],
    #     'G:5': ['E:50', 'F:1'],
    # }

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

    # node_val_sorted = pd.Series({
    #     'A:120': 120,
    #     'D:90': 90,
    #     'E:50': 50,
    #     'G:5': 5,
    #     'B:2': 2,
    #     'C:2': 2,
    #     'F:1': 1,
    # })
    # print("Counts sorted:\n{}".format(node_val_sorted))

    ### @@ data from UMI-tools
    graph_adj = {
        'A': ['B', 'C', 'D'],
        'B': ['A', 'C'],
        'C': ['A', 'B'],
        'D': ['A', 'E', 'F'],
        'E': ['D'],
        'F': ['D'],
    }
    print("An adjacency list of a graph:\n{}".format(graph_adj))

    node_val_sorted = pd.Series({
        'A': 456,
        'E': 90,
        'D': 72,
        'B': 2,
        'C': 2,
        'F': 1,
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
    # int_to_umi_dict = {
    #     'A:120': 'AGATCTCGCA',
    #     'B:2': 'AGATCCCGCA',
    #     'C:2': 'AGATCACGCA',
    #     'D:90': 'AGATCGCGCA',
    #     'E:50': 'AGATCGCGGA',
    #     'F:1': 'AGATCGCGTA',
    #     'G:5': 'TGATCGCGAA',
    # }
    int_to_umi_dict = {
        'A': 'ACGT',
        'B': 'TCGT',
        'C': 'CCGT',
        'D': 'ACAT',
        'E': 'ACAG',
        'F': 'AAAT',
    }
    df_mcl_ed = mcl.maxval_ed(
        df_mcl_ccs=df_mcl,
        df_umi_uniq_val_cnt=node_val_sorted,
        thres_fold=1,
        int_to_umi_dict=int_to_umi_dict,
    )
    dedup_res_mcl_ed_dc = mcl.decompose(list_nd=df_mcl_ed['clusters'].values)
    print("deduplicated clusters decomposed (mcl_ed):\n{}".format(dedup_res_mcl_ed_dc))
    dedup_res_mcl_ed_dc_full = mcl.get_full_subcc(ccs_dict=dedup_res_mcl_ed_dc, mcl_ccs_dict=dedup_res_mcl_dc)
    print("deduplicated clusters decomposed full list(mcl_ed):\n{}".format(dedup_res_mcl_ed_dc_full))

    fig, ax = plt.subplots(nrows=2, ncols=3, figsize=(14, 9.5))
    p = Graph(graph=graph_adj)
    p.color_list = ['cornflowerblue', 'lightcoral', 'mediumseagreen',]
    p.draw(ccs, ax=ax[0, 0], title='Cluster')
    p.draw(dedup_res_adj_dc, title='Adjacent', ax=ax[0, 1])
    p.draw(dedup_res_direc_dc, title='Directional', ax=ax[0, 2])
    p.draw(dedup_res_mcl_dc, title='MCL', ax=ax[1, 0])
    p.draw(dedup_res_mcl_val_dc_full, title='MCL-val', ax=ax[1, 1])
    p.draw(dedup_res_mcl_ed_dc_full, title='MCL-ed', ax=ax[1, 2])

    plt.subplots_adjust(
        top=0.92,
        bottom=0.04,
        left=0.04,
        right=0.98,
        hspace=0.10,
        wspace=0.15,
    )
    plt.show()
