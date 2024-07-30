__version__ = "v1.0"
__copyright__ = "Copyright 2024"
__license__ = "MIT"
__developer__ = "Jianfeng Sun"
__maintainer__ = "Jianfeng Sun"
__email__="jianfeng.sunmt@gmail.com"
__lab__ = "Cribbslab"

from typing import Dict

import pandas as pd

from umiche.deduplicate.method.Cluster import Cluster
from umiche.deduplicate.method.Adjacency import Adjacency
from umiche.deduplicate.method.Directional import Directional


def cluster(
        graph,
        method : str = 'deque'
):
    if method == 'deque':
        return Cluster().cc(
            graph_adj=graph,
        )
    if method == 'networkx':
        return Cluster().ccnx(
            edge_list=graph,
        )


def adjacency(
        connected_components : Dict,
        df_umi_uniq_val_cnt : pd.Series,
        graph_adj : Dict,
):
    return Adjacency().umi_tools(
        connected_components=connected_components,
        df_umi_uniq_val_cnt=df_umi_uniq_val_cnt,
        graph_adj=graph_adj
    )


def decompose(
        cc_sub_dict : Dict
):
    return Adjacency().decompose(
        cc_sub_dict=cc_sub_dict,
    )



if __name__ == "__main__":
    graph_adj = {
        'A': ['B', 'C', 'D'],
        'B': ['A', 'C'],
        'C': ['A', 'B'],
        'D': ['A', 'E', 'F'],
        'E': ['D', 'G'],
        'F': ['D', 'G'],
        'G': ['E', 'F'],
    }

    ccs = cluster(
        graph=graph_adj,
        method='deque',
    )
    # print(ccs)

    # edge_list = [('C', 'B'), ('C', 'A'), ('B', 'A'), ('G', 'E'), ('E', 'D'), ('G', 'F'), ('D', 'A'), ('F', 'D')]
    #
    # print(cluster(
    #     graph=edge_list,
    #     method='networkx',
    # ))

    node_val_sorted = pd.Series({
        'A': 120,
        'D': 90,
        'E': 50,
        'G': 5,
        'B': 2,
        'C': 2,
        'F': 1,
    })

    dedup_res = adjacency(
        connected_components=ccs,
        df_umi_uniq_val_cnt=node_val_sorted,
        graph_adj=graph_adj
    )
    dedup_count = dedup_res['count']
    dedup_clusters = dedup_res['clusters']
    print("deduplicated count:\n{}".format(dedup_count))
    print("deduplicated clusters:\n{}".format(dedup_clusters))

    dedup_clusters_dc = decompose(dedup_clusters)
    print("deduplicated clusters decomposed:\n{}".format(dedup_clusters_dc))