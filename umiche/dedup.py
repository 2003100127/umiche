__version__ = "v1.0"
__copyright__ = "Copyright 2025"
__license__ = "GPL-3.0"
__developer__ = "Jianfeng Sun"
__maintainer__ = "Jianfeng Sun"
__email__="jianfeng.sunmt@gmail.com"


from typing import List, Dict

import pandas as pd

from umiche.deduplicate.method.Cluster import Cluster
from umiche.deduplicate.method.Adjacency import Adjacency
from umiche.deduplicate.method.Directional import Directional
from umiche.deduplicate.method.MarkovClustering import MarkovClustering
from umiche.deduplicate.method.Clustering import Clustering
from umiche.deduplicate.method.trimer.MajorityVote import MajorityVote
from umiche.deduplicate.method.trimer.SetCover import SetCover


def cluster(
        graph,
        method : str = 'deque',
        verbose : bool = False,
):
    if method == 'deque':
        return Cluster(
            verbose=verbose,
        ).cc(
            graph_adj=graph,
        )
    if method == 'networkx':
        return Cluster(
            verbose=verbose,
        ).ccnx(
            edge_list=graph,
        )


def adjacency(
        connected_components : Dict,
        df_umi_uniq_val_cnt : pd.Series,
        graph_adj : Dict,
        verbose: bool = False,
):
    return Adjacency(
        verbose=verbose,
    ).umi_tools(
        connected_components=connected_components,
        df_umi_uniq_val_cnt=df_umi_uniq_val_cnt,
        graph_adj=graph_adj,
    )


def directional(
        connected_components : Dict,
        df_umi_uniq_val_cnt : pd.Series,
        graph_adj : Dict,
        verbose: bool = False,
):
    return Directional(
        verbose=verbose,
    ).umi_tools(
        connected_components=connected_components,
        df_umi_uniq_val_cnt=df_umi_uniq_val_cnt,
        graph_adj=graph_adj,
    )


def mcl(
        inflat_val : float,
        exp_val : int,
        iter_num : int,
        connected_components : Dict,
        graph_adj : Dict,
        verbose: bool = False,
):
    return MarkovClustering(
        inflat_val=inflat_val,
        exp_val=exp_val,
        iter_num=iter_num,
        verbose=verbose,
    ).dfclusters(
        connected_components=connected_components,
        graph_adj=graph_adj,
    )


def mcl_val(
        df_mcl_ccs : pd.DataFrame,
        df_umi_uniq_val_cnt : pd.Series,
        thres_fold : int,
        verbose: bool = False,
):
    return MarkovClustering(
        inflat_val=1.6,
        exp_val=2,
        iter_num=100,
        verbose=verbose,
    ).maxval_val(
        df_mcl_ccs=df_mcl_ccs,
        df_umi_uniq_val_cnt=df_umi_uniq_val_cnt,
        thres_fold=thres_fold,
    )


def mcl_ed(
        df_mcl_ccs : pd.DataFrame,
        df_umi_uniq_val_cnt : pd.Series,
        thres_fold : int,
        int_to_umi_dict : Dict,
        verbose: bool = False,
):
    return MarkovClustering(
        inflat_val=1.6,
        exp_val=2,
        iter_num=100,
        verbose=verbose,
    ).maxval_ed(
        df_mcl_ccs=df_mcl_ccs,
        df_umi_uniq_val_cnt=df_umi_uniq_val_cnt,
        thres_fold=thres_fold,
        int_to_umi_dict=int_to_umi_dict,
    )


def dbscan(
        connected_components : Dict,
        graph_adj : Dict,
        int_to_umi_dict : Dict,
        dbscan_eps=1.5,
        dbscan_min_spl=1,
        verbose: bool = False,
):
    return Clustering(
        clustering_method='dbscan',
        dbscan_eps=dbscan_eps,
        dbscan_min_spl=dbscan_min_spl,
        verbose=verbose,
    ).dfclusters(
        connected_components=connected_components,
        graph_adj=graph_adj,
        int_to_umi_dict=int_to_umi_dict,
    )


def birch(
        connected_components : Dict,
        graph_adj : Dict,
        int_to_umi_dict : Dict,
        birch_thres=1.8,
        birch_n_clusters=None,
        verbose: bool = False,
):
    return Clustering(
        clustering_method='birch',
        birch_thres=birch_thres,
        birch_n_clusters=birch_n_clusters,
        verbose=verbose,
    ).dfclusters(
        connected_components=connected_components,
        graph_adj=graph_adj,
        int_to_umi_dict=int_to_umi_dict,
    )


def affinity_propagation(
        connected_components : Dict,
        graph_adj : Dict,
        int_to_umi_dict : Dict,
        verbose: bool = False,
):
    return Clustering(
        clustering_method='aprop',
        verbose=verbose,
    ).dfclusters(
        connected_components=connected_components,
        graph_adj=graph_adj,
        int_to_umi_dict=int_to_umi_dict,
    )


def set_cover(
        multimer_list,
        recur_len=3,
        split_method='split_to_all',
        verbose=True,
):
    return SetCover(
        verbose=verbose,
    ).greedy(
        multimer_list=multimer_list,
        recur_len=recur_len,
        split_method=split_method,
    )


def majority_vote(
        multimer_list,
        recur_len=3,
        verbose=True,
):
    return MajorityVote(
        verbose=verbose,
    ).track(
        multimer_list=multimer_list,
        recur_len=recur_len,
    )


def decompose(
        cc_sub_dict : Dict,
        verbose: bool = False,
):
    return Adjacency(
        verbose=verbose,
    ).decompose(
        cc_sub_dict=cc_sub_dict,
    )


def decompose_mcl(
        list_nd : List,
        verbose: bool = False,
):
    return MarkovClustering(
        inflat_val=1.6,
        exp_val=2,
        iter_num=100,
        verbose=verbose,
    ).decompose(
        list_nd=list_nd,
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

    ### ++++++++++++++adjacency++++++++++++++++++++++
    # dedup_res = adjacency(
    #     connected_components=ccs,
    #     df_umi_uniq_val_cnt=node_val_sorted,
    #     graph_adj=graph_adj,
    #     verbose=True,
    # )
    # dedup_count = dedup_res['count']
    # dedup_clusters = dedup_res['clusters']
    # print("deduplicated count:\n{}".format(dedup_count))
    # print("deduplicated clusters:\n{}".format(dedup_clusters))
    # dedup_clusters_dc = decompose(dedup_clusters)
    # print("deduplicated clusters decomposed:\n{}".format(dedup_clusters_dc))

    ### ++++++++++++++directional++++++++++++++++++++++
    # dedup_res = directional(
    #     connected_components=ccs,
    #     df_umi_uniq_val_cnt=node_val_sorted,
    #     graph_adj=graph_adj,
    #     verbose=True,
    # )
    # dedup_count = dedup_res['count']
    # dedup_clusters = dedup_res['clusters']
    # print("deduplicated count:\n{}".format(dedup_count))
    # print("deduplicated clusters:\n{}".format(dedup_clusters))
    # dedup_clusters_dc = decompose(dedup_clusters)
    # print("deduplicated clusters decomposed:\n{}".format(dedup_clusters_dc))
    #
    # print(dedup_res['apv'])
    # print(dedup_res['disapv'])

    # ### ++++++++++++++MCL++++++++++++++++++++++
    # dedup_res = mcl(
    #     inflat_val=1.6,
    #     exp_val=2,
    #     iter_num=100,
    #     connected_components=ccs,
    #     graph_adj=graph_adj,
    #     verbose=True,
    # )
    # print(dedup_res.columns)
    # print("vertices: {}".format(dedup_res.loc[0, 'cc_vertices']))
    # print(dedup_res.loc[0, 'graph_cc_adj'])
    # print(dedup_res.loc[0, 'nt_to_int_map'])
    # print(dedup_res.loc[0, 'int_to_nt_map'])
    # print(dedup_res.loc[0, 'cc_adj_mat'])
    # print(dedup_res.loc[0, 'mcl_clusters'])
    # print(dedup_res.loc[0, 'clusters'])
    # print(dedup_res.loc[0, 'clust_num'])
    #
    # print(dedup_res['clusters'].values)
    # dedup_clusters_dc = decompose_mcl(list_nd=dedup_res['clusters'].values)
    # print("deduplicated clusters decomposed (mcl):\n{}".format(dedup_clusters_dc))

    ### ++++++++++++++MCL-val++++++++++++++++++++++
    # df_mcl_val = mcl_val(
    #     df_mcl_ccs=dedup_res,
    #     df_umi_uniq_val_cnt=node_val_sorted,
    #     thres_fold=1,
    # )
    # print(df_mcl_val)
    # dedup_count = df_mcl_val['count'].values[0]
    # dedup_clusters = df_mcl_val['clusters'].values[0]
    # print("deduplicated count (mcl_val):\n{}".format(dedup_count))
    # print("deduplicated clusters (mcl_val):\n{}".format(dedup_clusters))
    #
    # df_mcl_val = decompose_mcl(list_nd=df_mcl_val['clusters'].values)
    # print("deduplicated clusters decomposed (mcl_val):\n{}".format(df_mcl_val))

    ### ++++++++++++++MCL-ed++++++++++++++++++++++
    # int_to_umi_dict = {
    #     'A': 'AGATCTCGCA',
    #     'B': 'AGATCCCGCA',
    #     'C': 'AGATCACGCA',
    #     'D': 'AGATCGCGCA',
    #     'E': 'AGATCGCGGA',
    #     'F': 'AGATCGCGTA',
    #     'G': 'TGATCGCGAA',
    # }
    #
    # df_mcl_ed = mcl_ed(
    #     df_mcl_ccs=dedup_res,
    #     df_umi_uniq_val_cnt=node_val_sorted,
    #     thres_fold=1,
    #     int_to_umi_dict=int_to_umi_dict,
    # )
    # dedup_count = df_mcl_ed['count']
    # dedup_clusters = df_mcl_ed['clusters']
    # print('approval: {}'.format(df_mcl_ed['apv']))
    #
    # print("deduplicated count (mcl_ed):\n{}".format(dedup_count))
    # print("deduplicated clusters (mcl_ed):\n{}".format(dedup_clusters))
    #
    # df_mcl_ed = decompose_mcl(list_nd=df_mcl_ed['clusters'].values)
    # print("deduplicated clusters decomposed (mcl_ed):\n{}".format(df_mcl_ed))

    # ### ++++++++++++++dbscan++++++++++++++++++++++
    # int_to_umi_dict = {
    #     'A': 'AGATCTCGCA',
    #     'B': 'AGATCCCGCA',
    #     'C': 'AGATCACGCA',
    #     'D': 'AGATCGCGCA',
    #     'E': 'AGATCGCGGA',
    #     'F': 'AGATCGCGTA',
    #     'G': 'TGATCGCGAA',
    # }
    # df = dbscan(
    #     connected_components=ccs,
    #     graph_adj=graph_adj,
    #     # df_umi_uniq_val_cnt=node_val_sorted,
    #     int_to_umi_dict=int_to_umi_dict,
    #     dbscan_eps=1.5,
    #     dbscan_min_spl=1,
    # )
    # # print(df)
    # for i, col in enumerate(df.columns):
    #     print("{}.{}: \n{}".format(i+1, col, df[col].values[0]))
    # df_decomposed = decompose_mcl(list_nd=df['clusters'].values)
    # print("deduplicated clusters decomposed:\n{}".format(df_decomposed))

    # ### ++++++++++++++birch++++++++++++++++++++++
    # int_to_umi_dict = {
    #     'A': 'AGATCTCGCA',
    #     'B': 'AGATCCCGCA',
    #     'C': 'AGATCACGCA',
    #     'D': 'AGATCGCGCA',
    #     'E': 'AGATCGCGGA',
    #     'F': 'AGATCGCGTA',
    #     'G': 'TGATCGCGAA',
    # }
    # df = birch(
    #     connected_components=ccs,
    #     graph_adj=graph_adj,
    #     int_to_umi_dict=int_to_umi_dict,
    #     birch_thres=1.8,
    #     birch_n_clusters=None,
    # )
    # # print(df)
    # for i, col in enumerate(df.columns):
    #     print("{}.{}: \n{}".format(i + 1, col, df[col].values[0]))
    # df_decomposed = decompose_mcl(list_nd=df['clusters'].values)
    # print("deduplicated clusters decomposed:\n{}".format(df_decomposed))

    ### ++++++++++++++birch++++++++++++++++++++++
    # int_to_umi_dict = {
    #     'A': 'AGATCTCGCA',
    #     'B': 'AGATCCCGCA',
    #     'C': 'AGATCACGCA',
    #     'D': 'AGATCGCGCA',
    #     'E': 'AGATCGCGGA',
    #     'F': 'AGATCGCGTA',
    #     'G': 'TGATCGCGAA',
    # }
    # df = affinity_propagation(
    #     connected_components=ccs,
    #     graph_adj=graph_adj,
    #     int_to_umi_dict=int_to_umi_dict,
    # )
    # # print(df)
    # for i, col in enumerate(df.columns):
    #     print("{}.{}: \n{}".format(i + 1, col, df[col].values[0]))
    # df_decomposed = decompose_mcl(list_nd=df['clusters'].values)
    # print("deduplicated clusters decomposed:\n{}".format(df_decomposed))

    # ### ++++++++++++++set_cover++++++++++++++++++++++
    # from umiche.bam.Reader import Reader as alireader
    # alireader = alireader(
    #     bam_fpn="/mnt/d/Document/Programming/Python/umiche/umiche/data/simu/umi/trimer/seq_errs/permute_0/trimmed/seq_err_17.bam",
    #     verbose=True,
    # )
    # df_bam = alireader.todf(tags=['PO'])
    # print(df_bam.columns)
    # trimer_list = df_bam.query_name.apply(lambda x: x.split('_')[1]).values
    # print(trimer_list)
    # (
    #     dedup_cnt,
    #     multimer_umi_solved_by_sc,
    #     multimer_umi_not_solved,
    #     shortlisted_multimer_umi_list,
    #     monomer_umi_lens,
    #     multimer_umi_lens,
    # ) = set_cover(
    #     multimer_list=trimer_list,
    #     recur_len=3,
    #     split_method='split_to_all',
    #     verbose=True
    # )
    # print(dedup_cnt)

    ### ++++++++++++++set_cover++++++++++++++++++++++
    from umiche.bam.Reader import Reader as alireader
    alireader = alireader(
        bam_fpn="/mnt/d/Document/Programming/Python/umiche/umiche/data/simu/umi/trimer/seq_errs/permute_0/trimmed/seq_err_17.bam",
        verbose=True,
    )
    df_bam = alireader.todf(tags=['PO'])
    print(df_bam.columns)
    trimer_list = df_bam.query_name.apply(lambda x: x.split('_')[1]).values
    print(trimer_list)
    (
        dedup_cnt,
        uniq_multimer_cnt,
        shortlisted_multimer_umi_list,
    ) = majority_vote(
        multimer_list=trimer_list,
        recur_len=3,
        verbose=True
    )
    print(dedup_cnt)
    print(uniq_multimer_cnt)
    print(shortlisted_multimer_umi_list)