__version__ = "v1.0"
__copyright__ = "Copyright 2025"
__license__ = "GPL-3.0"
__developer__ = "Jianfeng Sun"
__maintainer__ = "Jianfeng Sun"
__email__="jianfeng.sunmt@gmail.com"


import numpy as np
import pandas as pd

from umiche.util.Hamming import Hamming
from umiche.util.Console import Console


class Gadgetry:

    def __init__(
            self,
            verbose=False,
    ):
        self.console = Console()
        self.console.verbose = verbose

    def num_removed_uniq_umis(
            self,
            df_row,
            by_col,
    ):
        """
        It calculates the number of the removed unique UMIs, observed on a given genomic position,
        after CC, Adjacency, Directional, MCL, MCL_val, or MCL_ed.

        Parameters
        ----------
        df_row
            object - a pandas-like df row
        by_col
            str - a column name in question

        Returns
        -------
            int - the sum of deduplicated unique UMI counts per position

        """
        return df_row['num_uniq_umis'] - len(df_row[by_col])

    def num_removed_reads(
            self,
            df_row,
            by_col,
    ):
        """
        It calculates the number of the removed reads in total.

        Parameters
        ----------
        df_row
            object - a pandas-like df row
        by_col
            str - a column name in question

        Returns
        -------
            int - the total counts of deduplicated reads per position

        """
        ### @@ set(df_row['uniq_repr_nodes'])
        # {0, 1, 2, ..., 1948}
        ### @@ len(set(df_row['uniq_repr_nodes']))
        # 1949
        ### @@ set(df_row[by_col])
        # {2}
        ### @@ len(set(df_row[by_col]))
        # 1
        diff_nodes = set(df_row['uniq_repr_nodes']) - set(df_row[by_col])
        ### @@ len(diff_nodes)
        # 1948

        if diff_nodes != set():
            # print(diff_nodes)
            umi_val_cnt_dict = df_row['vignette']['df_umi_uniq_val_cnt'].to_dict()
            ### @@ df_row['vignette']['df_umi_uniq_val_cnt'].to_dict()
            # {2: 55, 780: 51, 416: 47, ..., 1948: 1}
            return sum(umi_val_cnt_dict[node] for node in diff_nodes)
        else:
            return 0

    def umimax(
            self,
            df_row,
            by_col,
    ):
        """
        It returns ids of UMIs (i.e., representatives in their groupped UMIs) that has the
        highest count among all reads in their given genomic position.

        Examples
        --------
        {0: [0, 77, 81, ..., 1016], 1: [42, 46, 12], ..., 100: [2, 3, 5]}
        if no. 77 UMI has the highest count among all reads in 0, it will be added
        to umi_maxval_ids.

        Parameters
        ----------
        df_row
            a row of a pandas dataframe
        by_col
            a name of a pandas dataframe column

        Returns
        -------
        a list

        """
        umi_val_cnts = df_row['vignette']['df_umi_uniq_val_cnt']
        ### @@ umi_val_cnts
        # 2       55
        # 780     51
        # 416     47
        #         ..
        # 1948     1
        # Name: count, Length: 1949, dtype: int64
        umi_maxval_ids = []
        for k_c, nodes in df_row[by_col].items():
            self.console.print('key: {} nodes: {}'.format(k_c, nodes))
            ### @@ k_c, nodes
            # 0 [0, 77, 81, ..., 1016]
            ### @@ umi_val_cnts.loc[umi_val_cnts.index.isin(nodes)]
            # 2       55
            # 780     51
            # 416     47
            #         ..
            # 1948     1
            # Name: count, Length: 1949, dtype: int64
            ### @@ umi_val_cnts.loc[umi_val_cnts.index.isin(nodes)].idxmax()
            # 2
            ### @@ umi_val_cnts.loc[umi_val_cnts.index.isin(nodes)].max()
            # 55
            umi_max = umi_val_cnts.loc[umi_val_cnts.index.isin(nodes)].idxmax()
            umi_maxval_ids.append(umi_max)
        return umi_maxval_ids

    def length(
            self,
            df_val,
    ):
        """

        Parameters
        ----------
        df_val
            a pandas dataframe row

        Returns
        -------
            int - the length of the list

        """
        return len(df_val)

    def decompose(
            self,
            list_nd,
    ):
        """

        Parameters
        ----------
        x

        Returns
        -------

        """
        # print(list_nd)
        # print(len(list_nd))
        list_md = []
        for i in list_nd:
            list_md = list_md + i
        self.console.print('======># of the total reads left after deduplication: {}'.format(len(list_md)))
        # print(len(list_md))
        return list_md

    def bamids(
            self,
            df_row,
            by_col,
    ) -> list:
        """
        It outputs bamids of UMIs that are representative of all nodes in each group.

        Parameters
        ----------
        df_row
            a row of a pandas dataframe
        by_col
            a name of a pandas dataframe column

        Returns
        -------

        """
        uniq_umi_id_to_bam_id_dict = df_row['vignette']['uniq_umi_id_to_bam_id_dict']
        ### @@ uniq_umi_id_to_bam_id_dict
        # {0: 0, 1: 1, 2: 2, ..., 1948: 20609}
        ### @@ len(uniq_umi_id_to_bam_id_dict)
        # 1949
        list_1d = df_row[by_col]
        ### @@ list_1d
        # [2, 780, 416, ..., 1761]
        return [uniq_umi_id_to_bam_id_dict[node] for node in list_1d]

    def ed_ave(
            self,
            df_row,
            by_col,
    ):
        repr_nodes = df_row[by_col]
        # print(repr_nodes)
        node_len = len(repr_nodes)
        int_to_umi_dict = df_row['vignette']['int_to_umi_dict']
        if node_len != 1:
            ed_list = []
            for i in range(node_len):
                for j in range(i + 1, node_len):
                    ed_list.append(Hamming().general(
                        s1=int_to_umi_dict[repr_nodes[i]],
                        s2=int_to_umi_dict[repr_nodes[j]],
                    ))
            # print(pd.Series(ed_list).value_counts())
            return np.ceil(sum(ed_list) / (len(ed_list)))
        else:
            return -1

    def eds_(
            self,
            df_row,
            by_col,
    ):
        """"""
        # print(df_row.index)
        repr_nodes = df_row[by_col]
        int_to_umi_dict = df_row['vignette']['int_to_umi_dict']
        umi_val_cnts = df_row['vignette']['df_umi_uniq_val_cnt']
        # print(repr_nodes)
        # if len(repr_nodes) == len(np.unique(repr_nodes)):
        #     print(True)
        # else:
        #     print(False)
        node_len = len(repr_nodes)
        if node_len != 1:
            ed_list = []
            for i in range(node_len):
                for j in range(i + 1, node_len):
                    if repr_nodes[i] != repr_nodes[j]:
                        ed_list = ed_list + [Hamming().general(
                            int_to_umi_dict[repr_nodes[i]],
                            int_to_umi_dict[repr_nodes[j]])
                        ] * (umi_val_cnts.loc[repr_nodes[i]] * umi_val_cnts.loc[repr_nodes[j]])
            return round(sum(ed_list) / len(ed_list))
        else:
            return -1