__version__ = "v1.0"
__copyright__ = "Copyright 2024"
__license__ = "MIT"
__developer__ = "Jianfeng Sun"
__maintainer__ = "Jianfeng Sun"
__email__="jianfeng.sunmt@gmail.com"
__lab__ = "Cribbslab"

import time
import pandas as pd
from umiche.util.Console import Console


class Trace:

    def __init__(
            self,
            df_umi_uniq_val_cnt,
            umi_id_to_origin_id_dict,
            verbose=False,
    ):
        self.df_umi_uniq_val_cnt = df_umi_uniq_val_cnt
        self.umi_id_to_origin_id_dict = umi_id_to_origin_id_dict

        self.verbose = verbose
        self.console = Console()
        self.console.verbose = self.verbose

    def edge_class(
            self,
            series_2d_arr,
            sort='cnt',
    ):
        """

        Parameters
        ----------
        series_2d_arr
        sort

        Returns
        -------

        """
        stime = time.time()
        self.console.print('==================>edges to be traced is {}'.format(series_2d_arr.shape[0]))
        ### @@ series_2d_arr
        # cc_0     [[0, 65], [65, 55], [65, 188], [65, 190], [0, ...
        # cc_1     [[36, 1], [36, 46], [46, 183], [36, 103], [103...
        # cc_2     [[33, 2], [33, 30], [30, 215], [33, 73], [73, ...
        # cc_3                          [[3, 69], [3, 91], [3, 247]]
        # ...
        # cc_47    [[59, 106], [59, 133], [59, 180], [59, 194], [...
        # cc_48                    [[63, 101], [63, 219], [63, 238]]
        series_origin = series_2d_arr.apply(
            lambda x: self.is_same_origin(x),
        )
        ### @@ series_origin
        # cc_0        [1, 1, 1, 1, 1, 1, 1, 1]
        # cc_1     [1, 1, 1, 1, 1, 1, 1, 1, 1]
        # cc_2              [1, 1, 1, 1, 1, 1]
        # ...
        # cc_47                [1, 1, 1, 1, 1]
        # cc_48                      [1, 1, 1]
        if sort == 'cnt':
            total_len = series_origin.apply(lambda x: len(x))
            series_diff_cnt = series_origin.apply(lambda x: sum(1 for ele in x if ele == 0))
            series_same_cnt = series_origin.apply(lambda x: sum(1 for ele in x if ele == 1))
            cnt_total = total_len.sum()
            cnt_diff_origin = series_diff_cnt.sum()
            cnt_same_origin = series_same_cnt.sum()
            self.console.print('==================>trace edge cls time {time:.2f}s'.format(time=time.time()-stime))
            return {
                'diff_origin': cnt_diff_origin,
                'same_origin': cnt_same_origin,
                'total': cnt_total,
            }
        elif sort == 'pct':
            series_diff_pct = series_origin.apply(lambda x: sum(1 for ele in x if ele == 0) / len(x) if len(x) != 0 else 0)
            series_same_pct = series_origin.apply(lambda x: sum(1 for ele in x if ele == 1) / len(x) if len(x) != 0 else 0)
            # print(series_diff_pct)
            # print(series_same_pct)
            pct_diff_origin = series_diff_pct.mean()
            pct_same_origin = series_same_pct.mean()
            self.console.print('==================>trace edge cls time {time:.2f}s'.format(time=time.time()-stime))
            return {
                'diff_origin': pct_diff_origin,
                'same_origin': pct_same_origin,
                'total': 1,
            }

    def is_same_origin(self, arr_2d):
        """

        Parameters
        ----------
        arr_2d
            A 2d edge list, one with 2-sized vector consisting of 2 nodes.
            e.g. [[0, 65], [65, 55], [65, 188], [65, 190], [0, 76], [0, 162], [0, 237], [0, 256]]

        Returns
        -------

        """
        trace_marks = []
        for nodes in arr_2d:
            self.console.print('==================>2 nodes in the edge are {} and {}'.format(nodes[0], nodes[1]))
            umi_ori_node_1 = self.umi_id_to_origin_id_dict[nodes[0]]
            umi_ori_node_2 = self.umi_id_to_origin_id_dict[nodes[1]]
            intxn = set([umi_ori_node_1]) & set([umi_ori_node_2])
            self.console.print('==================>origins of the 2 nodes in the edge are {} and {}'.format(
                umi_ori_node_1,
                umi_ori_node_2,
            ))
            if len(intxn) != 0:
                trace_marks.append(1)
            else:
                trace_marks.append(0)
        return trace_marks

    def matchRepresentative(self, df):
        print('======>representative to be traced is {}'.format(df.shape[0]))
        stime = time.time()
        res = df.apply(
            lambda arr_2d: self.maxval(arr_2d),
        )
        # print(res)
        ttt = set()
        total1 = res.values
        for i in total1:
            ttt.update(i)
            # ttt = ttt + i
        total = res.apply(lambda x: len(x))
        # print(total)
        total_len = total.sum()
        # print(total_len)
        print('======>trace representative time {time:.2f}s'.format(time=time.time() - stime))
        return len(ttt)/50

    def maxval(self, arr_2d):
        repr_max_nodes = []
        # repr_max_nodes = set()
        for sub_cc_nodes in arr_2d:
            node_val_sorted = self.df_umi_uniq_val_cnt.loc[
                self.df_umi_uniq_val_cnt.index.isin(sub_cc_nodes)
            ].sort_values(ascending=False).to_dict()
            max_node = [*node_val_sorted.keys()][0]
            umi_ori_max_nodes = self.umi_id_to_origin_id_dict[max_node]
            if len(umi_ori_max_nodes) != 1:
                repr_max_nodes = repr_max_nodes + []
            else:
                repr_max_nodes = repr_max_nodes + umi_ori_max_nodes
        # print(len(repr_max_nodes))
        return repr_max_nodes

    def format_apv_disapv(
            self,
            cc_dict,
    ):
        """
        Format the apv or disapv dict, which has the following cc_dict format

        Parameters
        ----------
        cc_dict
            {'cc_0': {
              'node_0': [[0, 65], [65, 55], [65, 188], [65, 190], [0, 76], [0, 162], [0, 237], [0, 256]]
              },
            'cc_1': {
              'node_36': [[36, 1], [36, 46], [46, 183], [36, 103], [103, 216], [36, 108], [36, 121], [121, 142], [36, 197]]
              },
            'cc_47': {
            'node_59': [[59, 106], [59, 133], [59, 180], [59, 194], [59, 204]]
              },
            'cc_48': {
              'node_63': [[63, 101], [63, 219], [63, 238]]
              }
            }

        Returns
        -------

        """
        ### @@ cc_dict
        # {'cc_0': {
        #   'node_0': [[0, 65], [65, 55], [65, 188], [65, 190], [0, 76], [0, 162], [0, 237], [0, 256]]
        #   },
        # 'cc_1': {
        #   'node_36': [[36, 1], [36, 46], [46, 183], [36, 103], [103, 216], [36, 108], [36, 121], [121, 142], [36, 197]]
        #   },
        # 'cc_47': {
        # 'node_59': [[59, 106], [59, 133], [59, 180], [59, 194], [59, 204]]
        #   },
        # 'cc_48': {
        #   'node_63': [[63, 101], [63, 219], [63, 238]]
        #   }
        # }
        cc_series = pd.Series(cc_dict)
        ### @@ cc_series
        # cc_0     {'node_0': [[0, 65], [65, 55], [65, 188], [65,...
        # cc_1     {'node_36': [[36, 1], [36, 46], [46, 183], [36...
        # cc_2     {'node_33': [[33, 2], [33, 30], [30, 215], [33...
        # cc_3              {'node_3': [[3, 69], [3, 91], [3, 247]]}
        # ...
        # cc_47    {'node_59': [[59, 106], [59, 133], [59, 180], ...
        # cc_48       {'node_63': [[63, 101], [63, 219], [63, 238]]}
        cc_decomposed_2d_series  = cc_series.apply(lambda x: self.dictTo2d(x))
        ### @@ cc_decomposed_2d_series
        # cc_0     [[0, 65], [65, 55], [65, 188], [65, 190], [0, ...
        # cc_1     [[36, 1], [36, 46], [46, 183], [36, 103], [103...
        # cc_2     [[33, 2], [33, 30], [30, 215], [33, 73], [73, ...
        # ...
        # cc_47    [[59, 106], [59, 133], [59, 180], [59, 194], [...
        # cc_48                    [[63, 101], [63, 219], [63, 238]]
        return cc_decomposed_2d_series


    def dictTo2d(
            self,
            x,
    ):
        """

        Parameters
        ----------
        x
            {'node_0': [[0, 65], [65, 55], [65, 188], [65, 190], [0, 76], [0, 162], [0, 237], [0, 256]]}

        Returns
        -------
            It returns a pandas series with rows, each being represented by a 2D list.
            within which 2D lists of different nodes are decomposed into 1d list.

        """
        if isinstance(x , dict):
            node_list_3d = [*x.values()]
        else:
            node_list_3d = []
        ### @@ node_list_3d
        # [[[0, 65], [65, 55], [65, 188], [65, 190], [0, 76], [0, 162], [0, 237], [0, 256]]]
        # [[[36, 1], [36, 46], [46, 183], [36, 103], [103, 216], [36, 108], [36, 121], [121, 142], [36, 197]]]
        # [[[33, 2], [33, 30], [30, 215], [33, 73], [73, 119], [33, 127]]]
        # ...
        # [[[59, 106], [59, 133], [59, 180], [59, 194], [59, 204]]]
        # [[[63, 101], [63, 219], [63, 238]]]
        res_2d = []
        for i in node_list_3d:
            res_2d = res_2d + i
        return res_2d

    def formatCCS(
            self,
            cc_dict,
    ):
        """

        Parameters
        ----------
        cc_dict
            {'cc_0': {
              'node_0': [[0, 65], [65, 55], [65, 188], [65, 190], [0, 76], [0, 162], [0, 237], [0, 256]]
              },
            'cc_1': {
              'node_36': [[36, 1], [36, 46], [46, 183], [36, 103], [103, 216], [36, 108], [36, 121], [121, 142], [36, 197]]
              },
            'cc_47': {
            'node_59': [[59, 106], [59, 133], [59, 180], [59, 194], [59, 204]]
              },
            'cc_48': {
              'node_63': [[63, 101], [63, 219], [63, 238]]
              }
            }

        Returns
        -------
            It returns a pandas series with rows, each being represented by a 3D list.
            into which 2D lists of different nodes are merged.

        """
        cc_series = pd.Series(cc_dict)
        ### @@ cc_series
        # cc_0     {'node_0': [[0, 65], [65, 55], [65, 188], [65,...
        # cc_1     {'node_36': [[36, 1], [36, 46], [46, 183], [36...
        # cc_2     {'node_33': [[33, 2], [33, 30], [30, 215], [33...
        # cc_3              {'node_3': [[3, 69], [3, 91], [3, 247]]}
        # ...
        # cc_47    {'node_59': [[59, 106], [59, 133], [59, 180], ...
        # cc_48       {'node_63': [[63, 101], [63, 219], [63, 238]]}
        cc_decomposed_3d_series = cc_series.apply(lambda x: [*x.values()])
        ### @@ cc_series.apply(lambda x: [*x.values()])
        # cc_0     [[[0, 65], [65, 55], [65, 188], [65, 190], [0,...
        # cc_1     [[[36, 1], [36, 46], [46, 183], [36, 103], [10...
        # cc_2     [[[33, 2], [33, 30], [30, 215], [33, 73], [73,...
        # ...
        # cc_47    [[[59, 106], [59, 133], [59, 180], [59, 194], ...
        # cc_48                  [[[63, 101], [63, 219], [63, 238]]]
        return cc_decomposed_3d_series