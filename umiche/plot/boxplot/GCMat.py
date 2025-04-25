__version__ = "v1.0"
__copyright__ = "Copyright 2025"
__license__ = "GPL-3.0"


import numpy as np
import matplotlib.pyplot as plt
from umiche.plot.boxplot.Style import Style


class boxplot:

    def __init__(self, ):
        pass

    def draw1(
            self,
            data_dict,
    ):
        """
        ..  Description:
            ------------
            data_dict, which is a python dictionary, is subject to the data structure below:
            {
                {method_1: {val_dict_arr_1d_1: val_dict_containing_arr_1d_1}},
                {method_2: {val_dict_arr_1d_2: val_dict_containing_arr_1d_2}},
                {method_3: {val_dict_arr_1d_3: val_dict_containing_arr_1d_3}},
                ...
                {method_m: {val_dict_arr_1d_m: val_dict_containing_arr_1d_m}},
            }

            {
                'real': {
                    'val_dict_arr_1d':
                    {'CD8+ Cytotoxic T': array([31, 33, 17, ..., 30, 27, 23], dtype=int64),
                     'CD8+/CD45RA+ Naive Cytotoxic': array([31, 33, 17, ..., 27, 23, 21], dtype=int64),
                     'CD4+/CD45RO+ Memory': array([31, 33, 17, ..., 30, 23, 27], dtype=int64),
                     'CD19+ B': array([31, 33, 17, ..., 35, 27, 21], dtype=int64),
                     'CD4+/CD25 T Reg': array([31, 17, 33, ..., 30, 23, 27], dtype=int64),
                     'CD56+ NK': array([31, 17,  8, ..., 34,  8, 20], dtype=int64),
                     'CD4+ T Helper2': array([31, 33, 17, ..., 30, 27, 23], dtype=int64),
                     'CD4+/CD45RA+/CD25- Naive T': array([31, 33, 17, ..., 27, 21, 23], dtype=int64),
                     'CD34+': array([ 8,  0, 17, ..., 36,  0, 26], dtype=int64),
                     'Dendritic': array([31,  8,  3, ..., 26,  9, 35], dtype=int64),
                     'CD14+ Monocyte': array([ 3,  8, 31, ..., 35, 26, 12], dtype=int64)}
                }
                'method_2': {
                    ...
                }
            }

        ..  Example:
            --------

        :param data_dict: 2d dict
        :return: boxplot
        """
        key_1st = next(iter(data_dict))
        num_plot_pos = len(data_dict[key_1st]['val_dict_arr_1d'])
        num_plot_per_pos = len(data_dict)
        meanpointprops = {
            'marker': 'o',
            'markeredgecolor': 'grey',
            'markerfacecolor': 'grey',
            'markersize': 6
        }
        method_val_dict = {}
        method_pos_dict = {}
        offset_between_bplot_per_pos = 0.6
        for i_met, (met_key, v) in enumerate(data_dict.items()):
            method_pos_dict[met_key] = np.arange(num_plot_pos) * num_plot_per_pos - i_met * offset_between_bplot_per_pos
            method_val_dict[met_key] = [arr_1d for arr_1d in v['val_dict_arr_1d'].values()]
        print(method_val_dict)
        print(method_pos_dict)

        fig, ax = plt.subplots(nrows=1, ncols=1, figsize=(14, 7))
        bplot_handles = []
        # bplot_handles store methods with size equal to the number of the methods
        for val_key, pos_key in zip(method_val_dict, method_pos_dict):
             val_arr_2d = method_val_dict[val_key]
             pos_arr_1d = method_pos_dict[pos_key]
             bplot_handles.append(ax.boxplot(
                x=val_arr_2d,
                positions=pos_arr_1d,
                showmeans=True,
                meanprops=meanpointprops,
                patch_artist=True,
                showfliers=False, # hide outliers, the same as sym=''
            ))
        x_ticks = np.arange(0, num_plot_pos * num_plot_per_pos, num_plot_per_pos)
        x_tick_labels = [*data_dict[key_1st]['val_dict_arr_1d'].keys()]
        ax.set_xticks(x_ticks)
        ax.set_xticklabels(x_tick_labels, rotation=30, ha='right', fontsize=9)

        from matplotlib import rcParams
        rcParams['font.family'] = 'sans-serif'
        rcParams['font.sans-serif'] = ['Tahoma']
        ax.set_xlabel('Cell cluster', fontsize=18)
        ax.set_ylabel('Number of UMI counts', fontsize=18)
        ax.yaxis.grid(True)
        ax.spines['right'].set_color('none')
        ax.spines['top'].set_color('none')

        fig.subplots_adjust(
            bottom=0.26,
            # top=0.92, left=0.10, right=0.95, hspace=0.25, wspace=0.17
        )
        palette = [
            'magenta',
            'skyblue',
            'gray',
            'sienna',
            'crimson',
        ]
        Style().plain(boxplot_handles=bplot_handles, palette=palette)

        # add legend to the boxplot
        ax.legend(
            handles=[i["boxes"][0] for i in bplot_handles],
            labels=[*data_dict.keys()],
            bbox_to_anchor=(0.5, 1.15),
            loc=9,
            ncol=7,
        )
        # fig.tight_layout()
        plt.show()
        return 0

    def draw2(self, ):
        return

