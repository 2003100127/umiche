__version__ = "v1.0"
__copyright__ = "Copyright 2025"
__license__ = "GPL-3.0"
__developer__ = "Jianfeng Sun"
__maintainer__ = "Jianfeng Sun"
__email__="jianfeng.sunmt@gmail.com"


from typing import Dict

import numpy as np
import pandas as pd
import seaborn as sns
import matplotlib.pyplot as plt
from umiche.plot.gadget.Transmitter import Transmitter as transmitter


class TraceMultiple:

    def __init__(
            self,
            df_apv: pd.DataFrame,
            df_disapv: pd.DataFrame,
            scenarios,
            methods,
    ):
        self.df_apv = df_apv
        self.df_disapv = df_disapv
        self.scenarios = scenarios
        self.methods = methods

        sns.set(font="Helvetica")
        sns.set_style("ticks")

    def line(
            self,
            num_img_row=2,
            num_img_col=2,
    ):
        print(self.df_apv)
        print(self.df_disapv)
        fig, ax = plt.subplots(nrows=num_img_row, ncols=num_img_col, figsize=(12, 9), sharey=False, sharex=True)
        palette = sns.color_palette("Paired")

        dict_df = {
            "same origin (merged)": self.df_apv[['same_origin', 'metric', 'method']].rename(columns={"same_origin": 'val'}),
            "different origin (merged)": self.df_apv[['diff_origin', 'metric', 'method']].rename(columns={"diff_origin": 'val'}),
            "same origin (not merged)": self.df_disapv[['same_origin', 'metric', 'method']].rename(columns={"same_origin": 'val'}),
            "different origin (not merged)": self.df_disapv[['diff_origin', 'metric', 'method']].rename(columns={"diff_origin": 'val'}),
        }

        for n, (situation, df) in enumerate(dict_df.items()):

            for j, (method, method_formal) in enumerate(self.methods.items()):
                print(df)
                df_met = df[df['method'] == method_formal]

                if int(n / num_img_col) == num_img_row - 1:
                    xl_mark = True
                elif (int(n / num_img_col) == num_img_row - 2) and (
                        n % num_img_col >= (num_img_col - num_img_col * num_img_row + len([1]))):
                    xl_mark = True
                else:
                    xl_mark = False
                self.line_gadget(
                    ax=ax[int(n / num_img_col), n % num_img_col],
                    x=df_met.index,
                    y=df_met.val.values,
                    line_color=palette[j+1],
                    label=method_formal,  # " ".join(ds_key.split("_"))
                    linewidth=2,
                    marker="o",
                    marker_size=6,
                    marker_edge_width=1.2,
                    marker_face_color='none',
                    decoration_mark=True,
                    xl_mark=xl_mark,
                    yl_mark=True if n % num_img_col == 0 else False,
                    title='{}'.format(situation),
                    title_fs=16,
                    x_label='PCR cycle', # 'PCR cycle'
                    y_label="Percentage of UMI counts",
                    x_label_rotation=0,
                    x_label_rotation_align='right',
                    x_ticks=(np.arange(df_met.index.shape[0])).tolist(),
                    x_ticklabels=df_met['metric'].values.tolist(),
                    x_ticklabel_fs=12,
                    x_label_fs=14,
                    y_label_fs=18,
                    legend_fs=11,
                    legend_loc=None,
                )
        for i in range(num_img_col * num_img_row - (n + 1)):
            ax[num_img_row - 1, num_img_col - 1 - i].set_visible(False)
        plt.subplots_adjust(
            top=0.96,
            bottom=0.06,
            left=0.08,
            right=0.98,
            # hspace=0.40,
            # wspace=0.15,
        )
        plt.show()

    def line_apv(
            self,
            num_img_row=2,
            num_img_col=2,
    ):
        print(self.df_apv)
        fig, ax = plt.subplots(nrows=num_img_row, ncols=num_img_col, figsize=(12, 9), sharey=False, sharex=False)
        palette = sns.color_palette("Paired")

        dict_df = {
            "same origin (merged)": self.df_apv[['same_origin', 'metric', 'method']].rename(columns={"same_origin": 'val'}),
            "different origin (merged)": self.df_apv[['diff_origin', 'metric', 'method']].rename(columns={"diff_origin": 'val'}),
        }

        for n, (situation, df) in enumerate(dict_df.items()):

            for j, (method, method_formal) in enumerate(self.methods.items()):
                print(df)
                df_met = df[df['method'] == method_formal]

                if int(n / num_img_col) == num_img_row - 1:
                    xl_mark = True
                elif (int(n / num_img_col) == num_img_row - 2) and (
                        n % num_img_col >= (num_img_col - num_img_col * num_img_row + len([1]))):
                    xl_mark = True
                else:
                    xl_mark = False
                self.line_gadget(
                    ax=ax[int(n / num_img_col), n % num_img_col],
                    x=df_met.index,
                    y=df_met.val.values,
                    line_color=palette[j+1],
                    label=method_formal,  # " ".join(ds_key.split("_"))
                    linewidth=2,
                    marker="o",
                    marker_size=6,
                    marker_edge_width=1.2,
                    marker_face_color='none',
                    decoration_mark=True,
                    xl_mark=xl_mark,
                    yl_mark=True if n % num_img_col == 0 else False,
                    title='{}'.format(situation),
                    title_fs=16,
                    x_label='PCR cycle', # 'PCR cycle'
                    # y_label="Percentage of UMI counts",
                    y_label="UMI count",
                    x_label_rotation=0,
                    x_label_rotation_align='right',
                    x_ticks=(np.arange(df_met.index.shape[0])).tolist(),
                    x_ticklabels=df_met['metric'].values.tolist(),
                    x_ticklabel_fs=12,
                    x_label_fs=14,
                    y_label_fs=18,
                    legend_fs=11,
                    legend_loc=None,
                )
        for i in range(num_img_col * num_img_row - (n + 1)):
            ax[num_img_row - 1, num_img_col - 1 - i].set_visible(False)
        plt.subplots_adjust(
            top=0.96,
            bottom=0.06,
            left=0.08,
            right=0.98,
            # hspace=0.40,
            # wspace=0.15,
        )
        plt.show()

    @transmitter(type="line", task=None)
    def line_gadget(*args, **kwargs):
        return kwargs

    @transmitter(type="line_scatter", task=None)
    def line_scatter_gadget(*args, **kwargs):
        return kwargs


if __name__ == "__main__":
    from umiche.path import to

    from umiche.deduplicate.io.Stat import Stat as dedupstat

    scenarios = {
        'pcr_nums': 'PCR cycle',
        # 'pcr_errs': 'PCR error',
        # 'seq_errs': 'Sequencing error',
        # 'ampl_rates': 'Amplification rate',
        # 'umi_lens': 'UMI length',
        # 'seq_deps': 'Sequencing depth',
    }
    methods = {
        'unique': 'Unique',
        'cluster': 'Cluster',
        'adjacency': 'Adjacency',
        'directional': 'Directional',
        'dbscan_seq_onehot': 'DBSCAN',
        'birch_seq_onehot': 'Birch',
        'aprop_seq_onehot': 'Affinity Propagation',
        'mcl': 'MCL',
        'mcl_val': 'MCL-val',
        'mcl_ed': 'MCL-ed',
    }

    dedupstat = dedupstat(
        scenarios=scenarios,
        methods=methods,
        param_fpn=to('data/params.yml'),
    )

    p = TraceMultiple(
        df_apv=dedupstat.df_trace_cnt['apv'],
        df_disapv=dedupstat.df_trace_cnt['disapv'],
        scenarios=scenarios,
        methods=methods,
    )

    # p.line()
    p.line_apv()
