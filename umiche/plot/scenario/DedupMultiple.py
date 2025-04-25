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
from umiche.simu.Parameter import Parameter as params
from umiche.deduplicate.io.Stat import Stat as dedupstat
from umiche.plot.gadget.Transmitter import Transmitter as transmitter

from umiche.util.Reader import Reader as freader
from umiche.util.Console import Console


class DedupMultiple:

    def __init__(
            self,
            scenarios: Dict,
            methods: Dict,
            umi_gt_cnt: int = 50,
            param_fpn : str = None,
            verbose : bool = False,
    ):
        self.scenarios = scenarios
        self.methods = methods
        self.umi_gt_cnt = umi_gt_cnt
        self.param_fpn = param_fpn
        self.freader = freader()
        self.params = params(param_fpn=self.param_fpn)

        self.dedupstat = dedupstat(
            scenarios=self.scenarios,
            methods=self.methods,
            param_fpn=self.param_fpn,
        )

        self.df_dedup = self.dedupstat.df_dedup

        sns.set(font="Helvetica")
        sns.set_style("ticks")

        self.console = Console()
        self.console.verbose = verbose

    def line(
            self,
            num_img_row=2,
            num_img_col=3,
    ):
        print(self.df_dedup)
        fig, ax = plt.subplots(nrows=num_img_row, ncols=num_img_col, figsize=(17, 9), sharey=False, sharex=False)
        palette = sns.color_palette("Paired")
        for n, (scenario, scenario_formal) in enumerate(self.scenarios.items()):
            self.df_sce = self.df_dedup[self.df_dedup['scenario'] == scenario_formal]
            for j, (method, method_formal) in enumerate(self.methods.items()):
                self.df_sce_met = self.df_sce[self.df_sce['method'] == method_formal]

                if int(n / num_img_col) == num_img_row - 1:
                    xl_mark = True
                elif (int(n / num_img_col) == num_img_row - 2) and (
                        n % num_img_col >= (num_img_col - num_img_col * num_img_row + len([1]))):
                    xl_mark = True
                else:
                    xl_mark = False

                self.line_gadget(
                    ax=ax[int(n / num_img_col), n % num_img_col],
                    x=self.df_sce_met.index,
                    y=self.df_sce_met['mean'].values,
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
                    title='{}'.format(scenario_formal),
                    title_fs=16,
                    x_label=None, # 'Position'
                    y_label=r'($\frac{N_e-N_t}{N_t}$)',
                    x_label_rotation=20,
                    x_label_rotation_align='right',
                    x_ticks=(np.arange(self.df_sce_met.index.shape[0])).tolist(),
                    x_ticklabels=self.df_sce_met['metric'].values.tolist(),
                    x_ticklabel_fs=8,
                    x_label_fs=11,
                    y_label_fs=18,
                    legend_fs=11,
                    legend_loc=None,
                )
        for i in range(num_img_col * num_img_row - (n + 1)):
            ax[num_img_row - 1, num_img_col - 1 - i].set_visible(False)
        plt.subplots_adjust(
            top=0.96,
            bottom=0.06,
            left=0.06,
            right=0.98,
            # hspace=0.40,
            # wspace=0.15,
        )
        plt.show()
        return

    @transmitter(type="line", task=None)
    def line_gadget(*args, **kwargs):
        return kwargs

    @transmitter(type="line_scatter", task=None)
    def line_scatter_gadget(*args, **kwargs):
        return kwargs


if __name__ == "__main__":
    from umiche.path import to

    p = DedupMultiple(
        scenarios={
            'pcr_nums': 'PCR cycle',
            'pcr_errs': 'PCR error',
            'seq_errs': 'Sequencing error',
            'ampl_rates': 'Amplification rate',
            'umi_lens': 'UMI length',
            'seq_deps': 'Sequencing depth',
        },

        methods={
            # 'unique': 'Unique',
            # 'cluster': 'Cluster',
            # 'adjacency': 'Adjacency',
            'directional': 'Directional',
            'dbscan_seq_onehot': 'DBSCAN',
            'birch_seq_onehot': 'Birch',
            # 'aprop_seq_onehot': 'Affinity Propagation',
            'mcl': 'MCL',
            'mcl_val': 'MCL-val',
            'mcl_ed': 'MCL-ed',
        },

        param_fpn=to('data/params.yml'),
    )

    p.line()