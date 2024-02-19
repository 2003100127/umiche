__version__ = "v1.0"
__copyright__ = "Copyright 2024"
__license__ = "MIT"
__developer__ = "Jianfeng Sun"
__maintainer__ = "Jianfeng Sun"
__email__ = "jianfeng.sunmt@gmail.com"
__lab__ = "Cribbslab"

from typing import Dict

import numpy as np
import pandas as pd
import seaborn as sns
import matplotlib.pyplot as plt
from umiche.simu.Parameter import Parameter as params
from umiche.deduplicate.io.Stat import Stat as dedupstat
from umiche.plot.gadget.Transmitter import Transmitter as transmitter

from umiche.util.Reader import Reader as freader


class DedupMultipleTrimer:

    def __init__(
            self,
            scenarios: Dict,
            methods: Dict,
            umi_gt_cnt: int = 50,
            param_fpn: str = None,
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
            is_trans=False,
        )

        self.df_dedup = self.dedupstat.df_dedup

        print(self.df_dedup)
        sns.set(font="Helvetica")
        sns.set_style("ticks")

    def line(
            self,
            num_img_row=1,
            num_img_col=1,
    ):
        print(self.df_dedup)
        fig, ax = plt.subplots(nrows=num_img_row, ncols=num_img_col, figsize=(6.5, 5.5), sharey=False, sharex=False)
        palette = sns.color_palette("Paired")[2:]
        # palette = sns.color_palette("Set3")[3:]
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

                # ax[int(n / num_img_col), n % num_img_col].errorbar(
                ax.errorbar(
                    x=self.df_sce_met.index,
                    y=self.df_sce_met['mean'].values,
                    yerr=[self.df_sce_met['mean-min'], self.df_sce_met['max-mean']],
                    fmt='o',
                    # alpha=0.6,
                    # ecolor=palette[j],
                    color=palette[j],
                    linestyle='-',
                    linewidth=2,
                    elinewidth=1,
                    capsize=3,
                    markersize=6,
                    markeredgewidth=1.2,
                    markeredgecolor=palette[j],
                    markerfacecolor='none',
                    label=method_formal,
                )

                self.line_gadget(
                    ax=ax,
                    x=self.df_sce_met.index,
                    y=self.df_sce_met['mean'].values,
                    line_color=palette[j],
                    # label=method_formal,  # " ".join(ds_key.split("_"))
                    linewidth=2,
                    marker="o",
                    marker_size=6,
                    marker_edge_width=1.2,
                    marker_face_color='none',
                    decoration_mark=True,
                    xl_mark=xl_mark,
                    yl_mark=True if n % num_img_col == 0 else False,
                    # title='{}'.format(scenario_formal),
                    title_fs=16,
                    x_label='{}'.format(scenario_formal),  # 'Position'
                    y_label='Deduplicated homotrimer UMI count',
                    x_label_rotation=30,
                    x_label_rotation_align='right',
                    x_ticks=(np.arange(self.df_sce_met.index.shape[0])).tolist(),
                    x_ticklabels=self.df_sce_met['metric'].values.tolist(),
                    x_ticklabel_fs=12,
                    x_label_fs=18,
                    y_label_fs=18,
                    legend_fs=14,
                    legend_loc=None,
                )
                ax.spines['bottom'].set_linewidth(2)
                ax.spines['left'].set_linewidth(2)

        for i in range(num_img_col * num_img_row - (n + 1)):
            ax[num_img_row - 1, num_img_col - 1 - i].set_visible(False)
        plt.subplots_adjust(
            top=0.96,
            bottom=0.16,
            left=0.15,
            right=0.98,
            # hspace=0.40,
            # wspace=0.15,
        )
        plt.savefig('./solve_line.eps', format='eps')
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

    p = DedupMultipleTrimer(
        scenarios={
            # 'pcr_nums': 'PCR cycle',
            # 'pcr_errs': 'PCR error rate',
            'seq_errs': 'Sequencing error rate',
            # 'ampl_rates': 'Amplification rate',
            # 'umi_lens': 'UMI length',
            # 'seq_deps': 'Sequencing depth',
        },

        methods={
            'set_cover_solved_split_to_all': 'solved by spALL',
            'set_cover_solved_split_by_mv': 'solved by spMV',

            'set_cover_not_solved_split_to_all': 'not solved by spALL',
            'set_cover_not_solved_split_by_mv': 'not solved by spMV',
        },

        param_fpn=to('data/params_trimer.yml'),
    )

    p.line()
