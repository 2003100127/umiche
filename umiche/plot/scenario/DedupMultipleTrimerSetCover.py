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
import scipy.stats as st
from umiche.simu.Parameter import Parameter as params
from umiche.deduplicate.io.Stat import Stat as dedupstat
from umiche.plot.gadget.Transmitter import Transmitter as transmitter

from umiche.util.Reader import Reader as freader


class DedupMultipleTrimerSetCover:

    def __init__(
            self,
            scenarios: Dict,
            methods: Dict,
            by: str = 'pn0',
            umi_gt_cnt: int = 50,
            param_fpn: str = None,
    ):
        self.scenarios = scenarios
        self.methods = methods
        self.by = by
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

        self.df_dedup = self.dedupstat.df_dedup_set_cover_len

        # print(self.df_dedup)
        sns.set(font="Helvetica")
        sns.set_style("ticks")

    def violinplot(
            self,
            num_img_row=1,
            num_img_col=1,
            start=0,
            end=-1,
    ):
        y_dict = {}
        for n, (scenario, scenario_formal) in enumerate(self.scenarios.items()):
            self.df_sce = self.df_dedup[self.df_dedup['scenario'] == scenario_formal]
            for j, (method, method_formal) in enumerate(self.methods.items()):
                self.df_sce_met = self.df_sce[self.df_sce['method'] == method_formal]
                y_dict[method_formal] = []
                self.df_sce_met[self.by].apply(lambda x: y_dict[method_formal].append([int(i) for i in x.split(';')]))
        x_labels = self.df_sce_met['metric'].values[start:end]

        fig, ax = plt.subplots(nrows=num_img_row, ncols=num_img_col, figsize=(17, 5.5), sharey=False, sharex=False)
        # palette = sns.color_palette("Paired")[2:]
        palette = sns.color_palette("Set3")[3:]
        BLACK = "#282724"
        GREY_DARK = "#747473"
        displacement = [0.4, -0.4]

        for i, (method_formal, y_2d) in enumerate(y_dict.items()):
            y_2d = y_2d[start:end]

            position_init = [t for t in np.arange(len(y_2d) * 2)][::2]
            positions = [t - displacement[i] for t in position_init]
            print(positions)

            x_2d = [[t] * len(d) for t, d in zip(positions, y_2d)]
            x_jittered = [x_2d + st.t(df=6, scale=0.04).rvs(len(x_2d)) for x_2d in x_2d]

            violins = ax.violinplot(
                y_2d,
                positions=positions,
                widths=0.8,
                bw_method="silverman",
                showmeans=False,
                showmedians=False,
                showextrema=False
            )

            for pc in violins["bodies"]:
                pc.set_facecolor("none")
                pc.set_edgecolor(BLACK)
                pc.set_linewidth(1.4)
                pc.set_alpha(1)

            medianprops = dict(
                linewidth=4,
                color=GREY_DARK,
                solid_capstyle="butt"
            )
            boxprops = dict(
                linewidth=2,
                color=GREY_DARK
            )

            ax.boxplot(
                y_2d,
                widths=0.4,
                positions=positions,
                showfliers=False,
                showcaps=False,
                medianprops=medianprops,
                whiskerprops=boxprops,
                boxprops=boxprops,
                # showmeans=True,
            )

            for u, (x, y) in enumerate(zip(x_jittered, y_2d)):
                ax.scatter(
                    x,
                    y,
                    s=20,
                    color=palette[i],
                    label=method_formal if u == 0 else None,
                    # facecolor='none',
                    alpha=0.4,
                )

        ax.spines['right'].set_color('none')
        ax.spines['top'].set_color('none')
        ax.spines['bottom'].set_linewidth(2)
        ax.spines['left'].set_linewidth(2)
        ax.set_xticks(position_init)
        ax.set_xticklabels(x_labels, fontsize=12)
        ax.set_xlabel('Sequencing error rate', fontsize=18)
        ax.set_ylabel('Count', fontsize=18)
        plt.rcParams["legend.markerscale"] = 2
        ax.legend(frameon=False, fontsize=14)

        plt.subplots_adjust(
            top=0.96,
            bottom=0.12,
            left=0.05,
            right=0.98,
            # hspace=0.40,
            # wspace=0.15,
        )
        # plt.savefig('./solve_line.eps', format='eps')
        plt.show()
        return

    def violinplot_sgl(
            self,
            num_img_row=1,
            num_img_col=1,
            start=0,
            end=-1,
    ):
        y_dict = {}
        for n, (scenario, scenario_formal) in enumerate(self.scenarios.items()):
            self.df_sce = self.df_dedup[self.df_dedup['scenario'] == scenario_formal]
            for j, (method, method_formal) in enumerate(self.methods.items()):
                self.df_sce_met = self.df_sce[self.df_sce['method'] == method_formal]
                y_dict[method_formal] = []
                self.df_sce_met[self.by].apply(lambda x: y_dict[method_formal].append([int(i) for i in x.split(';')]))
        x_labels = self.df_sce_met['metric'].values[start:end]

        fig, ax = plt.subplots(nrows=num_img_row, ncols=num_img_col, figsize=(6, 5), sharey=False, sharex=False)
        # palette = sns.color_palette("Paired")[2:]
        palette = sns.color_palette("Set3")[3:]
        BLACK = "#282724"
        GREY_DARK = "#747473"
        displacement = [0.4, -0.4]

        for i, (method_formal, y_2d) in enumerate(y_dict.items()):
            y_2d = y_2d[start:end]

            position_init = [t for t in np.arange(len(y_2d) * 2)][::2]
            positions = [t - displacement[i] for t in position_init]
            print(positions)

            x_2d = [[t] * len(d) for t, d in zip(positions, y_2d)]
            x_jittered = [x_2d + st.t(df=6, scale=0.04).rvs(len(x_2d)) for x_2d in x_2d]

            violins = ax.violinplot(
                y_2d,
                positions=positions,
                widths=0.4,
                bw_method="silverman",
                showmeans=False,
                showmedians=False,
                showextrema=False
            )

            for pc in violins["bodies"]:
                pc.set_facecolor("none")
                pc.set_edgecolor(BLACK)
                pc.set_linewidth(1.4)
                pc.set_alpha(1)

            medianprops = dict(
                linewidth=4,
                color=GREY_DARK,
                solid_capstyle="butt"
            )
            boxprops = dict(
                linewidth=2,
                color=GREY_DARK
            )

            ax.boxplot(
                y_2d,
                widths=0.2,
                positions=positions,
                showfliers=False,
                showcaps=False,
                medianprops=medianprops,
                whiskerprops=boxprops,
                boxprops=boxprops,
                # showmeans=True,
            )

            for u, (x, y) in enumerate(zip(x_jittered, y_2d)):
                ax.scatter(
                    x,
                    y,
                    s=20,
                    color=palette[i],
                    label=method_formal if u == 0 else None,
                    # facecolor='none',
                    alpha=0.4,
                )

        ax.spines['right'].set_color('none')
        ax.spines['top'].set_color('none')
        ax.spines['bottom'].set_linewidth(2)
        ax.spines['left'].set_linewidth(2)
        ax.set_xticks(position_init)
        ax.set_xticklabels(x_labels, fontsize=12)
        ax.set_xlabel('Sequencing error rate', fontsize=18)
        ax.set_ylabel('Count', fontsize=18)
        plt.rcParams["legend.markerscale"] = 2
        ax.legend(frameon=False, fontsize=14, loc='lower center', bbox_to_anchor=(0.5, 1.0))

        plt.subplots_adjust(
            top=0.82,
            bottom=0.12,
            left=0.12,
            right=0.98,

            # top=0.82,
            # bottom=0.135,
            # left=0.15,
            # right=0.98
        )
        # plt.savefig('./solve_line.eps', format='eps')
        plt.show()
        return


if __name__ == "__main__":
    from umiche.path import to

    p = DedupMultipleTrimerSetCover(
        scenarios={
            # 'pcr_nums': 'PCR cycle',
            # 'pcr_errs': 'PCR error rate',
            'seq_errs': 'Sequencing error rate',
            # 'ampl_rates': 'Amplification rate',
            # 'umi_lens': 'UMI length',
            # 'seq_deps': 'Sequencing depth',
        },

        methods={
            # 'set_cover_mono_len_split_by_mv': 'monomer UMIs by spMV',
            # 'set_cover_mono_len_split_to_all': 'monomer UMIs by spALL',
            'set_cover_multi_len_split_by_mv': 'homotrimer UMIs by spMV',
            'set_cover_multi_len_split_to_all': 'homotrimer UMIs by spALL',
        },

        param_fpn=to('data/params_trimer.yml'),
    )

    # p.violinplot(
    #     start=0,
    #     end=-1,
    #
    #     # start=10,
    #     # end=-1,
    #
    #     # start=0,
    #     # end=10,
    # )

    p.violinplot_sgl(
        start=15,
        end=-1,
    )
