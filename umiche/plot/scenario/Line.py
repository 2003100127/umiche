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

from umiche.util.Reader import Reader as freader


class Line:

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

        self.params = params(param_fpn=self.param_fpn)

        self.freader = freader()

        self.dedupstat = dedupstat(
            scenarios=self.scenarios,
            methods=self.methods,
            param_fpn=self.param_fpn,
        )

        self.df_inflat, self.df_exp = self.dedupstat.df_inflat_exp

        # sns.set(font="Verdana")
        sns.set(font="Helvetica")
        sns.set_style("ticks")

    def draw(
            self,
            num_img_row=1,
            num_img_col=2,
    ):
        fig, ax = plt.subplots(nrows=num_img_row, ncols=num_img_col, figsize=(10, 5), sharey=False, sharex=False)
        palette = sns.color_palette("Set2") # Paired
        print(self.df_inflat)
        print(self.df_exp)
        dict_df = {
            'Inflation': self.df_inflat,
            'Expansion': self.df_exp,
        }
        scenarios = self.df_inflat.columns
        for i, (key, df) in enumerate(dict_df.items()):
            for j, scenario in enumerate(scenarios):
                ax[i].plot(
                    df.index.tolist(),
                    df[scenario].values.tolist(),
                    label=scenario,
                    linewidth=2.5,
                    marker='o',
                    markerfacecolor='white',
                    alpha=0.9,
                    color=palette[j],
                )
            ax[i].legend(title='', ncol=1, fontsize=12, )
            ax[i].set_xlabel(key, fontsize=14)
            ax[i].set_ylabel(r'($\frac{N_e-N_t}{N_t}$)', fontsize=18)
            ax[i].spines['right'].set_visible(False)
            ax[i].spines['top'].set_visible(False)
        ax[0].axvline(x=3.5, color='grey', linestyle='dotted')
        ax[1].axvline(x=4, color='grey', linestyle='dotted')

        plt.subplots_adjust(
            top=0.96,
            bottom=0.12,
            left=0.09,
            right=0.98,
            # hspace=0.40,
            # wspace=0.15,
        )
        plt.show()


if __name__ == "__main__":
    from umiche.path import to

    p = Line(
        scenarios={
            'pcr_nums': 'PCR cycle',
            'pcr_errs': 'PCR error',
            'seq_errs': 'Sequencing error',
            'ampl_rates': 'Amplification rate',
            'umi_lens': 'UMI length',
            'seq_deps': 'Sequencing depth',
        },

        methods={
            'mcl': 'MCL',
        },

        param_fpn=to('data/params.yml'),
    )

    p.draw()