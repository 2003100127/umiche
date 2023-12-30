__version__ = "v1.0"
__copyright__ = "Copyright 2024"
__license__ = "MIT"
__developer__ = "Jianfeng Sun"
__maintainer__ = "Jianfeng Sun"
__email__="jianfeng.sunmt@gmail.com"
__lab__ = "Cribbslab"

import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import seaborn as sns
from umiche.util.Reader import reader as greader
from umiche.simu.Parameter import Parameter as params


class DedupSingle:

    def __init__(
            self,
            df,
            df_melt,
            param_fpn,
    ):
        self.params = params(param_fpn=param_fpn)
        self.greader = greader()
        self.df = df
        self.df_melt = df_melt
        print(self.df)
        print(self.df_melt)

        self.df_direc = self.df[self.df['method'] == 'Directional']
        self.df_mcl_val = self.df[self.df['method'] == 'MCL-val']
        self.df_mcl_ed = self.df[self.df['method'] == 'MCL-ed']



    def jointplot(self, ):
        sns.set(font="Helvetica")
        sns.set_style("ticks")
        ppp = sns.jointplot(
            # x=self.df_mcl_ed['mean'].values,
            x=self.df_mcl_val['mean'].values,
            y=self.df_direc['mean'].values,
            kind="reg",
            color="crimson",
            label='asd',
        )
        ppp.ax_joint.plot([0, 8.5], [0, 8.5], 'grey', linewidth=2, alpha=1)
        # ppp.set_axis_labels('mcl_ed '+ r'($\frac{N_e-N_t}{N_t}$)' , 'directional ' + r'($\frac{N_e-N_t}{N_t}$)', fontsize=14)
        ppp.set_axis_labels('mcl_val '+ r'($\frac{N_e-N_t}{N_t}$)' , 'directional ' + r'($\frac{N_e-N_t}{N_t}$)', fontsize=14)
        ppp.ax_joint.text(8, 8.5, "confidence interval", horizontalalignment='right', size='medium', color='crimson',)
        ppp.ax_joint.text(3, 5.5, "regression", horizontalalignment='left', size='medium', color='crimson',)
        ppp.ax_joint.text(5.5, 5, "baseline", horizontalalignment='left', size='medium', color='black',)
        sns.despine(right=True, top=True)
        plt.tight_layout()
        plt.show()

    def jointgrid(self, ):
        sns.set(font="Helvetica")
        sns.set_style("ticks")
        # fig, ax = plt.subplots()
        # met = 'directional'
        # met = 'mcl_ed'
        met = 'mcl_val'

        ddd = pd.DataFrame()
        ddd['method'] = self.df_melt['method']
        ddd['Sequencing error'] = self.df_melt['Sequencing error']
        ddd[r'$\frac{N_e-N_t}{N_t}$'] = self.df_melt['value']

        g = sns.JointGrid(
            data=ddd[ddd['method'] == met],
            x="Sequencing error",
            y=r'$\frac{N_e-N_t}{N_t}$',
            marginal_ticks=True,
        )

        # Create an inset legend for the histogram colorbar
        cax = g.figure.add_axes([.20, .55, .02, .2])

        # Add the joint and marginal histogram plots
        g.plot_joint(
            sns.kdeplot, discrete=(True, False),
            cmap="light:#03012d", pmax=.8, cbar=True, cbar_ax=cax,
        )

        g.ax_joint.set_title(met, fontsize=14)
        g.ax_joint.set_xticks([int(x*100000) for x in self.seq_fix_errs[:]])
        g.ax_joint.set_xticklabels([1e-05, '', '', '', '', '', '', '', 0.001, 0.0025, 0.005, 0.0075, 0.01])
        plt.setp(g.ax_joint.get_xticklabels(), rotation=45)

        g.plot_marginals(
            sns.histplot,
            element="step",
            color="#03012d",
        )
        # ax.legend(ncol=2, loc="upper right", frameon=True)
        # g.set(ylabel="", xlabel="Automobile collisions per billion miles")
        # ax.set_ylabel('Sequencing error', fontsize=12)
        # ax.set_xlabel(r'$\frac{N_e-N_t}{N_t}$', fontsize=12)
        # sns.despine(left=True, bottom=True)
        plt.tight_layout()
        plt.show()
        return

    def strip(self, ):
        sns.set(font="Helvetica")
        sns.set_style("ticks")
        fig, ax = plt.subplots()
        sns.despine(bottom=True, left=True)
        cc = [
            'tab:green',
            'crimson',
            'tab:blue',
        ]
        # Show each observation with a scatterplot
        sns.stripplot(x="value", y="Sequencing error", hue='method', palette=cc,
                      data=self.df_melt, dodge=True, alpha=.25, zorder=1)
        # sns.pointplot(x="value", y="Sequencing error", hue='method',
        #               data=self.df_melt, dodge=.8 - .8 / 3,
        #               join=False,  palette=cc,
        #               markers="d", scale=.75, ci=None)

        # Improve the legend
        handles, labels = ax.get_legend_handles_labels()
        ax.legend(handles[3:], labels[3:], title='method',
                  handletextpad=0, columnspacing=1,
                  loc="upper right", ncol=3, frameon=True)
        ax.set_xlabel(r'$\frac{N_e-N_t}{N_t}$', fontsize=14)
        ax.set_ylabel('Sequencing error rate', fontsize=14)
        fig.subplots_adjust(
            top=0.98,
            bottom=0.15,
            left=0.15,
            right=0.98,
            # hspace=0.40,
            # wspace=0.15
        )
        plt.show()

    def stackedbar(self, ):
        fig, ax = plt.subplots(figsize=(4, 5))
        # sns.set(font="Verdana")
        sns.set(font="Helvetica")

        self.df_mcl_val["dmean"] = self.df_direc["mean"] - self.df_mcl_val["mean"]
        self.df_mcl_ed["dmean"] = self.df_direc["mean"] - self.df_mcl_ed["mean"]

        # self.df_mcl_val["dmean"] = np.exp(self.df_direc["mean"] - self.df_mcl_val["mean"])
        # self.df_mcl_ed["dmean"] = np.exp(self.df_direc["mean"] - self.df_mcl_ed["mean"])

        sns.set_color_codes("pastel")
        sns.barplot(
            x="dmean",
            y="metric",
            data=self.df_mcl_val,
            label="dFC_ed",
            color="b",
        )
        sns.set_color_codes("muted")
        sns.barplot(
            x="dmean",
            y="metric",
            data=self.df_mcl_ed,
            label="dFC_val",
            color="b",
        )

        # plt.fill_between(self.df_mcl_val["dmean"], self.df_mcl_val["metric"])
        ax.legend(ncol=2, loc="upper right", frameon=True)
        # ax.set(ylabel="", xlabel="Automobile collisions per billion miles")
        ax.set_ylabel('Sequencing error', fontsize=12)
        ax.set_xlabel(r'$\frac{N_e-N_t}{N_t}$', fontsize=12)
        sns.despine(left=True, bottom=True)
        fig.subplots_adjust(
            top=0.98,
            bottom=0.12,
            left=0.20,
            right=0.98,
            # hspace=0.40,
            # wspace=0.15
        )
        plt.show()

    def errorbar(self, ):
        sns.set(font="Helvetica")
        sns.set_style("ticks")
        # plt.rc('text', usetex=True)
        # plt.rc('font', **{'family': 'sans-serif', 'sans-serif': ['Helvetica']})
        fig, ax = plt.subplots()
        # ax.plot(
        #     file_n301.index,
        #     n301_means,
        #     label='Full TrainData dataset',
        #     alpha=0.9,
        #     linewidth=3.0,
        #     # s=,
        #     c='royalblue'
        # )
        cc = [
            'tab:green',
            'tab:blue',
            'crimson',
        ]

        df_gp = self.df.groupby(by=['method'])
        gp_keys = df_gp.groups.keys()

        for i, met in enumerate(gp_keys):
            df_met = df_gp.get_group(met)

            ax.errorbar(
                x=df_met.index,
                y=df_met['mean'],
                yerr=[df_met['mean-min'], df_met['max-mean']],
                fmt='o',
                alpha=0.6,
                ecolor=cc[i],
                color=cc[i],
                linewidth=1,
                elinewidth=1,
                capsize=2,
                markersize=6,
                label=met,
            )
            # print(df_met['metric'])
            ax.set_xticks(np.arange(df_met['metric'].shape[0]))
            ax.set_xticklabels(df_met['metric'], fontsize=8, rotation=45)
        # ax.set_xlabel('Sequencing error rate', fontsize=12)
        # ax.set_xlabel('Amplification rate', fontsize=12)
        # ax.set_xlabel('UMI length', fontsize=12)
        ax.set_xlabel('Polymerase error rate', fontsize=12)
        ax.set_ylabel(r'$\frac{N_e-N_t}{N_t}$', fontsize=16)
        # ax.set_title(DEFINE['title'], fontsize=12)
        sns.despine(right=True, top=True)
        fig.subplots_adjust(
            top=0.98,
            bottom=0.16,
            left=0.11,
            # left=0.12,
            right=0.95,
            # hspace=0.40,
            # wspace=0.15
        )
        plt.legend(fontsize=11, loc='upper left')
        plt.show()
        return

    def errorband(self, ):
        sns.set(font="Helvetica")
        sns.set_style("ticks")
        fig, ax = plt.subplots()
        cc = [
            'tab:green',
            'tab:blue',
            'crimson',
        ]
        df_gp = self.df.groupby(by=['method'])
        gp_keys = df_gp.groups.keys()

        for i, met in enumerate(gp_keys):
            df_met = df_gp.get_group(met)
            ax.errorbar(
                x=df_met.index,
                y=df_met['mean'],
                yerr=[df_met['mean-min'], df_met['max-mean']],
                fmt='o',
                alpha=0.7,
                ecolor=cc[i],
                color=cc[i],
                linestyle='-',
                linewidth=2,
                elinewidth=0.5,
                capsize=2,
                markersize=3,
                label=met,
            )
            ax.plot(df_met.index, df_met['mean'] - df_met['mean-min'], color=cc[i], linewidth=0.1,alpha=0.1)
            ax.plot(df_met.index, df_met['max-mean'] + df_met['mean'], color=cc[i], linewidth=0.1,alpha=0.1)
            ax.fill_between(
                df_met.index,
                df_met['mean'] - df_met['mean-min'],
                df_met['max-mean'] + df_met['mean'],
                alpha=0.1,
                color=cc[i],
            )
            ax.set_xticks(df_met['metric'].index)
            ax.set_xticklabels(df_met['metric'], fontsize=8, rotation=45)

        # ax.set_xlabel('Sequencing error rate', fontsize=12)
        ax.set_xlabel('Amplification rate', fontsize=12)
        ax.set_ylabel(r'$\frac{N_e-N_t}{N_t}$', fontsize=14)
        # ax.set_title(DEFINE['title'], fontsize=12)
        sns.despine(right=True, top=True)
        fig.subplots_adjust(
            top=0.98,
            bottom=0.16,
            left=0.13,
            right=0.95,
            # hspace=0.40,
            # wspace=0.15
        )
        plt.legend(fontsize=11)
        plt.show()
        return


if __name__ == "__main__":
    from umiche.path import to
    from umiche.deduplicate.io.Stat import Stat as dedupstat

    scenarios = {
        # 'pcr_nums': 'PCR cycle',
        # 'pcr_errs': 'PCR error',
        'seq_errs': 'Sequencing error',
        # 'ampl_rates': 'Amplification rate',
        # 'umi_lens': 'UMI length',
        # 'seq_deps': 'Sequencing depth',
    }
    methods = {
        # 'unique': 'Unique',
        # 'cluster': 'Cluster',
        # 'adjacency': 'Adjacency',
        'directional': 'Directional',
        # 'dbscan_seq_onehot': 'DBSCAN',
        # 'birch_seq_onehot': 'Birch',
        # 'aprop_seq_onehot': 'Affinity Propagation',
        # 'mcl': 'MCL',
        'mcl_val': 'MCL-val',
        'mcl_ed': 'MCL-ed',
    }
    dedupstat = dedupstat(
        scenarios=scenarios,
        methods=methods,
        param_fpn=to('data/params.yml'),
    )
    df_dedup = dedupstat.df_dedup
    df_dedup_perm_melt = dedupstat.df_dedup_perm_melt

    p = DedupSingle(
        df=df_dedup,
        df_melt=df_dedup_perm_melt,
        param_fpn=to('data/params.yml'),
    )
    # print(p.strip())
    # print(p.jointplot())
    print(p.jointgrid())
    # print(p.stackedbar())
    # print(p.errorbar())
    # print(p.errorband())