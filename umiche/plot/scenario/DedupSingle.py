__version__ = "v1.0"
__copyright__ = "Copyright 2025"
__license__ = "GPL-3.0"
__developer__ = "Jianfeng Sun"
__maintainer__ = "Jianfeng Sun"
__email__="jianfeng.sunmt@gmail.com"


import numpy as np
import pandas as pd
import seaborn as sns
import matplotlib.pyplot as plt
from umiche.util.Reader import Reader as freader


class DedupSingle:

    def __init__(
            self,
            df,
            df_melt,
    ):
        # sns.set(font="Verdana")
        sns.set(font="Helvetica")
        sns.set_style("ticks")

        self.freader = freader()
        self.df = df
        self.df_melt = df_melt
        print(self.df)
        print(self.df_melt)

    def jointplot(
            self,
            method,
    ):
        ppp = sns.jointplot(
            x=self.df[self.df['method'] == method]['mean'].values,
            y=self.df[self.df['method'] == 'Directional']['mean'].values,
            kind="reg",
            color="crimson",
            label='asd',
        )
        ppp.ax_joint.plot([0, 8.5], [0, 8.5], 'grey', linewidth=2, alpha=1)
        ppp.set_axis_labels(
            r'$\frac{N_e-N_t}{N_t}$' + ' (' + method + ')',
            r'$\frac{N_e-N_t}{N_t}$' + ' (Directional)',
            fontsize=14,
        )
        ppp.ax_joint.text(8, 8.5, "confidence interval", horizontalalignment='right', size='medium', color='crimson',)
        ppp.ax_joint.text(3, 5.5, "regression", horizontalalignment='left', size='medium', color='crimson',)
        ppp.ax_joint.text(5.5, 5, "baseline", horizontalalignment='left', size='medium', color='black',)
        sns.despine(right=True, top=True)
        plt.tight_layout()
        plt.show()

    def jointgrid(
            self,
            method='MCL-ed',
    ):
        # fig, ax = plt.subplots()
        self.df_melt = self.df_melt.rename(columns={'value': r'$\frac{N_e-N_t}{N_t}$'})

        self.df_melt['Sequencing error'] = self.df_melt['Sequencing error'].astype(float)
        self.df_melt = self.df_melt.loc[self.df_melt['Sequencing error'] <5.0e-03]
        print(self.df_melt)
        g = sns.JointGrid(
            data=self.df_melt[self.df_melt['method'] == method],
            x="Sequencing error",
            y=r'$\frac{N_e-N_t}{N_t}$',
            marginal_ticks=True,
        )

        # Create an inset legend for the histogram colorbar
        cax = g.figure.add_axes([.20, .55, .02, .2])
        g.ax_joint.set_xlabel('Sequencing error', fontsize=16)
        g.ax_joint.set_ylabel(r'$\frac{N_e-N_t}{N_t}$', fontsize=18)
        # Add the joint and marginal histogram plots
        g.plot_joint(
            sns.kdeplot,
            discrete=(True, False),
            cmap="light:#03012d", pmax=.8, cbar=True, cbar_ax=cax,
        )

        g.ax_joint.set_title(method, fontsize=14)
        xticks = self.df_melt['Sequencing error'].unique()
        print(xticks[xticks < 5.0e-03])
        g.ax_joint.set_xticks(np.linspace(xticks[xticks <5.0e-03][0], xticks[xticks <5.0e-03][-1], 4))
        # g.ax_joint.set_xticklabels(xticks[xticks <5.0e-03])
        # g.ax_joint.set_xticklabels([1e-05, '', '', '', '', '', '', '', 0.001, 0.0025, 0.005, 0.0075, 0.01, 0.025])
        g.ax_joint.set_xticklabels(np.linspace(xticks[xticks <5.0e-03][0], xticks[xticks <5.0e-03][-1], 4))
        plt.setp(g.ax_joint.get_xticklabels(), rotation=45)

        g.plot_marginals(
            sns.histplot,
            element="step",
            color="#03012d",
        )
        # sns.despine(left=True, bottom=True)
        plt.tight_layout()
        plt.show()

    def strip(self, ):
        fig, ax = plt.subplots()
        sns.despine(bottom=True, left=True)
        cc = [
            'tab:green',
            'tab:blue',
            'crimson',
        ]
        sns.stripplot(
            x="value",
            y="Sequencing error",
            hue='method',
            palette=cc,
            data=self.df_melt,
            dodge=True,
            alpha=.25,
            zorder=1,
        )
        sns.pointplot(
            x="value",
            y="Sequencing error",
            hue='method',
            data=self.df_melt,
            dodge=.8 - .8 / 3,
            join=False,
            palette=cc,
            markers="d",
            scale=.75,
            ci=None,
        )
        handles, labels = ax.get_legend_handles_labels()
        ax.legend(
            handles[len(handles)//2:],
            labels[len(handles)//2:],
            title='',
            handletextpad=0,
            columnspacing=1,
            loc="upper right",
            ncol=3,
            frameon=True,
            fontsize=12,
        )
        ax.set_xlabel(r'$\frac{N_e-N_t}{N_t}$', fontsize=16)
        ax.set_ylabel('Sequencing error rate', fontsize=16)
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

        self.df_direc = self.df[self.df['method'] == 'Directional']
        self.df_mcl_val = self.df[self.df['method'] == 'MCL-val']
        self.df_mcl_ed = self.df[self.df['method'] == 'MCL-ed']

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
        ax.legend(ncol=2, loc="upper right", frameon=True, fontsize=12)
        # ax.set(ylabel="", xlabel="Automobile collisions per billion miles")
        ax.set_ylabel('Sequencing error', fontsize=16)
        ax.set_xlabel(r'$\frac{N_e-N_t}{N_t}$', fontsize=16)
        sns.despine(left=True, bottom=True)
        fig.subplots_adjust(
            top=0.98,
            bottom=0.14,
            left=0.24,
            right=0.98,
            # hspace=0.40,
            # wspace=0.15
        )
        plt.show()

    def errorbar(self, ):
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

    def errorband(self, ):
        fig, ax = plt.subplots()
        cc = [
            'tab:green',
            'tab:blue',
            # 'dimgray',
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
            ax.plot(df_met.index, df_met['mean'] - df_met['mean-min'], color=cc[i], linewidth=0.1, alpha=0.1)
            ax.plot(df_met.index, df_met['max-mean'] + df_met['mean'], color=cc[i], linewidth=0.1, alpha=0.1)
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
        'mcl': 'MCL',
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
    )
    # print(p.strip())
    # print(p.jointplot(
    #     # method='MCL',
    #     method='MCL-val',
    #     # method='MCL-ed',
    # ))
    print(p.jointgrid(
        # method='Directional',
        # method='MCL',
        # method='MCL-val',
        method='MCL-ed',
    ))
    print(p.stackedbar())
    # print(p.errorbar())
    # print(p.errorband())