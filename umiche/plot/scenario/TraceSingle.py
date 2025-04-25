__version__ = "v1.0"
__copyright__ = "Copyright 2025"
__license__ = "GPL-3.0"
__developer__ = "Jianfeng Sun"
__maintainer__ = "Jianfeng Sun"
__email__="jianfeng.sunmt@gmail.com"


import seaborn as sns
import matplotlib.pyplot as plt


class TraceSingle:

    def __init__(
            self,
            df_apv,
            df_disapv,
    ):
        self.df_apv = df_apv
        self.df_disapv = df_disapv
        self.method = self.df_disapv['method'].unique()[0]
        sns.set(font="Helvetica")
        sns.set_style("ticks")

    def line_apv_disapv(self, ):
        fig, ax = plt.subplots(2, 1, figsize=(5, 6), sharex=True)
        colors = [
            'dimgray',  # black dimgray
            'palevioletred', # crimson mediumvioletred
        ]
        # colors = sns.color_palette("Set2")
        labels = [
            'different origin',
            'same origin',
        ]
        for i, col in enumerate(['diff_origin', 'same_origin']):
            ax[0].plot(
                self.df_disapv.index,
                self.df_disapv[col],
                label=labels[i],
                color=colors[i],
                lw=2.5,
                # alpha=0.8,
            )
            ax[1].plot(
                self.df_apv.index,
                self.df_apv[col],
                label=labels[i],
                color=colors[i],
                lw=2.5,
                # alpha=0.8,
            )

        ax[0].set_ylabel('UMI count', fontsize=14)
        ax[0].set_title('Not merged' + ' (' + self.method + ')', fontsize=12)
        ax[0].spines['right'].set_visible(False)
        ax[0].spines['top'].set_visible(False)

        ax[1].set_xlabel('PCR cycle', fontsize=12)
        ax[1].set_xticks(self.df_apv.index)
        # ax[1].set_xticklabels(df_apv['metric'].apply(lambda x: 'PCR #' + x), fontsize=7, rotation=30)
        ax[1].set_ylabel('UMI count', fontsize=14)
        ax[1].set_title('Merged' + ' (' + self.method + ')', fontsize=12)
        ax[1].spines['right'].set_visible(False)
        ax[1].spines['top'].set_visible(False)
        handles1, labels1 = ax[0].get_legend_handles_labels()
        ax[0].legend(
            handles1,
            labels1,
            fontsize=12,
        )
        handles2, labels2 = ax[1].get_legend_handles_labels()
        ax[1].legend(
            handles2,
            labels2,
            fontsize=12,
        )
        fig.subplots_adjust(
            top=0.95,
            bottom=0.1,
            left=0.18,
            right=0.95,
            hspace=0.20,
            # wspace=0.15
        )
        plt.show()

    def n2(self, df):
        fig, ax = plt.subplots(1, 1, figsize=(10, 6), sharex=True)
        df_gp = df.groupby(by=['method'])
        df_gp_keys = df_gp.groups.keys()
        method_map = {
            'ccs': 'cluster',
            'adj': 'adjacency',
            'direc': r'$directional$',
            'mcl_val': 'MCL-val',
            'mcl_ed': 'MCL-ed',
        }
        palette = {
            'ccs': 'red',
            'adj': 'black',
            'direc': 'steelblue',
            'mcl_val': 'chocolate',
            'mcl_ed': 'firebrick',
            # 'black',
            # 'chocolate',
            # 'saddlebrown',
            # 'darkgoldenrod',
            # 'firebrick',
        }
        for method in df_gp_keys:
            df_met = df_gp.get_group(method)
            print(df_met)
            # if method != 'adj':
            ax.plot(
                df_met['metric'],
                df_met['dedup_cnt'].apply(
                    # lambda x: x / 50
                    lambda x: (x - 50) / 50
                    # lambda x: np.exp((x - 50) / 50)
                ),
                label=method_map[method],
                color=palette[method],
                lw=2.5,
                alpha=0.7
            )
        # ax[0].set_xlabel('Time (ps)', fontsize=14)
        c = df.loc[df['method'] == 'mcl_ed']['metric']
        ax.set_xticks(c)
        ax.set_xticklabels(c, fontsize=8)
        # ax.set_xticklabels(c.astype(np.float).apply(lambda x: '{:.2e}'.format(x)), fontsize=8, rotation=30)
        # ax.set_xticklabels(c.astype(np.float).round(1), fontsize=8)

        # ax.set_xlabel('PCR cycle', fontsize=11)
        # ax.set_xlabel('Polymerase error', fontsize=11)
        # ax.set_xlabel('Sequencing error', fontsize=11)
        ax.set_xlabel('UMI length', fontsize=11)
        # ax.set_xlabel('Amplification rate', fontsize=11)
        ax.set_ylabel(r'$\frac{N_e-N_t}{N_t}$', fontsize=16)
        ax.spines['right'].set_visible(False)
        ax.spines['top'].set_visible(False)
        # sns.lineplot(data=data, palette="tab10", linewidth=2.5)
        handles1, labels1 = ax.get_legend_handles_labels()
        ax.legend(
            handles1,
            labels1,
            fontsize=10,
        )
        fig.subplots_adjust(
            # top=0.92,
            # bottom=0.13,
            # left=0.13,
            # right=0.95,
            hspace=0.40,
            # wspace=0.15
        )
        plt.show()

    def n2dist(self, df):
        sns.displot(data=df, x='dedup_cnt', hue='method', kind="kde", rug=True)
        plt.show()


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
        # 'unique': 'Unique',
        # 'cluster': 'Cluster',
        # 'adjacency': 'Adjacency',
        'directional': 'Directional',
        # 'dbscan_seq_onehot': 'DBSCAN',
        # 'birch_seq_onehot': 'Birch',
        # 'aprop_seq_onehot': 'Affinity Propagation',
        # 'mcl': 'MCL',
        # 'mcl_val': 'MCL-val',
        # 'mcl_ed': 'MCL-ed',
    }

    dedupstat = dedupstat(
        scenarios=scenarios,
        methods=methods,
        param_fpn=to('data/params.yml'),
    )

    p = TraceSingle(
        df_apv=dedupstat.df_trace_cnt['apv'],
        df_disapv=dedupstat.df_trace_cnt['disapv'],
    )

    p.line_apv_disapv()