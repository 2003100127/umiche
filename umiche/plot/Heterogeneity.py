__version__ = "v1.0"
__copyright__ = "Copyright 2024"
__license__ = "MIT"
__developer__ = "Jianfeng Sun"
__maintainer__ = "Jianfeng Sun"
__email__="jianfeng.sunmt@gmail.com"
__lab__ = "Cribbslab"

import pandas as pd
import seaborn as sns
import matplotlib.pyplot as plt


class Heterogeneity:

    def __init__(self, ):
        sns.set(font="Helvetica")
        sns.set_style("ticks")

    def n1(
            self,
            df_disapv,
            df_apv
    ):
        fig, ax = plt.subplots(2, 1, figsize=(8, 6), sharex=True)
        colors = [
            'dimgray', # black
            'crimson',
        ]
        labels = [
            'different',
            'same',
        ]
        for i, col in enumerate(['diff_origin', 'same_origin']):
            ax[0].plot(
                df_disapv.index,
                df_disapv[col],
                label=labels[i],
                color=colors[i],
                lw=3,
            )
            ax[1].plot(
                df_apv.index,
                df_apv[col],
                label=labels[i],
                color=colors[i],
                lw=3,
            )
        # ax[0].set_xlabel('Time (ps)', fontsize=14)

        ax[0].set_ylabel('UMI count', fontsize=14)
        ax[0].set_title('Not merged', fontsize=12)
        ax[0].spines['right'].set_visible(False)
        ax[0].spines['top'].set_visible(False)

        # ax[1].set_xlabel('Time (ps)', fontsize=14)
        ax[1].set_xticks(df_apv.index)
        # ax[1].set_xticklabels(df_apv['metric'].apply(lambda x: 'PCR #' + x), fontsize=7, rotation=30)
        ax[1].set_ylabel('UMI count', fontsize=14)
        ax[1].set_title('Merged', fontsize=12)
        ax[1].spines['right'].set_visible(False)
        ax[1].spines['top'].set_visible(False)
        # sns.lineplot(data=data, palette="tab10", linewidth=2.5)
        handles1, labels1 = ax[0].get_legend_handles_labels()
        ax[0].legend(
            handles1,
            labels1,
            fontsize=10,
        )
        handles2, labels2 = ax[1].get_legend_handles_labels()
        ax[1].legend(
            handles2,
            labels2,
            fontsize=10,
        )

        fig.subplots_adjust(
            top=0.92,
            bottom=0.13,
            left=0.09,
            right=0.95,
            hspace=0.40,
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
    from umiche.util.Reader import reader as freader
    from umiche.path import to

    p = Heterogeneity()

    scenario='pcr_nums'
    # scenario='pcr_errs'
    # scenario='seq_errs'
    # scenario='ampl_rates'
    # scenario='umi_lens'
    # scenario = 'seq_deps'

    # method='unique'
    # method='cluster'
    # method='adjacency'
    # method='directional'
    # method='mcl'
    # method='mcl_val'
    method='mcl_ed'
    # method='mcl_cc_all_node_umis'
    # method='dbscan_seq_onehot'
    # method='birch_seq_onehot'
    # method='aprop_seq_onehot'
    # method='hdbscan_seq_onehot'
    # method='set_cover'

    df_apv_cnt = freader().generic(
        df_fpn=to('data/simu/mclumi/') + scenario + '/' + method + '_apv_cnt.txt',
        header=0,
    )
    print(df_apv_cnt)
    df_disapv_cnt = freader().generic(
        df_fpn=to('data/simu/mclumi/') + scenario + '/' + method + '_disapv_cnt.txt',
        header=0,
    )
    print(df_disapv_cnt)

    df_apv_cnt = df_apv_cnt.groupby(by=['scenario']).agg({'diff_origin': 'mean', 'same_origin': 'mean'}).reset_index()
    df_disapv_cnt = df_disapv_cnt.groupby(by=['scenario']).agg({'diff_origin': 'mean', 'same_origin': 'mean'}).reset_index()
    print(df_apv_cnt)

    p.n1(
        df_apv=df_apv_cnt,
        df_disapv=df_disapv_cnt,
    )