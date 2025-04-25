__version__ = "v1.0"
__copyright__ = "Copyright 2025"
__license__ = "GPL-3.0"
__developer__ = "Jianfeng Sun"
__maintainer__ = "Jianfeng Sun"
__email__="jianfeng.sunmt@gmail.com"


import pandas as pd
import seaborn as sns
from matplotlib import pyplot as plt
import umiche.plot.bar.external.lca_standard_graphs as lsg
from umiche.util.Reader import Reader as freader
from umiche.path import to


class editdistance:

    def __init__(self, ):
        self.freader = freader()

        sns.set(font="Helvetica")
        sns.set_style("ticks")

    def umi_tools(self, ):
        tool = 'UMI-tools'
        df_raw = self.freader.generic(
            df_fpn=to('ext/data/seq/umi/example.bam/ed1/umi-tools/') + 'deduplicated_edit_distance.tsv',
            df_sep='\t',
            header=0,
            skiprows=1,
        )
        df_raw.columns = ['unique', 'unique_null', 'direc', 'direc_null', 'ed']
        unique = df_raw[['unique', 'ed']]
        unique = unique.rename(columns={'unique': 'cnt'}).T
        unique['met'] = 'unique'
        unique['tool'] = tool
        print(unique)

        direc = df_raw[['direc', 'ed']]
        direc = direc.rename(columns={'direc': 'cnt'}).T
        direc['met'] = 'directional'
        direc['tool'] = tool
        print(direc)
        df = pd.concat([unique, direc])
        return df

    def data(self, ):
        df_umi_tools = self.umi_tools()
        tool = 'Mclumi'
        umi_len = 12
        met_full = ['unique', 'cluster', 'adjacency', 'directional', 'mcl', 'mcl_val', 'mcl_ed']
        met_abbr = ['uniq', 'cc', 'adj', 'direc', 'mcl', 'mcl_val', 'mcl_ed']
        df = pd.DataFrame()
        for met_a, met_f in zip(met_abbr, met_full):
            asd = self.freader.generic(
                # df_fpn=to('ext/data/seq/umi/example.bam/ed1/mclumi/') + met_a + '_ave_ed_pos_bin.txt',
                df_fpn=to('ext/data/bulk/RM82CLK1_S3/umi/ed7/') + met_a + '_ave_ed_pos_bin.txt',
                # df_fpn=to('ext/data/sc/umi/hgmm100/ed6/') + met_a + '_ave_ed_pos_bin.txt',
                df_sep='\t',
                header=None,
                skiprows=1,
            )
            asd.columns = ['ed', 'cnt']
            t = dict(zip(asd.ed, asd.cnt))
            for i in range(umi_len + 1):
                if i not in [*t.keys()]:
                    t[i] = 0
            c = pd.DataFrame({'ed': t.keys(), 'cnt': t.values()}).sort_values(by='ed').reset_index(drop=True)
            c = c.T
            c['met'] = met_f
            c['tool'] = tool
            # print(c)
            df = pd.concat([df, c])
            print(df)
        # df = pd.concat([df, df_umi_tools])
        print(df)
        df = df.loc[df.index == 'cnt']
        # df = df.loc[(df.met == 'unique') | (df.met == 'directional')]
        df = df.rename(columns={i: 'ave edit distance ' + str(i) for i in range(umi_len + 1)})
        df1 = df.set_index(['met', 'tool'])
        print(df1)

        ax, fig = lsg.plot_grouped_stackedbars(
            df1, ix_categories='met',
            ix_entities_compared='tool',

            # df1, ix_categories='tool',
            # ix_entities_compared='met',

            xaxis_label='Relative (%) number of genomic positions observed',
            figsize=(10, 4)
        )
        plt.show()
        return

    def wk2(self, ):
        tool = 'Mclumi'
        umi_len = 10
        met_full = ['Cluster', 'Adjacency', 'Directional', 'MCL']
        met_full_s = ['cluster', 'adjacency', 'directional', 'mcl']
        met_abbr = ['cc', 'adj', 'direc', 'mcl']
        # df = pd.DataFrame()
        df_ex = self.freader.generic(
            df_fpn=to('data/simu/sc/hgmm100/ed6/uniq_dedup_sum.txt'),
            df_sep='\t',
            header=0,
            index_col=0,
            # skiprows=1,
        )
        df_ex['cell'] = df_ex.apply(lambda x: x.name.split(', ')[0].split('(')[1].split("'")[1], axis=1)
        df_ex['gene'] = df_ex.apply(lambda x: x.name.split(', ')[1].split(')')[0].split("'")[1], axis=1)
        cell_types = df_ex['cell'].unique()
        cell_maps = {v: i for i, v in enumerate(cell_types)}
        gene_types = df_ex['gene'].unique()
        gene_maps = {v: i for i, v in enumerate(gene_types)}
        gene_rev_maps = {i: v for i, v in enumerate(gene_types)}
        for met_a, met_f in zip(met_abbr, met_full):
            df = self.freader.generic(
                # df_fpn=to('ext/data/seq/umi/example.bam/ed1/mclumi/') + met_a + '_ave_ed_pos_bin.txt',
                # df_fpn=to('ext/data/bulk/RM82CLK1_S3/umi/') + met_a + '_ave_ed_pos_bin.txt',
                df_fpn=to('data/simu/sc/hgmm100/ed6/') + met_a + '_dedup_sum.txt',
                df_sep='\t',
                header=0,
                index_col=0,
                # skiprows=1,
            )
            df['cell'] = df.apply(lambda x: x.name.split(', ')[0].split('(')[1].split("'")[1], axis=1)
            df['gene'] = df.apply(lambda x: x.name.split(', ')[1].split(')')[0].split("'")[1], axis=1)
            df['cell_type'] = df['cell'].apply(lambda x: cell_maps[x])
            df['gene_type'] = df['gene'].apply(lambda x: gene_maps[x])
            # print(df.loc[df['cell_type'].isin(range(20))])
            # print(df.loc[df['cell_type'].isin(range(20))])
            from scipy.sparse import csr_matrix
            gc_mat = csr_matrix(
                (df[met_a + '_umi_len'].values, (df['cell_type'].values, df['gene_type'].values)),
                shape=(len(cell_maps), len(gene_maps))
            )
            print('===>gc matrix shape: {}'.format(gc_mat.shape))
            from umiche.plot.stat.gene.expression.Extreme import Extreme as ppp
            high_genes = ppp().geneIds(
                gc_mat=gc_mat.toarray(),
                oriented='maximal',
                n_top=10,
                scheme={'method': 'mean', 'cell_ratio': 0.9},
            )
            gene_ids = [gene_rev_maps[i] for i in high_genes]
            print(gene_ids)
            newdfd = df.loc[df['gene'].isin(gene_ids)]
            self.bplot(x_col_name=met_a + '_umi_len', y_col_name='gene', df=newdfd, xlabel=met_f)

            # for i in gene_ids:
            #     # print(df.loc[df['gene'].isin([gene_rev_maps[i]])][met_a + '_umi_len'].values)
            #     print(df.loc[df['gene'].isin([i])])
        return

    def asdas(self, ):
        import numpy as np
        np.random.seed(19680801)
        data = np.random.lognormal(size=(37, 4), mean=1.5, sigma=1.75)
        labels = list('ABCD')
        fs = 10  # fontsize
        fig, ax = plt.subplots(nrows=1, ncols=1, figsize=(6, 6), sharey=True)
        bplot = ax.boxplot(data, labels=labels)
        ax.set_title('Default', fontsize=fs)
        ax.spines['right'].set_color('none')
        ax.spines['top'].set_color('none')

        # for bplot_handle, color in zip(boxplot_handles, palette):
        for patch in bplot['boxes']:
            patch.set(color='black', linewidth=1.1)
            # patch.set(facecolor='white', alpha=0.3)
        for whisker in bplot['whiskers']:
            whisker.set(color='black', linewidth=1)
        for cap in bplot['caps']:
            cap.set(color='black', linewidth=1)
        for median in bplot['medians']:
            median.set(color='black', linewidth=1)
        for flier in bplot['fliers']:
            flier.set(marker='o', color='y', alpha=0.5)

        fig.subplots_adjust(hspace=0.4)
        plt.show()

    def bplot(self, x_col_name, y_col_name, xlabel, df):
        fig, ax = plt.subplots(figsize=(8, 6))
        # ax.set_xscale("log")

        # Plot the orbital period with horizontal boxes
        sns.boxplot(
            x=x_col_name,
            y=y_col_name,
            data=df,
            whis=[0, 100],
            width=0.6,
            fill=False,
            gap=0.1,
            color="grey",
            # palette="grey",
        )
        # Add in points to show each observation
        sns.stripplot(
            x=x_col_name,
            y=y_col_name,
            data=df,
            size=2,
            color="whitesmoke",
            linewidth=0.8,
            alpha=0.8,
        )
        # Tweak the visual presentation
        # ax.xaxis.grid(True)
        ax.set(ylabel="")
        ax.set(xlabel="")
        sns.despine(trim=True, left=True)

        ax.set_title(xlabel, fontsize=16)
        fig.subplots_adjust(
            # bottom=0.26,
            # top=0.92,
            left=0.25,
            # right=0.95, hspace=0.25, wspace=0.17
        )
        plt.show()


if __name__ == "__main__":

    p = editdistance()
    print(p.wk2())
    # print(p.asdas())

    # print(p.data())
# pd.options.display.multi_sparse = False

# # Define dataframes
# df1=pd.DataFrame(np.random.rand(4,2),index=["A","B","C","D"], columns=["zI","yJ"])
# df1.index.name='Criteria'

# df2=pd.DataFrame(np.random.rand(4,2),index=["A","B","C","D"], columns=["zI","yJ"])
# df2.index.name='Criteria'

# # Build comparison table
# comp = lsg.build_comparison_table([df1, df2], ['df1', 'df2'], fillna=0.0)

# ax, fig = lsg.plot_grouped_stackedbars(comp, ix_categories='Criteria', ix_entities_compared='Scenarios', norm=None)

# plt.show()