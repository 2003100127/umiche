__version__ = "v1.0"
__copyright__ = "Copyright 2025"
__license__ = "GPL-3.0"
__developer__ = "Jianfeng Sun"
__maintainer__ = "Jianfeng Sun"
__email__="jianfeng.sunmt@gmail.com"


import numpy as np
import pandas as pd


def read(
        gmat_fpn,
):
    df_gmat = pd.read_hdf(gmat_fpn, 'df')
    y = df_gmat['Y'].values
    if 'Y' in df_gmat.columns:
        df_gmat = df_gmat.drop(columns=['Y'])
    return df_gmat, y


df, y = read(
    gmat_fpn=to('mclumi/gmat_customized.h5'),
)

print(df)

condi_map ={
    'umi_nums_seq_err_0.1': 'Sequencing error 0.1',
    'umi_nums_seq_err_0.05': 'Sequencing error 0.05',
    'umi_nums_seq_err_0.01': 'Sequencing error 0.01',
    'umi_nums_seq_err_0.005': 'Sequencing error 0.005',
}

condi = 'umi_nums_seq_err_0.1'
# condi = 'umi_nums_seq_err_0.05'
# condi = 'umi_nums_seq_err_0.01'
# condi = 'umi_nums_seq_err_0.005'

direc_dedup = pd.read_csv(to('mclumi/') + condi + '/directional_dedup.txt', header=0)
mcl_dedup = pd.read_csv(to('mclumi/') + condi + '/mcl_dedup.txt', header=0)
# print(direc_dedup)
# print(mcl_dedup)

# print(np.unique(df.to_numpy().flatten()).tolist())
uniq = [1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 14, 15, 16, 17, 18, 19, 20, 21, 22, 23, 24, 25, 26, 27, 28, 29, 30, 31, 32, 33, 34, 35, 36, 37, 38, 39, 40, 41, 42, 43, 44, 45]

dedup_direc_map = {x:y for x, y in zip(uniq, direc_dedup['pn0'].values.tolist())}
dedup_mcl_map = {x:y for x, y in zip(uniq, mcl_dedup['pn0'].values.tolist())}
print(dedup_direc_map)
print(dedup_mcl_map)

x = df.values
x_direc = df.replace(uniq, direc_dedup['pn0'].values.tolist()).values
x_mcl = df.replace(uniq, mcl_dedup['pn0'].values.tolist()).values

from ext.sequelpy.plot.scatter.DimensionReduction import dimensionReduction as drplot # TSNE
# drplot().single(X=x, y=y, tech='TSNE', marker_size=3, cmap='Paired', title='Real')
# drplot().single(X=x_direc, y=y, tech='TSNE', marker_size=3, cmap='Paired', title='Directional' + ' (' +  condi_map[condi] + ')')
# drplot().single(X=x_mcl, y=y, tech='TSNE', marker_size=3, cmap='Paired', title='MCL' + ' (' + condi_map[condi] + ')')
# drplot().single(X=xs, y=ys, tech='PCA', marker_size=3, cmap='tab20b', title='CondiCVAE')
# drplot().single(X=xs, y=ys, tech='UMAP', marker_size=3, cmap='tab20b', title='CondiCVAE')

drplot().single(X=x_direc, y=y, tech='UMAP', marker_size=3, cmap='Paired', title='Directional' + ' (' +  condi_map[condi] + ')')
drplot().single(X=x_mcl, y=y, tech='UMAP', marker_size=3, cmap='Paired', title='MCL' + ' (' + condi_map[condi] + ')')