__version__ = "v1.0"
__copyright__ = "Copyright 2025"
__license__ = "GPL-3.0"
__developer__ = "Jianfeng Sun"
__maintainer__ = "Jianfeng Sun"
__email__="jianfeng.sunmt@gmail.com"


import numpy as np
import pandas as pd
from DimensionReduction import dimensionReduction as drplot


def read(
        gmat_fpn,
):
    df_gmat = pd.read_hdf(gmat_fpn, 'df')
    y = df_gmat['Y'].values
    if 'Y' in df_gmat.columns:
        df_gmat = df_gmat.drop(columns=['Y'])
    return df_gmat, y


title = '50 epoch 2000 cell'
# title = '50 epoch 1000 cell-new'
# title = '100 epoch 5000 cell'

df, y = read(
    gmat_fpn=to('tresor/') + title + '.h5',
)
x = df.values
print(123123123213, x.shape)
# print(np.unique(df.to_numpy().flatten()).tolist())
# print(x)
# print(y)

# condi = 'umi_nums_seq_err_0.1'
# condi = 'umi_nums_seq_err_0.05'
# condi = 'umi_nums_seq_err_0.01'
condi = 'umi_nums_seq_err_0.005'

direc_dedup = pd.read_csv(to('mclumi/') + condi + '/directional_dedup.txt', header=0)
# print(direc_dedup)

# print(np.unique(df.to_numpy().flatten()).tolist())
uniq = [1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 14, 15, 16, 17, 18, 19, 20, 21, 22, 23, 24, 25, 26, 27, 28, 29, 30, 31, 32, 33, 34, 35, 36, 37, 38, 39, 40, 41, 42, 43, 44, 45]

dedup_direc_map = {x:y for x, y in zip(uniq, direc_dedup['pn0'].values.tolist())}
# print(dedup_direc_map)
x_direc = df.replace(uniq, direc_dedup['pn0'].values.tolist()).values


# drplot().single(X=x, y=y, tech='PCA', marker_size=1, cmap='tab20b', title=title)
# drplot().single(X=x_direc, y=y, tech='PCA', marker_size=1, cmap='tab20b', title=title)

# drplot().single(X=x, y=y, tech='UMAP', marker_size=1, cmap='tab20b', title=title)
# drplot().single(X=x_direc, y=y, tech='UMAP', marker_size=1, cmap='tab20b', title=title)

drplot().single(X=x, y=y, tech='TSNE', marker_size=1, cmap='tab20b', title=title) # Set3 tab20b
drplot().single(X=x_direc, y=y, tech='TSNE', marker_size=1, cmap='tab20b', title=title) # Set3 tab20b

