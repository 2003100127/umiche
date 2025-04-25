__version__ = "v1.0"
__copyright__ = "Copyright 2025"
__license__ = "GPL-3.0"
__developer__ = "Jianfeng Sun"
__maintainer__ = "Jianfeng Sun"
__email__="jianfeng.sunmt@gmail.com"


import numpy as np
import pandas as pd
from path import to
from ext.sequelpy.plot.scatter.DimensionReduction import dimensionReduction as drplot


def read(
        gmat_fpn,
):
    df_gmat = pd.read_hdf(gmat_fpn, 'df')
    y = df_gmat['Y'].values
    if 'Y' in df_gmat.columns:
        df_gmat = df_gmat.drop(columns=['Y'])
    return df_gmat, y


# title = '10 epoch'
# title = '20 epoch'
# title = '50 epoch'
title = '100 epoch'

df, y = read(
    gmat_fpn=to('umiche/') + title + '.h5',
)
x = df.values
# print(x)
# print(y)

drplot().single(X=x, y=y, tech='PCA', marker_size=3, cmap='Set3', title=title)
drplot().single(X=x, y=y, tech='UMAP', marker_size=3, cmap='Set3', title=title)
drplot().single(X=x, y=y, tech='TSNE', marker_size=3, cmap='Set3', title=title) # Set3 tab20b
