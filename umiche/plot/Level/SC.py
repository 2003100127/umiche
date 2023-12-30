__version__ = "v1.0"
__copyright__ = "Copyright 2024"
__license__ = "MIT"
__developer__ = "Jianfeng Sun"
__maintainer__ = "Jianfeng Sun"
__email__="jianfeng.sunmt@gmail.com"
__lab__ = "Cribbslab"

import numpy as np
import pandas as pd


class SingleCell:

    def __init__(self, ):
        pass

    def read(
            self,
            gmat_fpn,
    ):
        df_gmat = pd.read_hdf(gmat_fpn, 'df')
        if 'Y' in df_gmat.columns:
            df_gmat = df_gmat.drop(columns=['Y'])
        return df_gmat


if __name__ == "__main__":
    from umiche.path import to

    p = SingleCell()

    df = p.read(
        gmat_fpn=to('data/gmat_customized.h5'),
    )

    print(df)

    print(np.unique(df.to_numpy().flatten()).tolist())