__version__ = "v1.0"
__copyright__ = "Copyright 2025"
__license__ = "GPL-3.0"
__developer__ = "Jianfeng Sun"
__maintainer__ = "Jianfeng Sun"
__email__="jianfeng.sunmt@gmail.com"


import pandas as pd
from functools import wraps


class Reader:



    def __call__(self, deal):
        generic = self.generic
        excel = self.excel
        @wraps(deal)
        def read(ph, *args, **kwargs):
            deal(ph, **kwargs)
            keys = [*kwargs.keys()]
            if kwargs['type'] == 'generic':
                return generic(
                    df_fpn=kwargs['df_fpn'],
                    df_sep='\t' if 'df_sep' not in keys else kwargs['df_sep'],
                    skiprows=False if 'skiprows' not in keys else kwargs['skiprows'],
                    header=None if 'header' not in keys else kwargs['header'],
                    is_utf8=False if 'is_utf8' not in keys else kwargs['is_utf8'],
                )
            elif kwargs['type'] == 'excel':
                return excel(
                    df_fpn=kwargs['df_fpn'],
                    sheet_name='Sheet1' if 'sheet_name' not in keys else kwargs['sheet_name'],
                    header=None if 'header' not in keys else kwargs['header'],
                    is_utf8=False if 'is_utf8' not in keys else kwargs['is_utf8'],
                )
        return read

    def generic(self, df_fpn, df_sep='\t', skiprows=None, header=None, index_col=None, is_utf8=False):
        """

        Parameters
        ----------
        df_fpn
        df_sep
        skiprows
        header
        is_utf8

        Returns
        -------

        """
        if is_utf8:
            return pd.read_csv(
                df_fpn,
                sep=df_sep,
                header=header,
                index_col=index_col,
                encoding='utf-8',
                skiprows=skiprows,
            )
        else:
            return pd.read_csv(
                df_fpn,
                sep=df_sep,
                header=header,
                index_col=index_col,
                skiprows=skiprows
            )

    def excel(self, df_fpn, sheet_name='Sheet1', header=None, is_utf8=False):
        """

        Parameters
        ----------
        df_fpn
        sheet_name
        header
        is_utf8

        Returns
        -------

        """
        if is_utf8:
            return pd.read_excel(
                df_fpn,
                sheet_name=sheet_name,
                header=header,
                encoding='utf-8',
                engine='openpyxl',
            )
        else:
            return pd.read_excel(
                df_fpn,
                sheet_name=sheet_name,
                header=header,
                engine='openpyxl',
            )