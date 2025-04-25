__version__ = "v1.0"
__copyright__ = "Copyright 2025"
__license__ = "GPL-3.0"
__developer__ = "Jianfeng Sun"
__maintainer__ = "Jianfeng Sun"
__email__="jianfeng.sunmt@gmail.com"


import time
import pandas as pd
from umiche.util.Console import Console


class Relation:

    def __init__(
            self,
            df,
            verbose=False,
    ):
        self.console = Console()
        self.console.verbose = verbose

        self.df = df
        self.df['umi#'] = self.df['query_name'].apply(lambda x: x.split('_')[0].split('-')[0])
        # print(self.df)
        ### @@ self.df
        #           id            query_name  flag  ...         umi  source  umi#
        # 0          0   23-pcr-8_GCATCAGTAG     0  ...  GCATCAGTAG       1    23
        # 1          1    6-pcr-8_CTGACTGCGC     0  ...  CTGACTGCGC       1     6
        # 2          2   13-pcr-7_CACACGATAC     0  ...  CACACGATAC       1    13
        # ...      ...                   ...   ...  ...         ...     ...   ...
        # 23813  23813   47-pcr-9_GCGAAGAGGA     0  ...  GCGAAGAGGA       1    47
        self.df['umi_pcr#'] = self.df['query_name'].apply(lambda x: self.pcrnum(x))
        # print(self.df)
        ### @@ self.df
        #           id            query_name  flag  ...  source  umi#  umi_pcr#
        # 0          0   23-pcr-8_GCATCAGTAG     0  ...       1    23         8
        # 1          1    6-pcr-8_CTGACTGCGC     0  ...       1     6         8
        # 2          2   13-pcr-7_CACACGATAC     0  ...       1    13         7
        # ...      ...                   ...   ...  ...     ...   ...       ...
        # 23813  23813   47-pcr-9_GCGAAGAGGA     0  ...       1    47         9
        ### @@@ note version 2 for resimpy
        # self.df['umi_src'] = self.df['query_name'].apply(lambda x: x.split('_')[0].split('-')[1])
        ### @@@ note version 2 for phylotres as well as resimpy
        self.df['umi_src'] = self.df['query_name'].apply(lambda x: x.split('_')[-2].split('-')[1])
        # print(self.df)
        ### @@ self.df
        #           id            query_name  flag  ...  umi#  umi_pcr#  umi_src
        # 0          0   23-pcr-8_GCATCAGTAG     0  ...    23         8      pcr
        # 1          1    6-pcr-8_CTGACTGCGC     0  ...     6         8      pcr
        # 2          2   13-pcr-7_CACACGATAC     0  ...    13         7      pcr
        # ...      ...                   ...   ...  ...   ...       ...      ...
        # 23813  23813   47-pcr-9_GCGAAGAGGA     0  ...    47         9      pcr
        umi_keymap_stime = time.time()
        self.df_umi_uniq = df.drop_duplicates(subset=['umi'], keep='first')
        ### @@ self.df_umi_uniq
        #           id            query_name  flag  ...  umi#  umi_pcr#  umi_src
        # 0          0   23-pcr-8_GCATCAGTAG     0  ...    23         8      pcr
        # 1          1    6-pcr-8_CTGACTGCGC     0  ...     6         8      pcr
        # 2          2   13-pcr-7_CACACGATAC     0  ...    13         7      pcr
        # ...      ...                   ...   ...  ...   ...       ...      ...
        # 23691  23691   23-pcr-8_GCATCAGGAG     0  ...    23         8      pcr
        # [257 rows x 17 columns]
        self.uniq_umis = self.df_umi_uniq['umi'].values
        self.uniq_umi_num = self.uniq_umis.shape[0]
        self.console.print('==================>unique UMI number: {}'.format(self.uniq_umi_num))

        self.umi_to_int_dict = {k: id for id, k in enumerate(self.uniq_umis)}
        ### @@ self.umi_to_int_dict
        # {'GCATCAGTAG': 0, 'CTGACTGCGC': 1, 'CACACGATAC': 2, 'TTAGATGATT': 3, ..., 'GAAGTATATT': 255, 'GCATCAGGAG': 256}
        self.int_to_umi_dict = {id: k for k, id in self.umi_to_int_dict.items()}
        self.df_umi_uniq_val_cnt = self.df['umi'].value_counts(ascending=False)
        ### @@ self.df_umi_uniq_val_cnt
        # umi
        # AGTACGCGAG    636
        # GGAGATCCGG    625
        # AGTAGAACCC    619
        #              ...
        # GCATCAGGAG      1
        # Name: count, Length: 257, dtype: int64
        df_umi_uniq_val_cnt_ids = self.df_umi_uniq_val_cnt.index
        self.df_umi_uniq_val_cnt.index = [self.umi_to_int_dict[i] for i in df_umi_uniq_val_cnt_ids]
        self.console.print('==================>umi keymap time: {:.3f}s'.format(time.time() - umi_keymap_stime))

        umi_trace_dict_stime = time.time()
        self.umi_id_to_origin_id_dict = {}
        umi_str_to_origin_dict = pd.Series(self.df_umi_uniq['umi#'].values, index=self.df_umi_uniq['umi'].values).to_dict()
        ### @@ umi_str_to_origin_dict
        # {'GCATCAGTAG': '23', 'CTGACTGCGC': '6', 'CACACGATAC': '13', 'TTAGATGATT': '10', ..., 'GAAGTATATT': '27', 'GCATCAGGAG': '23'}
        for uniq_umi in self.uniq_umis:
            self.umi_id_to_origin_id_dict[self.umi_to_int_dict[uniq_umi]] = int(umi_str_to_origin_dict[uniq_umi])
        ### @@ self.umi_id_to_origin_id_dict
        # {0: 23, 1: 6, 2: 13, 3: 10, ..., 255: 27, 256: 23}

        self.console.print('==================>umi trace dict time: {:.3f}s'.format(time.time() - umi_trace_dict_stime))

    def pcrnum(self, x):
        # print(x.split('_'))
        ### @@@ note version 1 for resimpy
        # c = x.split('_')[0].split('-')
        ### @@@ note version 2 for phylotres as well as resimpy
        c = x.split('_')[-2].split('-')
        # print(c)
        if c[1] == 'init':
            return -1
        else:
            return c[2]