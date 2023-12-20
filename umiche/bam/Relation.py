__version__ = "v1.0"
__copyright__ = "Copyright 2024"
__license__ = "MIT"
__developer__ = "Jianfeng Sun"
__maintainer__ = "Jianfeng Sun"
__email__="jianfeng.sunmt@gmail.com"
__lab__ = "Cribbslab"

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
        self.df['umi_pcr#'] = self.df['query_name'].apply(lambda x: self.pcrnum(x))
        self.df['umi_src'] = self.df['query_name'].apply(lambda x: x.split('_')[0].split('-')[1])
        # print(self.df)

        umi_keymap_stime = time.time()
        self.df_umi_uniq = df.drop_duplicates(subset=['umi'], keep='first')
        self.uniq_umis = self.df_umi_uniq['umi'].values
        self.uniq_umi_num = self.uniq_umis.shape[0]
        print('===>unique UMI number: {}'.format(self.uniq_umi_num))

        self.umi_to_int_dict = {k: id for id, k in enumerate(self.uniq_umis)}
        self.int_to_umi_dict = {id: k for k, id in self.umi_to_int_dict.items()}
        self.df_umi_uniq_val_cnt = self.df['umi'].value_counts(ascending=False)

        # print(self.df_umi_uniq_val_cnt)
        df_umi_uniq_val_cnt_ids = self.df_umi_uniq_val_cnt.index
        self.df_umi_uniq_val_cnt.index = [self.umi_to_int_dict[i] for i in df_umi_uniq_val_cnt_ids]
        self.console.print('=========>umi keymap time: {:.3f}s'.format(time.time() - umi_keymap_stime))

        umi_trace_dict_stime = time.time()
        self.umi_id_to_origin_id_dict = {}
        # self.umi_id_to_origin_id_dict1 = {}
        # self.df_umi_gp = self.df.groupby(['umi'])

        sad = pd.Series(self.df_umi_uniq['umi#'].values, index=self.df_umi_uniq['umi'].values).to_dict()

        for uniq_umi in self.uniq_umis:
            self.umi_id_to_origin_id_dict[self.umi_to_int_dict[uniq_umi]] = [int(sad[uniq_umi])]

            # print(self.df_fastq_umi_gp.get_group(uniq_umi))
            # self.umi_id_to_origin_id_dict1[self.umi_to_int_dict[uniq_umi]] = self.df_umi_gp.get_group(uniq_umi)['umi#'].unique().astype(int).tolist()
        # print(self.umi_id_to_origin_id_dict)
        # print(self.umi_id_to_origin_id_dict1)
        # print(len(self.umi_id_to_origin_id_dict))
        # print(len(self.umi_id_to_origin_id_dict1))
        # print(self.umi_id_to_origin_id_dict1 == self.umi_id_to_origin_id_dict)
        print('===>umi trace dict time: {:.3f}s'.format(time.time() - umi_trace_dict_stime))

    def pcrnum(self, x):
        c = x.split('_')[0].split('-')
        if c[1] == 'init':
            return -1
        else:
            return c[2]