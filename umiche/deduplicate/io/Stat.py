__version__ = "v1.0"
__copyright__ = "Copyright 2024"
__license__ = "MIT"
__developer__ = "Jianfeng Sun"
__maintainer__ = "Jianfeng Sun"
__email__="jianfeng.sunmt@gmail.com"
__lab__ = "Cribbslab"

from typing import Dict

import pandas as pd
from umiche.simu.Parameter import Parameter as params

from umiche.util.Reader import reader as freader


class Stat:

    def __init__(
            self,
            scenarios: Dict,
            methods: Dict,
            umi_gt_cnt: int = 50,
            param_fpn : str = None,
    ):
        self.scenarios = scenarios
        self.methods = methods
        self.umi_gt_cnt = umi_gt_cnt
        self.freader = freader()
        self.params = params(param_fpn=param_fpn)

    @property
    def df_dedup(self, ):
        df = pd.DataFrame()
        for scenario, scenario_formal in self.scenarios.items():
            for method, method_formal in self.methods.items():
                df_sce_met = self.freader.generic(
                    df_fpn=self.params.work_dir + scenario + '/' + method + '_dedup.txt',
                    header=0,
                )
                df_sce_met = (df_sce_met - 50) / 50
                df_sce_met['mean'] = df_sce_met.mean(axis=1)
                df_sce_met['max'] = df_sce_met.max(axis=1)
                df_sce_met['min'] = df_sce_met.min(axis=1)
                df_sce_met['std'] = df_sce_met.std(axis=1)
                # df_sce_met['mean-min'] = df_sce_met['std']
                # df_sce_met['max-mean'] = df_sce_met['std']
                df_sce_met['scenario'] = scenario_formal
                df_sce_met['method'] = method_formal
                df_sce_met['metric'] = [str(x) for x in self.params.varied[scenario]]
                df = pd.concat([df, df_sce_met], axis=0)
        return df

    @property
    def df_trace_cnt(self, ):
        df_apv = pd.DataFrame()
        df_disapv = pd.DataFrame()
        for scenario, scenario_formal in self.scenarios.items():
            for method, method_formal in self.methods.items():
                df_sce_met_apv_perm = self.freader.generic(
                    df_fpn=self.params.work_dir + scenario + '/' + method + '_apv_cnt.txt',
                    header=0,
                )
                df_sce_met_disapv_perm = self.freader.generic(
                    df_fpn=self.params.work_dir + scenario + '/' + method + '_disapv_cnt.txt',
                    header=0,
                )
                df_sce_met_apv_perm = df_sce_met_apv_perm.rename(columns={'scenario': 'metric'})
                df_sce_met_disapv_perm = df_sce_met_disapv_perm.rename(columns={'scenario': 'metric'})

                df_sce_met_apv = df_sce_met_apv_perm.groupby(by=['metric']).agg(
                    {'diff_origin': 'mean', 'same_origin': 'mean'}).reset_index()
                df_sce_met_apv['diff_origin_max'] = df_sce_met_apv_perm.groupby(by=['metric']).agg(
                    {'diff_origin': 'max'})['diff_origin'].values
                df_sce_met_apv['diff_origin_min'] = df_sce_met_apv_perm.groupby(by=['metric']).agg(
                    {'diff_origin': 'min'})['diff_origin'].values
                df_sce_met_apv['same_origin_max'] = df_sce_met_apv_perm.groupby(by=['metric']).agg(
                    {'same_origin': 'max'})['same_origin'].values
                df_sce_met_apv['same_origin_min'] = df_sce_met_apv_perm.groupby(by=['metric']).agg(
                    {'same_origin': 'min'})['same_origin'].values
                df_sce_met_apv['scenario'] = scenario_formal
                df_sce_met_apv['method'] = method_formal

                df_sce_met_disapv = df_sce_met_disapv_perm.groupby(by=['metric']).agg(
                    {'diff_origin': 'mean', 'same_origin': 'mean'}).reset_index()
                df_sce_met_disapv['diff_origin_max'] = df_sce_met_disapv_perm.groupby(by=['metric']).agg(
                    {'diff_origin': 'max'})['diff_origin'].values
                df_sce_met_disapv['diff_origin_min'] = df_sce_met_disapv_perm.groupby(by=['metric']).agg(
                    {'diff_origin': 'min'})['diff_origin'].values
                df_sce_met_disapv['same_origin_max'] = df_sce_met_disapv_perm.groupby(by=['metric']).agg(
                    {'same_origin': 'max'})['same_origin'].values
                df_sce_met_disapv['same_origin_min'] = df_sce_met_disapv_perm.groupby(by=['metric']).agg(
                    {'same_origin': 'min'})['same_origin'].values
                df_sce_met_disapv['scenario'] = scenario_formal
                df_sce_met_disapv['method'] = method_formal

                # print(df_sce_met_apv)
                # print(df_sce_met_disapv)

                df_apv = pd.concat([df_apv, df_sce_met_apv], axis=0)
                df_disapv = pd.concat([df_disapv, df_sce_met_disapv], axis=0)
        return {
            "apv": df_apv,
            "disapv": df_disapv,
        }

    @property
    def df_dedup_melt(self, ):
        df = pd.melt(
            frame=self.df_dedup[['method', 'metric', 'scenario', 'mean']],
            id_vars=['method', 'metric', 'scenario'],
            # value_vars=['mean',],
        )
        return df


if __name__ == "__main__":
    from umiche.path import to

    p = Stat(
        scenarios={
            'pcr_nums': 'PCR cycle',
            'pcr_errs': 'PCR error',
            # 'seq_errs': 'Sequencing error',
            # 'ampl_rates': 'Amplification rate',
            # 'umi_lens': 'UMI length',
            # 'seq_deps': 'Sequencing depth',
        },

        methods={
            # 'unique': 'Unique',
            # 'cluster': 'Cluster',
            # 'adjacency': 'Adjacency',
            'directional': 'Directional',
            # 'dbscan_seq_onehot': 'DBSCAN',
            # 'birch_seq_onehot': 'Birch',
            # 'aprop_seq_onehot': 'Affinity Propagation',
            # 'mcl': 'MCL',
            'mcl_val': 'MCL-val',
            'mcl_ed': 'MCL-ed',
        },

        param_fpn=to('data/params.yml'),
    )

    # print(p.df_dedup)
    # print(p.df_dedup_melt)
    print(p.df_trace_cnt)