__version__ = "v1.0"
__copyright__ = "Copyright 2025"
__license__ = "GPL-3.0"
__developer__ = "Jianfeng Sun"
__maintainer__ = "Jianfeng Sun"
__email__="jianfeng.sunmt@gmail.com"


from typing import Dict

import pandas as pd
from umiche.simu.Parameter import Parameter as params

from umiche.util.Reader import Reader as freader
from umiche.util.Console import Console


class Stat:

    def __init__(
            self,
            scenarios: Dict,
            methods: Dict,
            umi_gt_cnt: int = 50,
            param_fpn : str = None,
            is_trans : bool = True,
            verbose: bool = True,
    ):
        self.scenarios = scenarios
        self.methods = methods
        self.umi_gt_cnt = umi_gt_cnt
        self.freader = freader()
        self.params = params(param_fpn=param_fpn)
        self.is_trans = is_trans

        self.console = Console()
        self.console.verbose = verbose

    @property
    def df_dedup(self, ):
        df = pd.DataFrame()
        for scenario, scenario_formal in self.scenarios.items():
            self.console.print("======>scenario: {}".format(scenario_formal))
            for method, method_formal in self.methods.items():
                self.console.print("=========>method: {}".format(method_formal))
                ph = '/' if method_formal != 'Birch' else '/birch/'
                df_sce_met = self.freader.generic(
                    # df_fpn=self.params.work_dir + scenario + '/' + method + '_dedup.txt',
                    df_fpn=self.params.work_dir + scenario + ph + method + '_dedup.txt',
                    header=0,
                )
                if self.is_trans:
                    df_sce_met = (df_sce_met - self.umi_gt_cnt) / self.umi_gt_cnt
                df_sce_met['mean'] = df_sce_met.mean(axis=1)
                df_sce_met['max'] = df_sce_met.max(axis=1)
                df_sce_met['min'] = df_sce_met.min(axis=1)
                df_sce_met['std'] = df_sce_met.std(axis=1)
                df_sce_met['mean-min'] = df_sce_met['std']
                df_sce_met['max-mean'] = df_sce_met['std']
                df_sce_met['scenario'] = scenario_formal
                df_sce_met['method'] = method_formal
                df_sce_met['metric'] = [str(x) for x in self.params.varied[scenario]]
                df = pd.concat([df, df_sce_met], axis=0)
        return df

    @property
    def df_runtime(self, ):
        df = pd.DataFrame()
        for scenario, scenario_formal in self.scenarios.items():
            self.console.print("======>scenario: {}".format(scenario_formal))
            for method, method_formal in self.methods.items():
                self.console.print("=========>method: {}".format(method_formal))
                # print(self.params.work_dir + scenario + '/' + method + '.txt')
                df_sce_met = self.freader.generic(
                    df_fpn=self.params.work_dir + scenario + '/' + method + '.txt',
                    header=0,
                )
                # print(df_sce_met)
                if self.is_trans:
                    df_sce_met = (df_sce_met - self.umi_gt_cnt) / self.umi_gt_cnt
                df_sce_met['mean'] = df_sce_met.mean(axis=1)
                # print('1', df_sce_met)
                # print('asdsad', df_sce_met['mean'])
                df_sce_met['max'] = df_sce_met.max(axis=1)
                df_sce_met['min'] = df_sce_met.min(axis=1)
                df_sce_met['std'] = df_sce_met.std(axis=1)
                df_sce_met['mean-min'] = df_sce_met['std']
                df_sce_met['max-mean'] = df_sce_met['std']
                df_sce_met['scenario'] = scenario_formal
                df_sce_met['method'] = method_formal
                df_sce_met['metric'] = [str(x) for x in self.params.varied[scenario]]
                df = pd.concat([df, df_sce_met], axis=0)
        return df

    @property
    def df_dedup_set_cover_len(self, ):
        df = pd.DataFrame()
        for scenario, scenario_formal in self.scenarios.items():
            self.console.print("======>scenario: {}".format(scenario_formal))
            for method, method_formal in self.methods.items():
                self.console.print("=========>method: {}".format(method_formal))
                df_sce_met = self.freader.generic(
                    df_fpn=self.params.work_dir + scenario + '/' + method + '_dedup.txt',
                    header=0,
                )
                if self.is_trans:
                    df_sce_met = df_sce_met.applymap(lambda x: sum([int(x) for x in x.split(';')]))
                    # print(df_sce_met)
                    df_sce_met['mean'] = df_sce_met.mean(axis=1)
                    df_sce_met['max'] = df_sce_met.max(axis=1)
                    df_sce_met['min'] = df_sce_met.min(axis=1)
                    df_sce_met['std'] = df_sce_met.std(axis=1)
                    df_sce_met['mean-min'] = df_sce_met['std']
                    df_sce_met['max-mean'] = df_sce_met['std']
                df_sce_met['scenario'] = scenario_formal
                df_sce_met['method'] = method_formal
                df_sce_met['metric'] = [str(x) for x in self.params.varied[scenario]]
                df = pd.concat([df, df_sce_met], axis=0)
        return df

    @property
    def df_dedup_perm_melt(self, ):
        df = pd.DataFrame()
        for scenario, scenario_formal in self.scenarios.items():
            self.console.print("======>scenario: {}".format(scenario_formal))
            for method, method_formal in self.methods.items():
                self.console.print("=========>method: {}".format(method_formal))
                df_sce_met = self.freader.generic(
                    df_fpn=self.params.work_dir + scenario + '/' + method + '_dedup.txt',
                    header=0,
                )
                df_sce_met_T = df_sce_met.T
                df_sce_met_T = (df_sce_met_T - self.umi_gt_cnt) / self.umi_gt_cnt
                df_sce_met_T.columns = [str(x) for x in self.params.varied[scenario]]
                df_sce_met_T['method'] = method_formal
                df = pd.concat([df, df_sce_met_T], axis=0)
        df_melt = pd.melt(df, 'method', var_name="Sequencing error")
        return df_melt

    @property
    def df_dedup_melt(self, ):
        df = pd.melt(
            frame=self.df_dedup[['method', 'metric', 'scenario', 'mean']],
            id_vars=['method', 'metric', 'scenario'],
            # value_vars=['mean',],
        )
        return df

    @property
    def df_trace_cnt(self, ):
        df_apv = pd.DataFrame()
        df_disapv = pd.DataFrame()
        for scenario, scenario_formal in self.scenarios.items():
            self.console.print("======>scenario: {}".format(scenario_formal))
            for method, method_formal in self.methods.items():
                self.console.print("=========>method: {}".format(method_formal))
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

                df_sce_met_apv = df_sce_met_apv_perm.groupby(by=['metric']).agg({'diff_origin': 'mean', 'same_origin': 'mean'}).reset_index()
                df_sce_met_apv['diff_origin_max'] = df_sce_met_apv_perm.groupby(by=['metric']).agg({'diff_origin': 'max'})['diff_origin'].values
                df_sce_met_apv['diff_origin_min'] = df_sce_met_apv_perm.groupby(by=['metric']).agg({'diff_origin': 'min'})['diff_origin'].values
                df_sce_met_apv['same_origin_max'] = df_sce_met_apv_perm.groupby(by=['metric']).agg({'same_origin': 'max'})['same_origin'].values
                df_sce_met_apv['same_origin_min'] = df_sce_met_apv_perm.groupby(by=['metric']).agg({'same_origin': 'min'})['same_origin'].values
                df_sce_met_apv['scenario'] = scenario_formal
                df_sce_met_apv['method'] = method_formal

                df_sce_met_disapv = df_sce_met_disapv_perm.groupby(by=['metric']).agg({'diff_origin': 'mean', 'same_origin': 'mean'}).reset_index()
                df_sce_met_disapv['diff_origin_max'] = df_sce_met_disapv_perm.groupby(by=['metric']).agg({'diff_origin': 'max'})['diff_origin'].values
                df_sce_met_disapv['diff_origin_min'] = df_sce_met_disapv_perm.groupby(by=['metric']).agg({'diff_origin': 'min'})['diff_origin'].values
                df_sce_met_disapv['same_origin_max'] = df_sce_met_disapv_perm.groupby(by=['metric']).agg({'same_origin': 'max'})['same_origin'].values
                df_sce_met_disapv['same_origin_min'] = df_sce_met_disapv_perm.groupby(by=['metric']).agg({'same_origin': 'min'})['same_origin'].values
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
    def df_inflat_exp(self, ):
        df_inflat = pd.DataFrame()
        df_exp = pd.DataFrame()
        for scenario, scenario_formal in self.scenarios.items():
            self.console.print("======>scenario: {}".format(scenario_formal))
            df_inflat_sub = self.freader.generic(
                df_fpn=self.params.work_dir + scenario + '/inflat_val.txt',
                header=None,
            )
            df_inflat_sub = df_inflat_sub.rename(columns={0: 'id', 1: 'cnt'})
            df_exp_sub = self.freader.generic(
                df_fpn=self.params.work_dir + scenario + '/exp_val.txt',
                header=None,
            )
            df_exp_sub = df_exp_sub.rename(columns={0: 'id', 1: 'cnt'})
            # print(df_inflat_sub)
            # print(df_exp_sub)
            # print(scenario)
            df_inflat[scenario_formal] = df_inflat_sub['cnt']
            df_exp[scenario_formal] = df_exp_sub['cnt']
        df_inflat = (df_inflat - self.umi_gt_cnt) / self.umi_gt_cnt
        df_exp = (df_exp - self.umi_gt_cnt) / self.umi_gt_cnt
        df_inflat.index = df_inflat_sub['id'].values
        df_exp.index = df_exp_sub['id'].astype(int).values
        return df_inflat, df_exp


if __name__ == "__main__":
    from umiche.path import to

    p = Stat(
        scenarios={
            'pcr_nums': 'PCR cycle',
            'pcr_errs': 'PCR error',
            'seq_errs': 'Sequencing error',
            'ampl_rates': 'Amplification rate',
            'umi_lens': 'UMI length',
            'seq_deps': 'Sequencing depth',
        },

        methods={
            # 'unique': 'Unique',
            # 'cluster': 'Cluster',
            # 'adjacency': 'Adjacency',
            'directional': 'Directional',
            # 'dbscan_seq_onehot': 'DBSCAN',
            # 'birch_seq_onehot': 'Birch',
            # 'aprop_seq_onehot': 'Affinity Propagation',
            'mcl': 'MCL',
            # 'mcl_val': 'MCL-val',
            # 'mcl_ed': 'MCL-ed',
        },

        param_fpn=to('data/params.yml'),
        verbose=True
    )

    print(p.df_dedup)
    print(p.df_dedup_melt)
    print(p.df_trace_cnt)
    print(p.df_inflat_exp)