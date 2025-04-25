__version__ = "v1.0"
__copyright__ = "Copyright 2025"
__license__ = "GPL-3.0"
__developer__ = "Jianfeng Sun"
__maintainer__ = "Jianfeng Sun"
__email__="jianfeng.sunmt@gmail.com"


import warnings

import numpy as np

warnings.filterwarnings("ignore")

import json
import pandas as pd
from umiche.simu.Parameter import Parameter as params
from umiche.trim.Template import Template as trimmer
from umiche.fastq.Convert import Convert as fastqconverter

from umiche.deduplicate.OnePos import OnePos as dedupop
from umiche.plot.scenario.TraceSingle import TraceSingle as plothetero

from umiche.deduplicate.heterogeneity.Trace import Trace as umitrace
from umiche.bam.Relation import Relation as umirel

from umiche.util.Writer import Writer as fwriter
from umiche.util.Console import Console


class Heterogeneity:

    def __init__(
            self,
            scenario,
            method,
            param_fpn=None,
            is_trim=False,
            is_tobam=False,
            is_dedup=False,
            is_sv=False,
            verbose=False,
            **kwargs,
    ):
        self.scenario = scenario
        self.method = method
        self.is_sv = is_sv
        self.kwargs = kwargs
        self.verbose = verbose

        self.params = params(param_fpn=param_fpn)
        self.fwriter = fwriter()

        self.console = Console()
        self.console.verbose = self.verbose

        columns = ['diff_origin', 'same_origin', 'total', 'scenario', 'method', 'permutation']
        self.df_apv_cnt = pd.DataFrame(columns=columns)
        self.df_disapv_cnt = pd.DataFrame(columns=columns)
        self.df_apv_pct = pd.DataFrame(columns=columns)
        self.df_disapv_pct = pd.DataFrame(columns=columns)
        self.df_dedup = pd.DataFrame()
        self.node_repr_dict = {}
        print("===>dedup method: {}".format(self.method))
        for perm_num_i in range(self.params.fixed['permutation_num']):
            print("===>No.{} permutation".format(perm_num_i))
            self.console.print("===>permutation number {}".format(perm_num_i))
            dedup_arr = []
            self.node_repr_dict[perm_num_i] = {}
            for id, scenario_i in enumerate(self.params.varied[self.scenario]):
                if self.scenario == 'pcr_nums':
                    self.fn_mark = str(id) # scenario_i
                elif self.scenario == 'pcr_errs':
                    self.fn_mark = str(id)
                elif self.scenario == 'seq_errs':
                    self.fn_mark = str(id)
                elif self.scenario == 'ampl_rates':
                    self.fn_mark = str(id)
                elif self.scenario == 'seq_deps':
                    self.fn_mark = str(id)
                elif self.scenario == 'umi_lens':
                    self.fn_mark = str(id)
                elif self.scenario == 'umi_nums':
                    self.fn_mark = str(id)
                else:
                    self.fn_mark = str(scenario_i)
                self.fn_prefix = self.params.file_names[self.scenario] + self.fn_mark
                self.fastq_location = self.params.work_dir + self.scenario + '/permute_' + str(perm_num_i) + '/'
                if is_trim:
                    if self.scenario == 'umi_lens':
                        self.params.trimmed['umi_1']['len'] = scenario_i
                    self.console.print("======>fastq is being trimmed.")
                    self.params.trimmed['fastq']['fpn'] = self.fastq_location + self.fn_prefix + '.fastq.gz'
                    self.params.trimmed['fastq']['trimmed_fpn'] = self.fastq_location + 'trimmed/' + self.fn_prefix + '.fastq.gz'
                    umitrim_parser = trimmer(params=self.params.trimmed, verbose=self.verbose)
                    df = umitrim_parser.todf()
                    umitrim_parser.togz(df)
                if is_tobam:
                    self.console.print("======>fastq converter to bam is being used.")
                    fastqconverter(
                        fastq_fpn=self.fastq_location + 'trimmed/' + self.fn_prefix + '.fastq.gz',
                        bam_fpn=self.fastq_location + 'trimmed/bam/' + self.fn_prefix + '.bam',
                    ).tobam()
                if is_dedup:
                    self.console.print("======>reads are being deduplicated.")
                    # neither nor
                    if not (isinstance(self.params.dedup['inflat_val'], list) or isinstance(self.params.dedup['exp_val'], list)):
                        self.console.print("======>the working mode is not to optimise MCL inflation and expasion")
                        optim_arr = [self.params.dedup['inflat_val']]
                        inflat_flag = True
                    else:
                        self.console.print("======>the working mode is to optimise MCL inflation")
                        if isinstance(self.params.dedup['inflat_val'], list):
                            optim_arr = self.params.dedup['inflat_val']
                            inflat_flag = True
                        else:
                            self.console.print("======>the working mode is to optimise MCL expasion")
                            optim_arr = self.params.dedup['exp_val']
                            inflat_flag = False

                    for optim_val in optim_arr:
                        dedup_ob = dedupop(
                            # bam_fpn=self.fastq_location + 'trimmed/bam/pcr_num_12.bam',
                            # bam_fpn=self.fastq_location + 'trimmed/bam/pcr_err_10.bam',
                            # bam_fpn=self.fastq_location + 'trimmed/bam/seq_err_12.bam',
                            # bam_fpn=self.fastq_location + 'trimmed/bam/ampl_rate_8.bam',
                            # bam_fpn=self.fastq_location + 'trimmed/bam/umi_len_6.bam',
                            # bam_fpn=self.fastq_location + 'trimmed/bam/seq_dep_7.bam',
                            bam_fpn=self.fastq_location + 'trimmed/bam/' + self.fn_prefix + '.bam',
                            mcl_fold_thres=self.params.dedup['mcl_fold_thres'],
                            inflat_val=optim_val if inflat_flag else self.params.dedup['inflat_val'],
                            exp_val=optim_val if not inflat_flag else self.params.dedup['exp_val'],
                            # inflat_val=self.params.dedup['inflat_val'],
                            # exp_val=self.params.dedup['exp_val'],
                            iter_num=self.params.dedup['iter_num'],
                            ed_thres=self.params.dedup['ed_thres'],
                            dbscan_eps=self.params.dedup['dbscan_eps'],
                            dbscan_min_spl=self.params.dedup['dbscan_min_spl'],
                            birch_thres=self.params.dedup['birch_thres'],
                            birch_n_clusters=self.params.dedup['birch_n_clusters'],
                            work_dir=self.params.work_dir,
                            heterogeneity=True,
                            verbose=False,
                            **self.kwargs
                        )
                        df = self.tool(dedup_ob)[self.method]()
                        if self.method == "unique":
                            dedup_arr.append(df.num_uniq_umis.values[0])
                            print("No.{}->{} for {} dedup cnt: {}".format(
                                id,
                                scenario_i,
                                self.scenario,
                                df.num_uniq_umis.values[0]),
                            )
                        else:
                            dedup_arr.append(df.dedup_cnt.values[0])
                            print("No.{}->{} for {} dedup cnt: {}".format(
                                id,
                                scenario_i,
                                self.scenario,
                                df.dedup_cnt.values[0]),
                            )
                    # print(df.apv.values[0])
                    # print(len(df[self.method + '_repr_nodes'].loc[1]))

                        if self.method not in ['adjacency', 'cluster', 'unique']:
                            umiold = umirel(
                                df=dedup_ob.df_bam,
                                verbose=self.verbose,
                            )
                            umiidtrace = umitrace(
                                df_umi_uniq_val_cnt=umiold.df_umi_uniq_val_cnt,
                                umi_id_to_origin_id_dict=umiold.umi_id_to_origin_id_dict,
                            )
                            series_2d_arr_apv, series_2d_arr_disapv = umiidtrace.format(method=self.method, df=df)

                            series_dict_origin_apv = umiidtrace.match_representative(series_2d_arr=series_2d_arr_apv)
                            self.node_repr_dict[perm_num_i][scenario_i] = series_dict_origin_apv.to_dict()
                            # print(series_2d_arr_disapv)
                            if not series_2d_arr_apv.empty:
                                apv_cnt_dict = umiidtrace.edge_class(series_2d_arr=series_2d_arr_apv, sort='cnt')
                                # print(apv_cnt_dict)
                                apv_pct_dict = umiidtrace.edge_class(series_2d_arr=series_2d_arr_apv, sort='pct')
                                apv_cnt_dict['permutation'] = perm_num_i
                                apv_cnt_dict['method'] = self.method
                                apv_cnt_dict['scenario'] = scenario_i
                                apv_pct_dict['permutation'] = perm_num_i
                                apv_pct_dict['method'] = self.method
                                apv_pct_dict['scenario'] = scenario_i
                                self.df_apv_cnt = pd.concat([self.df_apv_cnt, pd.DataFrame.from_dict(apv_cnt_dict, orient='index').T]).reset_index(drop=True)
                                self.df_apv_pct = pd.concat([self.df_apv_pct, pd.DataFrame.from_dict(apv_pct_dict, orient='index').T]).reset_index(drop=True)
                            if not series_2d_arr_disapv.empty:
                                disapv_cnt_dict = umiidtrace.edge_class(series_2d_arr=series_2d_arr_disapv, sort='cnt')
                                disapv_pct_dict = umiidtrace.edge_class(series_2d_arr=series_2d_arr_disapv, sort='pct')
                                disapv_cnt_dict['permutation'] = perm_num_i
                                disapv_cnt_dict['method'] = self.method
                                disapv_cnt_dict['scenario'] = scenario_i
                                disapv_pct_dict['permutation'] = perm_num_i
                                disapv_pct_dict['method'] = self.method
                                disapv_pct_dict['scenario'] = scenario_i
                                self.df_disapv_cnt = pd.concat([self.df_disapv_cnt, pd.DataFrame.from_dict(disapv_cnt_dict, orient='index').T]).reset_index(drop=True)
                                self.df_disapv_pct = pd.concat([self.df_disapv_pct, pd.DataFrame.from_dict(disapv_pct_dict, orient='index').T]).reset_index(drop=True)

                            # print(self.df_apv_pct)
                            # plothetero(
                            #     df_apv=self.df_apv_cnt,
                            #     df_disapv=self.df_disapv_cnt,
                            # ).line_apv_disapv()
            self.df_dedup['pn' + str(perm_num_i)] = dedup_arr
            print(self.df_dedup)

        if self.is_sv:
            sv_dedup_fpn = self.params.work_dir + '/' + scenario + '/' + str(self.method) + '_dedup' + '.txt'
            sv_apv_cnt_fpn = self.params.work_dir + '/' + scenario + '/' + str(self.method) + '_apv_cnt' + '.txt'
            sv_disapv_cnt_fpn = self.params.work_dir + '/' + scenario + '/' + str(self.method) + '_disapv_cnt' + '.txt'
            sv_apv_pct_fpn = self.params.work_dir + '/' + scenario + '/' + str(self.method) + '_apv_pct' + '.txt'
            sv_disapv_pct_fpn = self.params.work_dir + '/' + scenario + '/' + str(self.method) + '_disapv_pct' + '.txt'
            sv_node_repr_fpn = self.params.work_dir + '/' + scenario + '/' + str(self.method) + '_node_repr' + '.json'
            self.fwriter.generic(df=self.df_dedup, sv_fpn=sv_dedup_fpn, header=True, )
            self.fwriter.generic(df=self.df_apv_cnt, sv_fpn=sv_apv_cnt_fpn, header=True, )
            self.fwriter.generic(df=self.df_disapv_cnt, sv_fpn=sv_disapv_cnt_fpn, header=True, )
            self.fwriter.generic(df=self.df_apv_pct, sv_fpn=sv_apv_pct_fpn, header=True, )
            self.fwriter.generic(df=self.df_disapv_pct, sv_fpn=sv_disapv_pct_fpn, header=True, )
            with open(sv_node_repr_fpn, 'w') as f:
                json.dump(self.node_repr_dict, f)

    def tool(self, dedup_ob):
        return {
            'unique': dedup_ob.unique,
            'cluster': dedup_ob.cluster,
            'adjacency': dedup_ob.adjacency,
            'directional': dedup_ob.directional,
            'mcl': dedup_ob.mcl,
            'mcl_val': dedup_ob.mcl_val,
            'mcl_ed': dedup_ob.mcl_ed,
            'mcl_cc_all_node_umis': dedup_ob.mcl_cc_all_node_umis,
            'dbscan_seq_onehot': dedup_ob.dbscan_seq_onehot,
            'birch_seq_onehot': dedup_ob.birch_seq_onehot,
            'hdbscan_seq_onehot': dedup_ob.hdbscan_seq_onehot,
            'aprop_seq_onehot': dedup_ob.aprop_seq_onehot,
            # 'set_cover': dedup_ob.set_cover,
        }

    @property
    def df_stat(self, ):
        return {
            i: pd.DataFrame(columns=['diff_origin', 'same_origin', 'total']) for i in [
                'df_apv_cnt',
                'df_disapv_cnt',
                'df_disapv_pct',
                'df_disapv_pct',
            ]}


if __name__ == "__main__":
    from umiche.path import to

    p = Heterogeneity(
        # scenario='pcr_nums',
        # scenario='pcr_errs',
        # scenario='seq_errs',
        # scenario='ampl_rates',
        # scenario='umi_lens',
        scenario='seq_deps',
        # scenario='umi_nums',

        # method='unique',
        # method='cluster',
        # method='adjacency',
        # method='directional',
        # method='mcl',
        # method='mcl_val',
        # method='mcl_ed',
        # method='mcl_cc_all_node_umis',
        # method='dbscan_seq_onehot',
        method='birch_seq_onehot',
        # method='aprop_seq_onehot',
        # method='hdbscan_seq_onehot',
        # method='set_cover',

        # is_trim=True,
        # is_tobam=False,
        # is_dedup=False,

        # is_trim=False,
        # is_tobam=True,
        # is_dedup=False,

        is_trim=False,
        is_tobam=False,
        is_dedup=True,

        is_sv=True,

        param_fpn=to('data/params.yml'),
    )