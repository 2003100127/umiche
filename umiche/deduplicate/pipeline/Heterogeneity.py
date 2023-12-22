__version__ = "v1.0"
__copyright__ = "Copyright 2024"
__license__ = "MIT"
__developer__ = "Jianfeng Sun"
__maintainer__ = "Jianfeng Sun"
__email__="jianfeng.sunmt@gmail.com"
__lab__ = "Cribbslab"

import pandas as pd
from umiche.simu.Parameter import Parameter as params
from umiche.trim.Template import Template as trimmer
from umiche.fastq.Convert import Convert as fastqconverter

from umiche.deduplicate.OnePos import OnePos as dedupop
from umiche.plot.Heterogeneity import Heterogeneity as plothetero

from umiche.deduplicate.heterogeneity.Trace import Trace as umitrace
from umiche.bam.Relation import Relation as umirel

from umiche.util.Writer import Writer as fwriter
from umiche.util.Console import Console


class Simulation:

    def __init__(
            self,
            scenario,
            method,
            param_fpn=None,
            is_trim=False,
            is_tobam=False,
            is_dedup=False,
            verbose=False,
            **kwargs,
    ):
        self.scenario = scenario
        self.method = method
        self.kwargs = kwargs

        self.params = params(param_fpn=param_fpn)
        self.fwriter = fwriter()
        self.plothetero = plothetero()

        self.verbose = verbose
        self.console = Console()
        self.console.verbose = self.verbose

        columns = ['diff_origin', 'same_origin', 'total', 'scenario', 'method', 'permutation']
        self.df_apv_cnt = pd.DataFrame(columns=columns)
        self.df_disapv_cnt = pd.DataFrame(columns=columns)
        self.df_apv_pct = pd.DataFrame(columns=columns)
        self.df_disapv_pct = pd.DataFrame(columns=columns)
        self.df_dedup = pd.DataFrame()
        for perm_num_i in range(self.params.fixed['permutation_num']):
            self.console.print("===>permutation number {}".format(perm_num_i))
            dedup_arr = []
            for id, scenario_i in enumerate(self.params.varied[self.scenario]):
                if self.scenario == 'pcr_nums':
                    self.fn_mark = str(scenario_i)
                elif self.scenario == 'pcr_errs':
                    self.fn_mark = str(id)
                elif self.scenario == 'seq_errs':
                    self.fn_mark = str(id)
                elif self.scenario == 'ampl_rates':
                    self.fn_mark = str(id)
                elif self.scenario == 'umi_lens':
                    self.fn_mark = str(scenario_i)
                else:
                    self.fn_mark = str(scenario_i)
                self.fn_prefix = self.params.file_names[self.scenario] + self.fn_mark
                self.fastq_location = self.params.work_dir + self.scenario + '/permute_' + str(perm_num_i) + '/'
                if is_trim:
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
                    dedup_ob = dedupop(
                        bam_fpn=self.fastq_location + 'trimmed/bam/' + self.fn_prefix + '.bam',
                        # pos_tag='PO',
                        mcl_fold_thres=self.params.dedup['mcl_fold_thres'],
                        inflat_val=self.params.dedup['inflat_val'],
                        exp_val=self.params.dedup['exp_val'],
                        iter_num=self.params.dedup['iter_num'],
                        ed_thres=self.params.dedup['ed_thres'],
                        work_dir=self.params.work_dir,
                        heterogeneity=True,
                        verbose=False,
                        # kwargs=self.kwargs
                    )
                    df = self.tool(dedup_ob)[self.method]()
                    print(df.dedup_cnt.values[0])
                    dedup_arr.append(df.dedup_cnt.values[0])
                    # print(df.apv.values[0])

                    # umiold = umirel(
                    #     df=dedup_ob.df_bam,
                    #     verbose=self.verbose,
                    # )
                    # umiidtrace = umitrace(
                    #     df_umi_uniq_val_cnt=umiold.df_umi_uniq_val_cnt,
                    #     umi_id_to_origin_id_dict=umiold.umi_id_to_origin_id_dict,
                    # )
                    # if self.method == 'directional':
                    #     series_2d_arr_apv = umiidtrace.format_apv_disapv(df.apv.values[0])
                    #     series_2d_arr_disapv = umiidtrace.format_apv_disapv(df.disapv.values[0])
                    # else:
                    #     series_2d_arr_apv = df.apv.values[0]
                    #     series_2d_arr_disapv = df.disapv.values[0]
                    # print(series_2d_arr_apv)
                    # # print(series_2d_arr_disapv)
                    # apv_cnt_dict = umiidtrace.edge_class(series_2d_arr=series_2d_arr_apv, sort='cnt')
                    # disapv_cnt_dict = umiidtrace.edge_class(series_2d_arr=series_2d_arr_disapv, sort='cnt')
                    # apv_pct_dict = umiidtrace.edge_class(series_2d_arr=series_2d_arr_apv, sort='pct')
                    # disapv_pct_dict = umiidtrace.edge_class(series_2d_arr=series_2d_arr_disapv, sort='pct')
                    # apv_cnt_dict['permutation'], apv_cnt_dict['method'], apv_cnt_dict['scenario'] = perm_num_i, self.method, scenario_i
                    # disapv_cnt_dict['permutation'], disapv_cnt_dict['method'], disapv_cnt_dict['scenario'] = perm_num_i, self.method, scenario_i
                    # apv_pct_dict['permutation'], apv_pct_dict['method'], apv_pct_dict['scenario'] = perm_num_i, self.method, scenario_i
                    # disapv_pct_dict['permutation'], disapv_pct_dict['method'], disapv_pct_dict['scenario'] = perm_num_i, self.method, scenario_i
                    #
                    # self.df_apv_cnt = pd.concat([self.df_apv_cnt, pd.DataFrame.from_dict(apv_cnt_dict, orient='index').T]).reset_index(drop=True)
                    # self.df_disapv_cnt = pd.concat([self.df_disapv_cnt, pd.DataFrame.from_dict(disapv_cnt_dict, orient='index').T]).reset_index(drop=True)
                    # self.df_apv_pct = pd.concat([self.df_apv_pct, pd.DataFrame.from_dict(apv_pct_dict, orient='index').T]).reset_index(drop=True)
                    # self.df_disapv_pct = pd.concat([self.df_disapv_pct, pd.DataFrame.from_dict(disapv_pct_dict, orient='index').T]).reset_index(drop=True)
                    # print(self.df_apv_pct)

                    # self.plothetero.n1(
                    #     df_apv=self.df_apv_pct,
                    #     df_disapv=self.df_disapv_pct,
                    # )
            self.df_dedup['pn' + str(perm_num_i)] = dedup_arr
            # print(df_dedup)
        sv_dedup_fpn = self.params.work_dir + '/' + scenario + '/' + str(self.method) + '_dedup' + '.txt'
        sv_apv_cnt_fpn = self.params.work_dir + '/' + scenario + '/' + str(self.method) + '_apv_cnt' + '.txt'
        sv_disapv_cnt_fpn = self.params.work_dir + '/' + scenario + '/' + str(self.method) + '_disapv_cnt' + '.txt'
        sv_apv_pct_fpn = self.params.work_dir + '/' + scenario + '/' + str(self.method) + '_apv_pct' + '.txt'
        sv_disapv_pct_fpn = self.params.work_dir + '/' + scenario + '/' + str(self.method) + '_disapv_pct' + '.txt'
        self.fwriter.generic(df=self.df_dedup, sv_fpn=sv_dedup_fpn, header=True, )
        self.fwriter.generic(df=self.df_apv_cnt, sv_fpn=sv_apv_cnt_fpn, header=True, )
        self.fwriter.generic(df=self.df_disapv_cnt, sv_fpn=sv_disapv_cnt_fpn, header=True, )
        self.fwriter.generic(df=self.df_apv_pct, sv_fpn=sv_apv_pct_fpn, header=True, )
        self.fwriter.generic(df=self.df_disapv_pct, sv_fpn=sv_disapv_pct_fpn, header=True, )

    def tool(self, dedup_ob):
        return {
            'unique': dedup_ob.unique,
            'cluster': dedup_ob.cluster,
            'adjacency': dedup_ob.adjacency,
            'directional': dedup_ob.directional,
            'mcl': dedup_ob.mcl,
            'mcl_val': dedup_ob.mcl_val,
            'mcl_ed': dedup_ob.mcl_ed,
            'dbscan_seq_onehot': dedup_ob.dbscan_seq_onehot,
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

    p = Simulation(
        # scenario='pcr_nums',
        # scenario='pcr_errs',
        scenario='seq_errs',
        # scenario='ampl_rates',
        # scenario='umi_lens',

        # method='unique',
        # method='cluster',
        # method='adjacency',
        # method='directional',
        # method='mcl',
        # method='mcl_val',
        # method='mcl_ed',
        method='dbscan_seq_onehot',
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
        dbscan_eps=1.5,
        dbscan_min_spl=1,

        param_fpn=to('data/seqerr_sl.yml')
    )
