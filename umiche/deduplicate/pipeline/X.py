__version__ = "v1.0"
__copyright__ = "Copyright 2024"
__license__ = "MIT"
__developer__ = "Jianfeng Sun"
__maintainer__ = "Jianfeng Sun"
__email__="jianfeng.sunmt@gmail.com"
__lab__ = "Cribbslab"

import pandas as pd
from umiche.fastq.Convert import convert as fas2bam
from umiche.trim.Template import template as umitrim
from umiche.graph.bfs.ConnectedComponent import ConnectedComponent as gbfscc
from umiche.simu.Parameter import Parameter as params

from umiche.deduplicate.MultiPos import MultiPos as deduppos
from umiche.plot.Valid import valid as plotv
from umiche.util.Writer import writer as gwriter
from umiche.util.Console import Console


class umi:

    def __init__(
            self,
            scenario,
            method,
            comp_cat,
            param_fpn=None,
            is_trim=False,
            is_tobam=False,
            is_dedup=False,
            verbose=False,
    ):
        self.scenario = scenario
        self.comp_cat = comp_cat
        self.method = method

        self.params = params(param_fpn=param_fpn)
        self.gbfscc = gbfscc()
        self.gwriter = gwriter()
        self.plotv = plotv()

        self.console = Console()
        self.console.verbose = verbose

        df_dedup = pd.DataFrame()
        for perm_num_i in range(self.params.fixed['permutation_num']):
            dedup_arr = []
            for id, scenario_i in enumerate(self.params.varied[self.scenario]):
                if self.scenario == 'pcr_nums':
                    self.console.print('===========>at PCR {}'.format(scenario_i))
                    fn_surf = str(scenario_i)
                    self.mcl_inflat = 2.3
                    self.mcl_exp = 2
                    self.mcl_fold_thres = 1
                elif self.scenario == 'pcr_errs':
                    self.mcl_inflat = scenario_i
                    self.console.print('===========>No.{} PCR error: {}'.format(id, scenario_i))
                    fn_surf = str(id)
                    # # /*** mcl_ed params ***/
                    # self.mcl_inflat = 1.1 if scenario_i > 0.005 else 1.7
                    # self.mcl_exp = 2
                    # self.mcl_fold_thres = 1

                    # # /*** mcl_val params ***/
                    self.mcl_inflat = 1.1 if scenario_i > 0.005 else 1.8
                    self.mcl_exp = 2
                    self.mcl_fold_thres = 2
                elif self.scenario == 'seq_errs':
                    self.console.print('===========>No.{} sequencing error: {}'.format(id, scenario_i))
                    # self.mcl_inflat = 1.1 if scenario_i > 0.005 else 2.7
                    # self.mcl_exp = 3
                    self.mcl_fold_thres = 1.6
                    self.mcl_inflat = 1.1 if scenario_i > 0.005 else 2.7
                    self.mcl_exp = 2
                    fn_surf = str(id)
                elif self.scenario == 'ampl_rates':
                    self.console.print('===========>No.{} amplification rate: {}'.format(id, scenario_i))
                    fn_surf = str(id)
                    # self.mcl_inflat = 1.3 if scenario_i > 0.5 else 2
                    # # /*** mcl_ed params ***/
                    # if scenario_i < 8:
                    #     self.mcl_inflat = 4
                    # if scenario_i >= 8 and scenario_i <= 11:
                    #     self.mcl_inflat = 2.3
                    # if scenario_i > 11:
                    #     self.mcl_inflat = 1.1
                    # self.mcl_exp = 3
                    # self.mcl_fold_thres = 1

                    # /*** mcl_val params ***/
                    if scenario_i < 8:
                        self.mcl_inflat = 2
                    if scenario_i >= 0.9:
                        self.mcl_inflat = 1.8
                    self.mcl_exp = 4
                    self.mcl_fold_thres = 11

                elif self.scenario == 'umi_lens':
                    self.console.print('===========>No.{} UMI length: {}'.format(id, scenario_i))
                    fn_surf = str(scenario_i)
                    # self.mcl_inflat = 1.1 if scenario_i > 11 else 2.3
                    # # /*** mcl_ed params ***/
                    # if scenario_i < 8:
                    #     self.mcl_inflat = 4
                    # if scenario_i >= 8 and scenario_i <= 11:
                    #     self.mcl_inflat = 2.3
                    # if scenario_i > 11:
                    #     self.mcl_inflat = 1.1
                    # self.mcl_exp = 3
                    # self.mcl_fold_thres = 1

                    # # # /*** mcl_val params ***/
                    # if scenario_i < 8:
                    #     self.mcl_inflat = 6
                    # if scenario_i >= 8 and scenario_i <= 11:
                    #     self.mcl_inflat = 2.3
                    # if scenario_i > 11:
                    #     self.mcl_inflat = 1.1
                    # self.mcl_exp = 4
                    # self.mcl_fold_thres = 11

                    # # /*** mcl_val params ***/
                    if scenario_i < 8:
                        self.mcl_inflat = 5.8
                        self.mcl_exp = 6
                    if scenario_i >= 8 and scenario_i <= 11:
                        self.mcl_inflat = 2.3
                        self.mcl_exp = 4
                    if scenario_i > 11:
                        self.mcl_inflat = 1.1
                        self.mcl_exp = 4
                    self.mcl_fold_thres = 11

                    self.umi_len = scenario_i
                else:
                    fn_surf = str(scenario_i)
                fn = self.fn_pref[self.scenario] + fn_surf
                if is_trim:
                    self.trim(
                        fastq_fpn=fastq_fp + self.scenario + '/permute_' + str(perm_num_i) + '/' + fn,
                        fastq_trimmed_fpn=fastq_fp + self.scenario + '/permute_' + str(perm_num_i) + '/trimmed/' + fn,
                        umi_len=self.umi_len,
                    )
                if is_tobam:
                    fas2bam(
                        fastq_fpn=fastq_fp + self.scenario + '/permute_' + str(perm_num_i) + '/trimmed/' + fn + '.fastq.gz',
                        # fastq_fpn=fastq_fp + self.scenario + '/permute_' + str(perm_num_i) + '/' + self.comp_cat + '/' + fn + '.fastq.gz',
                        bam_fpn=fastq_fp + self.scenario + '/permute_' + str(perm_num_i) + '/' + self.comp_cat + '/bam/' + fn,
                    ).tobam()
                if is_dedup:
                    # if self.scenario == 'seq_errs':
                    #     if scenario_i == 0.125 or scenario_i == 0.15:
                    #         continue
                    #     else:
                    dedup_ob = deduppos(
                        mode='internal',
                        method=self.method,
                        # bam_fpn=to('example/data/example.bam'),
                        bam_fpn=fastq_fp + self.scenario + '/permute_' + str(perm_num_i) + '/' + self.comp_cat + '/bam/' + fn + '.bam',
                        pos_tag='PO',
                        mcl_fold_thres=self.mcl_fold_thres,
                        inflat_val=self.mcl_inflat,
                        exp_val=self.mcl_exp,
                        iter_num=100,
                        verbose=False,
                        ed_thres=1,
                        is_sv=False,
                        sv_fpn=fastq_fp + self.scenario + '/permute_' + str(perm_num_i) + '/' + self.comp_cat + '/' + fn,
                    )
                    dedup_arr.append(dedup_ob.dedup_num)
            df_dedup['pn' + str(perm_num_i)] = dedup_arr
            print(df_dedup)
        self.gwriter.generic(
            df=df_dedup,
            sv_fpn=fastq_fp + self.scenario + '/' + str(self.method) + '_' + self.comp_cat + '.txt',
            header=True,
        )

    def trim(self, fastq_fpn, fastq_trimmed_fpn, umi_len):
        trim_params = {
            'read_struct': 'umi_1',
            'umi_1': {
                'len': umi_len,
            },
            'fastq': {
                'fpn': fastq_fpn + '.fastq.gz',
                'trimmed_fpn': fastq_trimmed_fpn + '.fastq.gz',
            },
        }
        umitrim_parser = umitrim(trim_params)
        df = umitrim_parser.todf()
        umitrim_parser.togz(df)
        return 0


if __name__ == "__main__":
    from umiche.path import to

    p = umi(
        # scenario='pcr_nums',
        # scenario='pcr_errs',
        scenario='seq_errs',
        # scenario='ampl_rates',
        # scenario='umi_lens',

        # method='unique',
        # method='cluster',
        # method='adjacency',
        method='directional',
        # method='mcl',
        # method='mcl_val',
        # method='mcl_ed',
        # method='set_cover',

        comp_cat='',
        # comp_cat='trimmed',
        # comp_cat='lmr',
        # comp_cat='ref',
        # comp_cat='bipartite',
        # comp_cat='double',
        # comp_cat='correct',

        # is_trim=True,
        # is_tobam=False,
        # is_dedup=False,

        # is_trim=False,
        # is_tobam=True,
        # is_dedup=False,

        is_trim=False,
        is_tobam=False,
        is_dedup=True,

        param_fpn=to('data/seqerr_sl.yml')
    )
