__version__ = "v1.0"
__copyright__ = "Copyright 2024"
__license__ = "MIT"
__developer__ = "Jianfeng Sun"
__maintainer__ = "Jianfeng Sun"
__email__="jianfeng.sunmt@gmail.com"
__lab__ = "Cribbslab"

import pandas as pd
from umiche.fastq.Convert import Convert as fastqconverter
from umiche.trim.Template import Template as trimmer
from umiche.graph.bfs.ConnectedComponent import ConnectedComponent as gbfscc
from umiche.simu.Parameter import Parameter as params

from umiche.deduplicate.MultiPos import MultiPos as deduppos
from umiche.plot.Heterogeneity import Heterogeneity as plotv
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
    ):
        self.scenario = scenario
        self.method = method

        self.params = params(param_fpn=param_fpn)
        self.gbfscc = gbfscc()
        self.fwriter = fwriter()
        self.plotv = plotv()

        self.verbose = verbose
        self.console = Console()
        self.console.verbose = self.verbose

        df_dedup = pd.DataFrame()
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
                    dedup_ob = deduppos(
                        bam_fpn=self.fastq_location + 'trimmed/bam/' + self.fn_prefix + '.bam',                        pos_tag='PO',
                        mcl_fold_thres=self.params.dedup['mcl_fold_thres'],
                        inflat_val=self.params.dedup['inflat_val'],
                        exp_val=self.params.dedup['exp_val'],
                        iter_num=self.params.dedup['iter_num'],
                        ed_thres=self.params.dedup['ed_thres'],
                        work_dir=self.params.work_dir,
                        verbose=False,
                    )
                    print(dedup_ob.directional().columns)
                    # dedup_arr.append(dedup_ob.dedup_num)
            # df_dedup['pn' + str(perm_num_i)] = dedup_arr
            # print(df_dedup)
        # self.fwriter.generic(
        #     df=df_dedup,
        #     sv_fpn=fastq_fp + self.scenario + '/' + str(self.method) + '_' + self.comp_cat + '.txt',
        #     header=True,
        # )


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
        method='directional',
        # method='mcl',
        # method='mcl_val',
        # method='mcl_ed',
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

        param_fpn=to('data/seqerr_sl.yml')
    )
