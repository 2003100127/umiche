__version__ = "v1.0"
__copyright__ = "Copyright 2025"
__license__ = "GPL-3.0"
__developer__ = "Jianfeng Sun"
__maintainer__ = "Jianfeng Sun"
__email__="jianfeng.sunmt@gmail.com"


import time
import textwrap
import pandas as pd
from umiche.deduplicate.method.trimer.Collapse import Collapse
from umiche.fastq.Reader import Reader as rfastq
from umiche.trim.Reader import Reader as trimreader
from umiche.util.Hamming import Hamming

from umiche.simu.Parameter import Parameter as params

from umiche.util.Reader import Reader as freader
from umiche.util.Writer import Writer as fwriter
from umiche.util.Console import Console


class selfHealing:
 
    def __init__(
            self,
            scenario,
            param_fpn=None,
            verbose=False,

            **kwargs,
    ):
        self.scenario = scenario
        self.kwargs = kwargs

        self.params = params(param_fpn=param_fpn)

        if 'fastq_fp' not in self.kwargs.keys():
            self.fastq_fp = self.params.fastq_fp
        else:
            self.fastq_fp = self.kwargs['fastq_fp']

        if 'umi_ref_fpn' not in self.kwargs.keys():
            self.umi_ref_fpn = self.params.umi_ref_fpn
        else:
            self.umi_ref_fpn = self.kwargs['umi_ref_fpn']

        self.rfastq = rfastq()
        self.collapse = Collapse()
        self.trimreader = trimreader()
        self.freader = freader()
        self.fwriter = fwriter()

        self.console = Console()
        self.console.verbose = verbose

        self.df_umi_ref_lib = self.freader.generic(df_fpn=self.umi_ref_fpn)
        self.df_umi_ref_lib = pd.DataFrame.from_dict({
            i: e for i, e in enumerate(self.df_umi_ref_lib[0].values)
        }, orient='index')
        self.console.print(self.df_umi_ref_lib)
        self.umi_monomer_ref_dict = self.df_umi_ref_lib[0].apply(lambda x: ''.join([i[0] for i in textwrap.wrap(x, 3)])).to_dict()
        self.console.print(self.umi_monomer_ref_dict)

    def rea(self, ) -> str:
        df_stat = pd.DataFrame()

        for perm_num_i in range(self.params.fixed['permutation_num']):
            for id, scenario_i in enumerate(self.params.varied[self.scenario]):

                read_stime = time.time()
                if self.scenario == 'pcr_nums':
                    self.fn_mark = str(id)  # scenario_i
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

                names, seqs, _, _ = self.rfastq.fromgz(
                    fastq_fpn=self.fastq_fp + self.fn_prefix + '.fastq.gz',
                )
                self.console.print('=========>read time: {:.3f}s'.format(time.time() - read_stime))

                df_fastq = self.trimreader.todf(names=names, seqs=seqs)

                mono_corr_stime = time.time()
                df_fastq['umi_mono_corr'] = df_fastq['umi'].apply(lambda x: self.collapse.majority_vote(x))
                # df_fastq['umi_mono_corr'] = df_fastq['umi'].apply(lambda x: self.correct(x))
                self.console.print('=========>mono_corr time: {:.3f}s'.format(time.time() - mono_corr_stime))

                hm_stime = time.time()
                df_stat['umi_hm' + str(id)] = df_fastq.apply(lambda x: Hamming().general(
                    x['umi_mono_corr'],
                    self.umi_monomer_ref_dict[x['umi#']],
                ), axis=1)
                self.console.print('=========>Hamming time: {:.3f}s'.format(time.time() - hm_stime))
                print(df_stat['umi_hm' + str(id)])

        self.fwriter.generic(
            df=df_stat,
            sv_fpn=to('data/simu/umi/seq_errs/trimer/trimmed/dasd1111.txt'),
        )
        return 'Finished.'

    def stat(self, x):
        return Hamming().general(x['umi_mono_corr'], self.umi_monomer_ref_dict[x['umi#']])


if __name__ == "__main__":
    from umiche.path import to

    p = selfHealing(
        # scenario='pcr_nums',
        # scenario='pcr_errs',
        scenario='seq_errs',
        # scenario='ampl_rates',
        # scenario='umi_lens',
        # scenario='seq_deps',
        # scenario='umi_nums',

        # umi_ref_fpn=to('data/simu/umi/trimer/seq_errs/umi.txt'),
        # fastq_fp=to('data/simu/umi/trimer/seq_errs/trimmed/'),

        param_fpn=to('data/trimer.yml'),

        verbose=True, # False True
    )

    print(p.rea())