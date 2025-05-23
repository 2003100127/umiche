__version__ = "v1.0"
__copyright__ = "Copyright 2025"
__license__ = "GPL-3.0"
__developer__ = "Jianfeng Sun"
__maintainer__ = "Jianfeng Sun"
__email__="jianfeng.sunmt@gmail.com"

import pandas as pd
from umiche.simu.Parameter import Parameter as params
from umiche.fastq.Reader import Reader as fastqreader
from umiche.util.Console import Console


class Anchor:

    def __init__(
            self,
            scenario,
            param_fpn,
            verbose=False,
    ):
        self.fastqreader = fastqreader()
        self.scenario = scenario
        self.param_fpn = param_fpn
        self.params = params(param_fpn=param_fpn)
        print(self.params.anchor)
        self.verbose = verbose
        self.console = Console()
        self.console.verbose = self.verbose

    def pct_read_captured(self, ):
        without_anchor = {}
        with_anchor = {}
        for permut_i in range(self.params.fixed['permutation_num']):
            self.console.print("===>permutation test {}".format(permut_i))
            without_anchor[permut_i] = []
            with_anchor[permut_i] = []
            for i, criterion in enumerate(self.params.varied['criteria']):
                arr_fastq = self.fastqreader.fromgz(fastq_fpn=self.params.work_dir + '/permute_' + str(permut_i) + "/" + self.scenario + "_" + str(i) + '.fastq.gz')
                df = pd.DataFrame(arr_fastq).T
                self.console.print('======>number of reads at criterion {}: {}'.format(criterion, df.shape[0]))
                df['len'] = df[1].apply(lambda x: len(x))
                without_anchor[permut_i].append(df.loc[df['len'] == 92].shape[0]/df.shape[0])
                df['is_anchor'] = df[1].apply(lambda x: x.__contains__(self.params.anchor['seq1']))
                with_anchor[permut_i].append(df[df['is_anchor']].shape[0]/df.shape[0])
        return with_anchor, without_anchor


if __name__ == "__main__":
    from umiche.path import to

    anchor = Anchor(
        scenario='pcr_del',
        param_fpn=to('data/params_anchor.yml'),
    )