__version__ = "v1.0"
__copyright__ = "Copyright 2025"
__license__ = "GPL-3.0"
__developer__ = "Jianfeng Sun"
__maintainer__ = "Jianfeng Sun"
__email__="jianfeng.sunmt@gmail.com"


import yaml
from umiche.util.Console import Console


class Parameter:

    def __init__(
            self,
            param_fpn,
            verbose=True,
    ):
        self.console = Console()
        self.console.verbose = verbose
        with open(param_fpn, "r") as f:
            self.params = yaml.safe_load(f)
            for i, (k, item) in enumerate(self.params.items()):
                self.console.print("======>key {}: {}".format(i+1, k))
                self.console.print("=========>value: {}".format(item))

    @property
    def fixed(self, ):
        return self.params['fixed']

    @property
    def varied(self, ):
        return self.params['varied']

    @property
    def trimmed(self, ):
        return self.params['trimmed']

    @property
    def dedup(self, ):
        return self.params['dedup']

    @property
    def work_dir(self, ):
        return self.params['work_dir']

    @property
    def anchor(self, ):
        return self.params['anchor']

    @property
    def fastq_fp(self, ):
        return self.params['fastq_fp']

    @property
    def umi_ref_fpn(self, ):
        return self.params['umi_ref_fpn']

    @property
    def file_names(self, ):
        return {
            'pcr_nums': 'pcr_num_',
            'pcr_errs': 'pcr_err_',
            'seq_errs': 'seq_err_',
            'ampl_rates': 'ampl_rate_',
            'umi_lens': 'umi_len_',
            'seq_deps': 'seq_dep_',
            'umi_nums': 'umi_num_',
        }


if __name__ == "__main__":
    p = Parameter(
        param_fpn='./params/param_fpn.txt'
    )
    print(p.file_names)