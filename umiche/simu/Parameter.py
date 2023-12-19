__version__ = "v1.0"
__copyright__ = "Copyright 2024"
__license__ = "MIT"
__developer__ = "Jianfeng Sun"
__maintainer__ = "Jianfeng Sun"
__email__="jianfeng.sunmt@gmail.com"
__lab__ = "Cribbslab"

import yaml
import numpy as np
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
        print(np.arange(self.params['varied']['sd']))
    @property
    def fixed(self, ):
        return self.params['fixed']

    @property
    def varied(self, ):
        return self.params['varied']

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