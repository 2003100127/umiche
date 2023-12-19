__version__ = "v1.0"
__copyright__ = "Copyright 2024"
__license__ = "MIT"
__developer__ = "Jianfeng Sun"
__maintainer__ = "Jianfeng Sun"
__email__="jianfeng.sunmt@gmail.com"
__lab__ = "Cribbslab"

import numpy as np


class Monomer:

    def __init__(
            self,
            param_fpn=None,
    ):
        if param_fpn:
            import yaml
            with open(param_fpn, "r") as f:
                configs = yaml.safe_load(f)
                for k, item in configs.items():
                    print(k, item)
            # len_params = configs['len_params']
            # kwargs['seq_params'] = configs['seq_params']
            # kwargs['material_params'] = configs['material_params']
            # condis = configs['condis']
            #
            # ampl_rate = configs['ampl_rate']
            # err_route = configs['err_route']
            # pcr_error = configs['pcr_error']
            # pcr_num = configs['pcr_num']
            # err_num_met = configs['err_num_met']
            # seq_errors = configs['seq_errors']
            # seq_sub_spl_number = configs['seq_sub_spl_number']
            # seq_sub_spl_rate = configs['seq_sub_spl_rate']
            # if "pcr_deletion" in configs.keys():
            #     pcr_deletion = configs['pcr_deletion']
            # if "pcr_insertion" in configs.keys():
            #     pcr_insertion = configs['pcr_insertion']
            # if "pcr_del_rate" in configs.keys():
            #     pcr_del_rate = configs['pcr_del_rate']
            # if "pcr_ins_rate" in configs.keys():
            #     pcr_ins_rate = configs['pcr_ins_rate']
            # if "seq_deletion" in configs.keys():
            #     seq_deletion = configs['seq_deletion']
            # if "seq_insertion" in configs.keys():
            #     seq_insertion = configs['seq_insertion']
            # if "seq_del_rate" in configs.keys():
            #     seq_del_rate = configs['seq_del_rate']
            # if "seq_ins_rate" in configs.keys():
            #     seq_ins_rate = configs['seq_ins_rate']

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

    @property
    def fixed(self, ):
        return {
            'pcr_num': 8,
            'pcr_err': 1e-5,
            'seq_err': 1e-3,
            'ampl_rate': 0.85,
            'umi_len': 36,
            'seq_dep': 400,
            'umi_num': 50,
            'permutation_num': 10,
            'umi_unit_pattern': 1,
            'umi_unit_len': 36,
            'seq_sub_spl_rate': 1,
            'sim_thres': 3,
        }

    @property
    def varied(self, ):
        return {
            'pcr_nums': np.arange(1, 20 + 1, 1),
            'pcr_errs': self.errors(),
            'seq_errs': self.errors(),
            'ampl_rates': np.linspace(0.1, 1, 10),
            'umi_lens': np.arange(6, 36 + 1, 1),
            'seq_deps': np.arange(100, 1000 + 100, 100),
            'umi_nums': np.arange(50, 1000 + 200, 200),
        }

    def errors(self, ):
        arr = []
        e = 1e-5
        while e < 3e-1:
            arr.append(e)
            if 5 * e < 3e-1:
                arr.append(2.5 * e)
                arr.append(5 * e)
                arr.append(7.5 * e)
            e = 10 * e
        arr.append(0.125)
        arr.append(0.15)
        arr.append(0.2)
        arr.append(0.225)
        arr.append(0.25)
        arr.append(0.3)
        return arr

    def errs(self, ):
        return [
            1e-05,
            2.5e-05,
            5e-05,
            7.5e-05,
            0.0001,
            0.00025,
            0.0005,
            0.00075,
            0.001,
            0.0025,
            0.005,
            0.0075,
            0.01,
        ]


if __name__ == "__main__":
    p = Monomer()

    print(p.file_names)

    print(p.varied)

    print(p.fixed)