__version__ = "v1.0"
__copyright__ = "Copyright 2025"
__license__ = "GPL-3.0"
__developer__ = "Jianfeng Sun"
__maintainer__ = "Jianfeng Sun"
__email__ = "jianfeng.sunmt@gmail.com"



from umiche.trim.Template import Template


def template(
        params,
):
    trimmer = Template(params)
    df = trimmer.todf()
    trimmer.togz(df)
    return df


if __name__ == "__main__":
    from umiche.path import to

    params = {
        'umi_1': {
            'len': 10,
        },
        'umi_2': {
            'len': 4,
        },
        'bc_1': {
            'len': 2,
        },
        'read_struct': 'umi_1',
        # 'read_struct': 'umi_1+seq_1',
        # 'read_struct': 'bc_1+umi_1+seq_1',
        'seq_1': {
            'len': 6,
        },
        'fastq': {
            'fpn': to('data/simu/mclumi/seq_errs/permute_0/seq_err_5.fastq.gz'),
            'trimmed_fpn': to('data/simu/mclumi/seq_errs/permute_0/seq_err_5_trimmed.fastq.gz'),
        },
    }
    print(template(params=params))