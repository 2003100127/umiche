__version__ = "v1.0"
__copyright__ = "Copyright 2025"
__license__ = "GPL-3.0"
__developer__ = "Jianfeng Sun"
__maintainer__ = "Jianfeng Sun"
__email__="jianfeng.sunmt@gmail.com"


from typing import Any

from umiche.deduplicate.io.Stat import Simulation as statsimu
from umiche.util.Reader import Reader as freader


def read(
        df_fpn,
        df_sep : str = '\t',
        skiprows : Any = None,
        header : Any = None,
        index_col : Any = None,
        is_utf8 : bool = False,
        sheet_name : str = 'Sheel1',
        engine : str = 'openpyxl',
        type : str = 'csv',
):
    m = {
        "csv": freader().generic,
        "tsv": freader().generic,
        "txt": freader().generic,
        "excel": freader().excel,
    }
    return m[type](
        df_fpn=df_fpn,
        df_sep=df_sep,
        skiprows=skiprows,
        header=header,
        index_col=index_col,
        is_utf8=is_utf8,
        sheet_name=sheet_name,
        engine=engine,
    )


def read_bam(
        bam_fpn : str,
        verbose : bool =True
):
    from umiche.bam.Reader import Reader as alireader
    return alireader(
        bam_fpn=bam_fpn,
        verbose=verbose,
    )


def stat(
        scenarios,
        methods,
        param_fpn,
        verbose : bool =True
):
    return statsimu(
        scenarios=scenarios,
        methods=methods,
        param_fpn=param_fpn,
        verbose=verbose,
    )


if __name__ == "__main__":
    from umiche.path import to

    # bam = read_bam(
    #     bam_fpn="/mnt/d/Document/Programming/Python/umiche/umiche/data/simu/umi/trimer/seq_errs/permute_0/trimmed/seq_err_17.bam",
    #     verbose=True,
    # )
    # print(bam.todf(tags=['PO']))

    stat_p = stat(
        scenarios={
            'pcr_nums': 'PCR cycle',
            'pcr_errs': 'PCR error',
            'seq_errs': 'Sequencing error',
            'ampl_rates': 'Amplification rate',
            'umi_lens': 'UMI length',
            'seq_deps': 'Sequencing depth',
        },
        methods={
            # 'unique': 'Unique',
            # 'cluster': 'Cluster',
            # 'adjacency': 'Adjacency',
            'directional': 'Directional',
            # 'dbscan_seq_onehot': 'DBSCAN',
            # 'birch_seq_onehot': 'Birch',
            # 'aprop_seq_onehot': 'Affinity Propagation',
            'mcl': 'MCL',
            # 'mcl_val': 'MCL-val',
            # 'mcl_ed': 'MCL-ed',
        },
        param_fpn=to('data/params.yml'),
        verbose=True,
    )

    print(stat_p.df_dedup)
    print(stat_p.df_trace_cnt)
    print(stat_p.df_inflat_exp)