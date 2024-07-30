__version__ = "v1.0"
__copyright__ = "Copyright 2024"
__license__ = "MIT"
__developer__ = "Jianfeng Sun"
__maintainer__ = "Jianfeng Sun"
__email__="jianfeng.sunmt@gmail.com"
__lab__ = "Cribbslab"

from typing import Dict

from umiche.bam.Reader import Reader as alireader

def read_bam(
        bam_fpn : str,
        verbose : bool =True
):
    return alireader(
        bam_fpn=bam_fpn,
        verbose=verbose,
    )


if __name__ == "__main__":
    bam = read_bam(
        bam_fpn="/mnt/d/Document/Programming/Python/umiche/umiche/data/simu/umi/trimer/seq_errs/permute_0/trimmed/seq_err_17.bam",
        verbose=True,
    )
    print(bam.todf(tags=['PO']))