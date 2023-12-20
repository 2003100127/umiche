__version__ = "v1.0"
__copyright__ = "Copyright 2024"
__license__ = "MIT"
__developer__ = "Jianfeng Sun"
__maintainer__ = "Jianfeng Sun"
__email__="jianfeng.sunmt@gmail.com"
__lab__ = "Cribbslab"

import pysam
from umiche.util.Console import Console


class Writer:

    def __init__(
            self,
             df,
            verbose=False,
    ):
        self.df = df
        self.console = Console()
        self.console.verbose = verbose

    def tobam(
            self,
            tobam_fpn,
            tmpl_bam_fpn,
            whitelist=[],
    ):
        tmpl_bam = pysam.AlignmentFile(tmpl_bam_fpn, "rb")
        write_to_bam = pysam.AlignmentFile(tobam_fpn, "wb", template=tmpl_bam)
        fs = self.df.loc[self.df['id'].isin(whitelist)]['read']
        for i in fs:
            # print(i)
            write_to_bam.write(i)
        write_to_bam.close()
        return write_to_bam