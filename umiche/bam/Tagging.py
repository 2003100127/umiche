__version__ = "v1.0"
__copyright__ = "Copyright 2025"
__license__ = "GPL-3.0"
__developer__ = "Jianfeng Sun"
__maintainer__ = "Jianfeng Sun"
__email__="jianfeng.sunmt@gmail.com"

import re

from umiche.util.Console import Console


class Tagging:

    def __init__(
            self,
            bam_fpn: str,
            sv_bam_fpn: str = 'tagged.bam',
            delimiter = r'[_:]',
            bc_idx: int = -2,
            umi_idx: int = -1,
            bc_tag_name: str = 'CB',
            umi_tag_name: str = 'MB',
            verbose: bool = False,
    ):
        self.bam_fpn = bam_fpn
        self.sv_bam_fpn = sv_bam_fpn

        self.delimiter = delimiter
        self.bc_idx = bc_idx
        self.umi_idx = umi_idx
        self.bc_tag_name = bc_tag_name
        self.umi_tag_name = umi_tag_name

        import pysam
        self.pysam = pysam

        self.console = Console()
        self.console.verbose = verbose

    @Console.vignette()
    def add_tags(
            self,
            strip_after: bool = False,
    ) -> str:
        with self.pysam.AlignmentFile(self.bam_fpn, "rb") as bam_in, \
                self.pysam.AlignmentFile(self.sv_bam_fpn, "wb", template=bam_in) as bam_out:

            for read in self.console._tqdm(
                    bam_in,
                    # total=bam_in.count(until_eof=True),
                    desc=f"[{self.bc_tag_name} & {self.umi_tag_name} tagging]",
                    unit="reads",
                    position=0,
                    leave=True,
                    dynamic_ncols=False,
            ):
                parts = re.split(self.delimiter, read.query_name)
                # print(parts)
                if len(parts) <= max(self.bc_idx, self.umi_idx):
                    # 当 read name is not qualified
                    bam_out.write(read)
                    continue

                bc = parts[self.bc_idx]
                umi = parts[self.umi_idx]
                # print(bc)
                # print(umi)

                read.set_tag(self.bc_tag_name, bc, value_type="Z")
                read.set_tag(self.umi_tag_name, umi, value_type="Z")
                if strip_after:
                    read.query_name = parts[0]

                bam_out.write(read)

        return 'Complete.'


if __name__ == "__main__":
    from umiche.path import to
    #  samtools view /mnt/d/Document/Programming/Python/umiche/umiche/data/r1/trumicount/10xn9k_10c_tagged.bam | head -n 5
    p = Tagging(
        bam_fpn=to('data/r1/trumicount/10xn9k_10c.bam'),
        sv_bam_fpn=to('data/r1/trumicount/10xn9k_10c_tagged.bam'),
        verbose=True,
    )
    print(p.add_tags())