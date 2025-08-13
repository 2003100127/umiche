__version__ = "v1.0"
__copyright__ = "Copyright 2025"
__license__ = "GPL-3.0"
__developer__ = "Jianfeng Sun"
__maintainer__ = "Jianfeng Sun"
__email__="jianfeng.sunmt@gmail.com"

from typing import Dict

import re

from umiche.bam.Gadgetry import Gadgetry
from umiche.util.Console import Console


class Tagging:

    def __init__(
            self,
            bam_fpn: str,
            sv_bam_fpn: str = 'tagged.bam',

            verbose: bool = False,
    ):
        self.bam_fpn = bam_fpn
        self.sv_bam_fpn = sv_bam_fpn

        import pysam
        self.pysam = pysam

        self.console = Console()
        self.console.verbose = verbose

    @Gadgetry().index_sort(do_sort=True, do_index=True)
    @Console.vignette()
    def add_from_header(
            self,
            delimiter=r'[_:]',
            tag_info: Dict={
                'CB' : {'idx': -2, 'val_type': 'Z'},
                'MB': {'idx': -1, 'val_type': 'Z'},
            },
            strip_after: bool = False,
    ) -> str:
        with self.pysam.AlignmentFile(self.bam_fpn, "rb") as bam_in, \
                self.pysam.AlignmentFile(self.sv_bam_fpn, "wb", template=bam_in) as bam_out:

            tag_running_title = ["&".join(i) for i in tag_info.keys()]
            tag_indices = [item['idx'] for item in tag_info.values()]
            for read in self.console._tqdm(
                    bam_in,
                    # total=bam_in.count(until_eof=True),
                    desc=f"[{tag_running_title} tagging]",
                    unit="reads",
                    position=0,
                    leave=True,
                    dynamic_ncols=False,
            ):
                parts = re.split(delimiter, read.query_name)
                # print(parts)
                if len(parts) <= max(tag_indices):
                    # 当 read name is not qualified
                    bam_out.write(read)
                    continue

                for tag_name, tag_val in tag_info.items():
                    target = parts[tag_val['idx']]
                    # print(read.get_tag("XF"))
                    read.set_tag(tag_name, target, value_type=tag_val['val_type'])
                if strip_after:
                    read.query_name = parts[0]

                bam_out.write(read)

        return self.sv_bam_fpn


if __name__ == "__main__":
    # from umiche.path import to
    # samtools view /mnt/d/Document/Programming/Python/umiche/umiche/data/r1/trumicount/10xn9k_10c/10xn9k_10c_tagged.bam | head -n 5
    # python Trumicount.py --input-bam /mnt/d/Document/Programming/Python/umiche/umiche/data/r1/trumicount/10xn9k_10c/10xn9k_10c_tagged.sorted.bam --cell-tag CB --umi-tag MB --gene-tag XF --threshold 2 --output-counts sample.counts.tsv --output-assigned sample.assigned.tsv
    p = Tagging(
        bam_fpn='/mnt/d/Document/Programming/Python/umiche/umiche/data/r1/trumicount/10xn9k_10c/10xn9k_10c.bam',
        sv_bam_fpn='/mnt/d/Document/Programming/Python/umiche/umiche/data/r1/trumicount/10xn9k_10c/10xn9k_10c_tagged.bam',

        # bam_fpn='/mnt/d/Document/Programming/Python/umiche/umiche/data/r1/umitools/umitools.test.RNA-seq.sorted.bam',
        # sv_bam_fpn='/mnt/d/Document/Programming/Python/umiche/umiche/data/r1/umitools/umitools.test.RNA-seq.sorted.tagged.bam',

        verbose=True,
    )
    print(p.add_from_header(
        tag_info={
            'CB': {'idx': -2, 'val_type': 'Z'},
            'MB': {'idx': -1, 'val_type': 'Z'},
        }
    ))