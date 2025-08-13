__version__ = "v1.0"
__copyright__ = "Copyright 2025"
__license__ = "GPL-3.0"
__developer__ = "Jianfeng Sun"
__maintainer__ = "Jianfeng Sun"
__email__="jianfeng.sunmt@gmail.com"

from typing import Dict, Union

import pandas as pd
from collections import OrderedDict
from umiche.bam.Gadgetry import Gadgetry
from umiche.util.Console import Console


class Generator:
    """
    Generate a minimal usable BAM for testing UMI + coordinate-based PCR
    duplicate removal, based on (label -> UMI) and (label -> count) mappings.

    Typical rule: Reads with the same coordinate + same UMI are
    considered duplicates.
    """

    def __init__(
            self,
            ref_name: str = "chr1",
            ref_len: int = 10_000_000,
            read_len: int = 50,
            mapq: int = 60,
            same_pos_for_all: bool = True,
            base_pos: int = 1000,
            base_pos_step: int = 500,
            umi_tag: str ="MB",
            cell_barcode: Union[str, None] = None,
            program_id: str = "umi_synth",
            program_version: str = "0.1",
            sv_bam_fpn: str = "umi_to_bam.bam",
            verbose: bool = True,
    ):
        """
        Parameters
        ----------
        ref_name/ref_len : Pseudo reference sequence name0 and length.
        read_len         : Length of pseudo reads.
        mapq             : Mapping quality (constant).
        same_pos_for_all : True=All UMIs use the same coordinates; False=Each UMI/label has staggered coordinates.
        base_pos         : Baseline starting coordinate (0-based).
        base_pos_step    : When same_pos_for_all=False, the coordinate step between different UMIs.
        umi_tag          : UMI tag, "MB".
        cell_barcode     : If not None, write the CB tag.
        """
        import pysam
        self.pysam = pysam
        self.ref_name = ref_name
        self.ref_len = int(ref_len)
        self.read_len = int(read_len)
        self.mapq = int(mapq)
        self.same_pos_for_all = bool(same_pos_for_all)
        self.base_pos = int(base_pos)
        self.base_pos_step = int(base_pos_step)
        self.umi_tag = umi_tag
        self.cell_barcode = cell_barcode
        self.program_id = program_id
        self.program_version = program_version

        self.verbose = verbose
        self.console = Console(verbose=self.verbose)
        self.sv_bam_fpn = sv_bam_fpn

    @Gadgetry().index_sort(do_sort=True, do_index=True)
    @Console.vignette_global()
    def mock_from_counts(
        self,
        umi_seq_dict: Dict[str, str],
        umi_cnt_dict: Union[Dict[str, int], pd.Series],
    ) -> str:
        """

        Parameters
        ----------


        Returns
        -------
            BAM path
        """
        ordered_counts = self._normalize_counts(umi_cnt_dict)

        missing = [k for k in ordered_counts if k not in umi_seq_dict]
        if missing:
            raise ValueError(f"label are missing in umi_seq_dict：{missing}")

        header = self._build_header()
        seq, qual = self._make_const_seq_and_qual(self.read_len)

        with self.pysam.AlignmentFile(self.sv_bam_fpn, "wb", header=header) as outf:
            for idx, (label, n_reads) in enumerate(ordered_counts.items()):
                umi = umi_seq_dict[label]
                ref_start = (
                    self.base_pos
                    if self.same_pos_for_all
                    else self.base_pos + idx * self.base_pos_step
                )

                for k in range(int(n_reads)):
                    a = self.pysam.AlignedSegment()
                    a.query_name = f"read|label:{label}|umi:{umi}|rep:{k+1}"
                    a.query_sequence = seq
                    a.flag = 0
                    a.reference_id = 0
                    a.reference_start = ref_start
                    a.mapping_quality = self.mapq
                    a.cigarstring = f"{self.read_len}M"
                    a.query_qualities = qual

                    a.set_tag(self.umi_tag, umi, value_type="Z")

                    if self.cell_barcode is not None:
                        a.set_tag("CB", self.cell_barcode, value_type="Z")

                    outf.write(a)

        return self.sv_bam_fpn

    def _build_header(self, ):
        return {
            "HD": {"VN": "1.6"},
            "SQ": [{
                "SN": self.ref_name,
                "LN": self.ref_len,
            }],
            "PG": [{
                "ID": self.program_id,
                "PN": self.program_id,
                "VN": self.program_version,
            }],
        }

    def _make_const_seq_and_qual(
            self,
            n: int,
    ):
        seq = "A" * n
        qual = self.pysam.qualitystring_to_array("I" * n)  # 'I'≈Q40
        return seq, qual

    @staticmethod
    def _normalize_counts(
            counts: Union[Dict[str, int], pd.Series],
    ) -> OrderedDict:
        if isinstance(counts, pd.Series):
            return OrderedDict((k, int(v)) for k, v in counts.items())
        elif isinstance(counts, dict):
            return OrderedDict((k, int(v)) for k, v in counts.items())
        else:
            raise TypeError("node_val_sorted must be a dict or a pandas.Series")


if __name__ == "__main__":
    # samtools view /mnt/d/Document/Programming/Python/umiche/umiche/data/r1/umi_to_bam.bam | head -n 5
    umi_seq_dict = {
        "A": "ACGT",
        "B": "TCGT",
        "C": "CCGT",
        "D": "ACAT",
        "E": "AAAT",
        "F": "ACAG",
    }

    umi_cnt_dict = pd.Series({
        'A': 456,
        'E': 90,
        'D': 72,
        'B': 2,
        'C': 2,
        'F': 1,
    })

    # umi_seq_dict = {
    #     "A": "AGATCTCGCA",
    #     "B": "AGATCCCGCA",
    #     "C": "AGATCACGCA",
    #     "D": "AGATCTGGCA",
    #     "E": "AGATCTGGGA",
    #     "F": "AGATCTGGCT",
    #     "G": "AGATCTGGGT",
    # }

    # umi_cnt_dict = pd.Series({
    #     'A': 120,
    #     'D': 90,
    #     'E': 50,
    #     'G': 5,
    #     'B': 2,
    #     'C': 2,
    #     'F': 1,
    # })

    p = Generator(
        sv_bam_fpn="/mnt/d/Document/Programming/Python/umiche/umiche/data/r1/mock/umitools_graph.bam",
    )
    bam_fpn = p.mock_from_counts(
        umi_seq_dict=umi_seq_dict,
        umi_cnt_dict=umi_cnt_dict,
    )
    print("saved in", bam_fpn)

    # gen2 = Generator(same_pos_for_all=False, base_pos=1000, base_pos_step=500)
    # bam2 = gen2.mock_from_counts(umi_seq_dict, umi_cnt_dict, "umi_sim_splitpos.bam")
    # print("Wrote:", bam2)
