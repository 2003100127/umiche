__version__ = "v1.0"
__copyright__ = "Copyright 2025"
__license__ = "GPL-3.0"
__developer__ = "Jianfeng Sun"
__maintainer__ = "Jianfeng Sun"
__email__="jianfeng.sunmt@gmail.com"


import os
from umiche.fastq.Reader import Reader as fastqread
from umiche.util.Folder import Folder as crtfolder
from umiche.util.Console import Console


class Convert:

    def __init__(
            self,
            fastq_fpn,
            bam_fpn,
            verbose=False,
    ):
        import pysam
        # import bamnostic as bs

        self.fastq_fpn = fastq_fpn
        self.bam_fpn = bam_fpn
        self.names, self.seqs, placeholders, qualities = fastqread().fromgz(fastq_fpn=fastq_fpn)
        # print(self.names)
        # print(self.seqs)

        self.console = Console()
        self.console.verbose = verbose

    def tobam(self, ):
        header = {
            'HD': {'VN': '1.0'},
            'SQ': [{'LN': 1575, 'SN': 'chr1'}, {'LN': 1584, 'SN': 'chr2'}]
        }
        crtfolder().osmkdir(DIRECTORY=os.path.dirname(self.bam_fpn))
        with pysam.AlignmentFile(self.bam_fpn, "wb", header=header) as outf:
            for name in self.names:
                a = pysam.AlignedSegment()
                a.query_name = name
                a.query_sequence = 'B'
                a.flag = 0
                a.reference_id = 0
                a.reference_start = 0
                a.mapping_quality = 0
                a.cigar = None
                a.next_reference_id = 0
                a.next_reference_start = 0
                a.template_length = 0
                a.query_qualities = None
                a.tags = (("PO", 1), )
                outf.write(a)
        return

    def tobamSimuBulk(self, ):
        header = {
            'HD': {'VN': '1.0'},
            'SQ': [{'LN': 1575, 'SN': 'chr1'}, {'LN': 1584, 'SN': 'chr2'}]
        }
        crtfolder().osmkdir(DIRECTORY=os.path.dirname(self.bam_fpn))
        with pysam.AlignmentFile(self.bam_fpn, "wb", header=header) as outf:
            for name in self.names:
                a = pysam.AlignedSegment()
                a.query_name = name
                a.query_sequence = 'B'
                a.flag = 0
                a.reference_id = 0
                a.reference_start = 0
                a.mapping_quality = 0
                a.cigar = None
                a.next_reference_id = 0
                a.next_reference_start = 0
                a.template_length = 0
                a.query_qualities = None
                a.tags = (
                    ("XT", int(name.split('-')[4])),
                    ("XS", "Assigned"),
                )
                outf.write(a)
        return

    def tobamsc(self, ):
        header = {
            'HD': {'VN': '1.0'},
            'SQ': [{'LN': 1575, 'SN': 'chr1'}, {'LN': 1584, 'SN': 'chr2'}]
        }
        crtfolder().osmkdir(DIRECTORY=os.path.dirname(self.bam_fpn))
        with pysam.AlignmentFile(self.bam_fpn, "wb", header=header) as outf:
            for name in self.names:
                asterisk_split = name.split('*')
                a = pysam.AlignedSegment()
                a.query_name = name
                a.query_sequence = 'B'
                a.flag = 0
                a.reference_id = 0
                a.reference_start = 0
                a.mapping_quality = 0
                a.cigar = None
                a.next_reference_id = 0
                a.next_reference_start = 0
                a.template_length = 0
                a.query_qualities = None
                a.tags = (
                    ("BC", int(asterisk_split[2])),
                    ("XT", int(asterisk_split[4])),
                    ("XS", "Assigned"),
                )
                outf.write(a)
        return

    def tobambulk(self, ):
        header = {
            'HD': {'VN': '1.0'},
            'SQ': [{'LN': 1575, 'SN': 'chr1'}, {'LN': 1584, 'SN': 'chr2'}]
        }
        crtfolder().osmkdir(DIRECTORY=os.path.dirname(self.bam_fpn))
        with pysam.AlignmentFile(self.bam_fpn, "wb", header=header) as outf:
            for name in self.names:
                asterisk_split = name.split('*')
                a = pysam.AlignedSegment()
                a.query_name = name
                a.query_sequence = 'B'
                a.flag = 0
                a.reference_id = 0
                a.reference_start = 0
                a.mapping_quality = 0
                a.cigar = None
                a.next_reference_id = 0
                a.next_reference_start = 0
                a.template_length = 0
                a.query_qualities = None
                a.tags = (
                    ("SP", int(asterisk_split[2])),
                    ("XT", int(asterisk_split[4])),
                    ("XS", "Assigned"),
                )
                outf.write(a)
        return


if __name__ == "__main__":
    from umiche.path import to
    p = Convert(
        # fastq_fpn=to('data/simu/transloc/trimer/bulk/seq_errs/permute_0/seq_err_0.fastq.gz'),
        # bam_fpn=to('data/simu/transloc/trimer/'),
        fastq_fpn=to('data/simu/trimer/pcr8/seq_errs/permute_0/trimmed/seq_err_18.fastq.gz'),
        bam_fpn=to('data/simu/trimer/pcr8/seq_errs/permute_0/trimmed/seq_err_18.bam'),
    )
    print(p.tobam())
    # print(p.tobamsc())
    # print(p.tobambulk())

    # print(p.tobambs())