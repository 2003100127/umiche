__version__ = "v1.0"
__copyright__ = "Copyright 2025"
__license__ = "GPL-3.0"
__developer__ = "Jianfeng Sun"
__maintainer__ = "Jianfeng Sun"
__email__="jianfeng.sunmt@gmail.com"


import pandas as pd
from umiche.util.Console import Console


class Writer:

    def __init__(
            self,
             df,
            verbose=False,
    ):
        import pysam

        self.pysam =pysam
        self.df = df
        self.console = Console()
        self.console.verbose = verbose

    def tobam(
            self,
            tobam_fpn,
            tmpl_bam_fpn,
            whitelist=[],
    ):
        tmpl_bam = self.pysam.AlignmentFile(tmpl_bam_fpn, "rb")
        write_to_bam = self.pysam.AlignmentFile(tobam_fpn, "wb", template=tmpl_bam)
        fs = self.df.loc[self.df['id'].isin(whitelist)]['read']
        for i in fs:
            # print(i)
            write_to_bam.write(i)
        write_to_bam.close()
        return write_to_bam

    def tobam_trust4(
            self,
            tobam_fpn,
            # tmpl_bam_fpn,
            # whitelist=[],
    ):
        # 'wb' indicates that it's a binary file opened for writing
        # print(self.pysam.AlignmentHeader())
        # with self.pysam.AlignmentFile(tobam_fpn, 'wb', header=self.pysam.AlignmentHeader()) as output_bam:
        #     # Create a new alignment object for your query sequence
        #     # Replace placeholders with actual values
        #     # fs = self.df.loc[self.df['id'].isin(whitelist)]['read']
        #     # for i in fs:
        #         new_alignment = self.pysam.AlignedSegment()
        #         new_alignment.query_name = "YourQueryName"
        #         new_alignment.query_sequence = "ATCG"  # Replace with the actual sequence
        #         new_alignment.flag = 0x0  # Replace with appropriate flag
        #         new_alignment.reference_id = -1  # No reference ID specified
        #         new_alignment.reference_start = 0  # Replace with the start position
        #         new_alignment.mapping_quality = 60  # Replace with the mapping quality
        #         new_alignment.cigar = [(0, 4)]  # Replace with the CIGAR string
        #         new_alignment.next_reference_id = -1  # No next reference ID specified
        #         new_alignment.next_reference_start = 0  # Replace with the next reference start position
        #         new_alignment.template_length = 0  # Replace with the template length
        #         new_alignment.query_qualities = [30, 30, 30, 30]  # Replace with the quality scores
        #
        #         # Add the new alignment to the BAM file
        #         output_bam.write(new_alignment)
        header = {'HD': {'VN': '1.0'},
                  'SQ': [{'LN': 1575, 'SN': 'chr1'},
                         {'LN': 1584, 'SN': 'chr2'}]}

        with self.pysam.AlignmentFile(tobam_fpn, "wb", header=header) as outf:
            a = self.pysam.AlignedSegment()
            a.query_name = "read_28833_29006_6945"
            a.query_sequence = "AGCTTAGCTAGCTACCTATATCTTGGTCTTGGCCG"
            a.flag = 99
            a.reference_id = 0
            a.reference_start = 32
            a.mapping_quality = 20
            a.cigar = ((0, 10), (2, 1), (0, 25))
            a.next_reference_id = 0
            a.next_reference_start = 199
            a.template_length = 167
            a.query_qualities = self.pysam.qualitystring_to_array("<<<<<<<<<<<<<<<<<<<<<:<9/,&,22;;<<<")
            a.tags = (("NM", 1),
                      ("RG", "L1"))
            outf.write(a)
        return 1

    def tobam_spikein(
            self,
            tobam_fpn,
            tmpl_bam_fpn,
    ):
        import sys, pysam, csv
        tags = {}
        with open(map_tsv, newline='') as f:
            r = csv.DictReader(f, delimiter='\t')
            for row in r:
                qn = row['qname']
                # 去掉空值
                tags[qn] = {k: v for k, v in row.items() if k not in ('qname',) and v}

        with pysam.AlignmentFile(tmpl_bam_fpn, "rb") as ib, \
                pysam.AlignmentFile(tobam_fpn, "wb", template=ib) as ob:
            for aln in ib:
                t = tags.get(aln.query_name)
                if t:
                    for tag, val in t.items():
                        aln.set_tag(tag, val, value_type='Z')  # 字符串型 TAG
                ob.write(aln)
        pysam.index(tobam_fpn)
        return


if __name__ == "__main__":
    p = Writer(df=pd.DataFrame())

    # p.tobam_trust4(
    #     tobam_fpn='./trimmed2.bam',
    #     # tmpl_bam_fpn=to('data/bone_marrow/merge_corrected_sorted_little.bam'),
    #     # whitelist=df.index,
    # )

    from umiche import io

    df = io.read(
        df_fpn="/mnt/d/Document/Programming/R/umiche/read2tags.tsv",
        df_sep="\t",
        header=0,
        type='tsv',
    )
    print(df)
    print(df.qname.unique().shape)
    # p.tobam()
