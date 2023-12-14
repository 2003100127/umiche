__version__ = "v1.0"
__copyright__ = "Copyright 2024"
__license__ = "MIT"
__developer__ = "Jianfeng Sun"
__lab__ = "cribbslab"

import pyfastx


class read():

    def __init__(self):
        pass

    def fromgz(self, fastq_fpn):
        """

        Note
        ----
        read and parse a Fastq file.

        Parameters
        ----------
        fastq_fpn
            the path of a fastq file
        Returns
        -------
            tuple consisting of names, seqs, placeholders, qualities
        """
        names = []
        seqs = []
        placeholders = []
        qualities = []
        fq = pyfastx.Fastx(fastq_fpn)
        print(seqs)
        for name, seq, qual in fq:
            seqs.append(''.join(seq))
            names.append(''.join(name))
        return names, seqs, placeholders, qualities
