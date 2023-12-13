__version__ = "v1.0"
__copyright__ = "Copyright 2023"
__license__ = "MIT"
__developer__ = "Jianfeng Sun"
__lab__ = "cribbslab"

import os
import gzip
from umiche.util.Folder import folder as crtfolder


class write:

    def __init__(self):
        pass

    def togz(self, list_2d, sv_fpn, symbol='_'):
        """

        Parameters
        ----------
        list_2d
        sv_fpn
        symbol

        Returns
        -------

        """
        crtfolder().osmkdir(DIRECTORY=os.path.dirname(sv_fpn))
        f = gzip.open(sv_fpn, 'wt')
        for i, read in enumerate(list_2d):
            seq = str(read[0])
            # print('No.{} saving in FASTQ format.'.format(i + 1))
            f.write('@' + symbol.join(read[1:]) + '\n')
            f.write(seq + '\n')
            f.write('+' + '\n')
            f.write('s' * len(seq) + '\n')
        f.close()