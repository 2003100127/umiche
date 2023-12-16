__version__ = "v1.0"
__copyright__ = "Copyright 2024"
__license__ = "MIT"
__developer__ = "Jianfeng Sun"
__maintainer__ = "Jianfeng Sun"
__email__="jianfeng.sunmt@gmail.com"
__lab__ = "Cribbslab"

import os


class folder:

    def osmkdir(self, DIRECTORY):
        """

        Parameters
        ----------
        DIRECTORY

        Returns
        -------

        """
        if not os.path.exists(DIRECTORY):
            os.makedirs(DIRECTORY)
        return 0
