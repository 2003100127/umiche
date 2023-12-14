__version__ = "v1.0"
__copyright__ = "Copyright 2024"
__license__ = "MIT"
__developer__ = "Jianfeng Sun"
__lab__ = "cribbslab"

import os


class folder:

    def __init__(self, ):
        pass

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
