__version__ = "v1.0"
__copyright__ = "Copyright 2024"
__license__ = "MIT"
__developer__ = "Jianfeng Sun"
__lab__ = "Cribbslab"


class hamming:

    def general(self, s1, s2):
        """

        Parameters
        ----------
        s1
        s2

        Returns
        -------

        """
        return sum(i != j for i, j in zip(s1, s2))