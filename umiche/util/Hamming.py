__version__ = "v1.0"
__copyright__ = "Copyright 2024"
__license__ = "MIT"
__developer__ = "Jianfeng Sun"
__maintainer__ = "Jianfeng Sun"
__email__="jianfeng.sunmt@gmail.com"
__lab__ = "Cribbslab"


class Hamming:

    def general(self, s1, s2):
        """

        Parameters
        ----------
        s1
            sequence 1
        s2
            sequence 2

        Returns
        -------
            int

        """
        return sum(i != j for i, j in zip(s1, s2))