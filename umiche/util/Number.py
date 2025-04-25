__version__ = "v1.0"
__copyright__ = "Copyright 2025"
__license__ = "GPL-3.0"
__developer__ = "Jianfeng Sun"
__maintainer__ = "Jianfeng Sun"
__email__="jianfeng.sunmt@gmail.com"


import numpy as np


class Number:

    def __init__(self, *args, **kwargs):
        self.args = args
        self.kwargs = kwargs

    def uniform(self, low, high, num, use_seed=True, seed=1):
        """

        Parameters
        ----------
        low
        high
        num
        use_seed
        seed

        Returns
        -------

        """
        if use_seed:
            state = np.random.RandomState(seed)
            return state.randint(
                low=low,
                high=high,
                size=num
            )
        else:
            return np.random.randint(
                low=low,
                high=high,
                size=num
            )