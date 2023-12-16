__version__ = "v1.0"
__copyright__ = "Copyright 2024"
__license__ = "MIT"
__developer__ = "Jianfeng Sun"
__maintainer__ = "Jianfeng Sun"
__email__="jianfeng.sunmt@gmail.com"
__lab__ = "Cribbslab"

import numpy as np


class number:

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