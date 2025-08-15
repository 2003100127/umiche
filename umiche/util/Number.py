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

    def binomial(self, n, p, use_seed=True, seed=1):
        if use_seed:
            state = np.random.RandomState(seed)
            return state.binomial(
                n,
                p,
            )
        else:
            return np.random.binomial(
                n,
                p,
            )

    def nbinomial(self, n, p, use_seed=True, seed=1):
        if use_seed:
            state = np.random.RandomState(seed)
            return state.negative_binomial(
                n,
                p,
            )
        else:
            return np.random.negative_binomial(
                n,
                p,
            )

    def uniform(self, low, high, num, use_seed=True, seed=1):
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

    def choice(self, high, num, use_seed=True, seed=1, replace=False):
        if use_seed:
            state = np.random.RandomState(seed)
            return state.choice(high, int(num), replace=replace)
        else:
            return np.random.choice(high, int(num), replace=replace)