__version__ = "v1.0"
__copyright__ = "Copyright 2025"
__license__ = "GPL-3.0"
__developer__ = "Jianfeng Sun"
__maintainer__ = "Jianfeng Sun"
__email__="jianfeng.sunmt@gmail.com"

from collections import Counter


class UMI:

    def __init__(
            self,
            umis,
    ):
        self.umis_cnt_dict = Counter(umis)
        # print(umis_cnt_dict)
        # @@ umis_cnt_dict
        # Counter({'AGATCTCGCA': 120, 'AGATCTGGCA': 90, 'AGATCTGGGA': 50, 'AGATCTGGGT': 5, 'AGATCCCGCA': 2, 'AGATCACGCA': 2, 'AGATCTGGCT': 1})
        self.umis_uniq_ordered = sorted(self.umis_cnt_dict, key=lambda u: (-self.umis_cnt_dict[u], u))
        # print(umis_uniq_ordered)
        # @@ umis_uniq_ordered
        # ['AGATCTCGCA', 'AGATCTGGCA', 'AGATCTGGGA', 'AGATCTGGGT', 'AGATCACGCA', 'AGATCCCGCA', 'AGATCTGGCT']
