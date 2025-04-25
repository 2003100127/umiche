__version__ = "v1.0"
__copyright__ = "Copyright 2025"
__license__ = "GPL-3.0"
__developer__ = "Jianfeng Sun"
__maintainer__ = "Jianfeng Sun"
__email__="jianfeng.sunmt@gmail.com"


from umiche.deduplicate.method.trimer.Collapse import Collapse
from umiche.util.Console import Console


class MajorityVote:

    def __init__(
            self,
            verbose=False,
    ):
        self.collapse = Collapse()

        self.console = Console()
        self.console.verbose = verbose

    def track(
            self,
            multimer_list,
            recur_len,
    ):
        """

        Parameters
        ----------
        multimer_list

        Returns
        -------

        """
        multimer_umi_to_mono_umi_map = {multimer_umi: self.collapse.majority_vote(
            umi=multimer_umi,
            recur_len=recur_len,
        ) for multimer_umi in multimer_list}
        mono_umi_to_multimer_umi_map = {self.collapse.majority_vote(
            umi=multimer_umi,
            recur_len=recur_len,
        ): multimer_umi for multimer_umi in multimer_list}

        uniq_multimer_cnt = len(multimer_umi_to_mono_umi_map)

        shortlisted_multimer_umi_list = [*mono_umi_to_multimer_umi_map.values()]
        dedup_cnt = len(shortlisted_multimer_umi_list)
        self.console.print('=========># of shortlisted multimer UMIs: {}'.format(len(shortlisted_multimer_umi_list)))
        self.console.print('=========>dedup cnt: {}'.format(dedup_cnt))
        return dedup_cnt, uniq_multimer_cnt, shortlisted_multimer_umi_list


if __name__ == "__main__":
    from umiche.path import to

    p = MajorityVote()