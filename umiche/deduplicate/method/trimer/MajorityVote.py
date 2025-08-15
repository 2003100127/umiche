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
            recur_len=3,
    ):
        """

        Parameters
        ----------
        multimer_list

        Returns
        -------

        """
        multimer_umi_to_id_map = {multimer_umi: i for i, multimer_umi in enumerate(multimer_list)}
        # multimer_id_to_umi_map = {i: multimer_umi for i, multimer_umi in enumerate(multimer_list)}
        multimer_umi_to_mono_umi_map = {multimer_umi: self.collapse.majority_vote(
            umi=multimer_umi,
            recur_len=recur_len,
        ) for multimer_umi in multimer_list}
        # print(len(multimer_umi_to_mono_umi_map))

        # print(multimer_umi_to_mono_umi_map)
        # @@ multimer_umi_to_mono_umi_map
        # {'GGGTTTGTGACCCCCTGTAAATTTCCCCGGAAAGTG': 'GTGCCTATCGAG',
        # 'GGGAAATTTTTTGTTCTCAAAGGGCAAGGGAAATTT': 'GATTTCAGAGAT',
        # 'TTTGGGAACAAAGGGTTTAGGTTTCGGAAAAAATTT': 'TGAAGTGTGAAT',
        # ... 'AAAGGGAAACCCAAATTTGGGTTTTCGTTTCCTTTT': 'AGACATGTTTCT'}
        # print(len(multimer_umi_to_mono_umi_map))
        # @@ len(multimer_umi_to_mono_umi_map)
        # 6892
        mono_umi_to_multimer_umi_map = {monomer_umi: multimer_umi for multimer_umi, monomer_umi in multimer_umi_to_mono_umi_map.items()}
        # print(mono_umi_to_multimer_umi_map)
        # @@ mono_umi_to_multimer_umi_map
        # {'GTGCCTATCGAG': 'GGGTTTGGGCCCCCCTTTGAATTTACCCGGAAAGGG',
        # 'GATTTCAGAGAT': 'AGGAAATTCTTTTCTCCCAAAGGGAAAGGGAAATTT',
        # 'TGAAGTGTGAAT': 'TTTGGGAAAAAAGGGTTAGGGTTTGGGAAAAAATTT',
        # ..., 'GCGGACTAGCGT': 'GGCCCCGGGTGGAAACACATTAAAGGTCCCGGGATT'}
        # print(len(mono_umi_to_multimer_umi_map))
        # @@ len(mono_umi_to_multimer_umi_map)
        # 1502

        mono_umi_to_id_map = {monomer_umi: i for i, monomer_umi in enumerate(mono_umi_to_multimer_umi_map.keys())}
        # print(mono_umi_to_id_map)
        # @@ mono_umi_to_id_map
        # {'GTGCCTATCGAG': 0, 'GATTTCAGAGAT': 1, 'TGAAGTGTGAAT': 2,
        # ... 'CCGTTAGGCTCA': 1506, 'AGACATGTTTCT': 1507}
        from collections import defaultdict
        clusters = defaultdict(list)
        for multimer_umi, monomer_umi in multimer_umi_to_mono_umi_map.items():
            mono_id = mono_umi_to_id_map[monomer_umi]
            # clusters[mono_id].append(multimer_umi)
            clusters[mono_id].append(multimer_umi_to_id_map[multimer_umi])
        # print(clusters)

        shortlisted_multimer_umi_list = [*mono_umi_to_multimer_umi_map.values()]
        dedup_cnt = len(shortlisted_multimer_umi_list)
        self.console.print('=========># of shortlisted multimer UMIs: {}'.format(len(shortlisted_multimer_umi_list)))
        self.console.print('=========>dedup cnt: {}'.format(dedup_cnt))
        return {
            'dedup_cnt': dedup_cnt,
            'clusters': dict(clusters), # defaultdict(<class 'list'> to normal dict
            'shortlisted_multimer_umi_list': shortlisted_multimer_umi_list,
        }


if __name__ == "__main__":
    from umiche.path import to

    p = MajorityVote(
        verbose=True,
    )

    from umiche.bam.Reader import Reader as alireader
    alireader = alireader(bam_fpn="/mnt/d/Document/Programming/Python/umiche/umiche/data/simu/umi/trimer/seq_errs/permute_0/trimmed/seq_err_17.bam", verbose=True)
    df_bam = alireader.todf(tags=['PO'])
    print(df_bam)
    print(df_bam.columns)
    print(df_bam.query_name.apply(lambda x: x.split('_')[1]).values.shape)
    print(df_bam.query_name.apply(lambda x: x.split('_')[1]).value_counts(ascending=False))
    print(df_bam.query_name.apply(lambda x: x.split('_')[1]).unique().shape)

    res = p.track(
        multimer_list=df_bam.query_name.apply(lambda x: x.split('_')[1]).values,
        recur_len=3,
    )
    print(res)
