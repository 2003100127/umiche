__version__ = "v1.0"
__copyright__ = "Copyright 2024"
__license__ = "MIT"
__developer__ = "Jianfeng Sun"
__maintainer__ = "Jianfeng Sun"
__email__="jianfeng.sunmt@gmail.com"
__lab__ = "Cribbslab"

import pysam
import pandas as pd
import collections
from umiche.deduplicate.method.trimer.Collapse import Collapse
from umiche.util.Console import Console


class SetCover:

    def __init__(
            self,
            verbose=False,
    ):
        self.collapse = Collapse()

        self.console = Console()
        self.console.verbose = verbose

    def umi_greedy_li(
            self,
            umi_dict,
    ):
        """

        Parameters
        ----------
        umi_arr

        Returns
        -------

        """
        # print(umi_arr)
        merged_umis = {}
        merged_umis_idx = {}
        trimer_umi_list = [*umi_dict.keys()]
        # @@ trimer_umi_list[0]
        # ['GGGTTTGTGACCCCCTGTAAATTTCCCCGGAAAGTG',
        # 'GGGAAATTTTTTGTTCTCAAAGGGCAAGGGAAATTT',
        # ...,
        # 'AAAGGGAAACCCAAATTTGGGTTTTCGTTTCCTTTT',]
        mono_umi_set_list = [*umi_dict.values()]
        # @@ mono_umi_set_list
        # [{'GTGCCTATCGAG', 'GTGACGATCGAT', ..., 'GTGCCGATCGAT'},
        # {'GATTGCAGAGAT', 'GATTTTAGCGAT', ..., 'GATTGCAGCGAT'},
        # ...,
        # {'AGACATGTGTTT', 'AGACATGTCTTT', ..., 'AGACATGTTTCT'}]
        n_merged = len(mono_umi_set_list)
        mono_umi_set_list_remaining = mono_umi_set_list
        n_steps = 0
        # print(left_umis)
        while n_merged > 1:
            mono_umi_to_trimer_id_dict = dict()
            for ii, mono_umi_set in enumerate(mono_umi_set_list_remaining):
                # print(len(left_umis))
                for mono_umi in mono_umi_set:
                    if mono_umi in mono_umi_to_trimer_id_dict:
                        mono_umi_to_trimer_id_dict[mono_umi].append(ii)
                    else:
                        mono_umi_to_trimer_id_dict[mono_umi] = [ii]
            # @@ mono_umi_to_trimer_id_dict
            # {'AGATGGACGGCA': [4523, 4989, 5607],
            # 'AGTAGGATGGCA': [4523, 5725],
            # ...,
            # 'GAAGCGAATTTC': [6521]}
            umi_counts = {kk: len(vv) for kk, vv in mono_umi_to_trimer_id_dict.items()}

            df_umi = pd.DataFrame({
                'mono_umi': umi_counts.keys(),
                'val_cnt': umi_counts.values(),
            }).sort_values(
                by='val_cnt',
                ascending=False,
            ).reset_index(drop=True)
            # print(df_umi)

            n_merged = df_umi.val_cnt[0]

            if n_merged > 1:
                merged_umis_idx[n_steps] = mono_umi_to_trimer_id_dict[df_umi.mono_umi[0]]
                merged_umis[n_steps] = df_umi.mono_umi[0]
                print(merged_umis_idx)
                # fix umis in path
                mono_umi_set_list_remaining_tmp = []
                for rr_idx, mono_left_umi_set in enumerate(mono_umi_set_list_remaining):
                    # print(rr_idx)
                    # print(mono_left_umi_set)
                    if rr_idx not in merged_umis_idx[n_steps]:
                        # print(rr_idx, merged_umis_idx[n_steps])
                        # break
                        mono_umi_set_list_remaining_tmp.append(mono_left_umi_set)
                # print(mono_umi_set_list_remaining_tmp)
                # print(len(mono_umi_set_list_remaining_tmp))
                # break
                mono_umi_set_list_remaining = mono_umi_set_list_remaining_tmp
                if len(mono_umi_set_list_remaining) == 0:
                    break
                n_steps += 1
            else:
                break
        print(df_umi.loc[df_umi['val_cnt'] == 1])
        print(df_umi.loc[df_umi['val_cnt'] == 1].shape)
        # print(mono_umi_set_list_remaining)
        # solutiom_unique_umis = set(merged_umis.values())
        # solution_umis = [ii for ii, cc in enumerate(mono_umi_set_list) if len(cc.intersection(solutiom_unique_umis)) > 0]
        shortest_path = len(mono_umi_set_list) - sum([len(ii) - 1 for ii in merged_umis_idx.values()])
        print('shortest_path: ', shortest_path)
        return (shortest_path)

    def count_li(
            self,
            inbam,
            tag,
            sep="_",
    ):
        tab = dict()
        # count = 0
        n_tag = 0
        with pysam.AlignmentFile(inbam) as bf:
            for i, r in enumerate(bf):
                if r.has_tag(tag) is True:

                    key = ('tag', r.get_tag(tag))

                    tab.setdefault(key, []).append(r.qname.split(sep)[1])
                    print(tab)
                    n_tag += 1
                else:
                    pass

            self.console.print("======>The total number of input reads is ".format(i + 1))
            self.console.print("======>The total number of input reads with XT tag is ".format(n_tag))

            trimer_with_mono = collections.defaultdict(list)

            n = 0
            for tag, trimer_list in tab.items():
                # print(str(tag))
                # print(len(trimer_list))
                if len(trimer_list) == 1:
                    count = 1
                    n += count

                    monocmi = self.collapse.cmi_li(trimer_list).pop()
                    # print(monocmi)
                    keys = str(tag) + '_' + str(trimer_list)
                    trimer_with_mono[keys].append(monocmi)

                else:
                    # corrected_cmis = set(tuple(self.collapse.ShuangLiCollapseCMI(uu)) for uu in trimer_list)
                    corrected_cmis = set(tuple(self.collapse.split_by_mv(uu)) for uu in trimer_list)
                    # print(corrected_cmis)

                    if len(corrected_cmis) == 1:
                        count = 1
                        n += count

                        for trimer in trimer_list:
                            keys = '_'.join(tag) + '_' + str(trimer)
                            # print(keys)
                            trimer_with_mono[keys].append(corrected_cmis)

                    else:
                        trimer_to_combos = {tri: self.collapse.split_to_all(tri) for tri in trimer_list}
                        # print(trimer_to_combos)

                        # count = self.umi_greedy([self.collapse.ShuangLiCollapseUMI(uu) for uu in ['TTTCCCTTTAAAGGGTTTGGGCCCCCC', 'TTGCCGTTTAAAGGGTTTGGGCCCCCC', 'TAGCAGTTGATAGGGTTTGGGCCCCCC']])
                        # count = self.umi_greedy([self.collapse.ShuangLiCollapseUMI(uu) for uu in trimer_list])
                        count = self.umi_greedy(trimer_to_combos)
                        # count = self.umi_greedy([self.collapse.split_to_all(uu) for uu in trimer_list])

                        # trimer_with_mono = {str('_'.join(tag)) + "_" + str(tri): mono for tri, mono in new_trimer_to_combos.items()}

                        n += count
        return n

    def greedy(
            self,
            multimer_list,
    ):
        """

        Parameters
        ----------
        umi_arr

        Returns
        -------

        """
        umi_dict = {multimer_umi: self.collapse.split_to_all(multimer_umi) for multimer_umi in multimer_list}

        merged_umis = {}
        merged_mono_umi_dict = {}
        merged_umis_idx = {}
        trimer_umi_to_id_map = {trimer_umi: k for k, trimer_umi in enumerate(umi_dict.keys())}
        trimer_id_to_umi_map = {k: trimer_umi for k, trimer_umi in enumerate(umi_dict.keys())}
        trimer_umi_list = [*umi_dict.keys()]
        # @@ trimer_umi_list[0]
        # ['GGGTTTGTGACCCCCTGTAAATTTCCCCGGAAAGTG',
        # 'GGGAAATTTTTTGTTCTCAAAGGGCAAGGGAAATTT',
        # ...,
        # 'AAAGGGAAACCCAAATTTGGGTTTTCGTTTCCTTTT',]
        mono_umi_set_list = [*umi_dict.values()]
        # @@ mono_umi_set_list
        # [{'GTGCCTATCGAG', 'GTGACGATCGAT', ..., 'GTGCCGATCGAT'},
        # {'GATTGCAGAGAT', 'GATTTTAGCGAT', ..., 'GATTGCAGCGAT'},
        # ...,
        # {'AGACATGTGTTT', 'AGACATGTCTTT', ..., 'AGACATGTTTCT'}]
        n_merged = len(mono_umi_set_list)
        mono_umi_set_list_remaining = umi_dict
        n_steps = 0
        # print(left_umis)
        while n_merged > 1:
            mono_umi_to_trimer_id_dict = dict()
            for multimer_umi, mono_umi_set in mono_umi_set_list_remaining.items():
                # print(len(left_umis))
                for mono_umi in mono_umi_set:
                    if mono_umi in mono_umi_to_trimer_id_dict:
                        mono_umi_to_trimer_id_dict[mono_umi].append(trimer_umi_to_id_map[multimer_umi])
                    else:
                        mono_umi_to_trimer_id_dict[mono_umi] = [trimer_umi_to_id_map[multimer_umi]]
            # print(mono_umi_to_trimer_id_dict)
            # @@ mono_umi_to_trimer_id_dict
            # {'AGATGGACGGCA': [4523, 4989, 5607],
            # 'AGTAGGATGGCA': [4523, 5725],
            # ...,
            # 'GAAGCGAATTTC': [6521]}
            umi_counts = {kk: len(vv) for kk, vv in mono_umi_to_trimer_id_dict.items()}

            df_umi = pd.DataFrame({
                'mono_umi': umi_counts.keys(),
                'val_cnt': umi_counts.values(),
            }).sort_values(
                by='val_cnt',
                ascending=False,
            ).reset_index(drop=True)
            # print(df_umi)

            n_merged = df_umi.val_cnt[0]

            if n_merged > 1:
                merged_umis_idx[n_steps] = mono_umi_to_trimer_id_dict[df_umi.mono_umi[0]]
                merged_umis[n_steps] = df_umi.mono_umi[0]
                merged_mono_umi_dict[df_umi.mono_umi[0]] = trimer_id_to_umi_map[mono_umi_to_trimer_id_dict[df_umi.mono_umi[0]][0]]
                print(merged_mono_umi_dict)
                # fix umis in path
                mono_umi_set_list_remaining_tmp = {}
                # print(len(mono_umi_set_list_remaining))
                for multimer_umi, mono_left_umi_set in mono_umi_set_list_remaining.items():
                    # print(multimer_umi)
                    # print(mono_left_umi_set)
                    if trimer_umi_to_id_map[multimer_umi] not in merged_umis_idx[n_steps]:
                        # print(multimer_umi, merged_umis_idx[n_steps])
                        # break
                        mono_umi_set_list_remaining_tmp[multimer_umi] = mono_left_umi_set
                # print(mono_umi_set_list_remaining_tmp)
                # print(len(mono_umi_set_list_remaining_tmp))
                # break
                mono_umi_set_list_remaining = mono_umi_set_list_remaining_tmp
                print(len(mono_umi_set_list_remaining))

                if len(mono_umi_set_list_remaining) == 0:
                    break
                n_steps += 1
            else:
                break
        multimer_umi_solved_by_sc = [*merged_mono_umi_dict.values()]
        multimer_umi_not_solved = [*mono_umi_set_list_remaining.keys()]
        print('=========>number of shortlisted multimer UMIs solved by set cover'.format(len(multimer_umi_solved_by_sc)))
        print('=========>number of shortlisted multimer UMIs cannot be solved by set cover'.format(len(multimer_umi_not_solved)))
        # solutiom_unique_umis = set(merged_umis.values())
        # solution_umis = [ii for ii, cc in enumerate(mono_umi_set_list) if len(cc.intersection(solutiom_unique_umis)) > 0]
        shortlisted_multimer_umi_list = multimer_umi_solved_by_sc + multimer_umi_not_solved
        dedup_cnt = len(mono_umi_set_list) - sum([len(ii) - 1 for ii in merged_umis_idx.values()])
        print('dedup_cnt: ', dedup_cnt)
        return dedup_cnt, shortlisted_multimer_umi_list


if __name__ == "__main__":
    from umiche.path import to

    p = SetCover()

    print(p.greedy(
        multimer_list=multimer_list,
    ))
    # print(p.count_li(
    #     inbam=to('data/simu/umi/trimer/seq_errs/permute_0/trimmed/seq_err_17.bam'),
    #     tag='PO',
    # ))