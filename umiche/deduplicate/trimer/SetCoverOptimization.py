__version__ = "v1.0"
__copyright__ = "Copyright 2024"
__license__ = "MIT"
__developer__ = "Jianfeng Sun"
__maintainer__ = "Jianfeng Sun"
__email__="jianfeng.sunmt@gmail.com"
__lab__ = "Cribbslab"

import pysam
import pandas as pd
from umiche.path import to
from umiche.deduplicate.trimer.Collapse import collapse


class setCoverOptimization:

    def __init__(self, ):
        self.collapse = collapse()

    def umi_greedy(self, input_umis):
        merged_umis = {}
        merged_umis_idx = {}
        n_merged = len(input_umis)
        left_umis = input_umis
        n_steps = 0
        while n_merged > 1:
            umi_ids = dict()
            for ii, uu_list in enumerate(left_umis):
                for uu in uu_list:
                    if uu in umi_ids:
                        umi_ids[uu].append(ii)
                    else:
                        umi_ids[uu] = [ii]
            umi_counts = {kk: len(vv) for kk, vv in umi_ids.items()}
            umi_df = pd.DataFrame(
                {'umi_id': umi_counts.keys(), 'umi_count': umi_counts.values()}
            ).sort_values('umi_count', ascending=False).reset_index(drop=True)

            n_merged = umi_df.umi_count[0]
            if n_merged > 1:
                merged_umis_idx[n_steps] = umi_ids[umi_df.umi_id[0]]
                merged_umis[n_steps] = umi_df.umi_id[0]
                # fix umis in path
                ll_umis = []
                for rr_idx, ll_uu in enumerate(left_umis):
                    if rr_idx not in merged_umis_idx[n_steps]:
                        ll_umis.append(ll_uu)
                left_umis = ll_umis
                if len(left_umis) == 0:
                    break
                n_steps += 1
            else:
                break
        shortest_path = len(input_umis) - sum([len(ii) - 1 for ii in merged_umis_idx.values()])
        return (shortest_path)

    def count(self, inbam, tag, sep="_"):
        tab = dict()
        # count = 0
        n_tag = 0
        with pysam.AlignmentFile(inbam) as bf:
            for i, r in enumerate(bf):
                if r.has_tag(tag) is True:
                    key = (r.get_tag(tag), "fake_cb")
                    tab.setdefault(key, []).append(r.qname.split(sep)[1])
                    n_tag += 1
                else:
                    pass
            print("The total number of input reads is ", i + 1)
            print("The total number of input reads with XT tag is ", n_tag)
            n = 0
            for kk, vv in tab.items():
                if len(vv) == 1:
                    count = 1
                    n += count
                else:
                    # corrected_cmis = set(tuple(self.collapse.ShuangLiCollapseCMI(uu)) for uu in vv)
                    corrected_cmis = set(tuple(self.collapse.splitByMV(uu)) for uu in vv)
                    if len(corrected_cmis) == 1:
                        count = 1
                        n += count
                    else:
                        # count = self.umi_greedy([self.collapse.ShuangLiCollapseUMI(uu) for uu in ['TTTCCCTTTAAAGGGTTTGGGCCCCCC', 'TTGCCGTTTAAAGGGTTTGGGCCCCCC', 'TAGCAGTTGATAGGGTTTGGGCCCCCC']])
                        # count = self.umi_greedy([self.collapse.ShuangLiCollapseUMI(uu) for uu in vv])
                        count = self.umi_greedy([self.collapse.splitToAll(uu) for uu in vv])
                        n += count
        return n


if __name__ == "__main__":
    p = setCoverOptimization()
    print(p.count(
        inbam=to('data/simu/trimer/pcr8/seq_errs/permute_0/trimmed/seq_err_18.bam'),
        tag='PO',
    ))
    # print(p.umi_greedy([collapse_umi(uu) for uu in vv]))
