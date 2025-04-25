__version__ = "v1.0"
__copyright__ = "Copyright 2025"
__license__ = "GPL-3.0"
__developer__ = "Jianfeng Sun"
__maintainer__ = "Jianfeng Sun"
__email__="jianfeng.sunmt@gmail.com"


from typing import Set

import itertools
import collections
from umiche.util.Number import Number as rannum
from umiche.util.Console import Console


class Collapse:

    def __init__(
            self,
            verbose=False,
    ):
        self.console = Console()
        self.console.verbose = verbose

    def majority_vote(
            self,
            umi,
            recur_len=3,
    ) -> str:
        """

        Notes
        --------
            elif len(s) == 2:
                sdict = {umi_trimer.count(i): i for i in s}
                t.append(sdict[2])

        Parameters
        ----------
        umi : str
            an umi CTTCCGCATTTTCCCTTTAAAGGGTTTGGGCCCCCC
        recur_len : int
            the length of a building block consisting of recurring nucleotides, e.g. trimer 3

        Returns
        -------
        str
            an umi, in which multipolymer nucleotides collapsed into a single nucleotide
        """
        scales = [i for i in range(len(umi)) if i%recur_len == 0]
        umi_blocks = [umi[v: v+recur_len] for v in scales]
        combos = []
        for umi_block in umi_blocks:
            s = set(umi_block)
            if len(s) == recur_len:
                rand_index = rannum().uniform(low=0, high=recur_len, num=1, use_seed=False, seed=3)[0]
                combos.append(umi_block[rand_index])
            elif len(s) == int(recur_len/2) + 1:
                sdict = {base: umi_block.count(base) for base in s}
                combos.append(max(sdict, key=sdict.get))
            else:
                combos.append(umi_block[0])
        return ''.join(combos)

    def take_by_order(
            self,
            umi,
            pos=0,
            recur_len=3,
    ):
        scales = [i for i in range(len(umi)) if i % recur_len == 0]
        umi_blocks = [umi[v: v + recur_len] for v in scales]
        combos = []
        for umi_block in umi_blocks:
            combos.append(umi_block[pos])
        return ''.join(combos)

    def majority_vote_most_common_deprecated(
            self,
            umi,
    ):
        from collections import Counter
        vernier = [i for i in range(36) if i % 3 == 0]
        umi_trimers = [umi[v: v+3] for v in vernier]
        # @@ textwrap is a bit slower than Python list.
        # import textwrap
        # umi_trimers = textwrap.wrap(umi, 3)
        t = []
        for umi_trimer in umi_trimers:
            s = Counter(umi_trimer).most_common()
            if len(s) == 3:
                rand_index = rannum().uniform(low=0, high=3, num=1, use_seed=False)[0]
                t.append(s[rand_index][0])
            elif len(s) == 2:
                t.append(s[0][0])
            else:
                t.append(umi_trimer[0])
        return ''.join(t)

    def split_to_all(
            self,
            umi,
            recur_len=3,
        ) -> Set:
        """

        Parameters
        ----------
        umi : str
            an umi CTTCCGCATTTTCCCTTTAAAGGGTTTGGGCCCCCC
        recur_len : int
            the length of a building block consisting of recurring nucleotides, e.g. trimer 3

        Returns
        -------
        Set
            a set containing umis, each of which is a nucleotide polymer (trinucleotides)
             collapsed into a single nucleotide.
            e.g. {'TCATCTAGTGCC', 'TGCTCTAGTGCC', 'CGTTCTAGTGCC', ..., 'TCCTCTAGTGCC'}

        """
        # @@ umi
        # CTTCCGCATTTTCCCTTTAAAGGGTTTGGGCCCCCC
        scales = [i for i in range(len(umi)) if i % recur_len == 0]
        # @@ scales
        # [0, 3, 6, 9, 12, 15, 18, 21, 24, 27, 30, 33]
        umi_blocks = [umi[v: v + recur_len] for v in scales]
        # @@ umi_blocks
        # ['CTT', 'CCG', 'CAT', 'TTT', 'CCC', 'TTT', 'AAA', 'GGG', 'TTT', 'GGG', 'CCC', 'CCC']
        combos = [umi_block for umi_block in set(umi_blocks[0])]
        # @@ combos
        # ['C', 'T']
        umi_blocks = umi_blocks[1:]
        # @@ umi_blocks
        # ['CCG', 'CAT', 'TTT', 'CCC', 'TTT', 'AAA', 'GGG', 'TTT', 'GGG', 'CCC', 'CCC']
        for id, umi_block in enumerate(umi_blocks):
            umi_block_bases = set(umi_block)
            c = []
            for i in range(len(combos)):
                for umi_block_base in umi_block_bases:
                    c.append(combos[i] + umi_block_base)
            combos = c
            # @@ combos
            # ['CC', 'CG', 'TC', 'TG']
            # ['CCC', 'CCT', 'CCA', 'CGC', 'CGT', 'CGA', 'TCC', 'TCT', 'TCA', 'TGC', 'TGT', 'TGA']
            # ['CCCT', 'CCTT', 'CCAT', 'CGCT', 'CGTT', 'CGAT', 'TCCT', 'TCTT', 'TCAT', 'TGCT', 'TGTT', 'TGAT']
            # ['CCCTC', 'CCTTC', 'CCATC', 'CGCTC', 'CGTTC', 'CGATC', 'TCCTC', 'TCTTC', 'TCATC', 'TGCTC', 'TGTTC', 'TGATC']
            # ['CCCTCT', 'CCTTCT', 'CCATCT', 'CGCTCT', 'CGTTCT', 'CGATCT', 'TCCTCT', 'TCTTCT', 'TCATCT', 'TGCTCT', 'TGTTCT', 'TGATCT']
            # ['CCCTCTA', 'CCTTCTA', 'CCATCTA', 'CGCTCTA', 'CGTTCTA', 'CGATCTA', 'TCCTCTA', 'TCTTCTA', 'TCATCTA', 'TGCTCTA', 'TGTTCTA', 'TGATCTA']
            # ['CCCTCTAG', 'CCTTCTAG', 'CCATCTAG', 'CGCTCTAG', 'CGTTCTAG', 'CGATCTAG', 'TCCTCTAG', 'TCTTCTAG', 'TCATCTAG', 'TGCTCTAG', 'TGTTCTAG', 'TGATCTAG']
            # ['CCCTCTAGT', 'CCTTCTAGT', 'CCATCTAGT', 'CGCTCTAGT', 'CGTTCTAGT', 'CGATCTAGT', 'TCCTCTAGT', 'TCTTCTAGT', 'TCATCTAGT', 'TGCTCTAGT', 'TGTTCTAGT', 'TGATCTAGT']
            # ['CCCTCTAGTG', 'CCTTCTAGTG', 'CCATCTAGTG', 'CGCTCTAGTG', 'CGTTCTAGTG', 'CGATCTAGTG', 'TCCTCTAGTG', 'TCTTCTAGTG', 'TCATCTAGTG', 'TGCTCTAGTG', 'TGTTCTAGTG', 'TGATCTAGTG']
            # ['CCCTCTAGTGC', 'CCTTCTAGTGC', 'CCATCTAGTGC', 'CGCTCTAGTGC', 'CGTTCTAGTGC', 'CGATCTAGTGC', 'TCCTCTAGTGC', 'TCTTCTAGTGC', 'TCATCTAGTGC', 'TGCTCTAGTGC', 'TGTTCTAGTGC', 'TGATCTAGTGC']
            # ['CCCTCTAGTGCC', 'CCTTCTAGTGCC', 'CCATCTAGTGCC', 'CGCTCTAGTGCC', 'CGTTCTAGTGCC', 'CGATCTAGTGCC', 'TCCTCTAGTGCC', 'TCTTCTAGTGCC', 'TCATCTAGTGCC', 'TGCTCTAGTGCC', 'TGTTCTAGTGCC', 'TGATCTAGTGCC']
        # @@ set(combos)
        # {'TCTTCTAGTGCC', 'CGATCTAGTGCC', 'CGTTCTAGTGCC', 'CCCTCTAGTGCC', 'TGATCTAGTGCC', 'CCATCTAGTGCC', 'CCTTCTAGTGCC', 'TGCTCTAGTGCC', 'TCATCTAGTGCC', 'TCCTCTAGTGCC', 'CGCTCTAGTGCC', 'TGTTCTAGTGCC'}
        return set(combos)

    def split_by_mv(
            self,
            umi,
            recur_len=3,
    ) -> Set:
        """

        Parameters
        ----------
        umi : str
            an umi CTTCCGCATTTTCCCTTTAAAGGGTTTGGGCCCCCC
        recur_len : int
            the length of a building block consisting of recurring nucleotides, e.g. trimer 3

        Returns
        -------
        Set
            a set containing umis, each of which is a nucleotide polymer (trinucleotides)
             collapsed into a single nucleotide.
            e.g. {'TCATCTAGTGCC', 'TCTTCTAGTGCC', 'TCCTCTAGTGCC'}

        """
        scales = [i for i in range(len(umi)) if i % recur_len == 0]
        # @@ scales
        # [0, 3, 6, 9, 12, 15, 18, 21, 24, 27, 30, 33]
        umi_blocks = [umi[v: v + recur_len] for v in scales]
        # @@ umi_blocks
        # ['CTT', 'CCG', 'CAT', 'TTT', 'CCC', 'TTT', 'AAA', 'GGG', 'TTT', 'GGG', 'CCC', 'CCC']
        umi_block_bases = self.vote(umi_blocks[0], recur_len)
        # @@ umi_block_bases
        # {'T'}
        combos = [umi_block for umi_block in umi_block_bases]
        # @@ combos
        # ['T']
        umi_blocks = umi_blocks[1:]
        for id, umi_block in enumerate(umi_blocks):
            umi_block_bases = self.vote(umi_block, recur_len)
            # @@ umi_block_bases
            # {'C'}
            # {'T', 'A', 'C'}
            # {'T'}
            # {'C'}
            # {'T'}
            # {'A'}
            # {'G'}
            # {'T'}
            # {'G'}
            # {'C'}
            # {'C'}
            c = []
            for i in range(len(combos)):
                for umi_block_base in umi_block_bases:
                    c.append(combos[i] + umi_block_base)
            combos = c
            # @@ combos
            # ['TC']
            # ['TCC', 'TCT', 'TCA']
            # ['TCCT', 'TCTT', 'TCAT']
            # ['TCCTC', 'TCTTC', 'TCATC']
            # ['TCCTCT', 'TCTTCT', 'TCATCT']
            # ['TCCTCTA', 'TCTTCTA', 'TCATCTA']
            # ['TCCTCTAG', 'TCTTCTAG', 'TCATCTAG']
            # ['TCCTCTAGT', 'TCTTCTAGT', 'TCATCTAGT']
            # ['TCCTCTAGTG', 'TCTTCTAGTG', 'TCATCTAGTG']
            # ['TCCTCTAGTGC', 'TCTTCTAGTGC', 'TCATCTAGTGC']
            # ['TCCTCTAGTGCC', 'TCTTCTAGTGCC', 'TCATCTAGTGCC']
        return set(combos)

    def vote(
            self,
            umi_block,
            recur_len,
    ) -> Set:
        """
        It makes a vote in a umi block by choosing the most common item.

        Parameters
        ----------
        umi_block : str
            a umi building block, e.g., homo-trinucleotide
        recur_len : int
            the length of a building block consisting of recurring nucleotides, e.g. trimer 3

        Returns
        -------
        Set
            unrepeated cleotides, e.g., {'T', 'A', 'C'}
        """
        if len(set(umi_block)) == recur_len:
            return set(umi_block)
        elif len(set(umi_block)) == int(recur_len / 2) + 1:
            sdict = {base: umi_block.count(base) for base in set(umi_block)}
            return set(max(sdict, key=sdict.get))
        else:
            return set(umi_block[0])

    def umi_li(
            self,
            umi,
    ) -> Set:
        """

        Parameters
        ----------
        umi : str
            an umi CTTCCGCATTTTCCCTTTAAAGGGTTTGGGCCCCCC

        Returns
        -------
        Set
            a set containing umis, each of which is a nucleotide polymer (trinucleotides)
             collapsed into a single nucleotide.
            e.g. {'TCATCTAGTGCC', 'TGCTCTAGTGCC', 'CGTTCTAGTGCC', ..., 'TCCTCTAGTGCC'}

        """
        trimers = list()
        umis = set()
        while len(umi) != 0:
            try:
                trimers.append(set(umi[0:3]))
                umi = umi[3:]
            except:
                self.console.print("umi existing indel or wrong")
        for val in itertools.product(*trimers):
            collapsed_umi = ''.join(val)
            umis.add(collapsed_umi)
        return umis

    def cmi_li(
            self,
            umi,
    ) -> Set:
        """
        This implements nu collapsing each trimer block into one
        by Canzar's lab (Stefan Canzar, Shuang Li, and Pablo Monteagudo-Mesas).

        Notes
        -----
        This implementation is the same as split_by_mv() implemented by Jianfeng Sun.

        Parameters
        ----------
        umi : str
            an umi CTTCCGCATTTTCCCTTTAAAGGGTTTGGGCCCCCC

        Returns
        -------
        Set
            a set containing umis, each of which is a nucleotide polymer (trinucleotides)
             collapsed into a single nucleotide.
            e.g. {'TCATCTAGTGCC', 'TCTTCTAGTGCC', 'TCCTCTAGTGCC'}

        """
        trimers = list()
        umis = set()
        while len(umi) != 0:
            try:
                if len(set(umi[0:3])) < 3:
                    base = collections.Counter(umi[0:3]).most_common(1)[0][0]
                    trimers.append(base)
                    umi = umi[3:]
                else:
                    trimers.append(set(umi[0:3]))
                    umi = umi[3:]
            except:
                self.console.print("umi existing indel or wrong")
        for val in itertools.product(*trimers):
            collapse = ''.join(val)
            umis.add(collapse)
        return umis


if __name__ == "__main__":
    p = Collapse()
    # set1 = p.umi_li('CTTCCGCATTTTCCCTTTAAAGGGTTTGGGCCCCCC')
    # set2 = p.split_to_all('CTTCCGCATTTTCCCTTTAAAGGGTTTGGGCCCCCC', recur_len=3)
    # print(set1)
    # print(set2)
    # print(set1 == set2)

    # set3 = p.cmi_li('CTTCCGCATTTTCCCTTTAAAGGGTTTGGGCCCCCC')
    # set4 = p.split_by_mv('CTTCCGCATTTTCCCTTTAAAGGGTTTGGGCCCCCC', recur_len=3)
    # print(set3)
    # print(set4)
    # print(set3 == set4)
    #
    # umi_collapsed = p.majority_vote('CTGCCGCATTTTCCCTTTAAAGGGTTTGGGCCCCCC') # T1TTCCGCATTTTCCCTTTAAAGGGTTTGGGCCCCCC
    # print(umi_collapsed)

    # print(p.split_by_mv('GGTTAAGGAAGGTTGGAAGTAATT', recur_len=2))
    # print(p.split_to_all('GGTTAAGGAAGGTTGGAAGTAATT', recur_len=2))
    # print(p.majority_vote('GGTTAAGGAAGGTTGGAAGTAATT', recur_len=2))


    # @@ UMIche work Figure 2
    print(p.split_by_mv('AAATCCGGATTTCGGAAATTTGGGCACCCC', recur_len=3))

    # print(p.split_to_all('AAATCCGGATTTCGGAAATTTGGGCACCCC', recur_len=3))
    # print(p.split_to_all('AAATCCGGATTTCGGAAATTTGGGCCCCCC', recur_len=3))
    # print(p.split_to_all('AAATCCGGATTTGGGAAATTTGGGCCCCCC', recur_len=3))
    # print(p.split_to_all('AAATCCGGGTTTGGGAAATTTGGGCCCCCC', recur_len=3))
    #
    # print(p.majority_vote('AAATCCGGGTTTGGGAAATTTGGGCCCCCC', recur_len=3))
    # a1 = p.split_to_all('AAATCCGGATTTCGGAAATTTGGGCACCCC', recur_len=3)
    # a2 = p.split_to_all('AAATCCGGATTTGGGAAACTTGGGCACCCC', recur_len=3)
    # a3 = p.split_to_all('AAATCCGGATTTCGGAAATTTGGGCCCCCC', recur_len=3)
    # print(a2)
    # print(a1 - a2)
    # print(a1 - a3)
    # print(a2 - a3)

    # set4 = p.split_by_mv('CTTCCGCATTTTCCCTTTAAAGGGTTTGGGCCCCCC', recur_len=3)
    # set2 = p.split_to_all('CTTCCGCATTTTCCCTTTAAAGGGTTTGGGCCCCCC', recur_len=3)
    # print(p.split_to_all('CTTCCGCATTTTCCCTTTAAAGGGTTTGGGCCCCCC', recur_len=3))
    # print(p.split_by_mv('CTTCCGCATTTTCCCTTTAAAGGGTTTGGGCCCCCC', recur_len=3))

    # print(p.split_to_all('TTTCCGCATTTTCCCTTTAAAGGGTTTGGGCCCCCC', recur_len=3))
    # print(p.split_by_mv('TTTCCGCATTTTCCCTTTAAAGGGTTTGGGCCCCCC', recur_len=3))
