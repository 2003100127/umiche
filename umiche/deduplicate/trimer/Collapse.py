__version__ = "v1.0"
__copyright__ = "Copyright 2024"
__license__ = "MIT"
__developer__ = "Jianfeng Sun"
__maintainer__ = "Jianfeng Sun"
__email__="jianfeng.sunmt@gmail.com"
__lab__ = "Cribbslab"

import os
import pysam
import pandas as pd
import numpy as np
import itertools
import gurobipy as gp
from gurobipy import GRB
import collections
import pdb
from umiche.util.Number import number as rannum
from umiche.path import to


class collapse:

    def majorityVoting(self, umi, recur_len=3):
        vernier = [i for i in range(len(umi)) if i%recur_len == 0]
        umi_oligomers = [umi[v: v+recur_len] for v in vernier]
        t = []
        for umi_oligomer in umi_oligomers:
            s = set(umi_oligomer)
            if len(s) == recur_len:
                rand_index = rannum().uniform(low=0, high=recur_len, num=1, use_seed=False, seed=3)[0]
                t.append(umi_oligomer[rand_index])
            elif len(s) == int(recur_len/2) + 1:
                sdict = {base: umi_oligomer.count(base) for base in s}
                t.append(max(sdict, key=sdict.get))
            else:
                t.append(umi_oligomer[0])
        return ''.join(t)

    def splitToAll(self, umi, recur_len=3):
        vernier = [i for i in range(len(umi)) if i % recur_len == 0]
        umi_oligomers = [umi[v: v + recur_len] for v in vernier]
        t = [umi_oligomer for umi_oligomer in set(umi_oligomers[0])]
        umi_oligomers = umi_oligomers[1:]
        for id, umi_oligomer in enumerate(umi_oligomers):
            umi_oligomer_bases = set(umi_oligomer)
            c = []
            for i in range(len(t)):
                for umi_oligomer_base in umi_oligomer_bases:
                    c.append(t[i] + umi_oligomer_base)
            t = c
        return set(t)

    def splitByMV(self, umi, recur_len=3):
        vernier = [i for i in range(len(umi)) if i % recur_len == 0]
        umi_oligomers = [umi[v: v + recur_len] for v in vernier]
        umi_oligomer_bases = self.vote(umi_oligomers[0], recur_len)
        t = [umi_oligomer for umi_oligomer in umi_oligomer_bases]
        umi_oligomers = umi_oligomers[1:]
        for id, umi_oligomer in enumerate(umi_oligomers):
            umi_oligomer_bases = self.vote(umi_oligomer, recur_len)
            c = []
            for i in range(len(t)):
                for umi_oligomer_base in umi_oligomer_bases:
                    c.append(t[i] + umi_oligomer_base)
            t = c
        return set(t)

    def vote(self, umi_oligomer, recur_len):
        if len(set(umi_oligomer)) == recur_len:
            return set(umi_oligomer)
        elif len(set(umi_oligomer)) == int(recur_len / 2) + 1:
            sdict = {base: umi_oligomer.count(base) for base in set(umi_oligomer)}
            return set(max(sdict, key=sdict.get))
        else:
            return set(umi_oligomer[0])

    def ShuangLiCollapseUMI(self, string):
        trimers = list()
        umis = set()
        while len(string) != 0:
            try:
                trimers.append(set(string[0:3]))
                string = string[3:]
            except:
                print("umi existing indel or wrong")
        for val in itertools.product(*trimers):
            collapse = ''.join(val)
            umis.add(collapse)
        return umis

    def allbasesSame(self, s):
        if len(set(s)) < 3:
            return True
        else:
            return False

    def ShuangLiCollapseCMI(self, string):
        trimers = list()
        umis = set()
        while len(string) != 0:
            try:
                if self.allbasesSame(string[0:3]):
                    base = collections.Counter(string[0:3]).most_common(1)[0][0]
                    trimers.append(base)
                    string = string[3:]
                else:
                    trimers.append(set(string[0:3]))
                    string = string[3:]
            except:
                print("umi existing indel or wrong")

        for val in itertools.product(*trimers):
            collapse = ''.join(val)
            umis.add(collapse)
        return umis


if __name__ == "__main__":
    p = collapse()
    # u1 = p.ShuangLiCollapseUMI('CTTCCGCATTTTCCCTTTAAAGGGTTTGGGCCCCCC')
    u2 = p.ShuangLiCollapseCMI('CTTCCGCATTTTCCCTTTAAAGGGTTTGGGCCCCCC')
    # u3 = p.majorityVoting('T1TTCCGCATTTTCCCTTTAAAGGGTTTGGGCCCCCC')
    # u4 = p.splitToAll('CTTCCGCATTTTCCCTTTAAAGGGTTTGGGCCCCCC')
    u5 = p.splitByMV('CTTCCGCATTTTCCCTTTAAAGGGTTTGGGCCCCCC', recur_len=3)
    # print(u1)
    print(u2)
    # print(u3)
    # print(u4)
    print(u5)
    # print(u1 == u4)
    print(u2 == u5)