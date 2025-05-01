__version__ = "v1.0"
__copyright__ = "Copyright 2025"
__license__ = "GPL-3.0"
__developer__ = "Jianfeng Sun"
__maintainer__ = "Jianfeng Sun"
__email__="jianfeng.sunmt@gmail.com"


import json
import networkx as nx
import seaborn as sns
import matplotlib.pyplot as plt
from umiche.plot.gadget.Element import Element as pele
from umiche.network.Adjacency import Adjacency as netadj


class Cluster:

    def __init__(
            self,
            json_fpn=None,
    ):
        self.pele = pele()
        self.palette = self.pele.color(which='tableau', is_random=True)

        sns.set(font="Helvetica")
        sns.set_style("ticks")

        if json_fpn:
            with open(json_fpn) as fp:
                self.cnt_dict = json.load(fp)
                print(self.cnt_dict['0']['0.025'])



if __name__ == "__main__":
    from umiche.path import to

    p = Cluster(
        json_fpn=to('data/simu/mclumi/seq_errs/mcl_node_repr.json')
    )
