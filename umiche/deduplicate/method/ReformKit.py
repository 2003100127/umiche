__version__ = "v1.0"
__copyright__ = "Copyright 2025"
__license__ = "GPL-3.0"
__developer__ = "Jianfeng Sun"
__maintainer__ = "Jianfeng Sun"
__email__="jianfeng.sunmt@gmail.com"



class ReformKit:

    def __init__(self, ):
        pass

    def maxid(
            self,
            df_umi_uniq_val_cnt,
            cc,
    ):
        """

        Parameters
        ----------
        df_umi_uniq_val_cnt
        cc

        Returns
        -------

        """
        xx = [*df_umi_uniq_val_cnt.loc[df_umi_uniq_val_cnt.index.isin(cc)].sort_values(ascending=False).to_dict().keys()]
        return xx[0]

    def sort_vals(
            self,
            df_umi_uniq_val_cnt,
            cc,
    ):
        """

        Parameters
        ----------
        df_umi_uniq_val_cnt
        cc

        Returns
        -------

        """
        return df_umi_uniq_val_cnt.loc[df_umi_uniq_val_cnt.index.isin(cc)].sort_values(ascending=False).to_dict()

    def breakpoint(self, x, connected_components):
        print('MCL clusters:\n {}'.format(x['clusters']))
        print('CCS:\n {}'.format(len(connected_components)))
        print('CC:\n {}'.format(x['cc_vertices']))
        print('Edge list of MCL clusters:\n {}'.format(x['edge_list']))
        print('Graph adj list of MCL clusters:\n {}'.format(x['graph_cc_adj']))
        print('\n')
        # df.to_csv('./sd.txt', sep='\t', index=False, header=True)
        return

    def keymap(
            self,
            graph_adj,
            reverse=False,
    ):
        """

        Parameters
        ----------
        graph_adj
        reverse

        Returns
        -------

        """
        keys = [*graph_adj.keys()]
        glen = len(keys)
        if reverse:
            return {k: keys[k] for k in range(glen)}
        else:
            return {keys[k]: k for k in range(glen)}

    def graph_cc_adj(self, cc, graph_adj):
        """

        Parameters
        ----------
        cc
            The first parameter.
        graph_adj
            The se parameter.

        Returns
        -------

        """
        return {node: graph_adj[node] for node in cc}

    def key2node(self, list_2d, keymap):
        """

        Parameters
        ----------
        list_2d
        keymap

        Returns
        -------

        """
        return [[keymap[i] for i in lis] for lis in list_2d]

    def onehot(
            self,
            umi,
    ):
        import numpy as np
        from umiche.util.Single import Single as dnasgl
        nt_to_int_dict = dnasgl().todict(nucleotides=dnasgl().get(universal=True), reverse=False)
        # int_to_nt_dict = dnasgl().todict(nucleotides=dnasgl().get(universal=True), reverse=True)
        umi_ltr2digi = [nt_to_int_dict[i] for i in list(umi)]
        ids_first_pos = [i*4 for i in range(len(umi))]
        ids_to_be_one = [i+j for i, j in zip(umi_ltr2digi, ids_first_pos)]
        # print(ids_to_be_one)
        one_hot = np.zeros(len(umi)*4)
        one_hot[ids_to_be_one] = 1
        return one_hot.astype(int)
