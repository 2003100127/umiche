import time
import numpy as np


class Extreme:

    def __init__(self, ):
        pass

    def pervector(self, gc_mat, axis=1):
        d_sum = np.sum(gc_mat, axis=axis)
        if axis == 1:
            return gc_mat / d_sum[:, np.newaxis]
        elif axis == 0:
            return (gc_mat.T / d_sum[:, np.newaxis]).T

    def geneIds(self, gc_mat, oriented='maximal', n_top=10, scheme={'method': 'mean',}, is_norm=False):
        """
        ..  Summary:
            --------
            find ids of extreme (highly-expressed or lowly-expressed) genes.

        ..  Description:
            ------------
            gene_ids is a list containing gene ids

        :param gc_mat:
        :param oriented: maximal or minimal
        :param n_top: the number of top-ranked maximal or minimal values
        :param scheme:
        :param is_norm:
        :return:
        """
        print('======>start to find extreme gene ids...')
        if is_norm:
            gc_mat = self.pervector(gc_mat)
        if scheme['method'] == 'mean':
            mean_stime = time.time()
            d = gc_mat.mean(axis=0)
            # print(d)
            print('======>time: {:.2f}s'.format(time.time() - mean_stime))
            return (-d).argsort()[:n_top].tolist() if oriented == 'maximal' else d.argsort()[:n_top].tolist()
        elif scheme['method'] == 'cell_ratio':
            cr_stime = time.time()
            import pandas as pd
            num_cells = gc_mat.shape[0]
            # print(gc_mat.argmax(axis=1))
            t = pd.Series(gc_mat.argmax(axis=1)).value_counts()
            c = t.rename_axis('unique_values').reset_index(name='counts')
            print('======>time for finding extreme gene ids: {:.2f}s'.format(time.time() - cr_stime))
            return c.loc[c.counts > num_cells * scheme['cell_ratio']].unique_values.values.tolist()

    def perCellClsGeneIds(self, gc_mat, cls_vec={}, oriented='maximal', n_top=10, scheme={'method': 'mean',}, is_norm=False):
        """
        ..  Summary:
            --------
            find ids of extreme (highly-expressed or lowly-expressed) genes in clusters.

        ..  Description:
            ------------
            gene_ids which is a python dictionary, is subject to the data structure below:
            {
                {cluster_1: gene_id_list_1},
                {cluster_2: gene_id_list_2},
                {cluster_3: gene_id_list_3},
                ...
                {cluster_n: gene_id_list_n},
            }

        :param gc_mat:
        :param cls_vec:
        :param oriented:
        :param n_top:
        :param scheme:
        :param is_norm:
        :return:
        """
        print('===>start to find extreme gene ids in clusters...')
        if not isinstance(cls_vec, dict):
            raise ValueError
        res = {}
        for k, cls_vec in cls_vec.items():
            cls_stime = time.time()
            # print(self.gc_mat[cls_vec, :])
            res[k] = self.geneIds(gc_mat[cls_vec, :], oriented=oriented, n_top=n_top, scheme=scheme, is_norm=is_norm)
            print('===>time: {}: {:.2f}s'.format(k, time.time() - cls_stime))
        return res


if __name__ == "__main__":
    gc_mat = np.array([
        [3., 3., 3., 6., 6.],
        [1., 1., 1., 2., 2.],
        [1., 22., 1., 2., 2.]
    ])
    p = Extreme()
    print(p.geneIds(gc_mat, oriented='maximal', n_top=3, scheme={'method': 'cell_ratio', 'cell_ratio': 0.2}))

    cls_vec = {0: [0, 2], 1: [1]}
    print(p.perCellClsGeneIds(
        gc_mat,
        cls_vec,
        oriented='maximal', n_top=3, scheme={'method': 'cell_ratio', 'cell_ratio': 0.2}
    ))