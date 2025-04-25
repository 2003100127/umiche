__version__ = "v1.0"
__copyright__ = "Copyright 2025"
__license__ = "GPL-3.0"


import time
import h5py
import json


class writer(object):

    def __init__(self, ):
        pass

    def gcmatLit(self, gc_dict, sv_fpn, compression_level=4):
        """
        ..  Description:
            ------------
            create a light version of a gene-by-cell matrix from a gc dict consisting of
            a gc matrix of a csr version (see link 1) with three vectors:
                1) gene expression values;
                2) row indices of cells;
                3) index pointers;
            and with three json encoded dicts of
                4) gene notations and
                5) cell names,
                6) cells clusters (using cell names)
            , respectively.

            gc_mat_vals: data of scipy csr matrix

            gc_mat_ids: indices of scipy csr matrix

            gc_mat_id_pointers: indptr of scipy csr matrix

            gene_notations:
            {
                gene_notation_1: 1, # int value
                gene_notation_2: 2, # int value
                gene_notation_3: 3, # int value
                ...
                gene_notation_n: n, # int value
            }

            cell_names:
            {
                cell_name_1: 1, # int value
                cell_name_2: 2, # int value
                cell_name_3: 3, # int value
                ...
                cell_name_m: m, # int value
            }

            dict_cell_cls (i.e., cell clusters) are structured as
            {
                cell_cluster_1: cell_indices_in_cluster_1, # python list
                cell_cluster_2: cell_indices_in_cluster_2, # python list
                cell_cluster_3: cell_indices_in_cluster_3, # python list
                ...
                cell_cluster_k: cell_indices_in_cluster_k, # python list
            }

        ..  See:
            ----
            1. https://docs.scipy.org/doc/scipy/reference/generated/scipy.sparse.csr_matrix.html;
            2. https://docs.h5py.org/en/stable/high/dataset.html?highlight=compression#lossless-compression-filters

        :param gc_dict:
        :param sv_fpn:
        :param compression_level: corresponds to compression_opts in hdf5 (see link 2).
        :return: a hdf5 ob
        """
        l2d_stime = time.time()
        with h5py.File(sv_fpn, 'w') as f:
            g1 = f.create_group('default')
            g1.create_dataset(
                'gc_mat_vals',
                data=gc_dict['gc_mat_vals'],
                compression='gzip',
                compression_opts=compression_level,
            )
            g1.create_dataset(
                'gc_mat_ids',
                data=gc_dict['gc_mat_ids'],
                compression='gzip',
                compression_opts=compression_level,
            )
            g1.create_dataset(
                'gc_mat_id_pointers',
                data=gc_dict['gc_mat_id_pointers'],
                compression='gzip',
                compression_opts=compression_level,
            )
            g1.create_dataset(
                'dict_genes',
                data=json.dumps(gc_dict['dict_genes'])
            )
            g1.create_dataset(
                'dict_cells', data=json.dumps(gc_dict['dict_cells'])
            )
            g1.create_dataset(
                'dict_cell_cls', data=json.dumps(gc_dict['dict_cell_cls'])
            )
        f.close()
        print('===>time for saving h5: {:.2f}'.format(time.time() - l2d_stime))
        # with h5py.File(sv_fpn, 'r') as f:
        #     print(f['default/gc_mat_vals'][()])
        #     print(json.loads(f['default/cell_names'][()]))
        return 0