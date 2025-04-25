__version__ = "v1.0"
__copyright__ = "Copyright 2025"
__license__ = "GPL-3.0"


import h5py
import json
import scanpy


class reader:

    def __init__(self, ):
        pass

    def anndata(self, h5_fpn):
        return scanpy.read_h5ad(h5_fpn)

    def gcmatLit(self, h5_fpn):
        from scipy.sparse import csr_matrix, coo_matrix, csc_matrix
        with h5py.File(h5_fpn, 'r') as f:
            dict_genes = json.loads(f['default/dict_genes'][()])
            dict_cells = json.loads(f['default/dict_cells'][()])
            dict_cell_cls = json.loads(f['default/dict_cell_cls'][()])
            dict_genes = {int(k): v for k, v in dict_genes.items()}
            dict_cells = {int(k): v for k, v in dict_cells.items()}
            # print(len(dict_genes))
            # print(f['default/gc_mat_vals'][()])
            # print(f['default/gc_mat_ids'][()])
            # print(f['default/gc_mat_id_pointers'][()])
            gc_mat = csr_matrix(
                (f['default/gc_mat_vals'][()], f['default/gc_mat_ids'][()], f['default/gc_mat_id_pointers'][()]),
                shape=(len(dict_cells), len(dict_genes))
            ).toarray()
            # print(gc_mat)
            # print(dict_genes)
            # print(dict_cells)
            # print(dict_cell_cls)
        return {
            'gc_mat': gc_mat,
            'dict_genes': dict_genes,
            'dict_cells': dict_cells,
            'dict_cell_cls': dict_cell_cls,
        }


if __name__ == "__main__":
    p = reader()
    print(p.anndata(h5_fpn='../../../../notebook/sc/Lechner/file.h5ad'))