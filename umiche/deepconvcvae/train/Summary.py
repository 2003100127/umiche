__version__ = "v1.0"
__copyright__ = "Copyright 2025"
__license__ = "GPL-3.0"
__developer__ = "Jianfeng Sun"
__maintainer__ = "Jianfeng Sun"
__email__="jianfeng.sunmt@gmail.com"


from umiche.deepconvcvae.Extreme import extreme as statextreme


class summary:

    def __init__(self, ):
        pass

    def extreme(self, gc_mat, cls_vec, n_top=10):
        print('===>train_sum_extreme')
        extrem_gene_ids = statextreme().perCellClsGeneIds(
            gc_mat=gc_mat,
            cls_vec=cls_vec,
            oriented='maximal',
            n_top=n_top,
            scheme={
                'method': 'mean',
            },
        )
        return {
            cluster_name: gc_mat[:, sel_gene_ids].flatten() for cluster_name, sel_gene_ids in extrem_gene_ids.items()
        }
