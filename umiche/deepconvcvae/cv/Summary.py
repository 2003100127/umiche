__version__ = "v1.0"
__copyright__ = "Copyright 2025"
__license__ = "GPL-3.0"
__developer__ = "Jianfeng Sun"
__maintainer__ = "Jianfeng Sun"
__email__="jianfeng.sunmt@gmail.com"


from umiche.deepconvcvae.Extreme import extreme as statextreme


class summary(object):

    def __init__(self, ):
        pass

    def extreme(self, cv_dict, gc_mat, type='train', n_top=10):
        std_dict = {}
        for cv_iter, v in cv_dict.items():
            print('===>cv_{}:'.format(cv_iter))
            std_dict[cv_iter] = {}
            extrem_gene_ids = statextreme().perCellClsGeneIds(
                gc_mat=gc_mat,
                cls_vec=v['train_cls'] if type == 'train' else v['test_cls'],
                oriented='maximal',
                n_top=n_top,
                scheme={
                    'method': 'mean',
                },
            )
            retrv_ids = v['train_index'] if type == 'train' else v['test_index']
            gc_mat_cv = gc_mat[retrv_ids, :]
            std_dict[cv_iter]['val_dict_arr_1d'] = {cluster_name: gc_mat_cv[:, sel_gene_ids].flatten() for cluster_name, sel_gene_ids in extrem_gene_ids.items()}
        return std_dict