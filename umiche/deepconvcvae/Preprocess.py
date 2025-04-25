__version__ = "v1.0"
__copyright__ = "Copyright 2025"
__license__ = "GPL-3.0"
__developer__ = "Jianfeng Sun"
__maintainer__ = "Jianfeng Sun"
__email__="jianfeng.sunmt@gmail.com"


import pandas as pd
from umiche.deepconvcvae.h5o.Reader import reader as h5reader
from umiche.deepconvcvae.cv.Summary import summary as cvsum
from umiche.deepconvcvae.train.Summary import summary as trainsum
from umiche.deepconvcvae.boxplot.GCMat import gcmat as splboxplot
from umiche.deepconvcvae.cv.Split import split as cvsplit


class preprocess:

    def __init__(self, data_fpn, train_ratio=3, test_ratio=2):
        self.cvsum = cvsum()
        self.trainsum = trainsum()
        self.splboxplot = splboxplot()
        self.h5reader = h5reader()
        self.cvsplit = cvsplit()
        self.data_dict = self.h5reader.gcmatLit(data_fpn)
        self.cv_dict = self.cvsplit.sss(data_dict=self.data_dict, train_ratio=train_ratio, test_ratio=test_ratio)
        # print(self.cv_dict)
        print([*self.data_dict['dict_cell_cls'].keys()])
        # print(len(self.data_dict['dict_genes']))
        # print(len(self.data_dict['dict_cells']))

    def plotSummary(self, ):
        # ### /* train sum plot */
        train_extreme_dict = self.trainsum.extreme(
            gc_mat=self.data_dict['gc_mat'],
            cls_vec=self.data_dict['dict_cell_cls']
        )
        train_std_dict = {'real': {'val_dict_arr_1d': train_extreme_dict}}
        print('asdasd',train_std_dict)
        self.splboxplot.boxplot(train_std_dict)

        # ### /* cv sum plot */
        cv_std_dict = self.cvsum.extreme(cv_dict=self.cv_dict, gc_mat=self.data_dict['gc_mat'], type='train')
        self.splboxplot.boxplot(cv_std_dict)

    def recipeData(self, cv=0):
        df_cell_cls = pd.DataFrame.from_dict(self.data_dict['dict_cells'], orient='index', columns=['cluster'])
        uniq = df_cell_cls.cluster.unique()
        clusters = {e: i for i, e in enumerate(uniq)}
        df_cell_cls['cls'] = df_cell_cls.cluster.apply(lambda x: clusters[x])
        # print(df_cell_cls)

        train_retrv_ids = self.cv_dict['cv' + str(cv)]['train_index']
        test_retrv_ids = self.cv_dict['cv' + str(cv)]['test_index']

        # from ext.dacube.dataset.gcmat.preprocess.Remove import remove as gcrv
        # self.data_dict['gc_mat'] = gcrv().lowlyexpressed(self.data_dict['gc_mat'], thres=1000)

        x_train = self.data_dict['gc_mat'][train_retrv_ids, :]
        y_train = df_cell_cls['cls'].values[train_retrv_ids]
        # print(x_train.shape)
        # print(y_train.shape)

        x_test = self.data_dict['gc_mat'][test_retrv_ids, :]
        y_test = df_cell_cls['cls'].values[test_retrv_ids]
        # print(x_test.shape)
        # print(y_test.shape)
        del self.data_dict
        del self.cv_dict
        return (x_train, y_train), (x_test, y_test)

    def augment(self, x_train, y_train, x_test, y_test, image_size):
        import numpy as np
        # print(x_train.shape)
        # print(y_train.shape)
        x_train_height = x_train.shape[0]
        x_test_height = x_test.shape[0]
        x_train_width = x_train.shape[1]

        img_col_offset_size = image_size * image_size - x_train_width
        train_col_offset_mat = np.zeros([x_train_height, img_col_offset_size])
        test_col_offset_mat = np.zeros([x_test_height, img_col_offset_size])
        x_train = np.concatenate([x_train, train_col_offset_mat], axis=1)
        x_test = np.concatenate([x_test, test_col_offset_mat], axis=1)
        x_train = np.reshape(x_train, newshape=(-1, image_size, image_size))
        x_test = np.reshape(x_test, newshape=(-1, image_size, image_size))
        # print(x_train)
        # print(x_test)

        x_train = np.reshape(x_train, [-1, image_size, image_size, 1])
        x_test = np.reshape(x_test, [-1, image_size, image_size, 1])
        print(x_train.shape)
        print(x_test.shape)
        std_max_mark = np.maximum(np.amax(x_train), np.amax(x_test))
        print(std_max_mark)
        x_train = x_train.astype('float32') / std_max_mark
        x_test = x_test.astype('float32') / std_max_mark
        print(y_train)
        print(y_test)
        return (x_train, y_train), (x_test, y_test)


if __name__ == "__main__":
    from umiche.path import to

    params = {
        'data_fpn': to('data/sc/gmat/pbmc68k/dA_flt/flt.h5'),
        'data_rv_zero_fpn': to('data/sc/gmat/pbmc68k/dA_flt/flt_rv_zero_col.h5'),
    }
    p = preprocess(params['data_rv_zero_fpn'])
    p.plotSummary()
    # (x_train, y_train), (x_test, y_test) = p.recipeData(cv=0)
    # print((x_train.shape, y_train.shape), (x_test.shape, y_test.shape))
    # print(len(p.data_dict))