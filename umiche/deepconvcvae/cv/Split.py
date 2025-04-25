__version__ = "v1.0"
__copyright__ = "Copyright 2025"
__license__ = "GPL-3.0"


import numpy as np
import pandas as pd
from sklearn.model_selection import StratifiedShuffleSplit


class split:

    def __init__(self,):
        pass

    def sss(self, data_dict, num_splits=5, train_ratio=1, test_ratio=1):
        """
        ..  Summary:
            --------
            split by the stratified shuffle split (SSS) method.

        ..  Description:
            ------------
            res is a dict
            {
                cv0: {
                    train: train_row_ids,
                    test: train_row_ids,
                    train_cls: train_cell_cls,
                    test_cls: test_cell_cls,
                },
                cv2: {
                    train: train_row_ids,
                    test: train_row_ids,
                    train_cls: train_cell_cls,
                    test_cls: test_cell_cls,
                }
                cv3: {
                    train: train_row_ids,
                    test: train_row_ids,
                    train_cls: train_cell_cls,
                    test_cls: test_cell_cls,
                },
                ...
                cvt: {
                    train: train_row_ids,
                    test: train_row_ids,
                    train_cls: train_cell_cls,
                    test_cls: test_cell_cls,
                }
            }

        :param data_dict:
        :param num_splits:
        :param train_ratio:
        :param test_ratio:
        :return:
        """
        df = pd.DataFrame.from_dict(data_dict['dict_cells'], orient='index', columns=['cluster'])
        uniq = df.cluster.unique()
        clusters = {e: i for i, e in enumerate(uniq)}
        clusters_ = {i: e for i, e in enumerate(uniq)}
        df['cls'] = df.cluster.apply(lambda x: clusters[x])
        # print(df)
        X = df.index.values
        y = df.cls.values
        sss = StratifiedShuffleSplit(
            n_splits=num_splits,
            test_size=test_ratio/(train_ratio + test_ratio),
            random_state=0,
        )
        res = {}
        for i_cv, (train_index, test_index) in enumerate(sss.split(X, y)):
            print("Train:", train_index.shape, "Test:", test_index.shape)
            df_cv_train = df.iloc[train_index]
            df_cv_test = df.iloc[test_index]
            df_cv_train_gp = df_cv_train.groupby(by=['cls'])
            df_cv_test_gp = df_cv_test.groupby(by=['cls'])
            cv_train_keys = df_cv_train_gp.groups.keys()
            cv_test_keys = df_cv_test_gp.groups.keys()
            # print(cv_train_keys)
            # print(cv_test_keys)
            cv_train_cls_dict = {clusters_[i]: df_cv_train_gp.get_group(i).index.values.astype(int) for i in cv_train_keys}
            cv_test_cls_dict = {clusters_[i]: df_cv_test_gp.get_group(i).index.values.astype(int) for i in cv_test_keys}
            res['cv' + str(i_cv)] = {
                'train_index': train_index,
                'test_index': test_index,
                'train_cls': cv_train_cls_dict,
                'test_cls': cv_test_cls_dict,
                'train_cls_all': df.iloc[train_index].cls.values,
                'test_cls_all': df.iloc[test_index].cls.values,
            }
        # print(res)
        return res

    def kfold(self, res):
        pass


if __name__ == "__main__":
    p = split()