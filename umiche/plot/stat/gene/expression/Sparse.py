
import numpy as np


class mean(object):

    def __init__(self, gc_mat):
        self.gc_mat = gc_mat

    def pergene(self, ):
        return self.gc_mat.mean(axis=0)

    def percell(self, ):
        return self.gc_mat.mean(axis=1)

    def percellcls(self, cls_vec={}):
        if not isinstance(cls_vec, dict):
            raise ValueError
        res = {}
        for k, v in cls_vec.items():
            res[k] = self.gc_mat[v, :].mean(axis=1).tolist()
        return res

    def allcellcls(self, cls_vec={}):
        if not isinstance(cls_vec, dict):
            raise ValueError
        res = {}
        for k, v in cls_vec.items():
            res[k] = self.gc_mat[v, :].mean().tolist()
        return res

    def pergenecls(self, cls):
        pass

    def allgenecls(self, cls):
        pass

    def all(self, ):
        return self.gc_mat.mean()


if __name__ == "__main__":
    gc_mat = np.array([
        [3., 3., 3., 6., 6.],
        [1., 1., 1., 2., 2.],
        [1., 22., 1., 2., 2.]
    ])
    p = mean(gc_mat)
    print(p.pergene())
    print(p.percell())
    print(p.all())

    cls_vec = {0: [0, 2], 1: [1]}
    print(p.percellcls(cls_vec))
    print(p.allcellcls(cls_vec))