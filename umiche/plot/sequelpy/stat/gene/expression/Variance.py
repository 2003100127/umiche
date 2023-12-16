
import numpy as np


class variance(object):

    def __init__(self, gc_mat):
        self.gc_mat = gc_mat

    def pergene(self, ):
        return self.gc_mat.var(axis=0)

    def percell(self, ):
        return self.gc_mat.var(axis=1)

    def all(self, ):
        return self.gc_mat.var()


if __name__ == "__main__":
    gc_mat = np.array([
        [3., 3., 3., 6., 6.],
        [1., 1., 1., 2., 2.],
        [1., 22., 1., 2., 2.]
    ])
    p = variance(gc_mat)
    print(p.pergene())
    print(p.percell())
    print(p.all())