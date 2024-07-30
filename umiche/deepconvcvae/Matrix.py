import numpy as np


class matrix(object):

    def __init__(self, gc_mat):
        """

        :param gc_mat: 2d array
        """
        self.gc_mat = gc_mat

    def pervector(self, axis=1):
        d_sum = np.sum(self.gc_mat, axis=axis)
        if axis == 1:
            return self.gc_mat / d_sum[:, np.newaxis]
        elif axis == 0:
            return (self.gc_mat.T / d_sum[:, np.newaxis]).T


if __name__ == "__main__":
    gc_mat = np.array([
        [3., 3., 3., 6., 6.],
        [1., 1., 1., 2., 2.],
        [1., 22., 1., 2., 2.]
    ])
    p = matrix(gc_mat)
    print(p.pervector(axis=0))