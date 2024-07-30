import numpy as np
import seaborn as sns
import matplotlib.pyplot as plt
from sklearn.decomposition import PCA
from sklearn.manifold import TSNE


class dimensionReduction:

    def __init__(self, ):
        pass

    def single(self, X, y, tech='TSNE', marker_size=2, cmap='tab20b', title=''):
        """https://matplotlib.org/stable/tutorials/colors/colormaps.html"""
        sns.set(font="Helvetica")
        sns.set_style("ticks")
        if tech == 'PCA':
            retech = PCA(n_components=2)
            embeddings = retech.fit_transform(np.array(X))

        elif tech == 'TSNE':
            retech = TSNE(n_components=2)
            embeddings = retech.fit_transform(np.array(X))

        elif tech == 'UMAP':
            import umap
            retech = umap.UMAP(n_components=2)
            embeddings = retech.fit_transform(np.array(X))

        else:
            embeddings = np.array(X)
        fig, ax = plt.subplots(1, 1, figsize=(10, 8))
        sp = ax.scatter(
            embeddings[:, 0], embeddings[:, 1],
            c=y,
            s=marker_size,
            cmap=plt.cm.get_cmap(cmap, 10),
        )
        ax.spines['right'].set_visible(False)
        ax.spines['top'].set_visible(False)
        ax.set_xlabel(tech + ' 1', fontsize=14)
        ax.set_ylabel(tech + ' 2', fontsize=14)
        ax.set_title(title, fontsize=14)
        fig.colorbar(sp)
        fig.subplots_adjust(
            top=0.95,
            bottom=0.1,
            left=0.1,
            right=0.98,
            # hspace=0.20,
            # wspace=0.15
        )
        plt.show()