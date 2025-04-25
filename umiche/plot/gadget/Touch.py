__author__ = "Jianfeng Sun"
__version__ = "v1.0"
__copyright__ = "Copyright 2025"
__license__ = "GPL v3.0"
__email__ = "jianfeng.sunmt@gmail.com"
__maintainer__ = "Jianfeng Sun"


class Touch:

    # def __init__(self, task):
    #     self.task = task

    def bar_height_label(
            self,
            canvas,
            ax,
            format,
            fontsize=8,
    ):
        """
        Attach a text label above each bar in *rects*, displaying its height.

        Parameters
        ----------
        canvas
        ax
        format
        fontsize

        Returns
        -------

        """
        heights = []
        for rect in canvas:
            height = rect.get_height()
            heights.append(height)
            # print(height)
            if format == 'pct':
                val_formatted = '{:.2%}'.format(height)
            elif format == 'int':
                val_formatted = '{:d}'.format(int(height))
            else:
                val_formatted = '{}'.format(height)
            ax.annotate(
                val_formatted,
                xy=(rect.get_x() + rect.get_width() / 2, height),
                xytext=(0, 3),  # 3 points vertical offset
                textcoords="offset points",
                ha='center',
                va='bottom',
                fontsize=fontsize,
            )
        return max(heights)