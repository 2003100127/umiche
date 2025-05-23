__version__ = "v1.0"
__copyright__ = "Copyright 2025"
__license__ = "GPL-3.0"
__developer__ = "Jianfeng Sun"
__maintainer__ = "Jianfeng Sun"
__email__="jianfeng.sunmt@gmail.com"

import numpy as np
import pandas as pd
import seaborn as sns
import matplotlib.pyplot as plt
import scipy.optimize as opt


class Anchor:

    def __init__(
            self,
            quant_captured,
            quant_anchor_captured,
            criteria,
            condition,
    ):
        self.quant_captured = quant_captured
        self.quant_anchor_captured = quant_anchor_captured

        self.df_quant_captured = pd.DataFrame.from_dict(self.quant_captured)
        self.df_quant_anchor_captured = pd.DataFrame.from_dict(self.quant_anchor_captured)

        self.df_quant_captured_mean = self.df_quant_captured.mean(axis=1).values
        self.df_quant_anchor_captured_mean = self.df_quant_anchor_captured.mean(axis=1).values

        self.err_captured = self.df_quant_captured.max(axis=1).values - self.df_quant_captured_mean
        self.err_anchor_captured = self.df_quant_anchor_captured.max(axis=1).values - self.df_quant_anchor_captured_mean

        # self.err_captured = pd.DataFrame.from_dict(quant_captured).std(axis=1).values
        # self.err_anchor_captured = pd.DataFrame.from_dict(quant_anchor_captured).std(axis=1).values

        self.criteria = criteria
        self.condition = condition

        self.palette = sns.color_palette("Set3")[3:5][::-1]

        self.meas = {
            'without anchor': self.df_quant_captured_mean,
            'with anchor (BAGC)': self.df_quant_anchor_captured_mean,
        }

        self.meas_errs = {
            'without anchor': self.err_captured,
            'with anchor (BAGC)': self.err_anchor_captured,
        }

        sns.set(font="Helvetica")
        sns.set_style("ticks")

    def fit_logistic(self, x, a, b, c, d):
        return a / (1. + np.exp(-c * (x - d))) + b

    def draw_line_broken(
            self,
            notion=False,
    ):
        x_notions = [str(criterion) for criterion in self.criteria]
        x = np.arange(len(x_notions))

        fig, ax = plt.subplots(nrows=2, ncols=1, figsize=(6, 5), sharex=True)
        for i in range(2):
            for j, (attr, measurement) in enumerate(self.meas.items()):
                # ax[i].scatter(
                #     x,
                #     measurement,
                #     color='grey',
                #     linewidths=2,
                #     s=60,
                #     facecolors='none',
                # )
                ax[i].errorbar(
                    x,
                    measurement,
                    yerr=self.meas_errs[attr],
                    elinewidth=2,
                    markersize=5,
                    fmt='o',
                    markerfacecolor='none',
                )

                def f(x, a, b, c, d):
                    return a / (1. + np.exp(-c * (x - d))) + b

                popt, pcov = opt.curve_fit(f, x, measurement, method="trf")
                print(popt)
                y_fit = f(x, *popt)

                # new_x = np.linspace(0, len(x_notions), 100)
                # import scipy.interpolate as spi
                # ipo3 = spi.splrep(x, measurement, k=3)+
                # iy3 = spi.splev(new_x, ipo3)
                ax[i].plot(
                    x,
                    y_fit,
                    # new_x,
                    # iy3,
                    color=self.palette[j],
                    label=attr,
                    linewidth=3,
                )
                ax[i].spines['right'].set_color('none')
                ax[i].spines['top'].set_color('none')
                ax[i].spines['bottom'].set_linewidth(2)
                ax[i].spines['left'].set_linewidth(2)

        ax[0].spines.bottom.set_visible(False)
        ax[0].xaxis.set_ticks_position('none')
        ax[1].xaxis.tick_bottom()
        ax[0].tick_params(axis='y', labelsize=14, width=2)
        ax[1].tick_params(axis='y', labelsize=14, width=2)

        if notion:
            # a / (1. + np.exp(-c * (x - d))) + b
            # import matplotlib.patches as mpatches
            ax[0].annotate(
                r'$\frac{-0.241}{1+e^{-0.709(x-11.128)}}$',
                # xy=(4, 0.9),
                # xytext=(0.38, 0.5),
                xy=(6, 0.89),
                xytext=(0.38, 0.5),
                textcoords='axes fraction',
                arrowprops=dict(
                    arrowstyle="fancy",
                    color="0.5",
                    # patchB=mpatches.Ellipse((5, 0.23), 0.3, 0.4, angle=30, alpha=0.2),
                    facecolor='dimgrey',
                    shrinkB=5,
                    connectionstyle="arc3,rad=0.3",
                ),
                fontsize=16,
                horizontalalignment='right', verticalalignment='top',
            )
            ax[1].annotate(
                r'$\frac{-0.364}{1+e^{-0.835(x-14.448)}}$',
                xy=(5, 0.23),
                xytext=(0.78, 1.2),
                textcoords='axes fraction',
                arrowprops=dict(
                    arrowstyle="fancy",
                    color="0.5",
                    # patchB=mpatches.Ellipse((5, 0.23), 0.3, 0.4, angle=30, alpha=0.2),
                    facecolor='dimgrey',
                    shrinkB=5,
                    connectionstyle="arc3,rad=0.3",
                ),
                fontsize=16,
                horizontalalignment='right', verticalalignment='top',
            )

        ax[1].set_xlabel('Insertion error rate', fontsize=18)

        ax[0].set_ylim(.59, 0.95)
        ax[1].set_ylim(-0.06, .25)

        d = .5  # proportion of vertical to horizontal extent of the slanted line
        kwargs = dict(marker=[(-1, -d), (1, d)], markersize=12,
                      linestyle="none", color='k', mec='k', mew=1, clip_on=False)
        ax[0].plot([0], [0], transform=ax[0].transAxes, **kwargs)
        ax[1].plot([0], [1], transform=ax[1].transAxes, **kwargs)

        def formatter(x, pos):
            del pos
            return '{}'.format(int(x * 100))

        ax[0].yaxis.set_major_formatter(formatter)
        ax[1].yaxis.set_major_formatter(formatter)

        # ax[1].arrow(2, 0.05, 4, 4, width=0.2)

        ax[1].set_xticks(np.linspace(0, 16, 5))
        # for i in [1e-05, 1e-04, 1e-03, 1e-02, 1e-01]
        ax[1].set_xticklabels(["$10^{-5}$", "$10^{-4}$", "$10^{-3}$", "$10^{-2}$", "$10^{-1}$"], fontsize=14)

        fig.text(
            0.04,
            0.5,
            'Percentage of reads captured (%)',
            va='center',
            rotation='vertical',
            fontsize=18,
        )

        plt.legend(ncol=1, fontsize=16)

        plt.subplots_adjust(
            top=0.94,
            bottom=0.129,
            left=0.14,
            right=0.98,
            # hspace=0.40,
            # wspace=0.15,
        )
        # plt.savefig('./data/Figure4d.eps', format='eps')
        plt.show()

    def draw(
            self,
            notion=False,
            xy=None,
            xytext=None,
    ):
        x_notions = [str(criterion) for criterion in self.criteria]
        x = np.arange(len(x_notions))
        fig, ax = plt.subplots(nrows=1, ncols=1, figsize=(6, 5), sharex=True)

        for j, (attr, measurement) in enumerate(self.meas.items()):
            # ax[i].scatter(
            #     x,
            #     measurement,
            #     color='grey',
            #     linewidths=2,
            #     s=60,
            #     facecolors='none',
            # )
            ax.errorbar(
                x,
                measurement,
                yerr=self.meas_errs[attr],
                elinewidth=2,
                markersize=5,
                fmt='o',
                markerfacecolor='none',
            )

            popt, pcov = opt.curve_fit(self.fit_logistic, x, measurement, method="trf")
            y_fit = self.fit_logistic(x, *popt)

            # new_x = np.linspace(0, len(x_notions), 100)
            # import scipy.interpolate as spi
            # ipo3 = spi.splrep(x, measurement, k=3)+
            # iy3 = spi.splev(new_x, ipo3)  # 根据观测点和样条参数，生成插值
            ax.plot(
                x,
                y_fit,
                # new_x,
                # iy3,
                color=self.palette[j],
                label=attr,
                linewidth=3,
            )
            ax.spines['right'].set_color('none')
            ax.spines['top'].set_color('none')
            ax.spines['bottom'].set_linewidth(2)
            ax.spines['left'].set_linewidth(2)

        ax.xaxis.set_ticks_position('none')
        ax.xaxis.tick_bottom()
        ax.tick_params(axis='y', labelsize=14, width=2)
        ax.set_xlabel(self.condition, fontsize=18)
        ax.set_ylim(-0.06, 1.1)
        ax.set_xticks(np.linspace(0, 16, 5))
        ax.set_xticklabels(["$10^{-5}$", "$10^{-4}$", "$10^{-3}$", "$10^{-2}$", "$10^{-1}$"], fontsize=14)

        if notion:
            # a / (1. + np.exp(-c * (x - d))) + b
            # import matplotlib.patches as mpatches
            ax.annotate(
                # pcr_ins
                # r'$\frac{-0.984}{1+e^{-0.834(x-8.747)}}$',
                # xy=(5.5, 0.88),
                # xytext=(0.38, 0.65),

                # # seq_ins
                # r'$\frac{-1.064}{1+e^{-0.706(x-11.176)}}$',
                # xy=(7, 0.92),
                # xytext=(0.48, 0.7),

                # seq_del
                # r'$\frac{-1.064}{1+e^{-0.704(x-11.165)}}$',
                # xy=(7, 0.92),
                # xytext=(0.48, 0.7),

                # pcr_del
                # r'$\frac{-0.979}{1+e^{-0.865(x-8.742)}}$',
                # xy=(5.5, 0.88),
                # xytext=(0.38, 0.65),

                r'$\frac{-0.979}{1+e^{-0.865(x-8.742)}}$',
                xy=xy,
                xytext=xytext,
                textcoords='axes fraction',
                arrowprops=dict(
                    arrowstyle="fancy",
                    color="0.5",
                    # patchB=mpatches.Ellipse((5, 0.23), 0.3, 0.4, angle=30, alpha=0.2),
                    facecolor='dimgrey',
                    shrinkB=5,
                    connectionstyle="arc3,rad=-0.3",
                ),
                fontsize=18,
                horizontalalignment='right', verticalalignment='top',
            )
            ax.annotate(
                # pcr_ins
                # r'$\frac{-0.689}{1+e^{-0.916(x-13.732)}}$',
                # xy=(14.2, 0.55),
                # xytext=(0.95, 0.35),

                # seq_ins
                # r'$\frac{-0.293}{1+e^{-0.919(x-14.337)}}$',
                # xy=(15.2, 0.8),
                # xytext=(1.0, 0.6),

                # seq_del
                # r'$\frac{-0.381}{1+e^{-0.895(x-14.318)}}$',
                # xy=(15, 0.75),
                # xytext=(1.0, 0.6),

                # pcr_del
                # r'$\frac{-0.832}{1+e^{-0.834(x-13.753)}}$',
                # xy=(14., 0.5),
                # xytext=(0.95, 0.35),

                r'$\frac{-0.832}{1+e^{-0.834(x-13.753)}}$',
                xy=xy,
                xytext=xytext,
                textcoords='axes fraction',
                arrowprops=dict(
                    arrowstyle="fancy",
                    color="0.5",
                    # patchB=mpatches.Ellipse((5, 0.23), 0.3, 0.4, angle=30, alpha=0.2),
                    facecolor='dimgrey',
                    shrinkB=5,
                    connectionstyle="arc3,rad=-0.3",
                ),
                fontsize=18,
                horizontalalignment='right', verticalalignment='top',
            )

        def formatter(x, pos):
            del pos
            return '{}'.format(int(x * 100))

        ax.yaxis.set_major_formatter(formatter)
        ax.yaxis.set_major_formatter(formatter)

        fig.text(
            0.04,
            0.5,
            'Percentage of reads recovered (%)',
            va='center',
            rotation='vertical',
            fontsize=18,
        )
        plt.legend(ncol=1, fontsize=16)
        plt.subplots_adjust(
            top=0.94,
            bottom=0.129,
            left=0.14,
            right=0.98,
            # hspace=0.40,
            # wspace=0.15,
        )
        # plt.savefig('../data/' + x_label + '.eps', format='eps')
        plt.show()

    def simple(self, ):
        x_notions = [str(criterion) for criterion in self.criteria]
        x = np.arange(len(x_notions))
        width = 0.35
        multiplier = 0

        fig, ax = plt.subplots(layout='constrained', figsize=(10, 6))

        for attr, measurement in self.meas.items():
            print(attr, measurement)
            offset = width * multiplier
            rects = ax.bar(x + offset, measurement, width, label=attr)
            ax.bar_label(rects, fmt="{:.1%}", padding=2, fontsize=10)
            multiplier += 1

        ax.set_ylabel('Successfully captured reads', fontsize=16)
        ax.set_title('scCOLOR-seq v2 bead design simulation', fontsize=14)
        ax.set_xticks(x + width, x_notions, ha='right', rotation=20)
        ax.set_xlabel(self.condition, fontsize=16)
        ax.spines['right'].set_visible(False)
        ax.spines['top'].set_visible(False)

        ax.legend(ncols=2)
        ax.set_ylim(0, 1)

        plt.show()