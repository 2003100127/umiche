__version__ = "v1.0"
__copyright__ = "Copyright 2025"
__license__ = "GPL-3.0"
__developer__ = "Jianfeng Sun"
__maintainer__ = "Jianfeng Sun"
__email__="jianfeng.sunmt@gmail.com"


from typing import Dict

import numpy as np
import matplotlib.pyplot as plt
import seaborn as sns


class TripletErrorCode:

    def __init__(self, ):
        sns.set(font="Helvetica")
        sns.set_style("ticks")

    def incorrect_all(
            self,
            num_nt=12,
    ):
        palette = sns.color_palette("Set2")[3:]

        fig, ax = plt.subplots(nrows=1, ncols=1, figsize=(7, 5), sharey=True, sharex=False)

        p = np.linspace(0.00001, 0.5, 500)
        # standard
        P_std_umi_error = 1 - (1 - p) ** num_nt
        # homotrimer
        # P_block_error = 3 * (p ** 2) * (1 - p) + p ** 3
        P_block_error = (7/3) * p**2 * (1 - p) + p**3
        P_block_error_upper_bound = 3 * p**2 - 2 * p**3
        # Homotrimer UMI 的错误率：num_nt 个三联体中至少一个出错
        P_homo_umi_error = 1 - (1 - P_block_error) ** num_nt
        # 绘图
        ax.grid(False)
        ax.plot(
            p,
            P_std_umi_error,
            label='Standard 12-bp UMI',
            color=palette[1],
            linewidth=2.5,
        )
        ax.plot(
            p,
            P_homo_umi_error,
            label='Homotrimer 36-bp UMI',
            color=palette[0],
            linewidth=2.5,
        )
        ax.spines['right'].set_color('none')
        ax.spines['top'].set_color('none')
        ax.tick_params(axis='both', which='major', labelsize=16)
        ax.tick_params(axis='both', which='minor', labelsize=16)
        ax.set_xlabel("Single-base error rate ($p$)", fontsize=16)
        ax.set_ylabel("Probability of incorrectly synthesising a UMI", fontsize=16)
        # ax.set_title("Comparison of UMI error probability vs. per-base error rate", fontsize=14)
        ax.legend(fontsize=12, frameon=False)
        plt.tight_layout()
        plt.show()
        return

    def correct_all(
            self,
            num_nt=12,
    ):
        palette = sns.color_palette("Set2")[3:]

        fig, ax = plt.subplots(nrows=1, ncols=1, figsize=(7, 5), sharey=True, sharex=False)

        # Define range of p values
        p = np.linspace(0.0001, 0.5, 500)

        # Standard UMI: probability all 12 bases are correct
        P_std_umi_correct = (1 - p) ** num_nt

        # Homotrimer UMI: block error probability
        # P_block_error = 3 * (p ** 2) * (1 - p) + p ** 3
        P_block_error = (7/3) * p**2 * (1 - p) + p**3
        # Homotrimer UMI: probability all 12 blocks are correct after majority voting
        P_homo_umi_correct = (1 - P_block_error) ** num_nt

        ax.grid(False)
        ax.plot(
            p,
            P_std_umi_correct,
            label='Standard 12-bp UMI',
            color=palette[1],
            linewidth=2.5,
        )
        ax.plot(
            p,
            P_homo_umi_correct,
            label='Homotrimer 36-bp UMI',
            color=palette[0],
            linewidth=2.5,
        )

        ax.spines['right'].set_color('none')
        ax.spines['top'].set_color('none')
        ax.tick_params(axis='both', which='major', labelsize=16)
        ax.tick_params(axis='both', which='minor', labelsize=16)
        ax.set_xlabel("Single-base error rate ($p$)", fontsize=16)
        ax.set_ylabel("Probability of correctly synthesising a UMI", fontsize=16)
        # ax.set_title("Comparison of UMI error probability vs. per-base error rate", fontsize=14)
        ax.legend(fontsize=12, frameon=False)
        plt.tight_layout()
        plt.show()
        return

    def incorrect_block(
            self,
            num_nt=12,
    ):
        palette = sns.color_palette("Set2")[3:]

        fig, ax = plt.subplots(nrows=1, ncols=1, figsize=(7, 5), sharey=True, sharex=False)

        # 定义 p 值范围：单碱基错误率从 0.0001 到 0.05
        p = np.linspace(0.00001, 0.5, 500)
        # 标准 12 碱基 UMI 的错误概率：至少有一个碱基出错
        P_std_umi_error = 1 - (1 - p) ** num_nt
        # Homotrimer 的单个三联体错误率（投票失败的概率）
        # P_block_error = 3 * (p ** 2) * (1 - p) + p ** 3
        P_block_error = (7/3) * p**2 * (1 - p) + p**3
        P_block_error_upper_bound = 3 * p**2 - 2 * p**3
        # Homotrimer UMI 的错误率：num_nt 个三联体中至少一个出错
        P_homo_umi_error = 1 - (1 - P_block_error) ** num_nt
        P_homo_umi_error_upper_bound = 1 - (1 - P_block_error_upper_bound) ** num_nt
        # 绘图
        ax.grid(False)
        ax.plot(
            p,
            P_std_umi_error,
            label='Standard 12-bp UMI',
            color="slategrey",
            linewidth=2.5,
        )
        ax.plot(
            p,
            P_homo_umi_error,
            label='Homotrimer 36-bp UMI',
            color="crimson",
            linewidth=2.5,
        )
        ax.plot(
            p,
            P_homo_umi_error_upper_bound,
            label='36 binary repetition codes (n=3)',
            # label='36-size binary codes with the 3-input majority gate',
            color="olivedrab",
            linewidth=2.5,
        )
        ax.spines['right'].set_color('none')
        ax.spines['top'].set_color('none')
        ax.tick_params(axis='both', which='major', labelsize=16)
        ax.tick_params(axis='both', which='minor', labelsize=16)
        ax.set_xlabel("Single-base error rate ($p$)", fontsize=16)
        ax.set_ylabel("Probability of incorrect synthesis", fontsize=16)
        # ax.set_title("Comparison of UMI error probability vs. per-base error rate", fontsize=14)
        ax.legend(fontsize=14, frameon=False)
        plt.tight_layout()
        plt.show()
        return

    def correct_block(
            self,
            num_nt=1,
    ):
        palette = sns.color_palette("Set2")[3:]

        fig, ax = plt.subplots(nrows=1, ncols=1, figsize=(7, 5), sharey=True, sharex=False)

        # Define range of p values
        p = np.linspace(0.00001, 0.5, 500)

        # Standard UMI: probability all 12 bases are correct
        P_std_umi_correct = (1 - p) ** num_nt

        # Homotrimer UMI: block error probability
        # P_block_error = 3 * (p ** 2) * (1 - p) + p ** 3
        P_block_error = (7 / 3) * p ** 2 * (1 - p) + p ** 3
        P_block_error_upper_bound = 3 * p**2 - 2 * p**3
        # Homotrimer UMI: probability all 12 blocks are correct after majority voting
        P_homo_umi_correct = (1 - P_block_error) ** num_nt
        P_homo_umi_correct_upper_bound = (1 - P_block_error_upper_bound) ** num_nt

        ax.grid(False)
        ax.plot(
            p,
            P_std_umi_correct,
            label='A single base',
            color="slategrey",
            linewidth=2.5,
        )
        ax.plot(
            p,
            P_homo_umi_correct,
            label='Homotrimer block',
            color="crimson",
            linewidth=2.5,
        )
        ax.plot(
            p,
            P_homo_umi_correct_upper_bound,
            label='Binary repetition code block (n=3)',
            # label='3-input majority gate',
            color="olivedrab",
            linewidth=2.5,
        )
        ax.spines['right'].set_color('none')
        ax.spines['top'].set_color('none')
        ax.tick_params(axis='both', which='major', labelsize=16)
        ax.tick_params(axis='both', which='minor', labelsize=16)
        ax.set_xlabel("Single-base error rate ($p$)", fontsize=16)
        ax.set_ylabel("Probability of correct synthesis", fontsize=16)
        # ax.set_title("building block", fontsize=14)
        ax.legend(fontsize=14, frameon=False)
        plt.tight_layout()
        plt.show()
        return


if __name__ == "__main__":
    p = TripletErrorCode()

    # p.correct_all(
    #     num_nt=1,
    # )
    # p.incorrect_all(
    #     num_nt=1,
    # )
    p.incorrect_block()
    p.correct_block()