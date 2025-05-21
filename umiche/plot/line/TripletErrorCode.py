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
from umiche.deduplicate.TMR import TMR


class TripletErrorCode:

    def __init__(self, error_rate):
        self.error_rate = error_rate
        self.tmr = TMR(error_rate=self.error_rate)

        sns.set(font="Helvetica")
        sns.set_style("ticks")

    def incorrect(
            self,
            num_nt=12,
    ):
        fig, ax = plt.subplots(nrows=1, ncols=1, figsize=(7, 5), sharey=True, sharex=False)
        self.tmr.set_num_nt(num_nt)
        P_std_umi_error = self.tmr.monomer_umi_error
        P_homo_umi_error = self.tmr.homotrimer_umi_error
        P_homo_umi_error_upper_bound = self.tmr.bitflip_code_error
        ax.grid(False)
        ax.plot(
            self.error_rate,
            P_std_umi_error,
            label='Standard 12-bp UMI',
            color="slategrey",
            linewidth=2.5,
        )
        ax.plot(
            self.error_rate,
            P_homo_umi_error,
            label='Homotrimer 36-bp UMI',
            color="crimson",
            linewidth=2.5,
        )
        ax.plot(
            self.error_rate,
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

    def correct(
            self,
            num_nt=1,
    ):
        fig, ax = plt.subplots(nrows=1, ncols=1, figsize=(7, 5), sharey=True, sharex=False)
        self.tmr.set_num_nt(num_nt)
        P_std_umi_correct = self.tmr.monomer_umi_error_free
        P_homo_umi_correct = self.tmr.homotrimer_umi_error_free
        P_homo_umi_correct_upper_bound = self.tmr.bitflip_code_error_free
        ax.grid(False)
        ax.plot(
            self.error_rate,
            P_std_umi_correct,
            label='A single base',
            color="slategrey",
            linewidth=2.5,
        )
        ax.plot(
            self.error_rate,
            P_homo_umi_correct,
            label='Homotrimer block',
            color="crimson",
            linewidth=2.5,
        )
        ax.plot(
            self.error_rate,
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
    e = np.linspace(0.00001, 0.5, 500)
    p = TripletErrorCode(error_rate=e)
    # p.incorrect(num_nt=12)
    p.correct(num_nt=1)