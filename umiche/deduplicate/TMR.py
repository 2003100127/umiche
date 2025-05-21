__version__ = "v1.0"
__copyright__ = "Copyright 2025"
__license__ = "GPL-3.0"
__developer__ = "Jianfeng Sun"
__maintainer__ = "Jianfeng Sun"
__email__ = "jianfeng.sunmt@gmail.com"

from typing import Union
import numpy as np

from umiche.util.Console import Console


class TMR:
    
    def __init__(
            self,
            error_rate: Union[float, np.ndarray] = 0.00001,
            num_nt: int = 12,
            verbose=True,
    ):
        self.error_rate = error_rate
        self._num_nt = num_nt
        # self.bitflip_block_error = 3 * self.error_rate**2 - 2 * self.error_rate**3
        # self.homotrimer_block_error = (7/3) * self.error_rate**2 * (1 - self.error_rate) + self.error_rate**3

        self.verbose = verbose
        self.console = Console()
        self.console.verbose = self.verbose

    def get_num_nt(self):
        self.console.print("=========>number of building blocks in a UMI is {}".format(self._num_nt))
        return self._num_nt

    def set_num_nt(self, value):
        # self.console.print("=========>number of building blocks in a UMI is set to {}:".format(value))
        self._num_nt = value

    @property
    def homotrimer_block_error(self, ):
        e = (7/3) * self.error_rate**2 * (1 - self.error_rate) + self.error_rate**3
        self.console.print("=========>homotrimer block error rate: {}".format(e))
        return e

    @property
    def bitflip_block_error(self, ):
        e = 3 * self.error_rate ** 2 - 2 * self.error_rate ** 3
        self.console.print("=========>binary repetition code block error rate is {}".format(e))
        return e

    @property
    def homotrimer_block_error_free(self, ):
        e = (7 / 3) * self.error_rate ** 2 * (1 - self.error_rate) + self.error_rate ** 3
        self.console.print("=========>homotrimer block error rate: {}".format(e))
        return 1 - e

    @property
    def bitflip_block_error_free(self, ):
        e = 3 * self.error_rate ** 2 - 2 * self.error_rate ** 3
        self.console.print("=========>binary repetition code block error rate is {}".format(e))
        return 1 - e

    @property
    def monomer_umi_error(self, ):
        self.console.print("======>monomer UMI error rate:")
        return 1 - (1 - self.error_rate) ** self.get_num_nt()

    @property
    def homotrimer_umi_error(self, ):
        self.console.print("======>homotrimer UMI error rate:")
        return 1 - (1 - self.homotrimer_block_error) ** self.get_num_nt()

    @property
    def bitflip_code_error(self, ):
        self.console.print("======>binary repetition code error rate:")
        return 1 - (1 - self.bitflip_block_error) ** self.get_num_nt()

    @property
    def monomer_umi_error_free(self, ):
        self.console.print("======>monomer UMI successful synthesis rate:")
        return (1 - self.error_rate) ** self.get_num_nt()

    @property
    def homotrimer_umi_error_free(self, ):
        self.console.print("======>homotrimer UMI successful synthesis rate:")
        return (1 - self.homotrimer_block_error) ** self.get_num_nt()

    @property
    def bitflip_code_error_free(self, ):
        self.console.print("======>binary repetition code successful transmission rate:")
        return (1 - self.bitflip_block_error) ** self.get_num_nt()