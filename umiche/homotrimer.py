__version__ = "v1.0"
__copyright__ = "Copyright 2025"
__license__ = "GPL-3.0"
__developer__ = "Jianfeng Sun"
__maintainer__ = "Jianfeng Sun"
__email__="jianfeng.sunmt@gmail.com"


from typing import List, Dict

from umiche.deduplicate.method.trimer.Collapse import Collapse



def collapse():
    return Collapse()


if __name__ == "__main__":
    splitter = collapse()
    # print(splitter.split_by_mv('AAACCCGGGTTTGGGAAATTTGGGCCCCCC', recur_len=3))

    # print(splitter.split_by_mv('AAATCCGGATTTCGGAAATTTGGGCACCCC', recur_len=3))
    # print(splitter.split_by_mv('AAATCCGGATTTCGGAAATTTGGGCCCCCC', recur_len=3))
    # print(splitter.split_by_mv('AAATCCGGATTTGGGAAATTTGGGCCCCCC', recur_len=3))
    # print(splitter.split_by_mv('AAATCCGGGTTTGGGAAATTTGGGCCCCCC', recur_len=3))

    # print(splitter.split_to_all('AAATCCGGATTTCGGAAATTTGGGCACCCC', recur_len=3))
    # print(splitter.split_to_all('AAATCCGGATTTCGGAAATTTGGGCCCCCC', recur_len=3))
    # print(splitter.split_to_all('AAATCCGGATTTGGGAAATTTGGGCCCCCC', recur_len=3))
    # print(splitter.split_to_all('AAATCCGGGTTTGGGAAATTTGGGCCCCCC', recur_len=3))

    # print(splitter.majority_vote('AAATCCGGGTTTGGGAAATTTGGGCCCCCC', recur_len=3))
    # print(splitter.take_by_order('AAATCCGGGTTTGGGAAATTTGGGCCCCCC', recur_len=3))

    print(splitter.vote('AAA', recur_len=3))
    print(splitter.vote('TAA', recur_len=3))
    print(splitter.vote('TGA', recur_len=3))