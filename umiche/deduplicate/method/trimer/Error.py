__version__ = "0.0.1"
__copyright__ = "Copyright 2025"
__license__ = "GPL-3.0"
__developer__ = "Jianfeng Sun"
__maintainer__ = "Jianfeng Sun"
__email__ = "jianfeng.sunmt@gmail.com"


from umiche.util.Single import Single as dnasgl
from umiche.util.Number import Number as rannum


class Error:

    def __init__(self, ):
        pass

    def deletion(
            self,
            read,
            del_rate,
            mode='normal',
    ):
        num_err_per_read = rannum().binomial(
            n=len(read), p=del_rate, use_seed=False, seed=False
        )
        if num_err_per_read != 0:
            mark = True
        else:
            mark = False
        pos_list = rannum().choice(
            high=len(read), num=num_err_per_read, use_seed=False, seed=False, replace=False,
        )
        for _, pos in enumerate(pos_list):
            read = read[:pos] + read[pos + 1:]
        if mode == 'bead_deletion':
            return read, {
                'mark': mark,
            }
        else:
            return read

    def insertion(
            self,
            read,
            ins_rate,
            mode='normal',
    ):
        num_err_per_read = rannum().binomial(
            n=len(read), p=ins_rate, use_seed=False, seed=False
        )
        # print(num_err_per_read)
        if num_err_per_read != 0:
            mark = True
        else:
            mark = False
        pos_list = rannum().choice(
            high=len(read), num=num_err_per_read, use_seed=False, seed=False, replace=False,
        )
        base_list = rannum().uniform(
            low=0, high=4, num=num_err_per_read, use_seed=False
        )
        for i, pos in enumerate(pos_list):
            dna_map = dnasgl().todict(
                nucleotides=dnasgl().get(
                    universal=True,
                ),
                reverse=True,
            )
            read = read[:pos] + dna_map[base_list[i]] + read[pos:]
            ### read
            ### pos, base_list[i], dna_map[base_list[i]], dna_map
            ### read
            # 5 3 G {0: 'A', 1: 'T', 2: 'C', 3: 'G'}
            # TTTTTTTTTGGGCCCGGGAAAAAACCCAAAGGGGGG
            # TTTTTGTTTTGGGCCCGGGAAAAAACCCAAAGGGGGG
            # 9 0 A {0: 'A', 1: 'T', 2: 'C', 3: 'G'}
            # CCCTTTCCCTTTGGGTTTGGGTTTCCCGGGAAACCC
            # CCCTTTCCCATTTGGGTTTGGGTTTCCCGGGAAACCC
            # 3 0 A {0: 'A', 1: 'T', 2: 'C', 3: 'G'}
            # AAATTTTTTAAACCCAAAAAAAAAAAATTTTTTCCC
            # AAAATTTTTTAAACCCAAAAAAAAAAAATTTTTTCCC
        if mode == 'bead_insertion':
            return read, {
                'mark': mark,
            }
        else:
            return read

    def mutated(
            self,
            read,
            mut_rate,
            mode='normal',
    ):
        num_err_per_read = rannum().binomial(
            n=len(read), p=mut_rate, use_seed=False, seed=False
        )
        if num_err_per_read != 0:
            mark = True
        else:
            mark = False
        # pos_list = rannum().uniform(
        #     low=0, high=len(read), num=num_err_per_read, use_seed=False, seed=False
        # )
        pos_list = rannum().choice(
            high=len(read), num=num_err_per_read, use_seed=False, seed=False, replace=False,
        )
        base_list = rannum().uniform(
            low=0, high=3, num=num_err_per_read, use_seed=False
        )
        read_l = list(read)
        dna_map = dnasgl().get(
            universal=True,
        )
        # print(dna_map)
        for i, pos in enumerate(pos_list):
            # print(read_l)
            base = read_l[pos]
            if base not in dna_map:
                base = 'A'
            dna_map = dnasgl().todict(
                nucleotides=dnasgl().getEleTrimmed(
                    ele_loo=base,
                    universal=True,
                ),
                reverse=True,
            )
            ### dna_map, base_list[i], dna_map[base_list[i]]
            # {0: 'A', 1: 'C', 2: 'G'} 0 A
            # {0: 'T', 1: 'C', 2: 'G'} 1 C
            # {0: 'A', 1: 'T', 2: 'C'} 1 T
            # {0: 'A', 1: 'T', 2: 'C'} 0 A
            read_l[pos] = dna_map[base_list[i]]
        if mode == 'bead_mutation':
            return ''.join(read_l), {
                'mark': mark,
            }
        else:
            return ''.join(read_l)

    def change(
            self,
            read,
            pos_list,
            base_list,
    ):
        read_l = list(read)
        for i, pos in enumerate(pos_list):
            dna_map = dnasgl().todict(
                nucleotides=dnasgl().getEleTrimmed(
                    ele_loo=read_l[pos],
                    universal=True,
                ),
                reverse=True,
            )
            read_l[pos] = dna_map[base_list[i]]
        return ''.join(read_l)


if __name__ == "__main__":
    p = Error()
    mut_umi = p.mutated("ACGTACGT", mut_rate=0.5, mode="normal")
    print(mut_umi)

    import pandas as pd
    df = pd.read_csv(
        '/mnt/d/Document/Programming/Python/umiche/umiche/data/r1/umicountr/df_res.txt',
        sep='\t',
        header=0,
    )

    from umiche.deduplicate.method.trimer.Expand import Expand
    df["htUMI"] = df["UX"].apply(Expand().homotrimer_umi)
    mut_rate = 0.5
    df["htUMI_" + str(mut_rate)] = df["htUMI"].apply(lambda umi: p.mutated(umi, mut_rate=mut_rate, mode="normal"))
    print(df)