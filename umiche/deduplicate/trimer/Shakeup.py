__version__ = "v1.0"
__copyright__ = "Copyright 2025"
__license__ = "GPL-3.0"
__developer__ = "Jianfeng Sun"
__maintainer__ = "Jianfeng Sun"
__email__="jianfeng.sunmt@gmail.com"


import time
import pandas as pd
from umiche.fastq.Reader import Reader as rfastq
from umiche.fastq.Writer import Writer as wfastq

from umiche.trim.Reader import Reader as trimreader
from umiche.simu.Parameter import Parameter as params

from umiche.util.Writer import Writer as fwriter
from umiche.util.Console import Console


class Shakeup:

    def __init__(
            self,
            scenario,
            umi_len=36,

            param_fpn=None,
            is_trim=False,
            is_split=False,

            verbose=False,

            **kwargs,
    ):
        self.scenario = scenario
        self.umi_len = umi_len
        self.is_trim = is_trim
        self.is_split = is_split
        self.kwargs = kwargs

        self.params = params(param_fpn=param_fpn)

        self.rfastq = rfastq()
        self.wfastq = wfastq()
        self.trimreader = trimreader()
        self.fwriter = fwriter()

        self.verbose = verbose
        self.console = Console()
        self.console.verbose = self.verbose


        for perm_num_i in range(self.params.fixed['permutation_num']):
            for id, scenario_i in enumerate(self.params.varied[self.scenario]):
                read_stime = time.time()
                print('=========>permutation {} under scenario No.{} {}'.format(perm_num_i, id, scenario_i))

                if self.scenario == 'pcr_nums':
                    self.fn_mark = str(id) # scenario_i
                elif self.scenario == 'pcr_errs':
                    self.fn_mark = str(id)
                elif self.scenario == 'seq_errs':
                    self.fn_mark = str(id)
                elif self.scenario == 'ampl_rates':
                    self.fn_mark = str(id)
                elif self.scenario == 'seq_deps':
                    self.fn_mark = str(id)
                elif self.scenario == 'umi_lens':
                    self.fn_mark = str(id)
                elif self.scenario == 'umi_nums':
                    self.fn_mark = str(id)
                else:
                    self.fn_mark = str(scenario_i)
                self.fn_prefix = self.params.file_names[self.scenario] + self.fn_mark
                self.fastq_location = self.params.work_dir + self.scenario + '/permute_' + str(perm_num_i) + '/'

                print(self.fn_prefix)
                print(self.fastq_location)

                if self.is_trim:
                    from umiche.trim.Template import Template as trimmer
                    if self.scenario == 'umi_lens':
                        self.params.trimmed['umi_1']['len'] = scenario_i
                    self.console.print("======>fastq is being trimmed.")
                    self.params.trimmed['fastq']['fpn'] = self.fastq_location + self.fn_prefix + '.fastq.gz'
                    self.params.trimmed['fastq']['trimmed_fpn'] = self.fastq_location + 'trimmed/' + self.fn_prefix + '.fastq.gz'
                    umitrim_parser = trimmer(params=self.params.trimmed, verbose=self.verbose)
                    df = umitrim_parser.todf()
                    umitrim_parser.togz(df)
                if self.is_split:
                    names, seqs, _, _ = self.rfastq.fromgz(
                        fastq_fpn=self.fastq_location + 'trimmed/' + self.fn_prefix + '.fastq.gz'
                    )
                    # print(seqs)
                    # print(names)

                    # df = self.trimreader.todf(names=names, seqs=seqs)
                    df = self.trimreader.todfFromTree(names=names, seqs=seqs)

                    # df['origin_info'] = df['name'].apply(lambda x: x.split('_')[0])
                    df['origin_info'] = df.apply(lambda x: str(x['umi#']) + '_' + x['umi_src'], axis=1)
                    mono_corr_stime = time.time()

                    df['umi_order'] = df['umi'].apply(lambda x: self.split_by_order(x, umi_len=self.umi_len))
                    df['umi_l_order'] = df['umi_order'].apply(lambda x: x.split(';')[0])
                    df['umi_m_order'] = df['umi_order'].apply(lambda x: x.split(';')[1])
                    df['umi_r_order'] = df['umi_order'].apply(lambda x: x.split(';')[2])
                    # print(df['umi_r_order'])
                    df['to_fas_l_order'] = df.apply(lambda x: x['origin_info'] + '_' + x['umi_l_order'], axis=1)
                    df['to_fas_m_order'] = df.apply(lambda x: x['origin_info'] + '_' + x['umi_m_order'], axis=1)
                    df['to_fas_r_order'] = df.apply(lambda x: x['origin_info'] + '_' + x['umi_r_order'], axis=1)
                    df_lmr = pd.concat([
                        pd.DataFrame({'seq_raw': df['seq_raw'].values, 'to_fas': df['to_fas_l_order'].values}),
                        pd.DataFrame({'seq_raw': df['seq_raw'].values, 'to_fas': df['to_fas_m_order'].values}),
                        pd.DataFrame({'seq_raw': df['seq_raw'].values, 'to_fas': df['to_fas_r_order'].values}),
                    ], axis=0).reset_index(drop=True)
                    # print(df_lmr)

                    df['umi_mv'] = df['umi'].apply(lambda x: self.split_by_mv(x, umi_len=self.umi_len))
                    df['umi_ref'] = df['umi_mv'].apply(lambda x: x.split(';')[0])
                    df['umi_l_mv'] = df['umi_mv'].apply(lambda x: x.split(';')[1])
                    df['umi_m_mv'] = df['umi_mv'].apply(lambda x: x.split(';')[2])
                    df['umi_r_mv'] = df['umi_mv'].apply(lambda x: x.split(';')[3])

                    df['to_fas'] = df.apply(lambda x: x['origin_info'] + '_' + x['umi_ref'], axis=1)
                    # print(df['to_fas'])
                    df['umi_mark'] = df['umi'].apply(lambda x: self.marker(x, umi_len=self.umi_len))
                    df_differ_3nt = df.loc[df['umi_mark'] == 1]


                    df_3differ_cp = df_differ_3bases.copy()

                    df_3differ = df_3differ.drop('umi_r_mv', 1)
                    df_3differ_cp = df_3differ_cp.drop('umi_m_mv', 1)
                    df_3differ = df_3differ.rename(columns={"umi_m_mv": "umi_bi"})
                    df_3differ_cp = df_3differ_cp.rename(columns={"umi_r_mv": "umi_bi"})

                    df_umi_3differ = pd.concat(
                        [df_3differ, df_3differ_cp],
                        axis=0,
                    ).reset_index(drop=True)
                    # print(df_umi_3differ)

                    if not df_umi_3differ.empty:
                        df_umi_3differ['to_fas'] = df_umi_3differ.apply(lambda x: x['origin_info'] + '_' + x['umi_bi'], axis=1)
                        print(df_umi_3differ)
                    # print(df_umi_3differ['umi_bi'])

                    df_not3differ = df.loc[df['umi_mark'] != 1]
                    df_not3differ['to_fas'] = df_not3differ.apply(lambda x: x['origin_info'] + '_' + x['umi_l_mv'], axis=1)
                    if df_umi_3differ.empty:
                        df_merge = df_not3differ[['seq_raw', 'to_fas']].reset_index(drop=True)
                    else:
                        df_merge = pd.concat(
                            [df_umi_3differ[['seq_raw', 'to_fas']], df_not3differ[['seq_raw', 'to_fas']]],
                            axis=0,
                        ).reset_index(drop=True)
                    print(df_merge)
                    self.wfastq.togz(
                        list_2d=df_lmr[['seq_raw', 'to_fas']].values,
                        sv_fpn=self.fastq_fp + 'seq_errs/permute_' + str(perm_num_i) + '/lmr/',
                        fn=self.cat + '_' + str(id),
                    )
                    self.wfastq.togz(
                        list_2d=df[['seq_raw', 'to_fas']].values,
                        sv_fp=self.fastq_fp + 'seq_errs/permute_' + str(perm_num_i) + '/ref/',
                        fn=self.cat + '_' + str(id),
                    )
                    self.wfastq.togz(
                        list_2d=df_merge.values,
                        sv_fp=self.fastq_fp + 'seq_errs/permute_' + str(perm_num_i) + '/bipartite/',
                        fn=self.cat + '_' + str(id),
                    )
                    print('===>getting it done with time: {:.3f}s'.format(time.time() - mono_corr_stime))

    def split_by_mv(
            self,
            umi,
            umi_len=36,
    ):
        """

        Parameters
        ----------
        umi

        umi_len

        Returns
        -------

        """
        vernier = [i for i in range(umi_len) if i % 3 == 0]
        umi_trimers = [umi[v: v+3] for v in vernier]
        ref = []
        l = []
        m = []
        r = []
        for umi_trimer in umi_trimers:
            ref.append(umi_trimer[0])
            s = set(umi_trimer)
            if len(s) == 3:
                l.append(umi_trimer[0])
                m.append(umi_trimer[1])
                r.append(umi_trimer[2])
            elif len(s) == 2:
                sdict = {umi_trimer.count(i): i for i in s}
                l.append(sdict[2])
                m.append(sdict[2])
                r.append(sdict[2])
            else:
                l.append(umi_trimer[0])
                m.append(umi_trimer[0])
                r.append(umi_trimer[0])
        ref = ''.join(ref)
        l = ''.join(l)
        m = ''.join(m)
        r = ''.join(r)
        return ref + ';' + l + ';' + m + ';' + r

    def split_by_order(
            self,
            umi,
            umi_len=36,
    ):
        vernier = [i for i in range(umi_len) if i % 3 == 0]
        umi_trimers = [umi[v: v + 3] for v in vernier]
        l = []
        m = []
        r = []
        for umi_trimer in umi_trimers:
            l.append(umi_trimer[0])
            m.append(umi_trimer[1])
            r.append(umi_trimer[2])
        l = ''.join(l)
        m = ''.join(m)
        r = ''.join(r)
        return l + ';' + m + ';' + r

    def marker(
            self,
            umi,
            umi_len=36,
    ):
        vernier = [i for i in range(umi_len) if i % 3 == 0]
        umi_trimers = [umi[v: v+3] for v in vernier]
        slens = [len(set(umi_trimer)) for umi_trimer in umi_trimers]
        if 3 in slens:
            return 1
        else:
            return 0


if __name__ == "__main__":
    from umiche.path import to

    p = Shakeup(
        # scenario='pcr_nums',
        # scenario='pcr_errs',
        scenario='seq_errs',
        # scenario='ampl_rates',
        # scenario='umi_lens',
        # scenario='seq_deps',
        # scenario='umi_nums',

        umi_len=30,

        is_trim=False, # False True
        is_split=True, # False True

        # umi_ref_fpn=to('data/simu/umi/trimer/seq_errs/umi.txt'),
        # fastq_fp=to('data/simu/umi/trimer/seq_errs/trimmed/'),

        param_fpn=to('data/trimer.yml'),

        verbose=True, # False True
    )
