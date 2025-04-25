__version__ = "v1.0"
__copyright__ = "Copyright 2025"
__license__ = "GPL-3.0"
__developer__ = "Jianfeng Sun"
__maintainer__ = "Jianfeng Sun"
__email__="jianfeng.sunmt@gmail.com"


import pandas as pd
from umiche.fastq.Reader import Reader as rfastq
from umiche.fastq.Writer import Writer as wfastq
from umiche.trim.Filter import Filter
from umiche.trim.BCRuleOut import bcRuleOut as bcro
from umiche.trim.UMIRuleOut import umiRuleOut as umiro
from umiche.trim.SeqRuleOut import seqRuleOut as seqro
from umiche.util.Console import Console


class Template:

    def __init__(
            self,
            params,
            verbose=True,
    ):
        self.params = params
        self.verbose = verbose
        self.filter = Filter()
        self.bcro = bcro(read_summary=self.params, verbose=self.verbose)
        self.umiro = umiro(read_summary=self.params, verbose=self.verbose)
        self.seqro = seqro(read_summary=self.params, verbose=self.verbose)
        self.rfastq = rfastq
        self.wfastq = wfastq

        self.console = Console()
        self.console.verbose = self.verbose

    def todf(self, ):
        """
        Notes
        -----
        Trimmed UMIs and genomic sequences.

        Examples
        --------
        >>> params = {
        ...  'read_struct': 'primer_1+umi_1+seq_1+seq_2+umi_2+primer_2',
        ...  'umi_1': {'len': 12},
        ...  'umi_2': {'len': 10},
        ...  'primer_1': {'len': 20},
        ...  'primer_2': {'len': 20},
        ...  'seq_1': {'len': 6},
        ...  'seq_2': {'len': 8},
        ...  'fastq': {
        ...      'path': to('example/data/'),
        ...      'name': 'pcr_1',
        ...      'trimmed_path': to('example/data/'),
        ...      'trimmed_name': 'pcr_1_trim',
        ...   },
        ... }
        >>> p = template(params)
        >>> p.todf()

        Returns
        -------
            pandas dataframe with column names: seq_raw, name, umi_1, umi_2, ..., umi_n, seq_2, seq_3, ..., seq_n

        """
        self.console.print('===>reading from fastq...')
        names, seqs, _, _ = self.rfastq().fromgz(
            fastq_fpn=self.params['fastq']['fpn'],
        )
        df = pd.DataFrame(seqs, columns=['seq_raw'])
        df['name'] = names
        print('before trimmed', df)
        self.console.print('===>umi structure: {}'.format(self.params['read_struct']))
        compo_struct = self.params['read_struct'].split('+')
        bc_pos_in_struct = [i for i, cstruct in enumerate(compo_struct) if 'bc' in cstruct.split('_')]
        # print(bc_pos_in_struct)
        umi_pos_in_struct = [i for i, cstruct in enumerate(compo_struct) if 'umi' in cstruct.split('_')]
        seq_pos_in_struct = [i for i, cstruct in enumerate(compo_struct) if 'seq' in cstruct.split('_')]
        self.console.print('===>bc positions in the read structure: {}'.format(', '.join(map(str, bc_pos_in_struct))))
        self.console.print('===>umi positions in the read structure: {}'.format(', '.join(map(str, umi_pos_in_struct))))
        self.console.print('===>seq positions in the read structure: {}'.format(', '.join(map(str, seq_pos_in_struct))))
        bc_rule_out_len_dict = self.bcro.sequential(
            compo_struct=compo_struct,
            bc_pos_in_struct=bc_pos_in_struct,
        )
        # print(bc_rule_out_len_dict)
        umi_rule_out_len_dict = self.umiro.sequential(
            compo_struct=compo_struct,
            umi_pos_in_struct=umi_pos_in_struct,
        )
        # print(df)
        seq_rule_out_len_dict = self.seqro.sequential(
            compo_struct=compo_struct,
            seq_pos_in_struct=seq_pos_in_struct,
        )
        # print(df)
        for key, val in bc_rule_out_len_dict.items():
            start = val
            end = val + self.params[key]['len']
            print('BC: start: {} end: {}'.format(start, end))
            df[key] = df.apply(lambda x: self.filter.singleStart(x['seq_raw'], start, end), axis=1)
            self.console.print('===>{} has been taken out'.format(key))
        # print(df)
        for key, val in umi_rule_out_len_dict.items():
            start = val
            end = val + self.params[key]['len']
            print('UMI: start: {} end: {}'.format(start, end))
            df[key] = df.apply(lambda x: self.filter.singleStart(x['seq_raw'], start, end), axis=1)
            self.console.print('===>{} has been taken out'.format(key))
        # print(df)
        for key, val in seq_rule_out_len_dict.items():
            start = val
            end = val + self.params[key]['len']
            df[key] = df.apply(lambda x: self.filter.singleStart(x['seq_raw'], start, end), axis=1)
            self.console.print('===>{} has been taken out'.format(key))
        print('after trimmed', df)
        return df

    def togz(self, df):
        """

        Note
        ----
            write to a Fastq file in gz format.

        Example
        -------
        >>>params = {
        ... 'read_struct': 'primer_1+umi_1+seq_1+seq_2+umi_2+primer_2',
        ... 'umi_1': {'len': 12},
        ... 'umi_2': {'len': 10},
        ... 'primer_1': {'len': 20},
        ... 'primer_2': {'len': 20},
        ... 'seq_1': {'len': 6},
        ... 'seq_2': {'len': 8},
        ... 'fastq': {
        ...     'path': to('example/data/'),
        ...     'name': 'pcr_1',
        ...     'trimmed_path': to('example/data/'),
        ...     'trimmed_name': 'pcr_1_trim',
        ... },
        ...}
        >>>p = Template(params)
        >>>df = p.todf()
        >>>p.togz(df)

        Parameters
        ----------
        df
            pandas Dataframe returned from todf().

        Returns
        -------
            pandas dataframe in the column orders: seq_1, name, umi_1, umi_2, ..., umi_n, seq_2, seq_3, ..., seq_n

        """
        self.console.print('===>start saving in gz format...')
        # # /* block. in case of only UMIs present in the reads. */
        if 'seq_1' not in df.columns.tolist():
            df['seq_1'] = 'B'
        df_col_names = df.columns.tolist()
        self.console.print(df_col_names)
        df_col_names.remove('seq_raw')
        rv_todo_list = [struct for i, struct in enumerate(df_col_names) if 'seq' in struct.split('_')]
        for i in rv_todo_list:
            df_col_names.remove(i)
        if len(rv_todo_list) != 1:
            df['seq'] = ''
            for i in rv_todo_list:
                df['seq'] = df.apply(lambda x: x['seq'] + x[i], axis=1)
        else:
            df['seq'] = df['seq_1']
        df_trimmed = df[['seq'] + df_col_names]
        self.wfastq().togz(
            list_2d=df_trimmed.values.tolist(),
            sv_fpn=self.params['fastq']['trimmed_fpn'],
        )
        self.console.print('===>trimmed UMIs have been saved in gz format.')
        # return df_trimmed


if __name__ == "__main__":
    from umiche.path import to

    params = {
        'umi_1': {
            'len': 12,
        },
        'umi_2': {
            'len': 4,
        },
        'umi_3': {
            'len': 5,
        },
        'bc_1': {
            'len': 2,
        },
        'read_struct': 'umi_1',
        # 'read_struct': 'umi_1+seq_1',
        # 'read_struct': 'bc_1+umi_1+seq_1',
        # 'read_struct': 'primer_1+umi_1+seq_1+seq_2+umi_2+umi_3',
        # 'read_struct': 'primer+umi_start+seq+umi_1+primer',
        # 'read_struct': 'primer+umi_1+seq_1+primer+umi_2+seq_2+umi_3+primer',
        'primer_1': {
            'len': 20,
        },
        'primer_2': {
            'len': 20,
        },
        'seq_1': {
            'len': 6,
        },
        'seq_2': {
            'len': 8,
        },
        'fastq': {
            'fpn': to('data/simu/mclumi/seq_errs/permute_0/seq_err_5.fastq.gz'),
            'trimmed_fpn': to('data/simu/mclumi/seq_errs/permute_0/seq_err_5_trimmed.fastq.gz'),
        },
    }
    p = Template(params)

    df = p.todf()

    print(df)
    print(p.togz(df))