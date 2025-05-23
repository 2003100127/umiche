__version__ = "v1.0"
__copyright__ = "Copyright 2025"
__license__ = "GPL-3.0"
__developer__ = "Jianfeng Sun"
__maintainer__ = "Jianfeng Sun"
__email__="jianfeng.sunmt@gmail.com"

import pandas as pd
from umiche.fastq.Convert import Convert as fastqconverter
from umiche.trim.Template import Template as trimmer
from umiche.network.bfs.ConnectedComponent import ConnectedComponent as gbfscc
from umiche.simu.Parameter import Parameter as params

from umiche.deduplicate.MultiPos import MultiPos as deduppos
from umiche.util.Writer import Writer as fwriter
from umiche.util.Console import Console


class Standard:

    def __init__(
            self,
            scenario,
            method,
            param_fpn=None,
            is_trim=False,
            is_tobam=False,
            is_dedup=False,
            is_sv=False,
            verbose=False,
            **kwargs,
    ):
        self.scenario = scenario
        self.method = method
        self.is_sv = is_sv
        self.kwargs = kwargs
        if 'is_collapse_block' in self.kwargs.keys():
            self.is_collapse_block = self.kwargs['is_collapse_block']
        else:
            self.is_collapse_block = False
        if 'deduped_method' in self.kwargs.keys():
            self.deduped_method = self.kwargs['deduped_method']
        else:
            self.deduped_method = ''
        if 'split_method' in self.kwargs.keys():
            self.split_method = self.kwargs['split_method']
        else:
            self.split_method = ''
        if 'collapse_block_method' in self.kwargs.keys():
            self.collapse_block_method = self.kwargs['collapse_block_method']
        else:
            self.collapse_block_method = ''

        self.params = params(param_fpn=param_fpn)
        self.gbfscc = gbfscc()
        self.fwriter = fwriter()
        print('UMI homopolymer recurring pattern: {}'.format(self.params.fixed['umi_unit_pattern']))

        self.verbose = verbose
        self.console = Console()
        self.console.verbose = self.verbose

        df_dedup = pd.DataFrame()
        df_slv = pd.DataFrame()
        df_not_slv = pd.DataFrame()
        df_mono_len = pd.DataFrame()
        df_multi_len = pd.DataFrame()
        for perm_num_i in range(self.params.fixed['permutation_num']):
            print("===>Permutation number: {}".format(perm_num_i))
            dedup_arr = []
            slv_arr = []
            not_slv_arr = []
            mono_len_arr = []
            multi_len_arr = []
            for id, scenario_i in enumerate(self.params.varied[self.scenario]):
                self.console.print("======>No.{} scenario: {}".format(id+1, scenario_i))
                if self.scenario == 'pcr_nums':
                    self.fn_mark = str(scenario_i)
                elif self.scenario == 'pcr_errs':
                    self.fn_mark = str(id)
                elif self.scenario == 'seq_errs':
                    self.fn_mark = str(id)
                elif self.scenario == 'ampl_rates':
                    self.fn_mark = str(id)
                elif self.scenario == 'umi_lens':
                    self.fn_mark = str(scenario_i)
                else:
                    self.fn_mark = str(scenario_i)
                self.fn_prefix = self.params.file_names[self.scenario] + self.fn_mark
                self.fastq_location = self.params.work_dir + self.scenario + '/permute_' + str(perm_num_i) + '/'
                if is_trim:
                    self.console.print("======>fastq is being trimmed.")
                    self.params.trimmed['fastq']['fpn'] = self.fastq_location + self.fn_prefix + '.fastq.gz'
                    self.params.trimmed['fastq']['trimmed_fpn'] = self.fastq_location + 'trimmed/' + self.fn_prefix + '.fastq.gz'
                    umitrim_parser = trimmer(params=self.params.trimmed, verbose=self.verbose)
                    df = umitrim_parser.todf()
                    umitrim_parser.togz(df)
                if is_tobam:
                    self.console.print("======>fastq converter to bam is being used.")
                    fastqconverter(
                        fastq_fpn=self.fastq_location + 'trimmed/' + self.fn_prefix + '.fastq.gz',
                        bam_fpn=self.fastq_location + 'trimmed/bam/' + self.fn_prefix + '.bam',
                    ).tobam()
                if is_dedup:
                    self.console.print("======>reads are being deduplicated.")
                    if self.is_collapse_block:
                        if self.deduped_method == '' and self.split_method == '':
                            plus = ''
                        else:
                            plus = self.deduped_method + '_' + self.split_method
                        bam_fpn = self.fastq_location + 'trimmed/bam/' + plus + '/' + self.fn_prefix + '.bam'
                    else:
                        bam_fpn = self.fastq_location + 'trimmed/bam/' + self.fn_prefix + '.bam'
                    if self.method == 'set_cover' or self.method == 'majority_vote':
                        self.is_build_graph = False
                    else:
                        self.is_build_graph = True
                    # print(self.is_collapse_block)
                    dedup_ob = deduppos(
                        bam_fpn=bam_fpn,
                        pos_tag='PO',
                        mcl_fold_thres=self.params.dedup['mcl_fold_thres'],
                        inflat_val=self.params.dedup['inflat_val'],
                        exp_val=self.params.dedup['exp_val'],
                        iter_num=self.params.dedup['iter_num'],
                        ed_thres=self.params.dedup['ed_thres'],
                        work_dir=self.params.work_dir,
                        sv_interm_bam_fpn=self.fastq_location + 'trimmed/bam/' + self.method + '_' + self.split_method + '/' + self.fn_prefix + '.bam',
                        heterogeneity=False, # False True
                        is_build_graph=self.is_build_graph,
                        is_collapse_block=self.is_collapse_block,
                        umi_unit_pattern=self.params.fixed['umi_unit_pattern'],
                        split_method=self.split_method,
                        collapse_block_method=self.collapse_block_method,
                        verbose=self.verbose,
                        # **self.kwargs
                    )
                    df = self.tool(dedup_ob)[self.method]()
                    if self.method == "unique":
                        print("============>No.{}, dedup cnt: {}".format(id, df.num_uniq_umis.values[0]))
                        dedup_arr.append(df.num_uniq_umis.values[0])
                    elif self.method == "set_cover":
                        print("============>No.{}, dedup cnt: {}".format(id, df.dedup_cnt.values[0]))
                        # print("============>No.{}, num solved: {}".format(id, df.num_solved.values[0]))
                        # print("============>No.{}, num not solved: {}".format(id, df.num_not_solved.values[0]))
                        dedup_arr.append(df.dedup_cnt.values[0])
                        slv_arr.append(df.num_solved.values[0])
                        not_slv_arr.append(df.num_not_solved.values[0])
                        mono_len_arr.append(df.monomer_umi_len.values[0])
                        multi_len_arr.append(df.multimer_umi_len.values[0])
                    else:
                        print("============>No.{}, dedup cnt: {}".format(id, df.dedup_cnt.values[0]))
                        dedup_arr.append(df.dedup_cnt.values[0])
            df_dedup['pn' + str(perm_num_i)] = dedup_arr
            if self.method == 'set_cover':
                df_slv['pn' + str(perm_num_i)] = slv_arr
                df_not_slv['pn' + str(perm_num_i)] = not_slv_arr
                df_mono_len['pn' + str(perm_num_i)] = mono_len_arr
                df_multi_len['pn' + str(perm_num_i)] = multi_len_arr
            # print(df_dedup)
        if not self.is_collapse_block:
            remark = 'collapse_block'
        else:
            remark = 'no_collapse_block'
        self.sv_dedup_dir = self.params.work_dir + self.scenario + '/' + remark + '/'
        from umiche.util.Folder import Folder as crtfolder
        crtfolder().osmkdir(DIRECTORY=self.sv_dedup_dir)
        # self.fwriter.generic(
        #     df=df_dedup,
        #     sv_fpn=self.sv_dedup_dir + self.method + '_dedupby_' + self.deduped_method + '_splitby_' + self.split_method + '_collblockby_' + self.collapse_block_method + '_dedup.txt',
        #     header=True,
        # )
        # if self.method == 'set_cover':
            # self.fwriter.generic(
            #     df=df_slv,
            #     sv_fpn=self.sv_dedup_dir + self.method + '_solved_' + self.split_method + '_dedup.txt',
            #     header=True,
            # )
            # self.fwriter.generic(
            #     df=df_not_slv,
            #     sv_fpn=self.sv_dedup_dir + self.method + '_not_solved_' + self.split_method + '_dedup.txt',
            #     header=True,
            # )
            # self.fwriter.generic(
            #     df=df_mono_len,
            #     sv_fpn=self.sv_dedup_dir + self.method + '_mono_len_' + self.split_method + '_dedup.txt',
            #     header=True,
            # )
            # self.fwriter.generic(
            #     df=df_multi_len,
            #     sv_fpn=self.sv_dedup_dir + self.method + '_multi_len_' + self.split_method + '_dedup.txt',
            #     header=True,
            # )

    def tool(self, dedup_ob):
        return {
            'unique': dedup_ob.unique,
            'cluster': dedup_ob.cluster,
            'adjacency': dedup_ob.adjacency,
            'directional': dedup_ob.directional,
            'mcl': dedup_ob.mcl,
            'mcl_val': dedup_ob.mcl_val,
            'mcl_ed': dedup_ob.mcl_ed,
            'mcl_cc_all_node_umis': dedup_ob.mcl_cc_all_node_umis,
            # 'dbscan_seq_onehot': dedup_ob.dbscan_seq_onehot,
            # 'birch_seq_onehot': dedup_ob.birch_seq_onehot,
            # 'hdbscan_seq_onehot': dedup_ob.hdbscan_seq_onehot,
            # 'aprop_seq_onehot': dedup_ob.aprop_seq_onehot,
            'set_cover': dedup_ob.set_cover,
            'majority_vote': dedup_ob.majority_vote,
        }


if __name__ == "__main__":
    from umiche.path import to

    p = Standard(
        # scenario='pcr_nums',
        # scenario='pcr_errs',
        scenario='seq_errs',
        # scenario='ampl_rates',
        # scenario='umi_lens',

        # method='unique',
        # method='cluster',
        # method='adjacency',
        # method='directional',
        # method='mcl',
        # method='mcl_val',
        # method='mcl_ed',
        method='set_cover',
        # method='majority_vote',

        # is_trim=True,
        # is_tobam=False,
        # is_dedup=False,

        # is_trim=False,
        # is_tobam=True,
        # is_dedup=False,

        is_trim=False,
        is_tobam=False,
        is_dedup=True,

        # @@ for directional on multimer umis deduplicated by set_cover
        is_collapse_block=False,
        deduped_method='set_cover',
        split_method='split_to_all', # split_to_all split_by_mv

        # @@ for directional on multimer umis deduplicated by majority_vote
        # is_collapse_block=False,
        # deduped_method='majority_vote',
        # split_method='',

        # @@ for directional but on monomer umis of set_cover or majority_vote
        # is_collapse_block=True, # True False
        # collapse_block_method='take_by_order', # majority_vote take_by_order
        # deduped_method='set_cover', # majority_vote set_cover
        # split_method='split_by_mv', # split_to_all split_by_mv

        # @@ for directional but on monomer umis without other methods.
        # is_collapse_block=True,  # True False
        # collapse_block_method='majority_vote',  # majority_vote take_by_order
        # deduped_method='',  # majority_vote set_cover
        # split_method='',  # split_to_all split_by_mv

        # param_fpn=to('data/params_dimer.yml'),
        param_fpn=to('data/params_trimer.yml'),
        # param_fpn=to('data/params.yml'),

        verbose=False, # True False
    )