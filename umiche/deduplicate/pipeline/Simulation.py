__version__ = "v1.0"
__copyright__ = "Copyright 2024"
__license__ = "MIT"
__developer__ = "Jianfeng Sun"
__maintainer__ = "Jianfeng Sun"
__email__="jianfeng.sunmt@gmail.com"
__lab__ = "Cribbslab"

import pandas as pd
from umiche.fastq.Convert import Convert as fastqconverter
from umiche.trim.Template import Template as trimmer
from umiche.graph.bfs.ConnectedComponent import ConnectedComponent as gbfscc
from umiche.simu.Parameter import Parameter as params

from umiche.deduplicate.MultiPos import MultiPos as deduppos
from umiche.util.Writer import Writer as fwriter
from umiche.util.Console import Console


class Simulation:

    def __init__(
            self,
            scenario,
            method,
            param_fpn=None,
            is_trim=False,
            is_tobam=False,
            is_dedup=False,
            verbose=False,
            **kwargs,
    ):
        self.scenario = scenario
        self.method = method
        self.kwargs = kwargs

        self.params = params(param_fpn=param_fpn)
        self.gbfscc = gbfscc()
        self.fwriter = fwriter()

        self.verbose = verbose
        self.console = Console()
        self.console.verbose = self.verbose

        df_dedup = pd.DataFrame()
        for perm_num_i in range(self.params.fixed['permutation_num']):
            self.console.print("===>permutation number {}".format(perm_num_i))
            dedup_arr = []
            for id, scenario_i in enumerate(self.params.varied[self.scenario]):
                self.console.print("===>No.{} scenario: {}".format(id+1, scenario_i))
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
                    dedup_ob = deduppos(
                        bam_fpn=self.fastq_location + 'trimmed/bam/set_cover/' + self.fn_prefix + '.bam',
                        pos_tag='PO',
                        mcl_fold_thres=self.params.dedup['mcl_fold_thres'],
                        inflat_val=self.params.dedup['inflat_val'],
                        exp_val=self.params.dedup['exp_val'],
                        iter_num=self.params.dedup['iter_num'],
                        ed_thres=self.params.dedup['ed_thres'],
                        work_dir=self.params.work_dir,
                        sv_setcover_bam_fpn=self.fastq_location + 'trimmed/bam/set_cover/' + self.fn_prefix + '.bam',
                        heterogeneity=False, # False True
                        is_build_graph=True, # False True
                        verbose=False,
                        **self.kwargs
                    )
                    df = self.tool(dedup_ob)[self.method]()

                    print(df.dedup_cnt)
                    # if self.method == "unique":
                    #     print("No.{}, dedup cnt: {}".format(id, df.num_uniq_umis.values[0]))
                    #     dedup_arr.append(df.num_uniq_umis.values[0])
                    # else:
                    #     print("No.{}, dedup cnt: {:.2f}\t{}".format(id, c, df.dedup_cnt.values[0]))
                    #     dedup_arr.append(df.dedup_cnt.values[0])
                    # dedup_arr.append(dedup_ob.dedup_num)
            # df_dedup['pn' + str(perm_num_i)] = dedup_arr
            # print(df_dedup)
        # self.fwriter.generic(
        #     df=df_dedup,
        #     sv_fpn=fastq_fp + self.scenario + '/' + str(self.method) + '_' + self.comp_cat + '.txt',
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
        }


if __name__ == "__main__":
    from umiche.path import to

    p = Simulation(
        # scenario='pcr_nums',
        # scenario='pcr_errs',
        scenario='seq_errs',
        # scenario='ampl_rates',
        # scenario='umi_lens',

        # method='unique',
        # method='cluster',
        # method='adjacency',
        method='directional',
        # method='mcl',
        # method='mcl_val',
        # method='mcl_ed',
        # method='set_cover',

        # is_trim=True,
        # is_tobam=False,
        # is_dedup=False,

        # is_trim=False,
        # is_tobam=True,
        # is_dedup=False,

        is_trim=False,
        is_tobam=False,
        is_dedup=True,

        param_fpn=to('data/params_trimer.yml'),
        # param_fpn=to('data/params.yml'),

        verbose=True
    )