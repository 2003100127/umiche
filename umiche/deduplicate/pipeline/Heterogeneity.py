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

from umiche.deduplicate.OnePos import OnePos as dedupop
from umiche.plot.Valid import valid as plotv
from umiche.util.Writer import Writer as gwriter
from umiche.util.Console import Console

from umiche.deduplicate.heterogeneity.Trace import Trace as umitrace

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
    ):
        self.scenario = scenario
        self.method = method

        self.params = params(param_fpn=param_fpn)
        self.gbfscc = gbfscc()
        self.gwriter = gwriter()
        self.plotv = plotv()

        self.verbose = verbose
        self.console = Console()
        self.console.verbose = self.verbose

        df_dedup = pd.DataFrame()
        for perm_num_i in range(self.params.fixed['permutation_num']):
            self.console.print("===>permutation number {}".format(perm_num_i))
            dedup_arr = []
            for id, scenario_i in enumerate(self.params.varied[self.scenario]):
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
                    dedup_ob = dedupop(
                        bam_fpn=self.fastq_location + 'trimmed/bam/' + self.fn_prefix + '.bam',
                        # pos_tag='PO',
                        mcl_fold_thres=self.params.dedup['mcl_fold_thres'],
                        inflat_val=self.params.dedup['inflat_val'],
                        exp_val=self.params.dedup['exp_val'],
                        iter_num=self.params.dedup['iter_num'],
                        ed_thres=self.params.dedup['ed_thres'],
                        work_dir=self.params.work_dir,
                        heterogeneity=True,
                        verbose=False,
                    )
                    df = self.tool(dedup_ob)[self.method]()
                    # print(dedup_ob.directional())
                    # print(df.columns)
                    # print(df.dedup_cnt.values[0])

                    # print(df.apv.values)
                    from umiche.deduplicate.method.Directional import Directional as umitoolmonodirec

                    umitooldirec = umitoolmonodirec()

                    df_direc_apv = umitooldirec.formatApvsDisapv(df.apv.values[0])
                    df_direc_disapv = umitooldirec.formatApvsDisapv(df.disapv.values[0])
                    # print(df_direc_apv)

                    from umiche.bam.Relation import Relation as umirel
                    umiold = umirel(
                        df=dedup_ob.df_bam
                    )
                    # print('jsun', umiold.umi_trace_dict)
                    umiidtrace = umitrace(
                        df_umi_uniq_val_cnt=umiold.df_umi_uniq_val_cnt,
                        umi_trace_dict=umiold.umi_id_to_origin_id_dict,
                    )
                    print(dedup_ob.df_bam.columns)
                    print(list(umiidtrace.edgecls(df_list_2d=df_direc_apv, sort='cnt')) + [str(scenario_i)] + ['direc'])
                    print(list(umiidtrace.edgecls(df_list_2d=df_direc_disapv, sort='cnt')) + [str(scenario_i)] + ['direc'])
                    # print(df.disapv.values)
                    # print(df.direc.values)
                    # dedup_arr.append(dedup_ob.dedup_num)
            # df_dedup['pn' + str(perm_num_i)] = dedup_arr
            # print(df_dedup)
        # self.gwriter.generic(
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
            # 'set_cover': dedup_ob.set_cover,
        }


    def statistics(self, ):
        return {
            'ccs': {
                'df_apv_cnt': pd.DataFrame(columns=['0', '1', 'all', 'scenario', 'method']),
                'df_disapv_cnt': pd.DataFrame(columns=['0', '1', 'all', 'scenario', 'method']),
                'df_apv_pct': pd.DataFrame(columns=['0', '1', 'scenario', 'method']),
                'df_disapv_pct': pd.DataFrame(columns=['0', '1', 'scenario', 'method']),
                'df_dedup_cnt': pd.DataFrame(columns=['dedup_cnt', 'scenario', 'method']),
            },
            'adj': {
                'df_apv_cnt': pd.DataFrame(columns=['0', '1', 'all', 'scenario', 'method']),
                'df_disapv_cnt': pd.DataFrame(columns=['0', '1', 'all', 'scenario', 'method']),
                'df_apv_pct': pd.DataFrame(columns=['0', '1', 'scenario', 'method']),
                'df_disapv_pct': pd.DataFrame(columns=['0', '1', 'scenario', 'method']),
                'df_dedup_cnt': pd.DataFrame(columns=['dedup_cnt', 'scenario', 'method']),
            },
            'direc': {
                'df_apv_cnt': pd.DataFrame(columns=['0', '1', 'all', 'scenario', 'method']),
                'df_disapv_cnt': pd.DataFrame(columns=['0', '1', 'all', 'scenario', 'method']),
                'df_apv_pct': pd.DataFrame(columns=['0', '1', 'scenario', 'method']),
                'df_disapv_pct': pd.DataFrame(columns=['0', '1', 'scenario', 'method']),
                'df_dedup_cnt': pd.DataFrame(columns=['dedup_cnt', 'scenario', 'method']),
            },
            'mcl_val': {
                'df_apv_cnt': pd.DataFrame(columns=['0', '1', 'all', 'scenario', 'method']),
                'df_disapv_cnt': pd.DataFrame(columns=['0', '1', 'all', 'scenario', 'method']),
                'df_apv_pct': pd.DataFrame(columns=['0', '1', 'scenario', 'method']),
                'df_disapv_pct': pd.DataFrame(columns=['0', '1', 'scenario', 'method']),
                'df_dedup_cnt': pd.DataFrame(columns=['dedup_cnt', 'scenario', 'method']),
            },
            'mcl_ed': {
                'df_apv_cnt': pd.DataFrame(columns=['0', '1', 'all', 'scenario', 'method']),
                'df_disapv_cnt': pd.DataFrame(columns=['0', '1', 'all', 'scenario', 'method']),
                'df_apv_pct': pd.DataFrame(columns=['0', '1', 'scenario', 'method']),
                'df_disapv_pct': pd.DataFrame(columns=['0', '1', 'scenario', 'method']),
                'df_dedup_cnt': pd.DataFrame(columns=['dedup_cnt', 'scenario', 'method']),
            },
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

        param_fpn=to('data/seqerr_sl.yml')
    )
