__version__ = "v1.0"
__copyright__ = "Copyright 2025"
__license__ = "GPL-3.0"
__developer__ = "Jianfeng Sun"
__maintainer__ = "Jianfeng Sun"
__email__="jianfeng.sunmt@gmail.com"

from typing import Union, List, Dict
import numpy as np

from umiche.plot.graph.Cluster import Cluster
from umiche.plot.scenario.DedupMultiple import DedupMultiple
from umiche.plot.scenario.DedupSingle import DedupSingle
from umiche.plot.scenario.DedupMultipleTrimer import DedupMultipleTrimer
from umiche.plot.line.TripletErrorCode import TripletErrorCode
from umiche.plot.line.Anchor import Anchor
from umiche.plot.scenario.DedupMultipleTrimerSetCover import DedupMultipleTrimerSetCover
from umiche.plot.scenario.TraceSingle import TraceSingle
from umiche.plot.scenario.TraceMultiple import TraceMultiple
from umiche.plot.scenario.Line import Line


def graph_cluster(
        dedup_cluster,
        ccs,
        edge_list,
        umi_uniq_val_cnt_fpn,
):
    return Cluster(
        dedup_cluster=dedup_cluster,
        ccs=ccs,
        edge_list=edge_list,
        umi_uniq_val_cnt_fpn=umi_uniq_val_cnt_fpn,
    )


def dedup_multiple(
        scenarios,
        methods,
        param_fpn,
):
    return DedupMultiple(
        scenarios=scenarios,
        methods=methods,
        param_fpn=param_fpn,
    )


def dedup_single(
        df_dedup,
        df_dedup_perm_melt,
):
    return DedupSingle(
        df=df_dedup,
        df_melt=df_dedup_perm_melt,
    )


def dedup_multiple_trimer(
        scenarios,
        methods,
        param_fpn,
):
    return DedupMultipleTrimer(
        scenarios=scenarios,
        methods=methods,
        param_fpn=param_fpn,
    )


def trace_single(
        df_apv,
        df_disapv,
):
    return TraceSingle(
        df_apv=df_apv,
        df_disapv=df_disapv,
    )


def trace_multiple(
        df_apv,
        df_disapv,
        scenarios,
        methods,
):
    return TraceMultiple(
        df_apv=df_apv,
        df_disapv=df_disapv,
        scenarios=scenarios,
        methods=methods,
    )


def inflat_exp(
        scenarios,
        methods,
        param_fpn,
):
    return Line(
        scenarios=scenarios,
        methods=methods,
        param_fpn=param_fpn,
    )


def prob_correct(
        error_rate: Union[float, np.ndarray]=0.00001,
        num_nt=12,
):
    TripletErrorCode(error_rate=error_rate).correct(num_nt=num_nt)


def prob_incorrect(
        error_rate: Union[float, np.ndarray]=0.00001,
        num_nt = 12,
):
    TripletErrorCode(error_rate=error_rate).incorrect(num_nt=num_nt)


def anchor_efficiency(
        quant_captured: Dict,
        quant_anchor_captured: Dict,
        criteria : List,
        condition : str,
):
    Anchor(
        quant_captured=quant_captured,
        quant_anchor_captured=quant_anchor_captured,
        criteria=criteria,
        condition=condition,
    ).draw()


def anchor_efficiency_broken(
        quant_captured: Dict,
        quant_anchor_captured: Dict,
        criteria : List,
        condition : str,
):
    Anchor(
        quant_captured=quant_captured,
        quant_anchor_captured=quant_anchor_captured,
        criteria=criteria,
        condition=condition,
    ).draw_line_broken()


def anchor_efficiency_simple(
        quant_captured: Dict,
        quant_anchor_captured: Dict,
        criteria : List,
        condition : str,
):
    Anchor(
        quant_captured=quant_captured,
        quant_anchor_captured=quant_anchor_captured,
        criteria=criteria,
        condition=condition,
    ).simple()


if __name__ == "__main__":
    from umiche.path import to

    plot_dm = dedup_multiple(
        scenarios={
            'pcr_nums': 'PCR cycle',
            'pcr_errs': 'PCR error',
            'seq_errs': 'Sequencing error',
            'ampl_rates': 'Amplification rate',
            'umi_lens': 'UMI length',
            'seq_deps': 'Sequencing depth',
        },
        methods={
            'unique': 'Unique',
            'cluster': 'Cluster',
            'adjacency': 'Adjacency',
            'directional': 'Directional',
            'dbscan_seq_onehot': 'DBSCAN',
            'birch_seq_onehot': 'Birch',
            'aprop_seq_onehot': 'Affinity Propagation',
            'mcl': 'MCL',
            'mcl_val': 'MCL-val',
            'mcl_ed': 'MCL-ed',
        },
        param_fpn=to('data/params.yml'),
    )
    # print(plot_dm.line())

    # print(graph_cluster(
    #     umi_uniq_val_cnt_fpn=to('data/simu/mclumi/seq_errs/umi_uniq_val_cnt-simu.txt'),
    #     # umi_uniq_val_cnt_fpn=to('data/simu/mclumi/seq_errs/umi_uniq_val_cnt-iclip.txt'),
    #     dedup_cluster=[
    #         [69, 72, 838, 1221, 97, 210, 249, 315, 324, 374, 457, 658, 727, 760, 771, 933, 1073, 1126, 1198, 1260, 1271,
    #          1307, 1498, 1505, 1541, 1563, 914, 946, 1083, 684, 1288, 1543, 822],
    #         [1174, 119, 290, 303, 204, 218, 289, 302, 404, 545, 586, 633, 674, 709, 720, 802, 884, 943, 980, 1355, 1436,
    #          1488, 1553, 786, 537, 867, 1649, 1255, 701, 1080, 347, 251], [1315, 1549]],
    #     ccs={69: [72, 838, 1221],
    #          72: [69, 97, 210, 249, 315, 324, 374, 457, 658, 727, 760, 771, 838, 933, 1073, 1126, 1198, 1221, 1260,
    #               1271, 1307, 1498, 1505, 1541, 1563], 838: [69, 72, 1221], 1221: [69, 72, 838, 914, 946],
    #          97: [72, 658, 1505], 210: [72, 1083, 1271], 249: [72, 684, 727, 1073], 315: [72, 374, 946, 1126],
    #          324: [72, 771], 374: [72, 315, 1126], 457: [72, 1288, 1307], 658: [72, 97, 1505], 727: [72, 249, 1073],
    #          760: [72, 933, 1543], 771: [72, 324], 933: [72, 760], 1073: [72, 249, 727], 1126: [72, 315, 374],
    #          1198: [72, 1260, 1541], 1260: [72, 1083, 1198, 1541, 1543], 1271: [72, 210], 1307: [72, 457],
    #          1498: [72, 684, 1563], 1505: [72, 97, 658, 822], 1541: [72, 1198, 1260], 1563: [72, 914, 1498],
    #          914: [1221, 1563], 946: [315, 1221], 1083: [210, 1260], 684: [249, 1498], 1288: [457],
    #          1543: [760, 1174, 1260], 822: [1505], 1174: [119, 290, 303, 1315, 1543],
    #          119: [204, 218, 289, 290, 302, 303, 404, 545, 586, 633, 674, 709, 720, 802, 884, 943, 980, 1174, 1355,
    #                1436, 1488, 1553], 290: [119, 303, 786, 1174], 303: [119, 290, 1174], 1315: [1174, 1549],
    #          204: [119, 537, 720, 1488], 218: [119, 302, 674, 867, 1649], 289: [119, 980],
    #          302: [119, 218, 674, 1255], 404: [119, 802, 943], 545: [119, 701, 1436], 586: [119, 709, 884],
    #          633: [119, 1080, 1553], 674: [119, 218, 302], 709: [119, 586, 884], 720: [119, 204, 347, 1488],
    #          802: [119, 404, 537, 943], 884: [119, 586, 701, 709], 943: [119, 404, 802], 980: [119, 251, 289],
    #          1355: [119, 251, 867], 1436: [119, 545], 1488: [119, 204, 720], 1553: [119, 633, 786], 786: [290, 1553],
    #          1549: [1315], 537: [204, 802], 867: [218, 1355, 1649], 1649: [218, 867], 1255: [302], 701: [545, 884],
    #          1080: [633], 347: [720], 251: [980, 1355]},
    #     edge_list=[(97, 72), (658, 97), (1355, 119), (204, 119), (1649, 218), (1541, 72), (933, 72), (914, 1563),
    #                (1505, 72), (1221, 69), (404, 119), (1649, 867), (1307, 72), (302, 119), (802, 119), (1543, 760),
    #                (1126, 72), (709, 586), (786, 1553), (374, 315), (946, 315), (1073, 72), (720, 204), (933, 760),
    #                (1174, 1543), (210, 72), (771, 72), (537, 204), (727, 249), (1563, 72), (684, 1498), (822, 1505),
    #                (1549, 1315), (1505, 97), (545, 119), (374, 72), (1553, 119), (1436, 119), (867, 218), (457, 72),
    #                (1563, 1498), (709, 119), (720, 119), (1260, 1198), (1543, 1260), (218, 119), (251, 1355),
    #                (701, 545), (1541, 1260), (1488, 204), (1553, 633), (1083, 210), (980, 289), (119, 1174),
    #                (727, 72), (290, 119), (1073, 727), (1541, 1198), (884, 586), (802, 404), (303, 119), (347, 720),
    #                (674, 218), (1221, 838), (1221, 72), (914, 1221), (943, 119), (586, 119), (943, 802), (838, 72),
    #                (1271, 210), (980, 119), (674, 302), (786, 290), (251, 980), (1488, 720), (1488, 119),
    #                (946, 1221), (1073, 249), (884, 709), (1288, 457), (1126, 315), (943, 404), (537, 802),
    #                (302, 218), (249, 72), (315, 72), (303, 290), (290, 1174), (884, 119), (633, 119), (771, 324),
    #                (658, 72), (1505, 658), (1307, 457), (303, 1174), (72, 69), (289, 119), (838, 69), (684, 249),
    #                (1498, 72), (1315, 1174), (1126, 374), (701, 884), (1255, 302), (674, 119), (1080, 633),
    #                (1083, 1260), (1436, 545), (1260, 72), (1271, 72), (867, 1355), (760, 72), (1198, 72), (324, 72)],
    #
    # ).draw())

    # print(inflat_exp(
    #     scenarios={
    #         'pcr_nums': 'PCR cycle',
    #         'pcr_errs': 'PCR error',
    #         'seq_errs': 'Sequencing error',
    #         'ampl_rates': 'Amplification rate',
    #         'umi_lens': 'UMI length',
    #         'seq_deps': 'Sequencing depth',
    #     },
    #     methods={
    #         'mcl': 'MCL',
    #     },
    #     param_fpn=to('data/params.yml'),
    # ).draw())

    # from umiche.io import stat
    # scenarios = {
    #     # 'pcr_nums': 'PCR cycle',
    #     # 'pcr_errs': 'PCR error',
    #     'seq_errs': 'Sequencing error',
    #     # 'ampl_rates': 'Amplification rate',
    #     # 'umi_lens': 'UMI length',
    #     # 'seq_deps': 'Sequencing depth',
    # }
    # methods = {
    #     # 'unique': 'Unique',
    #     # 'cluster': 'Cluster',
    #     # 'adjacency': 'Adjacency',
    #     'directional': 'Directional',
    #     # 'dbscan_seq_onehot': 'DBSCAN',
    #     # 'birch_seq_onehot': 'Birch',
    #     # 'aprop_seq_onehot': 'Affinity Propagation',
    #     'mcl': 'MCL',
    #     'mcl_val': 'MCL-val',
    #     'mcl_ed': 'MCL-ed',
    # }
    # dedupstat = stat(
    #     scenarios=scenarios,
    #     methods=methods,
    #     param_fpn=to('data/params.yml'),
    # )
    # df_dedup = dedupstat.df_dedup
    # df_dedup_perm_melt = dedupstat.df_dedup_perm_melt
    # t = dedup_single(
    #     df_dedup=df_dedup,
    #     df_dedup_perm_melt=df_dedup_perm_melt,
    # )
    # # print(t.jointgrid(
    # #     # method='Directional',
    # #     # method='MCL',
    # #     # method='MCL-val',
    # #     method='MCL-ed',
    # # ))
    # # print(t.stackedbar())
    # # print(t.strip())

    # print(dedup_multiple_trimer(
    #     scenarios={
    #         'seq_errs': 'Sequencing error rate',
    #     },
    #     methods={
    #         'directional_dedupby_majority_vote_splitby__collblockby_take_by_order': 'UMI-tools+drMV+cbRAN',
    #         'directional_dedupby_majority_vote_splitby__collblockby_majority_vote': 'UMI-tools+drMV+cbMV',
    #         'directional_dedupby_set_cover_splitby_split_by_mv_collblockby_take_by_order': 'UMI-tools+drSC+spMV+cbRAN',
    #         'directional_dedupby_set_cover_splitby_split_by_mv_collblockby_majority_vote': 'UMI-tools+drSC+spMV+cbMV',
    #         'directional_dedupby_set_cover_splitby_split_to_all_collblockby_take_by_order': 'UMI-tools+drSC+spALL+cbRAN',
    #         'directional_dedupby_set_cover_splitby_split_to_all_collblockby_majority_vote': 'UMI-tools+drSC+spALL+cbMV',
    #     },
    #     param_fpn=to('data/params_trimer.yml'),
    # ).line())


    ### ++++++++++++++++++++trace++++++++++++++++++++++++++
    from umiche.deduplicate.io.Stat import Stat as dedupstat
    scenarios = {
        'pcr_nums': 'PCR cycle',
        # 'pcr_errs': 'PCR error',
        # 'seq_errs': 'Sequencing error',
        # 'ampl_rates': 'Amplification rate',
        # 'umi_lens': 'UMI length',
        # 'seq_deps': 'Sequencing depth',
    }
    methods = {
        # 'unique': 'Unique',
        # 'cluster': 'Cluster',
        # 'adjacency': 'Adjacency',
        'directional': 'Directional',
        # 'dbscan_seq_onehot': 'DBSCAN',
        # 'birch_seq_onehot': 'Birch',
        # 'aprop_seq_onehot': 'Affinity Propagation',
        # 'mcl': 'MCL',
        # 'mcl_val': 'MCL-val',
        # 'mcl_ed': 'MCL-ed',
    }
    dedupstat11 = dedupstat(
        scenarios=scenarios,
        methods=methods,
        param_fpn=to('data/params.yml'),
    )
    # trace_single(
    #     df_apv=dedupstat11.df_trace_cnt['apv'],
    #     df_disapv=dedupstat11.df_trace_cnt['disapv'],
    # ).line_apv_disapv()

    scenarios = {
        'pcr_nums': 'PCR cycle',
        # 'pcr_errs': 'PCR error',
        # 'seq_errs': 'Sequencing error',
        # 'ampl_rates': 'Amplification rate',
        # 'umi_lens': 'UMI length',
        # 'seq_deps': 'Sequencing depth',
    }
    methods = {
        'unique': 'Unique',
        'cluster': 'Cluster',
        'adjacency': 'Adjacency',
        'directional': 'Directional',
        'dbscan_seq_onehot': 'DBSCAN',
        'birch_seq_onehot': 'Birch',
        'aprop_seq_onehot': 'Affinity Propagation',
        'mcl': 'MCL',
        'mcl_val': 'MCL-val',
        'mcl_ed': 'MCL-ed',
    }
    dedupstat22 = dedupstat(
        scenarios=scenarios,
        methods=methods,
        param_fpn=to('data/params.yml'),
    )

    # trace_multiple(
    #     df_apv=dedupstat22.df_trace_cnt['apv'],
    #     df_disapv=dedupstat22.df_trace_cnt['disapv'],
    #     scenarios=scenarios,
    #     methods=methods,
    # ).line_apv()

    error_rates= np.linspace(0.00001, 0.5, 500)
    prob_correct(error_rate=error_rates, num_nt=1)
    prob_incorrect(error_rate=error_rates, num_nt=12)