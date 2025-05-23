__version__ = "v1.0"
__copyright__ = "Copyright 2025"
__license__ = "GPL-3.0"
__developer__ = "Jianfeng Sun"
__maintainer__ = "Jianfeng Sun"
__email__="jianfeng.sunmt@gmail.com"

from umiche.deduplicate.pipeline.Heterogeneity import Heterogeneity as pipeheter
from umiche.deduplicate.pipeline.Standard import Standard as pipestd
from umiche.recovery.pipeline.Anchor import Anchor


def standard(
        param_fpn : str,
        scenario : str,
        method : str,
        deduped_method : str,
        split_method : str,
        is_collapse_block: bool,
        is_trim : bool = False,
        is_tobam : bool = False,
        is_dedup : bool = False,
        is_sv : bool = False,
        verbose : bool = True,
):
    return pipestd(
        scenario=scenario,
        method=method,
        is_trim=is_trim,
        is_tobam=is_tobam,
        is_dedup=is_dedup,
        is_collapse_block=is_collapse_block,
        deduped_method=deduped_method,
        split_method=split_method,
        param_fpn=param_fpn,
        is_sv=is_sv,
        verbose=verbose,
    )


def heterogeneity(
        param_fpn: str,
        scenario : str,
        method : str,
        is_trim : bool = False,
        is_tobam : bool = False,
        is_dedup : bool = False,
        is_sv : bool = False,
        verbose : bool = False,
):
    return pipeheter(
        param_fpn=param_fpn,
        scenario=scenario,
        method=method,
        is_trim=is_trim,
        is_tobam=is_tobam,
        is_dedup=is_dedup,
        is_sv=is_sv,
        verbose=verbose,
    )


def anchor(
        scenario: str,
        param_fpn: str,
):
    return Anchor(
        scenario=scenario,
        param_fpn=param_fpn,
    ).pct_read_captured()


if __name__ == "__main__":
    from umiche.path import to

    # print(heterogeneity(
    #     # scenario='pcr_nums',
    #     # scenario='pcr_errs',
    #     scenario='seq_errs',
    #     # scenario='ampl_rates',
    #     # scenario='umi_lens',
    #     # scenario='seq_deps',
    #     # scenario='umi_nums',
    #
    #     # method='unique',
    #     # method='cluster',
    #     # method='adjacency',
    #     method='directional',
    #     # method='mcl',
    #     # method='mcl_val',
    #     # method='mcl_ed',
    #     # method='mcl_cc_all_node_umis',
    #     # method='dbscan_seq_onehot',
    #     # method='birch_seq_onehot',
    #     # method='aprop_seq_onehot',
    #     # method='hdbscan_seq_onehot',
    #     # method='set_cover',
    #
    #     # is_trim=True,
    #     # is_tobam=False,
    #     # is_dedup=False,
    #
    #     # is_trim=False,
    #     # is_tobam=True,
    #     # is_dedup=False,
    #
    #     is_trim=False,
    #     is_tobam=False,
    #     is_dedup=True,
    #     is_sv=True,
    #
    #     param_fpn=to('data/params.yml'),
    # ))

    print(standard(
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
    ))