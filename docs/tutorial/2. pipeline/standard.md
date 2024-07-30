The standard pipeline is built to perform UMI deduplication with 11 methods on a regualr basis. However, this pipeline still runs based on the output of the Tresor tool or other simulated reads.  

!!! abstract "Feature"

    The standard pipeline can give you the hints about how experimental researchers can optimise their sequencing libraries or how computational scientists can design more effective tools for UMI deduplication. Therefore, it will provide you with the statistics during the deduplication as well the deduplicated UMI count under multiple conditions.

We show an example of using the `set_cover` method to deduplicate trimer UMIs. Please be sure of the working directory set to be the root of a batch of simulated reads by `Tresor` software, where you should be able to see different permutation folders. We can use the below code to do the batch deduplication.

:material-language-python: `Python`
``` py linenums="1"
import umiche as uc

uc.pipeline.standard(
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

```

:material-console: `console`
``` shell
30/07/2024 21:25:34 logger: ======>key 1: work_dir
30/07/2024 21:25:34 logger: =========>value: /mnt/d/Document/Programming/Python/umiche/umiche/data/simu/umiche/trimer/
30/07/2024 21:25:34 logger: ======>key 2: trimmed
30/07/2024 21:25:34 logger: =========>value: {'fastq': {'fpn': 'None', 'trimmed_fpn': 'None'}, 'umi_1': {'len': 36}, 'seq': {'len': 100}, 'read_struct': 'umi_1'}
30/07/2024 21:25:34 logger: ======>key 3: fixed
30/07/2024 21:25:34 logger: =========>value: {'pcr_num': 8, 'pcr_err': 1e-05, 'seq_err': 0.001, 'ampl_rate': 0.85, 'seq_dep': 400, 'umi_num': 50, 'permutation_num': 10, 'umi_unit_pattern': 3, 'umi_unit_len': 12, 'seq_sub_spl_rate': 0.333, 'sim_thres': 3}
30/07/2024 21:25:34 logger: ======>key 4: varied
30/07/2024 21:25:34 logger: =========>value: {'pcr_nums': [1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 14, 15, 16], 'pcr_errs': [1e-05, 2.5e-05, 5e-05, 7.5e-05, 0.0001, 0.00025, 0.0005, 0.00075, 0.001, 0.0025, 0.005, 0.0075, 0.01, 0.05], 'seq_errs': [1e-05, 2.5e-05, 5e-05, 7.5e-05, 0.0001, 0.00025, 0.0005, 0.00075, 0.001, 0.0025, 0.005, 0.0075, 0.01, 0.025, 0.05, 0.075, 0.1, 0.2], 'ampl_rates': [0.1, 0.2, 0.3, 0.4, 0.5, 0.6, 0.7, 0.8, 0.9, 1.0], 'umi_lens': [6, 7, 8, 9, 10, 11, 12, 13, 14, 15, 16, 17, 18], 'umi_nums': [1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 14, 15, 16, 17, 18, 19, 20, 21, 22, 23, 24, 25, 26, 27, 28, 29, 30, 31, 32, 33, 34, 35, 36, 37, 38, 39, 40, 41, 42, 43, 44, 45], 'seq_deps': [100, 200, 500, 600, 800, 1000, 2000, 3000, 5000]}
30/07/2024 21:25:34 logger: ======>key 5: dedup
30/07/2024 21:25:34 logger: =========>value: {'dbscan_eps': 1.5, 'dbscan_min_spl': 1, 'birch_thres': 1.8, 'birch_n_clusters': 'None', 'hdbscan_min_spl': 3, 'aprop_preference': 'None', 'aprop_random_state': 0, 'ed_thres': 1, 'mcl_fold_thres': 1.6, 'inflat_val': 2.7, 'exp_val': 2, 'iter_num': 100}
UMI homopolymer recurring pattern: 3
===>Permutation number: 0
============>No.0, dedup cnt: 50.0
============>No.1, dedup cnt: 50.0
============>No.2, dedup cnt: 50.0
============>No.3, dedup cnt: 50.0
============>No.4, dedup cnt: 50.0
============>No.5, dedup cnt: 50.0
============>No.6, dedup cnt: 50.0
============>No.7, dedup cnt: 50.0
============>No.8, dedup cnt: 50.0
============>No.9, dedup cnt: 50.0
============>No.10, dedup cnt: 50.0
============>No.11, dedup cnt: 50.0
============>No.12, dedup cnt: 50.0
============>No.13, dedup cnt: 51.0
============>No.14, dedup cnt: 55.0
============>No.15, dedup cnt: 77.0
============>No.16, dedup cnt: 115.0
============>No.17, dedup cnt: 338.0
===>Permutation number: 1
============>No.0, dedup cnt: 50.0
============>No.1, dedup cnt: 50.0
============>No.2, dedup cnt: 50.0
============>No.3, dedup cnt: 50.0
============>No.4, dedup cnt: 50.0
============>No.5, dedup cnt: 50.0
============>No.6, dedup cnt: 50.0
============>No.7, dedup cnt: 50.0
============>No.8, dedup cnt: 50.0
============>No.9, dedup cnt: 50.0
============>No.10, dedup cnt: 50.0
============>No.11, dedup cnt: 50.0
============>No.12, dedup cnt: 50.0
============>No.13, dedup cnt: 51.0
============>No.14, dedup cnt: 55.0
============>No.15, dedup cnt: 76.0
============>No.16, dedup cnt: 112.0
...
```


!!! note 

    The :material-text-box-multiple: `params_trimer.yml` file configures the parameters for simulated reads and those for deduplication. Please see the relevant page for details.