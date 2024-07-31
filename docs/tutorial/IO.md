# Read BAM files

A [BAM](https://support.illumina.com/help/BS_App_RNASeq_Alignment_OLH_1000000006112/Content/Source/Informatics/BAM-Format.htm) file can be fed into UMIche as shown below.

:material-language-python: `Python`
``` py linenums="1"
import umiche as uc

bam = uc.io.read_bam(
    bam_fpn="/mnt/d/Document/Programming/Python/umiche/umiche/data/simu/umi/trimer/seq_errs/permute_0/trimmed/seq_err_17.bam",
    verbose=True,
)
```

If you know where each read is related or originates with a certain tag, you can extract only this proportion of reads from the BAM file by `tags=['PO']`.

:material-language-python: `Python`
``` py linenums="1"
print(bam.todf(tags=['PO']))
```

:material-console: `console`
``` shell
31/07/2024 01:36:42 logger: ===>reading the bam file... /mnt/d/Document/Programming/Python/umiche/umiche/data/simu/umi/trimer/seq_errs/permute_0/trimmed/seq_err_17.bam
31/07/2024 01:36:42 logger: ===>reading BAM time: 0.03s
31/07/2024 01:36:42 logger: =========>start converting bam to df...
31/07/2024 01:36:42 logger: =========>time to df: 0.033s
        id  ... PO
0        0  ...  1
1        1  ...  1
2        2  ...  1
3        3  ...  1
4        4  ...  1
...    ...  ... ..
6944  6944  ...  1
6945  6945  ...  1
6946  6946  ...  1
6947  6947  ...  1
6948  6948  ...  1

[6949 rows x 13 columns]
```

# Read deduplication files

We build a powerful module `uc.io.stat` for processing deduplicated UMI counts stored in file `{method}_dedup.txt`. This file can be obtained by running the UMIche pipelines or solely the deduplication methods. It can handle files of multiple sequencing conditions and multiple methods.

:material-language-python: `Python`
``` py linenums="1"
import umiche as uc

stat_data = uc.io.stat(
    scenarios={
        'pcr_nums': 'PCR cycle',
        'pcr_errs': 'PCR error',
        'seq_errs': 'Sequencing error',
        'ampl_rates': 'Amplification rate',
        'umi_lens': 'UMI length',
        'seq_deps': 'Sequencing depth',
    },
    methods={
        # 'unique': 'Unique',
        # 'cluster': 'Cluster',
        # 'adjacency': 'Adjacency',
        'directional': 'Directional',
        # 'dbscan_seq_onehot': 'DBSCAN',
        # 'birch_seq_onehot': 'Birch',
        # 'aprop_seq_onehot': 'Affinity Propagation',
        'mcl': 'MCL',
        # 'mcl_val': 'MCL-val',
        # 'mcl_ed': 'MCL-ed',
    },
    param_fpn=to('data/params.yml'),
    verbose=True,
)
print(stat_data)
```

:material-console: `console`
``` shell
31/07/2024 01:31:24 logger: ======>key 1: work_dir
31/07/2024 01:31:24 logger: =========>value: /mnt/d/Document/Programming/Python/umiche/umiche/data/simu/mclumi/
31/07/2024 01:31:24 logger: ======>key 2: trimmed
31/07/2024 01:31:24 logger: =========>value: {'fastq': {'fpn': 'None', 'trimmed_fpn': 'None'}, 'umi_1': {'len': 10}, 'seq': {'len': 100}, 'read_struct': 'umi_1'}
31/07/2024 01:31:24 logger: ======>key 3: fixed
31/07/2024 01:31:24 logger: =========>value: {'pcr_num': 8, 'pcr_err': 1e-05, 'seq_err': 0.001, 'ampl_rate': 0.85, 'seq_dep': 400, 'umi_num': 50, 'permutation_num': 2, 'umi_unit_pattern': 1, 'umi_unit_len': 10, 'seq_sub_spl_rate': 0.333, 'sim_thres': 3}
31/07/2024 01:31:24 logger: ======>key 4: varied
31/07/2024 01:31:24 logger: =========>value: {'pcr_nums': [1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 14, 15, 16], 'pcr_errs': [1e-05, 2.5e-05, 5e-05, 7.5e-05, 0.0001, 0.00025, 0.0005, 0.00075, 0.001, 0.0025, 0.005, 0.0075, 0.01, 0.05], 'seq_errs': [1e-05, 2.5e-05, 5e-05, 7.5e-05, 0.0001, 0.00025, 0.0005, 0.00075, 0.001, 0.0025, 0.005, 0.0075, 0.01, 0.025, 0.05, 0.075, 0.1], 'ampl_rates': [0.1, 0.2, 0.3, 0.4, 0.5, 0.6, 0.7, 0.8, 0.9, 1.0], 'umi_lens': [6, 7, 8, 9, 10, 11, 12, 13, 14, 15, 16, 17, 18], 'umi_nums': [1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 14, 15, 16, 17, 18, 19, 20, 21, 22, 23, 24, 25, 26, 27, 28, 29, 30, 31, 32, 33, 34, 35, 36, 37, 38, 39, 40, 41, 42, 43, 44, 45], 'seq_deps': [100, 200, 500, 600, 800, 1000, 2000, 3000, 5000]}
31/07/2024 01:31:24 logger: ======>key 5: dedup
31/07/2024 01:31:24 logger: =========>value: {'dbscan_eps': 1.5, 'dbscan_min_spl': 1, 'birch_thres': 1.8, 'birch_n_clusters': 'None', 'hdbscan_min_spl': 3, 'aprop_preference': 'None', 'aprop_random_state': 0, 'ed_thres': 1, 'mcl_fold_thres': 1.6, 'iter_num': 100, 'inflat_val': [1.1, 2.7, 3.6], 'exp_val': 2}
31/07/2024 01:31:24 logger: ======>scenario: PCR cycle
31/07/2024 01:31:24 logger: =========>method: Directional
31/07/2024 01:31:24 logger: =========>method: MCL
31/07/2024 01:31:24 logger: ======>scenario: PCR error
31/07/2024 01:31:24 logger: =========>method: Directional
31/07/2024 01:31:24 logger: =========>method: MCL
31/07/2024 01:31:24 logger: ======>scenario: Sequencing error
31/07/2024 01:31:24 logger: =========>method: Directional
31/07/2024 01:31:24 logger: =========>method: MCL
31/07/2024 01:31:24 logger: ======>scenario: Amplification rate
31/07/2024 01:31:24 logger: =========>method: Directional
31/07/2024 01:31:24 logger: =========>method: MCL
31/07/2024 01:31:24 logger: ======>scenario: UMI length
31/07/2024 01:31:24 logger: =========>method: Directional
31/07/2024 01:31:24 logger: =========>method: MCL
31/07/2024 01:31:24 logger: ======>scenario: Sequencing depth
31/07/2024 01:31:24 logger: =========>method: Directional
31/07/2024 01:31:24 logger: =========>method: MCL
    pn0  pn1  pn2  pn3  ...  max-mean          scenario       method  metric
0   0.0  0.0  0.0  0.0  ...       0.0         PCR cycle  Directional       1
1   0.0  0.0  0.0  0.0  ...       0.0         PCR cycle  Directional       2
2   0.0  0.0  0.0  0.0  ...       0.0         PCR cycle  Directional       3
3   0.0  0.0  0.0  0.0  ...       0.0         PCR cycle  Directional       4
4   0.0  0.0  0.0  0.0  ...       0.0         PCR cycle  Directional       5
..  ...  ...  ...  ...  ...       ...               ...          ...     ...
4   0.0  0.0  0.0  0.0  ...       0.0  Sequencing depth          MCL     800
5   0.0  0.0  0.0  0.0  ...       0.0  Sequencing depth          MCL    1000
6   0.0  0.0  0.0  0.0  ...       0.0  Sequencing depth          MCL    2000
7   0.0  0.0  0.0  0.0  ...       0.0  Sequencing depth          MCL    3000
8   0.0  0.0  0.0  0.0  ...       0.0  Sequencing depth          MCL    5000

[158 rows x 19 columns]
```


# Read Inflation and expansion files

Fold change of deduplication with respective to different inflation and expansion values.

:material-language-python: `Python`
``` py linenums="1"
print(stat_data.df_inflat_exp)
```

:material-console: `console`
=== "`inflation`"

    ``` shell
    PCR cycle  PCR error  ...  UMI length  Sequencing depth
    1.10       0.00       0.14  ...        0.02               0.0
    1.45       0.04       0.16  ...        0.02               0.0
    1.80       0.04       0.26  ...        0.02               0.0
    2.15       0.04       0.46  ...        0.02               0.0
    2.50       0.06       0.50  ...        0.02               0.0
    2.85       0.12       0.58  ...        0.02               0.0
    3.20       0.16       0.78  ...        0.02               0.0
    3.55       0.18       0.90  ...        0.02               0.0
    3.90       0.20       0.90  ...        0.02               0.0
    4.25       0.20       0.90  ...        0.02               0.0
    4.60       0.22       0.90  ...        0.02               0.0
    4.95       0.22       0.90  ...        0.02               0.0
    5.30       0.22       0.92  ...        0.02               0.0
    5.65       0.22       1.02  ...        0.02               0.0
    6.00       0.22       1.04  ...        0.02               0.0
    
    [15 rows x 6 columns]
    ```

=== "`expansion`"

    ``` shell
    PCR cycle  PCR error  ...  UMI length  Sequencing depth
    2       0.06       0.56  ...        0.02               0.0
    3       0.04       0.22  ...        0.02               0.0
    4       0.04       0.16  ...        0.02               0.0
    5       0.04       0.16  ...        0.02               0.0
    6       0.04       0.16  ...        0.02               0.0
    7       0.04       0.16  ...        0.02               0.0
    8       0.04       0.16  ...        0.02               0.0
    9       0.04       0.16  ...        0.02               0.0
    
    [8 rows x 6 columns]
    ```

# UMI trajectory files

Merged UMI nodes (`apv`)



# Read UMI trajectory files

:material-language-python: `Python`
``` py linenums="1"
print(stat_data.df_trace_cnt)
```

:material-console: `console`
=== "`approved`"

    ``` shell
    metric  diff_origin  ...          scenario       method
    0      1.0          0.0  ...         PCR cycle  Directional
    1      2.0          0.0  ...         PCR cycle  Directional
    2      3.0          0.0  ...         PCR cycle  Directional
    3      4.0          0.0  ...         PCR cycle  Directional
    4      5.0          0.0  ...         PCR cycle  Directional
    ..     ...          ...  ...               ...          ...
    4    800.0          0.0  ...  Sequencing depth          MCL
    5   1000.0          0.0  ...  Sequencing depth          MCL
    6   2000.0          0.0  ...  Sequencing depth          MCL
    7   3000.0          0.0  ...  Sequencing depth          MCL
    8   5000.0          0.0  ...  Sequencing depth          MCL
    
    [158 rows x 9 columns]
    ```

=== "`disapproved`"

    ``` shell
    metric  diff_origin  ...          scenario       method
    0       1          0.0  ...         PCR cycle  Directional
    1       2          0.0  ...         PCR cycle  Directional
    2       3          0.0  ...         PCR cycle  Directional
    3       4          0.0  ...         PCR cycle  Directional
    4       5          0.0  ...         PCR cycle  Directional
    ..    ...          ...  ...               ...          ...
    4     800          0.0  ...  Sequencing depth  Directional
    5    1000          0.0  ...  Sequencing depth  Directional
    6    2000          0.0  ...  Sequencing depth  Directional
    7    3000          0.0  ...  Sequencing depth  Directional
    8    5000          0.0  ...  Sequencing depth  Directional
    
    [79 rows x 9 columns]
    ```