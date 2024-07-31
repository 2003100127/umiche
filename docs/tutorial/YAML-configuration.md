In UMIche, all calculation methods rely on a `YAML` file for parameter initialisation.

It defines 5 sections

- [x] `work_dir` - working directory
- [x] `trimmed` - trimming FastQ reads and extracting barcodes and UMIs
- [x] `fixed` - single-valued simulation parameters
- [x] `varied` - varying-valued simulation parameters
- [x] `dedup` - UMI deduplication parameters

# Monomer UMI pipeline settings

``` shell
work_dir: /mnt/d/Document/Programming/Python/umiche/umiche/data/simu/mclumi/
#work_dir: D:/Document/Programming/Python/umiche/umiche/data/simu/general/seq_errs/
# work_dir data/simu/tree/trimer/
# work_dir data/simu/monomer/pcr8/
# work_dir data/simu/trimer/pcr8/
# work_dir data/simu/dimer/pcr8/
# work_dir data/simu/dimer/treepcr22_250/
# work_dir data/simu/dimer/pcr8_mono24/

trimmed:
  fastq:
    fpn: None
    trimmed_fpn: None

  umi_1:
    len: 10

  seq:
    len: 100

  read_struct: 'umi_1'


fixed:
  pcr_num: 8
  pcr_err: 0.00001
  seq_err: 0.001
  ampl_rate: 0.85
  seq_dep: 400
  umi_num: 50
  permutation_num: 2
  umi_unit_pattern: 1
  umi_unit_len: 10
  seq_sub_spl_rate: 0.333
  sim_thres: 3


varied:
  pcr_nums: [ # pcr_nums_err_2d_spl0.33
    1, 2, 3,
    4, 5, 6, 7,
    8, 9, 10, 11, 12,
    13, 14, 15, 16,
#    17, 18
#    17, 18, 19, 20,
  ]
  pcr_errs: [
    0.00001,
    0.000025,
    0.00005,
    0.000075,
    0.0001,
    0.00025,
    0.0005,
    0.00075,
    0.001,
    0.0025,
    0.005,
    0.0075,
    0.01,
#    0.025,
    0.05,
#    0.075,
#    0.1,
#    0.2,
#    0.3,
  ]
  seq_errs: [
    0.00001,
    0.000025,
    0.00005,
    0.000075,
    0.0001,
    0.00025,
    0.0005,
    0.00075,
    0.001,
    0.0025,
    0.005,
    0.0075,
    0.01,
    0.025,
    0.05,
    0.075,
    0.1,
#    0.2,
#    0.3,
  ]
  ampl_rates: [0.1, 0.2, 0.3, 0.4, 0.5, 0.6, 0.7, 0.8, 0.9, 1.0]
  umi_lens: [6, 7, 8, 9, 10, 11, 12, 13, 14, 15, 16, 17, 18]
#  umi_lens: [6, 7, 8, 9, 10, 11, 12, 13, 14, 15, 16, 17, 18, 19, 20, 21, 22, 23, 24, 25, 26, 27, 28, 29, 30, 31, 32, 33, 34, 35, 36]
#  umi_nums: [50, 250, 450, 650, 850, 1050]
  umi_nums: [1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 14, 15, 16, 17, 18, 19, 20, 21, 22, 23, 24, 25, 26, 27, 28, 29, 30, 31, 32, 33, 34, 35, 36, 37, 38, 39, 40, 41, 42, 43, 44, 45]
  seq_deps: [100, 200, 500, 600, 800, 1000, 2000, 3000, 5000 ]
#  seq_deps: [100, 200, 300, 400, 500, 600, 700, 800, 900, 1000]


dedup:
  dbscan_eps: 1.5 # 1.5
  dbscan_min_spl: 1
  birch_thres: 1.8 # 1.8
  birch_n_clusters: None
  hdbscan_min_spl: 3
  aprop_preference: None
  aprop_random_state: 0

  ed_thres: 1
  mcl_fold_thres: 1.6 # 1.6
  iter_num: 100

#  inflat_val: 2.7 # 1.1 2.7
#  exp_val: 2 # 2 3
#
  inflat_val: [1.1, 2.7, 3.6]
  exp_val: 2

#  exp_val: [2, 3, 4]

  # mcl_ed trace!!!
#  ed_thres: 1
#  mcl_fold_thres: 2 # 1.6
#  inflat_val: 2.7 # 1.1 2.7
#  exp_val: 2 # 2 3
#  iter_num: 100

  # pcr_nums
  # mcl_inflat: 2.3
  # mcl_exp: 2
  # mcl_fold_thres: 1
```


# Homotrimer UMI pipeline settings

``` shell
#work_dir: d:/Document/Programming/Python/umiche/umiche/data/simu/umiche/trimer/
work_dir: /mnt/d/Document/Programming/Python/umiche/umiche/data/simu/umiche/trimer/

trimmed:
  fastq:
    fpn: None
    trimmed_fpn: None

  umi_1:
    len: 36

  seq:
    len: 100

  read_struct: 'umi_1'


fixed:
  pcr_num: 8
  pcr_err: 0.00001
  seq_err: 0.001
  ampl_rate: 0.85
  seq_dep: 400
  umi_num: 50
  permutation_num: 10
  umi_unit_pattern: 3
  umi_unit_len: 12
  seq_sub_spl_rate: 0.333
  sim_thres: 3


varied:
  pcr_nums: [ # pcr_nums_err_2d_spl0.33
    1, 2, 3,
    4, 5, 6, 7,
    8, 9, 10, 11, 12,
    13, 14, 15, 16,
#    17, 18
#    17, 18, 19, 20,
  ]
  pcr_errs: [
    0.00001,
    0.000025,
    0.00005,
    0.000075,
    0.0001,
    0.00025,
    0.0005,
    0.00075,
    0.001,
    0.0025,
    0.005,
    0.0075,
    0.01,
#    0.025,
    0.05,
#    0.075,
#    0.1,
#    0.2,
#    0.3,
  ]
  seq_errs: [
    0.00001,
    0.000025,
    0.00005,
    0.000075,
    0.0001,
    0.00025,
    0.0005,
    0.00075,
    0.001,
    0.0025,
    0.005,
    0.0075,
    0.01,
    0.025,
    0.05,
    0.075,
    0.1,
    0.2,
#    0.3,
  ]
  ampl_rates: [0.1, 0.2, 0.3, 0.4, 0.5, 0.6, 0.7, 0.8, 0.9, 1.0]
  umi_lens: [6, 7, 8, 9, 10, 11, 12, 13, 14, 15, 16, 17, 18]
#  umi_lens: [6, 7, 8, 9, 10, 11, 12, 13, 14, 15, 16, 17, 18, 19, 20, 21, 22, 23, 24, 25, 26, 27, 28, 29, 30, 31, 32, 33, 34, 35, 36]
#  umi_nums: [50, 250, 450, 650, 850, 1050]
  umi_nums: [1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 14, 15, 16, 17, 18, 19, 20, 21, 22, 23, 24, 25, 26, 27, 28, 29, 30, 31, 32, 33, 34, 35, 36, 37, 38, 39, 40, 41, 42, 43, 44, 45]
  seq_deps: [100, 200, 500, 600, 800, 1000, 2000, 3000, 5000 ]
#  seq_deps: [100, 200, 300, 400, 500, 600, 700, 800, 900, 1000]


dedup:
  dbscan_eps: 1.5 # 1.5
  dbscan_min_spl: 1
  birch_thres: 1.8 # 1.8
  birch_n_clusters: None
  hdbscan_min_spl: 3
  aprop_preference: None
  aprop_random_state: 0

  ed_thres: 1
  mcl_fold_thres: 1.6 # 1.6
  inflat_val: 2.7 # 1.1 2.7
  exp_val: 2 # 2 3
  iter_num: 100

  # mcl_ed trace!!!
#  ed_thres: 1
#  mcl_fold_thres: 2 # 1.6
#  inflat_val: 2.7 # 1.1 2.7
#  exp_val: 2 # 2 3
#  iter_num: 100

  # pcr_nums
  # mcl_inflat: 2.3
  # mcl_exp: 2
  # mcl_fold_thres: 1
```