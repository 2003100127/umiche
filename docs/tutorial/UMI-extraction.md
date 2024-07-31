Most of sequencing short-read technologies produce sequencing reads in a fixed form. UMI segment can be counted from the beginning of read 1. For these fix-length reads, we designed a trimming method for extracting UMIs, barcodes, or other read components in bulk.

Given a FastQ file, We can specifically extract UMIs with the code below.

:material-language-python: `Python`
``` py linenums="1"
import umiche as uc

params = {
    'umi_1': {
        'len': 10,
    },
    'umi_2': {
        'len': 4,
    },
    'bc_1': {
        'len': 2,
    },
    'read_struct': 'umi_1',
    # 'read_struct': 'umi_1+seq_1',
    # 'read_struct': 'bc_1+umi_1+seq_1',
    'seq_1': {
        'len': 6,
    },
    'fastq': {
        'fpn': to('data/simu/mclumi/seq_errs/permute_0/seq_err_5.fastq.gz'),
        'trimmed_fpn': to('data/simu/mclumi/seq_errs/permute_0/seq_err_5_trimmed.fastq.gz'),
    },
}

df = uc.trim.template(params=params)
print(df)
```

The following output contains three dataframes. The first dataframe is converted from raw FastQ reads directly. The second one is with UMIs trimmed from reads. The third one is with the rest of the reads for genomic sequences.

:material-console: `console`
``` shell
31/07/2024 02:48:51 logger: ===>reading from fastq...
before trimmed          seq_raw                name
0     TCGAATTCAG      17_2_4_5-pcr-5
1     TCAGAGTGAC     8_1_5_6_8-pcr-8
2     GTACCCGATT   1_3_5_6_7_8-pcr-8
3     ATGGACTTCG  21_2_4_5_6_8-pcr-8
4     CAAGCAGCTG  47_1_2_4_6_7-pcr-7
...          ...                 ...
6944  GCACTTCGAC      44_2_5_7-pcr-7
6945  TCAGAGTGAC           8_8-pcr-8
6946  GGCGACGGCA    27_2_3_7_8-pcr-8
6947  GCACGATCAC      10_1_2_7-pcr-7
6948  AACGTAAAGG    41_1_3_6_8-pcr-8

[6949 rows x 2 columns]
31/07/2024 02:48:51 logger: ===>umi structure: umi_1
31/07/2024 02:48:51 logger: ===>bc positions in the read structure: 
31/07/2024 02:48:51 logger: ===>umi positions in the read structure: 0
31/07/2024 02:48:51 logger: ===>seq positions in the read structure: 
31/07/2024 02:48:51 logger: ======>finding the starting positions of all UMIs...
31/07/2024 02:48:51 logger: ======>finding the starting positions of all UMIs...
31/07/2024 02:48:51 logger: =========>umi_1 starting position: 0
31/07/2024 02:48:51 logger: ======>finding the starting positions of all genomic sequence...
UMI: start: 0 end: 10
31/07/2024 02:48:51 logger: ===>umi_1 has been taken out
after trimmed          seq_raw                name       umi_1
0     TCGAATTCAG      17_2_4_5-pcr-5  TCGAATTCAG
1     TCAGAGTGAC     8_1_5_6_8-pcr-8  TCAGAGTGAC
2     GTACCCGATT   1_3_5_6_7_8-pcr-8  GTACCCGATT
3     ATGGACTTCG  21_2_4_5_6_8-pcr-8  ATGGACTTCG
4     CAAGCAGCTG  47_1_2_4_6_7-pcr-7  CAAGCAGCTG
...          ...                 ...         ...
6944  GCACTTCGAC      44_2_5_7-pcr-7  GCACTTCGAC
6945  TCAGAGTGAC           8_8-pcr-8  TCAGAGTGAC
6946  GGCGACGGCA    27_2_3_7_8-pcr-8  GGCGACGGCA
6947  GCACGATCAC      10_1_2_7-pcr-7  GCACGATCAC
6948  AACGTAAAGG    41_1_3_6_8-pcr-8  AACGTAAAGG

[6949 rows x 3 columns]
31/07/2024 02:48:51 logger: ===>start saving in gz format...
31/07/2024 02:48:51 logger: ['seq_raw', 'name', 'umi_1', 'seq_1']
31/07/2024 02:48:52 logger: ===>trimmed UMIs have been saved in gz format.
         seq_raw                name       umi_1 seq_1 seq
0     TCGAATTCAG      17_2_4_5-pcr-5  TCGAATTCAG     B   B
1     TCAGAGTGAC     8_1_5_6_8-pcr-8  TCAGAGTGAC     B   B
2     GTACCCGATT   1_3_5_6_7_8-pcr-8  GTACCCGATT     B   B
3     ATGGACTTCG  21_2_4_5_6_8-pcr-8  ATGGACTTCG     B   B
4     CAAGCAGCTG  47_1_2_4_6_7-pcr-7  CAAGCAGCTG     B   B
...          ...                 ...         ...   ...  ..
6944  GCACTTCGAC      44_2_5_7-pcr-7  GCACTTCGAC     B   B
6945  TCAGAGTGAC           8_8-pcr-8  TCAGAGTGAC     B   B
6946  GGCGACGGCA    27_2_3_7_8-pcr-8  GGCGACGGCA     B   B
6947  GCACGATCAC      10_1_2_7-pcr-7  GCACGATCAC     B   B
6948  AACGTAAAGG    41_1_3_6_8-pcr-8  AACGTAAAGG     B   B

[6949 rows x 5 columns]
```

!!! info "Note"

    In the above example, our sequencing reads are simulated with only UMIs and `read_struct` is the read structure where all read components are present in reads sequentially. If `umi_1` is specified in the structure, then its length has to be specified.