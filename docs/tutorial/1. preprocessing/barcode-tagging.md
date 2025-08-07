
## Analogy with other tools

### 1) UMI-tools

#### Read groupping

The `--per-cell` and `--per-gene` attributes are used to tell UMI-tools to group reads within their respective cells and genes. But what they eactly used for read groupping is needed to be specified separately via `--cell-tag` and `--gene-tag`, which use a 2-character tag (e.g., `XT`) found in the `tag` field of each read in the `.bam` file.

1. For `--per-gene`, it must be used in combination with `--gene-tag`.

2. For `--per-cell`, it can optionally be used in combination with `--cell-tag`.

``` shell
--per-gene
--gene-tag
```

``` shell
--gene-tag
--cell-tag
--umi-tag
```

!!! note "gene"

    If you have used `--gene-tag`, then you also need to use `--assigned-status-tag` too, which is a bool parameter, specified in `tag` region as well and first introduced by FeatureCounts.

    ``` shell
    --assigned-status-tag=XS
    ```

#### UMI retrieval

It uses the `--umi-tag` attribute to extract a UMI per read, with its respective tag stored in a `.bam` file. For example, to extract UMIs tagged with the **`MB`** tag, it is used in combination with the `--extract-umi-method` attribute. 

``` shell
--extract-umi-method=tag --umi-tag=MB
```

When this is not available, it uses `--umi-separator` to get UMIs from read names. For example, to extract UMIs separated with `_`, it does

``` shell
--umi-separator=_
```

