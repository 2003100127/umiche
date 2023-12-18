from umiche.align import BundlePos as bp


def umi_tools(
        bam_fpn,
        options,
        sv_fpn,
):
    bp.convert(options=options, in_fpn=bam_fpn, out_fpn=sv_fpn)
    return


if __name__ == "__main__":
    from umiche.path import to

    options = {
        'stats': 'deduplicated',
        'get_umi_method': 'read_id',
        'umi_sep': '_',
        'umi_tag': 'RX',
        'umi_tag_split': None,
        'umi_tag_delim': None,
        'cell_tag': None,
        'cell_tag_split': '-',
        'cell_tag_delim': None,
        'filter_umi': None,
        'umi_whitelist': None,
        'umi_whitelist_paired': None,
        'method': 'directional',
        'threshold': 1,
        'spliced': False,
        'soft_clip_threshold': 4,
        'read_length': False,
        'per_gene': False,
        'gene_tag': None,
        'assigned_tag': None,
        'skip_regex': '^(__|Unassigned)',
        'per_contig': False,
        'gene_transcript_map': None,
        'per_cell': False,
        'whole_contig': False,
        'detection_method': None,
        'mapping_quality': 0,
        'output_unmapped': False,
        'unmapped_reads': 'discard',
        'chimeric_pairs': 'use',
        'unpaired_reads': 'use',
        'ignore_umi': False,
        'ignore_tlen': False,
        'chrom': None,
        'subset': None,
        'in_sam': False,
        'paired': False,
        'out_sam': False,
        'no_sort_output': False,
        'stdin': "<_io.TextIOWrapper name='example.bam' mode='r' encoding='UTF-8'>",
        'stdlog': "<_io.TextIOWrapper name='<stdout>' mode='w' encoding='UTF-8'>", 'log2stderr': False,
        'compresslevel': 6,
        'timeit_file': None,
        'timeit_name': 'all',
        'timeit_header': None,
        'loglevel': 1,
        'short_help': None,
        'random_seed': None
    }

    print(umi_tools(
        bam_fpn=to('data/example.bam'),
        options=options,
        sv_fpn=to('data/example_bundle1.bam'),
    ))
