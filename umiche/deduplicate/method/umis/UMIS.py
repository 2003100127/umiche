

import pysam, re
import pandas as pd
from collections import defaultdict

class UmiCounter:

    def __init__(self, tag_map=None, min_evidence=1.0, weighted=True, positional=False):
        self.tag_map = {'cell': 'CR', 'umi': 'UM', 'gene': 'GX'}
        if tag_map: self.tag_map.update(tag_map)
        self.min_evidence = float(min_evidence)
        self.weighted = bool(weighted)
        self.positional = bool(positional)
        self._qname_re = re.compile(r'.*:CELL_(?P<CB>.*):UMI_(?P<MB>.*)')

    def count(self, bam_path, whitelist=None):
        bam = pysam.AlignmentFile(bam_path, 'rb' if bam_path.endswith('.bam') else 'r')
        ev = defaultdict(float)
        cells, genes = set(), set()

        for aln in bam.fetch(until_eof=True):
            if aln.is_unmapped: continue

            cb, um = None, None
            try:
                cb = aln.get_tag(self.tag_map['cell'])
                um = aln.get_tag(self.tag_map['umi'])
            except KeyError:
                m = self._qname_re.match(aln.query_name)
                if not m: continue
                cb, um = m.group('CB'), m.group('MB')

            if whitelist and cb not in whitelist: continue

            gene = None
            try:
                gene = aln.get_tag(self.tag_map['gene']).split(',')[0]
            except KeyError:
                gene = bam.get_reference_name(aln.reference_id)
            if gene is None: continue

            key = (cb, gene, um, aln.pos) if self.positional else (cb, gene, um)
            w = 1.0
            if self.weighted:
                nh = dict(aln.get_tags()).get('NH', 1) or 1
                w = 1.0 / nh
            ev[key] += w
            cells.add(cb); genes.add(gene)

        collapsed = defaultdict(int)
        for (cb, gene, um, *rest), evid in ev.items():
            if evid >= self.min_evidence:
                collapsed[(cb, gene)] += 1

        df = pd.DataFrame(0, index=sorted(genes), columns=sorted(cells), dtype=int)
        for (cb, gene), cnt in collapsed.items():
            df.at[gene, cb] = cnt
        return df


if __name__ == "__main__":
    counter = UmiCounter(
        tag_map={'cell': 'CB', 'umi': 'MB', 'gene': 'XF'},
        min_evidence=1.0,
        weighted=True,
    )
    counts = counter.count(
        "/mnt/d/Document/Programming/Python/umiche/umiche/data/r1/umis/tagcount.bam"
    )
    print(counts)
