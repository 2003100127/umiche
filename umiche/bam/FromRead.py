__version__ = "v1.0"
__copyright__ = "Copyright 2025"
__license__ = "GPL-3.0"
__developer__ = "Jianfeng Sun"
__maintainer__ = "Jianfeng Sun"
__email__="jianfeng.sunmt@gmail.com"

from typing import Optional, Tuple

import os
import gzip


class FromRead:

    def __init__(
            self,
    ):
        import pysam

        self.pysam = pysam

    def asarusim(
            self,
            ground_truth_path: str,
            fastq_path: str,
            bam_out_path: str,
            cb_len: int = 16,
            umi_len: int = 12,
            adapter_seq: Optional[str] = None,
            sample_name: str = "sample1",
    ):
        """
        Convert AsaruSim simulated FASTQ + ground_truth to an unmapped BAM with CB/MB/TB/XT tags.

        Parameters
        ----------
        ground_truth_path : str
            Path to results_demo/ground_truth.tsv(.gz)  (ONLY used to fetch TRUE UMI -> TB)
        fastq_path : str
            Path to results_demo/simulated.fastq(.gz)
        bam_out_path : str
            Output BAM path.
        cb_len : int
            Observed CB length at read start.
        umi_len : int
            Observed UMI length immediately following CB.
        adapter_seq : Optional[str]
            If provided and present, observed [CB][UMI] are taken from the prefix before ADAPTER.
        sample_name : str
            RG:SM value for header (optional metadata).
        """

        def open_text(path: str):
            return gzip.open(path, "rt") if path.endswith(".gz") else open(path, "r", encoding="utf-8")

        def _extract_cb_umi(seq: str, cb_len: int, umi_len: int, adapter_seq: Optional[str]) -> Tuple[str, str]:
            """Extract observed CB/UMI from a read sequence."""
            if adapter_seq:
                idx = seq.find(adapter_seq)
                if idx >= 0:
                    prefix = seq[:idx]
                    if len(prefix) >= cb_len + umi_len:
                        return prefix[:cb_len], prefix[cb_len:cb_len + umi_len]
                    # if adapter found but prefix too short, fall back to start-based parsing
            cb = seq[:cb_len] if len(seq) >= cb_len else seq
            umi = seq[cb_len:cb_len + umi_len] if len(seq) >= cb_len + umi_len else ""
            return cb, umi

        def _parse_tx_from_read_id(read_id: str) -> str:
            """Parse transcript ID from simulated read header (composite id)."""
            toks = read_id.rsplit("-", 3)
            if len(toks) == 4:
                return toks[2]
            return ""  # fallback

        def _load_true_umi_map(ground_truth_path: str) -> dict:
            """
            Build dict: read_id -> true_umi from ground_truth.tsv(.gz)
            ground truth col1 is "<CB>-<UMI>-<TranscriptID>-<ReadIndex>"
            """
            m = {}
            with open_text(ground_truth_path) as f:
                for line in f:
                    line = line.strip()
                    if not line:
                        continue
                    first = line.split("\t", 1)[0]
                    toks = first.rsplit("-", 3)  # robust from right
                    if len(toks) != 4:
                        raise ValueError(f"Unexpected ground truth ID format: {first}")
                    cb_true, umi_true, _tx, _rid = toks
                    m[first] = umi_true
            return m

        # 1) map read_id -> TRUE UMI (TB)
        true_umi = _load_true_umi_map(ground_truth_path)

        # 2) BAM header (minimal, unmapped)
        header = {
            "HD": {"VN": "1.6", "SO": "unknown"},
            "SQ": [{"SN": "DUMMY", "LN": 1}],
            "RG": [{"ID": "RG1", "SM": sample_name, "PL": "ONT"}],
        }
        os.makedirs(os.path.dirname(bam_out_path) or ".", exist_ok=True)
        out_bam = self.pysam.AlignmentFile(bam_out_path, "wb", header=header)

        # 3) iterate FASTQ
        with open_text(fastq_path) as f:
            line_no = 0
            rid = None
            tx_from_header = ""
            seq = ""
            for line in f:
                line_no += 1
                s = line.rstrip("\n")
                mod = line_no % 4
                if mod == 1:
                    # header line: @<read_id> (take up to first whitespace)
                    rid = s[1:].split()[0] if s.startswith("@") else s.split()[0]
                    tx_from_header = _parse_tx_from_read_id(rid)  # XT from simulated read header
                elif mod == 2:
                    seq = s
                elif mod == 3:
                    continue
                elif mod == 0:
                    qual = s

                    # observed CB/UMI from SEQUENCE (CB -> tag CB, UMI -> tag MB)
                    cb_obs, umi_obs = _extract_cb_umi(seq, cb_len, umi_len, adapter_seq)

                    # true UMI from ground truth mapping (TB)
                    tb_true = true_umi.get(rid, "")

                    # build unmapped record
                    a = self.pysam.AlignedSegment(out_bam.header)
                    a.query_name = rid
                    a.flag = 4  # unmapped
                    a.reference_id = -1
                    a.reference_start = -1
                    a.mapping_quality = 255
                    a.cigarstring = None
                    a.next_reference_id = -1
                    a.next_reference_start = -1
                    a.template_length = 0
                    a.query_sequence = seq
                    a.query_qualities = self.pysam.qualitystring_to_array(qual)

                    a.set_tag("RG", "RG1")
                    a.set_tag("CB", cb_obs, value_type="Z")  # observed cell barcode
                    a.set_tag("MB", umi_obs, value_type="Z")  # observed UMI
                    a.set_tag("TB", tb_true, value_type="Z")  # TRUE UMI
                    a.set_tag("XT", tx_from_header, value_type="Z")  # transcript from simulated read header

                    out_bam.write(a)

                    # reset
                    rid = ""
                    tx_from_header = ""
                    seq = ""

        out_bam.close()
        return bam_out_path

    def screadsim(
            self,
            fq_gz_path: str,
            out_bam_path: str,
            cb_len: int = 16,
            umi_len: int = 10,
    ) -> None:
        """
        Convert 10x-style R1 (CB+UMI) FASTQ (.gz) into an unmapped BAM and attach the following tags:
        CB: cell barcode (the first cb_len bases of R1)
        MB: observed UMI (bases cb_len+1 through cb_len+umi_len of R1)
        TB: true UMI (taken from the 26-mer <CB+UMI> in the FASTQ header; slice out the UMI segment)
        XT: transcript/gene name recorded in the FASTQ header during simulation (prefer fields matching GENE*)

        Supported FASTQ header format (generated by scReadSim), e.g.:
        @<26mer>:<cell_id>:<region>#<idx>
        e.g. @GTGATGATGTAGAGGTAAAAGTCCGC:CellNo1:chr1_3340000_3350000_GENE3_0_+#0214

        Produce an unmapped BAM (FLAG=4) that contains only the sequence, quality, and the tags above.
        """
        import gzip
        import pysam

        def open_text(path: str):
            return gzip.open(path, "rt") if path.endswith(".gz") else open(path, "r")

        header = {
            "HD": {"VN": "1.6", "SO": "unknown"},
            "SQ": [{"SN": "DUMMY", "LN": 1}],
        }

        n_written = 0
        with pysam.AlignmentFile(out_bam_path, "wb", header=header) as bam_out, open_text(fq_gz_path) as fq:
            while True:
                h = fq.readline().rstrip("\n")
                if not h:
                    break
                seq = fq.readline().rstrip("\n")
                _ = fq.readline()  # '+'
                qual = fq.readline().rstrip("\n")
                if not qual:
                    break

                # @<26mer>:<cell_id>:<region>#<idx>
                head = h[1:] if h.startswith("@") else h
                parts = head.split(":")
                bc26_truth = parts[0] if len(parts) > 0 else ""  # CB+UMI
                region = parts[2] if len(parts) > 2 else ""  # GENE
                # gt UMI from 26mer
                true_umi = bc26_truth[cb_len:cb_len + umi_len] if len(bc26_truth) >= cb_len + umi_len else ""

                # XT
                xt = "NA"
                if region:
                    region_core = region.split("#", 1)[0]
                    toks = region_core.split("_")
                    xt = next((t for t in toks if t.startswith("GENE")), None) or \
                         next((t for t in toks if t.startswith("ENS")), None) or \
                         region_core

                cb_measured = seq[:cb_len] if len(seq) >= cb_len else ""
                umi_measured = seq[cb_len:cb_len + umi_len] if len(seq) >= cb_len + umi_len else ""

                seg = pysam.AlignedSegment(bam_out.header)
                seg.query_name = head
                seg.flag = 4
                seg.reference_id = -1
                seg.reference_start = 0
                seg.mapping_quality = 0
                seg.cigar = None
                seg.next_reference_id = -1
                seg.next_reference_start = 0
                seg.template_length = 0
                seg.query_sequence = seq
                seg.query_qualities = pysam.qualitystring_to_array(qual)

                tags = []
                if cb_measured: tags.append(("CB", cb_measured))
                if umi_measured: tags.append(("MB", umi_measured))
                if true_umi:    tags.append(("TB", true_umi))
                if xt:          tags.append(("XT", xt))
                seg.tags = tags

                bam_out.write(seg)
                n_written += 1

        print(f"Finished with {n_written} reads to {out_bam_path}")

    def tksm(
            self,
            fq_path: str,
            bam_out: str,
            tsb_mdf: str,
            pcr_mdf: str,
            cb_len: int = 16,
            umi_len: int = 10
    ):
        """
        Convert `reads.fq(.gz)` to an **unmapped BAM** and write the following tags:

        CB: cell barcode
        MB: sequenced UMI (the last `umi_len` bases at the read's 3′ end)
        TB: ground-truth UMI (prefer `info(TB=...)` from `tsb.mdf`; otherwise from
        `pcr.mdf`; otherwise fall back to MB)
        XT: simulated transcript name (prefer the contig of the first interval in
        `pcr.mdf`; otherwise from `tsb.mdf`)

        Molecule matching:
        Parse `molecule_id=XXXX` from the read header (example: `M0_23.1.2.3`).
        Look up this molecule in `pcr.mdf` to obtain its information.
        Ancestor molecule ID = the first segment of `molecule_id` (e.g., `M0_23`);
        then query `tsb.mdf` to retrieve the true CB/TB/XT.

        Return: the number of reads written (`int`).

        """
        import re

        def open_text(path):
            return gzip.open(path, "rt") if path.endswith(".gz") else open(path, "rt", encoding="utf-8")

        def _parse_mdf(path, umi_len):
            mol2info, mol_first_tx = {}, {}
            with open_text(path) as fh:
                cur = None
                seen_polyA = False
                first_interval_done = False
                for line in fh:
                    if not line.strip():
                        continue
                    if line.startswith('+'):
                        cur = None
                        seen_polyA = False
                        first_interval_done = False
                        parts = line[1:].strip().split()
                        if parts:
                            cur = parts[0]
                            mol2info.setdefault(cur, {})
                        m = re.search(r"info\((.*?)\)", line)
                        if m:
                            for k, v in re.findall(r"([A-Za-z0-9_]+)=([^;]+)", m.group(1)):
                                if k == "tid":
                                    mol2info[cur]["tid"] = v;
                                    continue
                                if k in ("CB", "TB", "XT", "TX", "TR", "Transcript", "transcript"):
                                    if k in ("XT", "TX", "TR", "Transcript", "transcript"):
                                        mol2info[cur]["XT"] = v
                                    else:
                                        mol2info[cur][k] = v
                        header_kv = dict(re.findall(r"\b([A-Za-z0-9_]+)=([^\s;]+)", line.split("info(")[0]))
                        if "tid" in header_kv: mol2info[cur]["tid"] = header_kv["tid"]
                        if "TB" in header_kv and "TB" not in mol2info[cur]: mol2info[cur]["TB"] = header_kv["TB"]
                        if "CB" in header_kv and "CB" not in mol2info[cur]: mol2info[cur]["CB"] = header_kv["CB"]
                        continue

                    if not cur:
                        continue
                    cols = line.strip().split()
                    if not cols:
                        continue

                    if re.fullmatch(r"[ACGTNacgtn]+", cols[0]):
                        seqtok = cols[0].upper()
                        if set(seqtok) <= {"A"}:
                            seen_polyA = True
                            continue
                        if seen_polyA and len(seqtok) == int(umi_len) and "TB" not in mol2info[cur]:
                            mol2info[cur]["TB"] = seqtok
                        continue

                    if not first_interval_done:
                        mol_first_tx[cur] = cols[0]
                        first_interval_done = True

            for mid, tx in mol_first_tx.items():
                mol2info.setdefault(mid, {})
                mol2info[mid].setdefault("XT", tx)
            return mol2info

        def _read_fastq_iter(path):
            with open_text(path) as fh:
                while True:
                    name = fh.readline()
                    if not name:
                        return
                    seq = fh.readline().rstrip("\n")
                    plus = fh.readline()
                    qual = fh.readline().rstrip("\n")
                    if not qual:
                        return
                    yield name.strip(), seq, qual

        def _root_molecule(mid: str):
            return mid.split('.', 1)[0] if (mid and '.' in mid) else mid

        # mdf
        tsb_info = _parse_mdf(tsb_mdf, umi_len) if tsb_mdf else {}
        pcr_info = _parse_mdf(pcr_mdf, umi_len) if pcr_mdf else {}

        header = {
            "HD": {"VN": "1.6", "SO": "unknown"},
            "SQ": [{"SN": "DUMMY", "LN": 1}],
        }
        out = self.pysam.AlignmentFile(bam_out, "wb", header=header)

        n = 0
        for name, seq, qual in _read_fastq_iter(fq_path):
            rid = name[1:] if name.startswith('@') else name
            # parse molecule_id
            m = re.search(r"molecule_id=([^\s;]+)", rid)
            mol_id = m.group(1) if m else None
            root_id = mol_id.split('.', 1)[0] if mol_id and '.' in mol_id else mol_id
            pinfo = pcr_info.get(mol_id, {}) if mol_id else {}
            tinfo = tsb_info.get(root_id, {}) if root_id else {}
            # CB: TSB -> PCR ->
            CB = tinfo.get("CB") or pinfo.get("CB")
            if not CB:
                CB = seq[:cb_len] if cb_len > 0 and len(seq) >= cb_len else "NA"
            # MB
            MB = seq[-umi_len:] if umi_len > 0 and len(seq) >= umi_len else "NA"
            # TB
            TB = pinfo.get("TB") or tinfo.get("TB") or MB
            # XT: XT/TX/TR
            XT = (
                    pinfo.get("tid") or tinfo.get("tid") or
                    pinfo.get("XT") or tinfo.get("XT") or
                    pinfo.get("TX") or tinfo.get("TX") or
                    pinfo.get("TR") or tinfo.get("TR") or
                    "NA"
            )
            a = self.pysam.AlignedSegment()
            a.query_name = rid.split()[0]
            a.query_sequence = seq
            a.flag = 0x4
            a.mapping_quality = 0
            a.reference_id = -1
            a.reference_start = -1
            a.cigar = None
            a.next_reference_id = -1
            a.next_reference_start = -1
            a.template_length = 0
            a.query_qualities = self.pysam.qualitystring_to_array(qual)
            # MI molecule_id
            # mi = mol_id if mol_id else "NA"
            a.tags = [("CB", CB), ("MB", MB), ("TB", TB), ("XT", XT),
                      ("MI", mol_id if mol_id else "NA")]
            out.write(a)
            n += 1
        out.close()
        return n

    def minnow(self, r1_fq: str, out_bam: str, cb_len: int = 16, umi_len: int = 12):
        """
        Convert MINNOW R1 FASTQ into an unmapped BAM with tags:
          CB: cell barcode  (R1[0:cb_len])
          MB: observed UMI  (R1[cb_len:cb_len+umi_len])
          TB: true UMI      (no GT provided -> set equal to MB for now)
          XT: transcript id (from read name like '@<...>:tx10:...')

        Only R1 is used. Supports .gz.
        """
        from pathlib import Path

        def _open(path):
            return gzip.open(path, "rt") if str(path).endswith(".gz") else open(path, "rt")

        # prepare BAM (no reference; unmapped)
        header = {
            "HD": {"VN": "1.6", "SO": "unknown"},
            "SQ": [{"SN": "DUMMY", "LN": 1}],
            "PG": [{"ID": "r1_fastq_to_bam", "PN": "r1_fastq_to_bam", "VN": "0.1"}]
        }
        Path(out_bam).parent.mkdir(parents=True, exist_ok=True)
        bam = self.pysam.AlignmentFile(out_bam, "wb", header=header)

        with _open(r1_fq) as f1:
            n = 0
            while True:
                h = f1.readline()
                if not h:
                    break
                seq = f1.readline().rstrip()
                plus = f1.readline()
                qual = f1.readline().rstrip()
                if not qual:
                    raise RuntimeError("Malformed FASTQ: missing quality line.")

                qname = h.strip()
                if qname.startswith("@"): qname = qname[1:]
                # XT: transcript from second colon-delimited token (e.g., 'tx10')
                xt = ""
                parts = qname.split(":")
                if len(parts) >= 2:
                    xt = parts[1]

                if len(seq) < cb_len + umi_len:
                    print(seq)
                    raise RuntimeError(
                        f"R1 too short for CB({cb_len})+UMI({umi_len}): len={len(seq)}; QNAME={qname}"
                    )
                cb = seq[:cb_len]
                mb = seq[cb_len:cb_len + umi_len]
                tb = mb  # no ground-truth supplied -> use observed UMI

                seg = self.pysam.AlignedSegment()
                seg.query_name = qname
                seg.query_sequence = seq
                seg.query_qualities = self.pysam.qualitystring_to_array(qual)
                seg.is_unmapped = True
                seg.mapping_quality = 0

                seg.set_tag("CB", cb, value_type="Z")
                seg.set_tag("MB", mb, value_type="Z")
                seg.set_tag("TB", tb, value_type="Z")
                seg.set_tag("XT", xt, value_type="Z")

                bam.write(seg)
                n += 1

        bam.close()
        print(f"Finished with {n} records to {out_bam}")

    def tresor(
            self,
            fastq_fpn: str,
            bam_fpn: str,
            umi_start: int,
            umi_len: int,
    ) -> None:
        """
        Convert a given FASTQ file to an unmapped BAM file and attach the following tags to each read:

        CB: cell barcode (extracted from the FASTQ header using s*<number>)
        XF: gene ID (extracted from the FASTQ header using g*<number>)
        MB: UMI (subsequence taken from the read sequence at a given start and length)

        Applicable FASTQ header example
        @29*s*0*g*5*_1_4_5-pcr-5
        Here, the number after s (0) is recorded as CB, and the number after g (5) is recorded as XF.

        Parameters
        ----------
        fastq_path : str
        Path to the input FASTQ file. Supports both .gz and plain text.
        bam_path : str
        Path to the output BAM file (will be created/overwritten).
        umi_start : int
        Zero-based start index of the UMI within the read sequence (inclusive).
        If negative, count from the end of the sequence (e.g., -10 means the 10th base
        from the end is the start).
        umi_len : int
        Length of the UMI (must be positive).

        Notes
        -----
        The output is an unmapped BAM: flag = 4, with no reference coordinates.
        Only the CB, XF, and MB tags are parsed and written; the original read sequence
        and qualities are unchanged.
        If a read is too short to extract the specified UMI segment, a ValueError is raised.
        The s and g values are captured with the regexes r's\*(\d+)' and r'g\*(\d+)'; if not
        found, an empty string is written.

        """
        import re
        import io

        if umi_len <= 0:
            raise ValueError("umi_len needs to be greater than 0")

        def _open_text_maybe_gz(path: str) -> io.TextIOBase:
            if path.endswith(".gz"):
                return gzip.open(path, "rt")
            return open(path, "rt")

        header = {
            "HD": {"VN": "1.6"},
            "SQ": [{"SN": "DUMMY", "LN": 1}],
        }

        re_s = re.compile(r"s\*(\d+)")
        re_g = re.compile(r"g\*(\d+)")

        with _open_text_maybe_gz(fastq_fpn) as fq, self.pysam.AlignmentFile(bam_fpn, "wb", header=header) as out_bam:
            lineno = 0
            while True:
                name = fq.readline()
                if not name:
                    break  # EOF
                seq = fq.readline()
                plus = fq.readline()
                qual = fq.readline()
                lineno += 4

                if not qual:
                    raise ValueError(f"FASTQ structure is imcomplete: {lineno}")

                name = name.rstrip("\n\r")
                seq = seq.rstrip("\n\r")
                qual = qual.rstrip("\n\r")

                if not name.startswith("@"):
                    raise ValueError(f"FASTQ header line format error (should start with “@”): line {lineno - 3} -> {name}")

                # gene and cell no.
                m_s = re_s.search(name)
                m_g = re_g.search(name)
                cb = m_s.group(1) if m_s else ""
                xf = m_g.group(1) if m_g else ""

                # UMI pos
                if umi_start < 0:
                    start = len(seq) + umi_start
                else:
                    start = umi_start
                end = start + umi_len

                if start < 0 or end > len(seq):
                    raise ValueError(
                        f"读取 UMI 越界：read='{name[1:]}' 序列长度={len(seq)}，请求区间=[{start}, {end})。"
                    )
                umi = seq[start:end]

                a = self.pysam.AlignedSegment()
                a.query_name = name[1:]  # del '@'
                a.query_sequence = seq
                a.flag = 4  # unmapped
                a.mapping_quality = 0
                a.is_unmapped = True
                a.reference_id = -1
                a.reference_start = -1
                a.next_reference_id = -1
                a.next_reference_start = -1
                a.template_length = 0
                a.cigar = None

                a.query_qualities = self.pysam.qualitystring_to_array(qual)

                a.set_tag("CB", cb)
                a.set_tag("XF", xf)
                a.set_tag("MB", umi)
                out_bam.write(a)


if __name__ == "__main__":
    # Example usage — adjust paths as needed
    p = FromRead()
    # out = p.asarusim(
    #     ground_truth_path="/mnt/d/Document/Programming/Python/umiche/umiche/data/r1/cm/asarusim/results_demo/ground_truth.tsv.gz",
    #     fastq_path="/mnt/d/Document/Programming/Python/umiche/umiche/data/r1/cm/asarusim/results_demo/simulatedp3-0.005.fastq.gz",
    #     bam_out_path="/mnt/d/Document/Programming/Python/umiche/umiche/data/r1/cm/asarusim/results_demo/simulated_with_tags.bam",
    #     cb_len=16,
    #     umi_len=12,
    #     adapter_seq=None,  # or e.g. "ACTAAAGGCCATTACGGCCTACACGACGCTCTTCCGATCT"
    #     sample_name="demo",
    # )
    # print(f"[done] wrote: {out}")

    # p.screadsim(
    #     fq_gz_path= "/mnt/d/Document/Programming/Python/umiche/umiche/data/r1/cm/screadsim/outputs/10X_RNA_chr1_3073253_4526737.syntheticBAM.gene10.read1.bed2fa.sorted.fq.gz",
    #     out_bam_path= "/mnt/d/Document/Programming/Python/umiche/umiche/data/r1/cm/screadsim/outputs/reads_with_tags.bam",
    # )

    # count = p.tksm(
    #     fq_path="/mnt/d/Document/Programming/Python/umiche/umiche/data/r1/cm/tksm/reads.fq.gz",
    #     bam_out="/mnt/d/Document/Programming/Python/umiche/umiche/data/r1/cm/tksm/reads_with_tags.bam",
    #     tsb_mdf="/mnt/d/Document/Programming/Python/umiche/umiche/data/r1/cm/tksm/tsb.mdf",
    #     pcr_mdf="/mnt/d/Document/Programming/Python/umiche/umiche/data/r1/cm/tksm/pcr.mdf",
    #     cb_len=16,
    #     umi_len=10,
    # )

    # count = p.minnow(
    #     r1_fq='/mnt/d/Document/Programming/Python/umiche/umiche/data/r1/cm/minnow/out/sim_pcr_0p001/hg_100_S1_L001_R1_001.fastq.gz',
    #     out_bam='/mnt/d/Document/Programming/Python/umiche/umiche/data/r1/cm/minnow/out/sim_pcr_0p001/reads_with_tags.bam',
    #     cb_len=16,
    #     umi_len=10,
    # )

    p.tresor(
        # @@ tresor
        fastq_fpn="/mnt/d/Document/Programming/Python/umiche/umiche/data/r1/simuread/tresor/pcr_err_0.fastq.gz",
        bam_fpn = "/mnt/d/Document/Programming/Python/umiche/umiche/data/r1/simuread/tresor/tresor.bam",
        umi_start=0,
        umi_len=12,

        # @@ scifiseq
        # fastq_fpn="/mnt/d/Document/Programming/Python/umiche/umiche/data/r1/seqtech/simu/scifiseq/pcr_err_0.fastq.gz",
        # bam_fpn = "/mnt/d/Document/Programming/Python/umiche/umiche/data/r1/seqtech/simu/scifiseq/scifiseq.bam",
        # umi_start=0,
        # umi_len=8,

        # @@ scrbseq
        # fastq_fpn="/mnt/d/Document/Programming/Python/umiche/umiche/data/r1/seqtech/simu/scrbseq/pcr_err_0.fastq.gz",
        # bam_fpn="/mnt/d/Document/Programming/Python/umiche/umiche/data/r1/seqtech/simu/scrbseq/scrbseq.bam",
        # umi_start=6,
        # umi_len=10,

        # @@ shareseq
        # fastq_fpn="/mnt/d/Document/Programming/Python/umiche/umiche/data/r1/seqtech/simu/shareseq/pcr_err_0.fastq.gz",
        # bam_fpn="/mnt/d/Document/Programming/Python/umiche/umiche/data/r1/seqtech/simu/shareseq/shareseq.bam",
        # umi_start=0,
        # umi_len=10,

        # @@ snareseq
        # fastq_fpn="/mnt/d/Document/Programming/Python/umiche/umiche/data/r1/seqtech/simu/snareseq/pcr_err_0.fastq.gz",
        # bam_fpn="/mnt/d/Document/Programming/Python/umiche/umiche/data/r1/seqtech/simu/snareseq/snareseq.bam",
        # umi_start=12,
        # umi_len=8,

        # @@ splitseq
        # fastq_fpn="/mnt/d/Document/Programming/Python/umiche/umiche/data/r1/seqtech/simu/splitseq/pcr_err_0.fastq.gz",
        # bam_fpn="/mnt/d/Document/Programming/Python/umiche/umiche/data/r1/seqtech/simu/splitseq/splitseq.bam",
        # umi_start=0,
        # umi_len=10,

        # # @@ strtseqc1
        # fastq_fpn="/mnt/d/Document/Programming/Python/umiche/umiche/data/r1/seqtech/simu/strtseqc1/pcr_err_0.fastq.gz",
        # bam_fpn="/mnt/d/Document/Programming/Python/umiche/umiche/data/r1/seqtech/simu/strtseqc1/strtseqc1.bam",
        # umi_start=0,
        # umi_len=5,

        # @@ strtseq2i
        # fastq_fpn="/mnt/d/Document/Programming/Python/umiche/umiche/data/r1/seqtech/simu/strtseq2i/pcr_err_0.fastq.gz",
        # bam_fpn="/mnt/d/Document/Programming/Python/umiche/umiche/data/r1/seqtech/simu/strtseq2i/strtseq2i.bam",
        # umi_start=0,
        # umi_len=6,

        # # @@ vasaseq
        # fastq_fpn="/mnt/d/Document/Programming/Python/umiche/umiche/data/r1/seqtech/simu/vasaseq/pcr_err_0.fastq.gz",
        # bam_fpn="/mnt/d/Document/Programming/Python/umiche/umiche/data/r1/seqtech/simu/vasaseq/vasaseq.bam",
        # umi_start=0,
        # umi_len=6,

        # # @@ quartzseq2
        # fastq_fpn="/mnt/d/Document/Programming/Python/umiche/umiche/data/r1/seqtech/simu/quartzseq2/pcr_err_0.fastq.gz",
        # bam_fpn="/mnt/d/Document/Programming/Python/umiche/umiche/data/r1/seqtech/simu/quartzseq2/quartzseq2.bam",
        # umi_start=15,
        # umi_len=8,

        # @@ petipseq
        # fastq_fpn="/mnt/d/Document/Programming/Python/umiche/umiche/data/r1/seqtech/simu/petipseq/pcr_err_0.fastq.gz",
        # bam_fpn="/mnt/d/Document/Programming/Python/umiche/umiche/data/r1/seqtech/simu/petipseq/petipseq.bam",
        # umi_start=0,
        # umi_len=7,

        # @@ pipseq
        # fastq_fpn="/mnt/d/Document/Programming/Python/umiche/umiche/data/r1/seqtech/simu/pipseq/pcr_err_0.fastq.gz",
        # bam_fpn="/mnt/d/Document/Programming/Python/umiche/umiche/data/r1/seqtech/simu/pipseq/pipseq.bam",
        # umi_start=28,
        # umi_len=12,

        # # @@ pairedseq
        # fastq_fpn="/mnt/d/Document/Programming/Python/umiche/umiche/data/r1/seqtech/simu/pairedseq/pcr_err_0.fastq.gz",
        # bam_fpn="/mnt/d/Document/Programming/Python/umiche/umiche/data/r1/seqtech/simu/pairedseq/pairedseq.bam",
        # umi_start=22,
        # umi_len=10,

        # # @@ marsseq2
        # fastq_fpn="/mnt/d/Document/Programming/Python/umiche/umiche/data/r1/seqtech/simu/marsseq2/pcr_err_0.fastq.gz",
        # bam_fpn="/mnt/d/Document/Programming/Python/umiche/umiche/data/r1/seqtech/simu/marsseq2/marsseq2.bam",
        # umi_start=7,
        # umi_len=8,

        # @@ issaacseq
        # fastq_fpn="/mnt/d/Document/Programming/Python/umiche/umiche/data/r1/seqtech/simu/issaacseq/pcr_err_0.fastq.gz",
        # bam_fpn="/mnt/d/Document/Programming/Python/umiche/umiche/data/r1/seqtech/simu/issaacseq/issaacseq.bam",
        # umi_start=0,
        # umi_len=10,

        # @@ indrop
        # fastq_fpn="/mnt/d/Document/Programming/Python/umiche/umiche/data/r1/seqtech/simu/indrop/pcr_err_0.fastq.gz",
        # bam_fpn="/mnt/d/Document/Programming/Python/umiche/umiche/data/r1/seqtech/simu/indrop/indrop.bam",
        # umi_start=8,
        # umi_len=6,

        # @@ dropseq
        # fastq_fpn="/mnt/d/Document/Programming/Python/umiche/umiche/data/r1/seqtech/simu/dropseq/pcr_err_0.fastq.gz",
        # bam_fpn="/mnt/d/Document/Programming/Python/umiche/umiche/data/r1/seqtech/simu/dropseq/dropseq.bam",
        # umi_start=12,
        # umi_len=8,

        # # @@ flashseq
        # fastq_fpn="/mnt/d/Document/Programming/Python/umiche/umiche/data/r1/seqtech/simu/flashseq/pcr_err_0.fastq.gz",
        # bam_fpn="/mnt/d/Document/Programming/Python/umiche/umiche/data/r1/seqtech/simu/flashseq/flashseq.bam",
        # umi_start=0,
        # umi_len=8,

        # # @@ celseq2
        # fastq_fpn="/mnt/d/Document/Programming/Python/umiche/umiche/data/r1/seqtech/simu/celseq2/pcr_err_0.fastq.gz",
        # bam_fpn="/mnt/d/Document/Programming/Python/umiche/umiche/data/r1/seqtech/simu/celseq2/celseq2.bam",
        # umi_start=0,
        # umi_len=6,

        # @@ 10xv2
        # fastq_fpn="/mnt/d/Document/Programming/Python/umiche/umiche/data/r1/seqtech/simu/10xv2/pcr_err_0.fastq.gz",
        # bam_fpn="/mnt/d/Document/Programming/Python/umiche/umiche/data/r1/seqtech/simu/10xv2/10xv2.bam",
        # umi_start=16,
        # umi_len=10,

        # # @@ 10xv3
        # fastq_fpn="/mnt/d/Document/Programming/Python/umiche/umiche/data/r1/seqtech/simu/10xv3/pcr_err_0.fastq.gz",
        # bam_fpn="/mnt/d/Document/Programming/Python/umiche/umiche/data/r1/seqtech/simu/10xv3/10xv3.bam",
        # umi_start=16,
        # umi_len=12,
    )