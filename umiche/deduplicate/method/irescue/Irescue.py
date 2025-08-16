__version__ = "v1.0"
__copyright__ = "Copyright 2025"
__license__ = "GPL-3.0"
__developer__ = "Jianfeng Sun"
__maintainer__ = "Jianfeng Sun"
__email__="jianfeng.sunmt@gmail.com"

from typing import Dict, Iterable, List, Optional, Sequence, Set, Tuple

import numpy as np
import pandas as pd
import networkx as nx

from dataclasses import dataclass



@dataclass(frozen=True)
class EquivalenceClass:
    """
    An equivalence class (EC) as defined in IRescue: reads carrying the same
    UMI *and* mapping to the same *set* of features (e.g., TE subfamilies).

    Attributes
    ----------
    index : int
        Unique index of this EC within the cell (1-based for readability).
    umi : str
        UMI sequence.
    features : frozenset[str]
        The set of feature IDs (names) this EC maps to.
    count : int
        Read count for this EC (number of reads supporting this UMI–feature-set).
    read_names : tuple[str, ...] | None
        Optional read names included in this EC (kept for ec_dump provenance).
    """
    index: int
    umi: str
    features: frozenset[str]
    count: int
    read_names: Optional[Tuple[str, ...]] = None

    def hdist(self, other_umi: str) -> int:
        """(directional adjacency, 1MM criterion)"""
        return sum(1 for a, b in zip(self.umi, other_umi) if a != b)

    def connects_to(self, other: "EquivalenceClass", max_hd: int = 1) -> bool:
        """Return True if a directed edge should go self -> other (template->dup).

        Conditions follow IRescue / UMI-tools directional adjacency:
          • hamming(self.umi, other.umi) ≤ max_hd
          • self.count ≥ 2×other.count − 1
          • |self.features ∩ other.features| ≥ 1
        """
        return (
            self.hdist(other.umi) <= max_hd
            and self.count >= (2 * other.count) - 1
            and len(self.features.intersection(other.features)) > 0
        )


def _build_substr_index(eqcs: Sequence[EquivalenceClass]) -> Dict[Tuple[int, int], Dict[str, Set[EquivalenceClass]]]:
    """
    near-linear neighbor generation for 1-mismatch UMI pairing

    Index UMIs by all leave-one-position substrings (L positions → L buckets).
    This lets us retrieve all ECs within Hamming distance ≤1 in ~O(L·N).
    Returns a dict: (start, stop) -> { substring -> set(ECs) }.
    """
    if not eqcs:
        return {}
    L = len(eqcs[0].umi)
    index: Dict[Tuple[int, int], Dict[str, Set[EquivalenceClass]]] = {}
    for i in range(L):
        sl = (0, i)
        sr = (i + 1, L)
        key = (i,)
        # store map from substring (left+right) to ECs
        m: Dict[str, Set[EquivalenceClass]] = {}
        for ec in eqcs:
            sub = ec.umi[: i] + ec.umi[i + 1 :]
            m.setdefault(sub, set()).add(ec)
        index[(0, i, i + 1, L)] = m
    return index


def _candidate_pairs(eqcs: Sequence[EquivalenceClass]) -> Iterable[Tuple[EquivalenceClass, EquivalenceClass]]:
    """
    Yield potentially-1MM pairs using the substring index. Falls back to
    O(N^2) for tiny N.
    """
    n = len(eqcs)
    if n <= 25:
        for i in range(n):
            for j in range(i + 1, n):
                yield eqcs[i], eqcs[j]
        return

    idx = _build_substr_index(eqcs)
    seen: Set[Tuple[int, int]] = set()
    for _, bucket in idx.items():
        for group in bucket.values():
            group_list = list(group)
            for i in range(len(group_list)):
                for j in range(i + 1, len(group_list)):
                    a, b = group_list[i], group_list[j]
                    key = (a.index, b.index) if a.index < b.index else (b.index, a.index)
                    if key not in seen:
                        seen.add(key)
                        yield a, b


def _pathfinder(g: nx.DiGraph, node: int, path: Optional[List[int]] = None, features: Optional[Set[str]] = None) -> List[int]:
    """Pathfinder: recursively follow successors that share ≥1 feature"""
    if path is None:
        path = []
    if features is None:
        features = set(g.nodes[node]["ft"])  # seed with parent's features
    path = path + [node]
    for nxt in g.successors(node):
        if nxt in path:
            continue
        if len(features.intersection(g.nodes[nxt]["ft"])) > 0:
            path = _pathfinder(g, nxt, path, features)
    return path


class Irescue:
    """
    Minimal, class-based re-implementation of **IRescue**'s UMI deduplication
    (directional adjacency with 1MM + EC graph) and multi-mapping EM.

    Designed to run **inside Python** on your DataFrame (no CLI), mirroring the
    logic described in the NAR methods.

    Parameters
    ----------
    cell_col, umi_col, feature_col : str
        Column names in the input DataFrame. `feature_col` can be TE subfamily
        (IRescue’s default) or gene/any feature – your choice.
    read_id_col : str | None, default='read_id'
        Optional column of read names. If absent, each row is treated as a
        single read; you can still deduplicate by (UMI, feature-set) but true
        per-read ECs cannot be reconstructed.
    max_hd : int, default=1
        Max Hamming distance for UMI adjacency.
    em_cycles : int, default=100
        Max iterations of EM redistribution.
    em_tol : float, default=1e-5
        Convergence tolerance of EM.
    dump_ec : bool, default=True
        Whether to build an EC dump DataFrame compatible with IRescue’s
        `ec_dump.tsv.gz` schema.
    no_umi : bool, default=False
        UMI-less mode (e.g., SMART-seq); every read contributes 1 to its
        feature or to EM if multi-mapping.
    """

    def __init__(
        self,
        cell_col: str = "cell",
        umi_col: str = "umi",
        feature_col: str = "feature",
        read_id_col: Optional[str] = "read_id",
        *,
        max_hd: int = 1,
        em_cycles: int = 100,
        em_tol: float = 1e-5,
        dump_ec: bool = True,
        no_umi: bool = False,
    ) -> None:
        self.cell_col = cell_col
        self.umi_col = umi_col
        self.feature_col = feature_col
        self.read_id_col = read_id_col
        self.max_hd = max_hd
        self.em_cycles = em_cycles
        self.em_tol = em_tol
        self.dump_ec = dump_ec
        self.no_umi = no_umi

    def fit_transform(
        self,
        df: pd.DataFrame,
    ) -> Tuple[pd.DataFrame, Optional[pd.DataFrame]]:
        """
        Run UMI deduplication + EM per cell.

        Parameters
        ----------
        df : DataFrame
            Must contain columns for cell / feature. If `no_umi=False`, must
            also contain the UMI column. If `read_id_col` is present, we use it
            to reconstruct true read-level ECs; otherwise we derive ECs directly
            from (UMI, feature) groupings.

        Returns
        -------
        counts_long : DataFrame
            Long-form counts with columns [cell, feature, count] (float).
        ec_dump : DataFrame | None
            If `dump_ec=True`, returns the EC dump with columns matching
            IRescue’s `ec_dump.tsv.gz`. Otherwise `None`.
        """
        required = {self.cell_col, self.feature_col}
        if not self.no_umi:
            required.add(self.umi_col)
        missing = required.difference(df.columns)
        if missing:
            raise ValueError(f"Missing required columns: {sorted(missing)}")

        # Feature indexing (1-based, like IRescue)
        feat_vals = df[self.feature_col].astype(str).unique().tolist()
        feat_index: Dict[str, int] = {f: i for i, f in enumerate(sorted(feat_vals), start=1)}

        counts_accum: Dict[Tuple[str, str], float] = {}
        dump_rows: List[Tuple] = []

        # Process per cell
        for cell, dfc in df.groupby(self.cell_col, sort=False):
            cell_counts, cell_dump = self._process_cell(dfc, feat_index)
            # accumulate counts
            for feat_idx, val in cell_counts.items():
                if val <= 0:
                    continue
                feat_name = self._revindex(feat_index, feat_idx)
                counts_accum[(cell, feat_name)] = counts_accum.get((cell, feat_name), 0.0) + float(val)
            # collect dump
            if self.dump_ec and cell_dump is not None:
                dump_rows.extend([
                    (cell,) + row for row in cell_dump
                ])

        # Build outputs
        counts_long = (
            pd.DataFrame(
                [
                    {self.cell_col: c, self.feature_col: f, "count": v}
                    for (c, f), v in counts_accum.items()
                ]
            )
            .sort_values([self.cell_col, self.feature_col])
            .reset_index(drop=True)
        )

        ec_dump_df: Optional[pd.DataFrame] = None
        if self.dump_ec:
            # Match IRescue column names (Barcode_id as 1..N per cell; we also add cell string)
            # Here we have only the cell string, so we set Barcode_id to NaN; caller can map if needed.
            ec_dump_df = pd.DataFrame(
                dump_rows,
                columns=[
                    self.cell_col,  # not in original ec_dump, but very handy
                    "Barcode_id",  # optional numeric id (set None here)
                    "Barcode",     # the barcode string
                    "EqClass",     # EC index (within cell)
                    "UMI" if not self.no_umi else "Read_name",
                    "Features",
                    "Read_count",
                    "Dedup_UMI" if not self.no_umi else None,
                    "Dedup_feature" if not self.no_umi else None,
                ],
            )
            # Drop Nones in column names for no-UMI mode
            ec_dump_df = ec_dump_df[[c for c in ec_dump_df.columns if c is not None]]

        ec_dump_df = self.summarize_ec_dump(
            ec_dump_df,
            tuple_index=True,
            ave_ed_mode="canonical_pairwise",
        )

        return counts_long, ec_dump_df

    def _revindex(self, idx: Dict[str, int], i: int) -> str:
        # inverse lookup (small dicts per run; safe)
        return next(k for k, v in idx.items() if v == i)

    def _build_ecs(self, dfc: pd.DataFrame) -> List[EquivalenceClass]:
        """Build ECs for a *single cell* DataFrame slice.

        If `read_id_col` exists, we first collate per-read feature sets, then
        group by (UMI, feature_set). Otherwise we use (UMI, feature) directly.
        """
        if self.no_umi:
            # UMI-less: each read is an EC. If read_id missing, treat each row as a read.
            if self.read_id_col and self.read_id_col in dfc.columns:
                g = dfc.groupby(self.read_id_col)
                per_read = (
                    g[self.feature_col].agg(lambda s: frozenset(map(str, s)))
                    .to_frame("features")
                    .assign(count=g.size())
                    .reset_index()
                )
            else:
                per_read = (
                    dfc.assign(_rowid=np.arange(len(dfc)))
                    .groupby("_rowid")[self.feature_col]
                    .agg(lambda s: frozenset(map(str, s)))
                    .to_frame("features")
                    .assign(count=1)
                    .reset_index(drop=True)
                )
            eqcs: List[EquivalenceClass] = []
            for i, row in enumerate(per_read.itertuples(index=False), start=1):
                eqcs.append(EquivalenceClass(i, umi="", features=row.features, count=int(row.count)))
            return eqcs

        # UMI mode
        if self.read_id_col and self.read_id_col in dfc.columns:
            # Build per-read feature set
            rg = dfc.groupby(self.read_id_col)
            read_feats = rg[self.feature_col].agg(lambda s: frozenset(map(str, s)))
            read_meta = rg[[self.umi_col]].first()
            per_read = pd.concat([read_meta, read_feats.rename("features")], axis=1).reset_index()
            # ECs: group reads by (UMI, feature_set)
            eg = per_read.groupby([self.umi_col, "features"], sort=False)
            ecs = eg.size().rename("count").reset_index()
            # collect read names per EC (for dump)
            read_lists = eg[self.read_id_col].apply(tuple).rename("read_names").reset_index(drop=True)
            ecs = pd.concat([ecs, read_lists], axis=1)
        else:
            # Fallback: ECs by (UMI, single feature) – best-effort when read ids missing
            eg = dfc.groupby([self.umi_col, self.feature_col], sort=False)
            ecs = (
                eg.size().rename("count").reset_index()
                .rename(columns={self.feature_col: "features"})
            )
            ecs["features"] = ecs["features"].map(lambda x: frozenset([str(x)]))
            ecs["read_names"] = None

        out: List[EquivalenceClass] = []
        for i, row in enumerate(ecs.itertuples(index=False), start=1):
            out.append(
                EquivalenceClass(
                    index=i,
                    umi=str(row[0]),
                    features=row[1],
                    count=int(row[2]),
                    read_names=row[3] if len(row) > 3 else None,
                )
            )
        return out

    def _process_cell(
        self,
        dfc: pd.DataFrame,
        feat_index: Dict[str, int],
    ) -> Tuple[Dict[int, float], Optional[List[Tuple]]]:
        """Compute counts for a single cell and (optionally) ec_dump rows.

        Returns
        -------
        counts : dict[int -> float]
            Feature-indexed counts (1-based feature ids).
        dump_rows : list[tuple] | None
            When `dump_ec=True`, rows for the ec_dump with columns:
            (Barcode_id, Barcode, EqClass, UMI/Read, FeaturesCSV, Read_count,
             Dedup_UMI, Dedup_feature)
        """
        # Map features to 1..N for EM matrix columns
        all_feats_sorted = sorted(feat_index.keys())
        nfeat = len(all_feats_sorted)
        feat_to_idx = feat_index  # alias

        # Build ECs
        eqcs = self._build_ecs(dfc)

        counts: Dict[int, float] = {i: 0.0 for i in range(1, nfeat + 1)}
        em_rows: List[List[int]] = []
        dump: Dict[int, Tuple] = {}

        # Helper to stringify features
        def feats_to_idx_list(feats: Iterable[str]) -> List[int]:
            return [feat_to_idx[f] for f in feats]

        def feats_to_csv(feats: Iterable[str]) -> str:
            return ",".join(sorted(map(str, feats)))

        # UMI-less path: trivial (each EC contributes 1; multimappers go to EM)
        if self.no_umi:
            for ec in eqcs:
                feats = list(ec.features)
                if len(feats) == 1:
                    counts[feat_to_idx[feats[0]]] += 1.0
                else:
                    row = [1 if (f in feats) else 0 for f in all_feats_sorted]
                    em_rows.append(row)
                if self.dump_ec:
                    dump[ec.index] = (
                        None,  # Barcode_id (optional, fill later in caller)
                        dfc[self.cell_col].iloc[0],  # Barcode (string)
                        ec.index,  # EqClass
                        dfc[self.read_id_col].iloc[0] if self.read_id_col and self.read_id_col in dfc.columns else "",
                        feats_to_csv(ec.features),
                        ec.count,
                    )
        else:
            # Build directional UMI graph over ECs
            g = nx.DiGraph()
            # Node attributes: ft (features set), c (count), umi
            for ec in eqcs:
                g.add_node(ec.index, ft=set(ec.features), c=ec.count, umi=ec.umi)
                if self.dump_ec:
                    dump[ec.index] = (
                        None,  # Barcode_id placeholder
                        dfc[self.cell_col].iloc[0],
                        ec.index,
                        ec.umi,
                        feats_to_csv(ec.features),
                        ec.count,
                    )

            # Edges by directional adjacency with feature intersection
            for a, b in _candidate_pairs(eqcs):
                if a.connects_to(b, self.max_hd):
                    g.add_edge(a.index, b.index)
                elif b.connects_to(a, self.max_hd):
                    g.add_edge(b.index, a.index)

            # Traverse connected components
            subgraphs = [g.subgraph(n).copy() for n in nx.connected_components(g.to_undirected())]
            for subg in subgraphs:
                parents = [n for n in subg.nodes if subg.in_degree(n) == 0]
                if not parents:
                    # bidirected / tie case: treat *all* as parents; union features
                    parents = list(subg.nodes)
                    parent_features: List[List[str]] = [sorted(set().union(*[set(subg.nodes[x]["ft"]) for x in subg.nodes]))]
                else:
                    parent_features = None  # will be taken from the winning path config

                # For each parent, repeatedly find feature-compatible paths, removing nodes
                paths_by_parent: Dict[int, List[List[int]]] = {}
                for p in parents:
                    sub_copy = subg.copy()
                    paths_by_parent[p] = []
                    blacklist: Set[int] = set()
                    while True:
                        # choose next root: first node w/o predecessors that isn't used yet
                        roots = [n for n in sub_copy.nodes if sub_copy.in_degree(n) == 0 and n not in blacklist]
                        if not roots:
                            # if none, take any remaining node (fallback)
                            remaining = [n for n in sub_copy.nodes if n not in blacklist]
                            if not remaining:
                                break
                            root = remaining[0]
                        else:
                            root = roots[0]
                        path = _pathfinder(sub_copy, root, path=[], features=None)
                        for x in path:
                            blacklist.add(x)
                            if sub_copy.has_node(x):
                                sub_copy.remove_node(x)
                        paths_by_parent[p].append(path)

                # Select the path configuration with the *fewest* paths (minimum dedup UMIs)
                min_paths = min(len(v) for v in paths_by_parent.values()) if paths_by_parent else 0
                chosen_parent = next(k for k, v in paths_by_parent.items() if len(v) == min_paths)
                path_config = paths_by_parent[chosen_parent]

                # Determine features per *path* (parent EC features)
                if parent_features is None:
                    parent_features = [sorted(subg.nodes[path[0]]["ft"]) for path in path_config]
                else:
                    parent_features = parent_features * len(path_config)

                # Assign 1 to unique-feature paths, otherwise send row to EM
                for feats in parent_features:
                    if len(feats) == 1:
                        counts[feat_to_idx[feats[0]]] += 1.0
                    elif len(feats) > 1:
                        row = [1 if (f in feats) else 0 for f in all_feats_sorted]
                        em_rows.append(row)

                # For ec_dump: record parent->child relationships (children get parent UMI & parent features)
                if self.dump_ec:
                    for path, feats in zip(path_config, parent_features):
                        # Mark parent with empty dedup fields (root of dedup)
                        parent = path[0]
                        # children inherit parent UMI and features
                        for child in path[1:]:
                            # Append Dedup_UMI and Dedup_feature
                            parent_umi = subg.nodes[parent]["umi"]
                            dedup_feat_csv = feats_to_csv(feats)
                            # Rebuild the stored tuple with dedup fields appended
                            base = dump[child]
                            dump[child] = base + (parent_umi, dedup_feat_csv)
                        # ensure parent has placeholders if not set yet
                        if len(dump[parent]) == 6:  # no dedup fields yet
                            dump[parent] = dump[parent] + ("", "")

        # Run EM if needed
        if em_rows:
            M = np.array(em_rows, dtype=float)
            # drop all-zero columns and remember kept feature indices
            col_any = M.any(axis=0)
            keep_idx = np.where(col_any)[0]
            M = M[:, keep_idx]
            # EM on binary compatibility matrix – return normalized feature weights
            weights, (cycles, converged, loglik, dloglik) = self._run_em(M)
            # Scale by #multimapped UMIs (rows)
            weights = weights * M.shape[0]
            # map back to 1-based feature indices
            for j_in_M, w in enumerate(weights):
                if w <= 0:
                    continue
                feat_name = all_feats_sorted[int(keep_idx[j_in_M])]
                counts[feat_to_idx[feat_name]] += float(w)

        # Build dump rows (compatible with IRescue ec_dump)
        dump_rows: Optional[List[Tuple]] = None
        if self.dump_ec:
            dump_rows = []
            barcode_str = str(dfc[self.cell_col].iloc[0])
            for ecid in sorted(dump.keys()):
                row = dump[ecid]
                # row layout set above; ensure it has dedup fields (parent has "", "")
                if len(row) == 6 and not self.no_umi:
                    row = row + ("", "")
                dump_rows.append((None, barcode_str) + row[2:])  # (Barcode_id, Barcode, EqClass, ...)

        return counts, dump_rows

    def _run_em(self, M: np.ndarray) -> Tuple[np.ndarray, Tuple[int, bool, float, float]]:
        """
        Simple EM identical in spirit to IRescue's
        Run EM on a binary compatibility matrix (rows=UMIs, cols=features).
        Returns feature weights (sum to 1)."""
        # initialize feature counts uniformly
        counts = np.ones(M.shape[1], dtype=float) / M.shape[1]
        prev_loglik = -np.inf
        converged = False
        cycles = 0
        for cycles in range(1, self.em_cycles + 1):
            # E-step: assign each UMI proportionally to current feature counts
            denom = (M * counts).sum(axis=1, keepdims=True)
            # Avoid division by zero
            denom[denom == 0] = 1e-12
            R = (M / denom)  # each row sums to 1 over its compatible features
            # M-step: update feature counts from expected assignments
            counts = R.sum(axis=0)
            counts = counts / counts.sum()  # normalize
            # log-likelihood
            ll = np.log((M * counts).sum(axis=1) + 1e-300).sum()
            dll = ll - prev_loglik
            if abs(dll) < self.em_tol:
                converged = True
                break
            prev_loglik = ll
        return counts, (cycles, converged, float(prev_loglik), float(dll))

    def summarize_ec_dump(
            self,
            df,
            cell_col="BC",
            gene_col="Features",
            umi_col="UMI",
            dedup_umi_col="Dedup_UMI",
            tuple_index=True,
            expand_multifeature=True,
            ave_ed_mode="canonical_pairwise",  # or "raw_to_canon"
    ):
        """
        Attributes
        ----------
        ec : str | pd.DataFrame
            Path to tsv or DataFrame.
        cell_col, gene_col, umi_col, dedup_umi_col : str
            Corresponding column names, which can be modified according to your actual file names (e.g., using 'Barcode' for cell).
        tuple_index : bool
            True: Use (cell, gene) tuple as index (consistent with your example); False: Return explicit two columns.
        expand_multifeature : bool
            If a row's 'Features' contains comma-separated multiple genes, expand into multiple rows and count each gene separately.
        ave_ed_mode : str
            - "canonical_pairwise": Use the pairwise average Hamming distance of the "deduplicated canonical UMI set" as ave_ed;
                                   returns -1.0 if the number of canonical UMIs is less than 2.
            - "raw_to_canon": Only calculate the average Hamming distance between raw UMI and canonical UMI on rows where replacements occurred;
                              returns -1.0 if there are no replacement rows.
        """
        import itertools

        for c in [cell_col, gene_col, umi_col, dedup_umi_col]:
            if c in df.columns:
                df[c] = df[c].astype(str)

        # # Calculate the canonical UMI per row (if Dedup_UMI is empty or missing, take its own UMI)
        df[dedup_umi_col] = df[dedup_umi_col].replace({"nan": np.nan})
        df["_canon_umi"] = df[dedup_umi_col].where(
            df[dedup_umi_col].notna() & (df[dedup_umi_col] != ""),
            df[umi_col]
        ).astype(str)

        # If there are multiple genes (comma-separated), expand as needed
        if expand_multifeature:
            df[gene_col] = df[gene_col].str.split(r"\s*,\s*")
            df = df.explode(gene_col, ignore_index=True)

        def ave_ham_pairwise(umis):
            uniq = list(dict.fromkeys(umis))
            if len(uniq) < 2:
                return -1.0
            # Use the most frequent length to avoid mixing UMIs of different lengths
            from collections import Counter
            L = Counter([len(u) for u in uniq]).most_common(1)[0][0]
            uniq = [u for u in uniq if len(u) == L]
            if len(uniq) < 2:
                return -1.0
            dists = [sum(a != b for a, b in zip(x, y)) for x, y in itertools.combinations(uniq, 2)]
            return float(np.mean(dists)) if dists else -1.0

        def ave_ham_raw_to_canon(sub):
            mask = sub[umi_col] != sub["_canon_umi"]
            if not mask.any():
                return -1.0
            d = []
            for a, b in zip(sub.loc[mask, umi_col], sub.loc[mask, "_canon_umi"]):
                if len(a) == len(b):
                    d.append(sum(x != y for x, y in zip(a, b)))
            return float(np.mean(d)) if d else -1.0

        rows = []
        for (cell, gene), sub in df.groupby([cell_col, gene_col], sort=False):
            raw_uniq = set(sub[umi_col].unique())
            canon_uniq = set(sub["_canon_umi"].unique())
            dedup_cnt = len(canon_uniq)
            num_uniq_umis = len(raw_uniq)
            num_diff_dedup_uniq_umis = len(raw_uniq - canon_uniq)
            num_diff_dedup_reads = int(((sub[dedup_umi_col].notna()) & (sub[dedup_umi_col] != "") & (
                        sub[umi_col] != sub[dedup_umi_col])).sum())
            ave_ed = ave_ham_pairwise(
                list(canon_uniq)) if ave_ed_mode == "canonical_pairwise" else ave_ham_raw_to_canon(sub)
            rows.append((cell, gene, dedup_cnt, ave_ed, num_uniq_umis, num_diff_dedup_uniq_umis, num_diff_dedup_reads))

        out = pd.DataFrame(
            rows,
            columns=[
                "cell",
                "gene",
                "dedup_cnt",
                "ave_ed",
                "num_uniq_umis",
                "num_diff_dedup_uniq_umis",
                "num_diff_dedup_reads",
            ]
        )

        if tuple_index:
            out.index = list(zip(out["cell"], out["gene"]))
            out = out.drop(columns=["cell", "gene"])

        return out


def write_dedup_bam(
    df: pd.DataFrame,
    src_bam_fpn: str,
    out_bam_fpn: str,
    *,
    cell_col: str = "cell",
    umi_col: str = "umi",
    feature_col: str = "feature",
    read_col: str = "read",
    read_id_col: Optional[str] = None,
    max_hd: int = 1,
    no_umi: bool = False,
    mark_only: bool = False,
    duplicate_tag: Optional[str] = None,
    threads: int = 2,
) -> Tuple[int, int]:
    """
    Write a deduplicated BAM file based on IRescue-style UMI collapsing decisions.

    This function reuses the existing IRescue logic to reconstruct, per cell,
    which ECs are canonical (parents) and which are duplicates (children),
    then writes either only representatives (default) or all reads with proper
    duplicate flags/tags.

    Parameters
    ----------
    df : pd.DataFrame
        The same DataFrame you pass into IRescue (must include `cell_col`,
        `feature_col`, and if `no_umi=False`, also `umi_col`). For efficient
        writing, it should contain a column `read_col` with pysam
        AlignedSegment objects (e.g., produced by your ReaderChunk with
        bam_fields including "read").
    src_bam_fpn : str
        Path to the source BAM whose header will be reused.
    out_bam_fpn : str
        Path to the output BAM to be written.
    cell_col, umi_col, feature_col : str
        Column names for cell/UMI/feature in `df`.
    read_col : str
        Column name holding pysam AlignedSegment objects.
    read_id_col : Optional[str]
        If you originally ran IRescue with per-read grouping, provide the
        read id column name here to be consistent with EC building; otherwise
        leave as None (the fallback EC definition (UMI, feature) is used).
    max_hd : int
        Max Hamming distance used by directional adjacency (should match your
        IRescue instance).
    no_umi : bool
        UMI-less mode. When True, representatives are chosen per (cell, feature).
    mark_only : bool
        When False (default): write **only** representative reads (hard
        deduplicated BAM). When True: write **all** reads but mark duplicates
        via the duplicate flag and optional tag.
    duplicate_tag : Optional[str]
        Optional auxiliary tag to set on reads (0=not duplicate, 1=duplicate).
        If provided (e.g., "PD"), the tag will be stored as an int.
    threads : int
        Number of threads for BAM writing.

    Returns
    -------
    kept_count, duplicate_count : Tuple[int, int]
        Numbers of representatives written and duplicates handled
        (either dropped or marked, depending on `mark_only`).

    Notes
    -----
    - Canonical ECs (parents; i.e., rows with empty `Dedup_UMI`) define which
      (cell, UMI, feature) triplets are kept. For multi-feature ECs, this
      function keeps one representative per feature present in `df` for that
      canonical UMI in the given cell (typical single-feature per row is fully
      supported).
    - If `df` does not contain a `read_col` with AlignedSegment objects, this
      routine cannot write reads efficiently and will raise a ValueError.
    """
    # Local imports to avoid modifying global imports in your file
    import pysam
    from typing import Optional as _Optional, Tuple as _Tuple, Set as _Set

    # Basic sanity checks
    if read_col not in df.columns:
        raise ValueError(
            f"`{read_col}` column with pysam.AlignedSegment objects is required."
        )
    for col in [cell_col, feature_col] + ([] if no_umi else [umi_col]):
        if col not in df.columns:
            raise ValueError(f"Missing required column: {col}")

    # Open source BAM for header, and allocate writer
    src = pysam.AlignmentFile(src_bam_fpn, "rb")
    out = pysam.AlignmentFile(out_bam_fpn, "wb", header=src.header, threads=threads)

    # Instantiate your existing IRescue to reuse its EC logic
    ir = Irescue(
        cell_col=cell_col,
        umi_col=umi_col,
        feature_col=feature_col,
        read_id_col=read_id_col,
        max_hd=max_hd,
        dump_ec=True,
        no_umi=no_umi,
    )

    kept = 0
    dups = 0

    # Process per cell to build canonical EC keys
    for cell, dfc in df.groupby(cell_col, sort=False):
        # Build per-cell feature index as in the core logic
        feats = sorted(dfc[feature_col].astype(str).unique())
        feat_index = {f: i for i, f in enumerate(feats, start=1)}

        # Use the existing internal method to get EC dump rows for this cell
        _, dump_rows = ir._process_cell(dfc, feat_index)

        # Collect canonical (parent) triplets as (cell, umi, feature)
        canonical_keys: _Set[_Tuple[str, str, str]] = set()
        if dump_rows:
            for row in dump_rows:
                # Expected layout in UMI mode:
                # (Barcode_id, Barcode, EqClass, UMI, FeaturesCSV, Read_count, Dedup_UMI, Dedup_feature)
                # Expected layout in no-UMI mode:
                # (Barcode_id, Barcode, EqClass, Read_name, FeaturesCSV, Read_count)
                barcode = str(row[1])
                features_csv = str(row[4])
                if not no_umi:
                    ec_umi = str(row[3])
                    dedup_umi = str(row[6]) if len(row) > 6 else ""
                    # A parent EC is identified by an empty Dedup_UMI
                    if dedup_umi == "":
                        feats_here = [f.strip() for f in features_csv.split(",") if f.strip() != ""]
                        for f in feats_here:
                            canonical_keys.add((barcode, ec_umi, f))
                else:
                    # UMI-less mode: keep one representative per (cell, feature)
                    feats_here = [f.strip() for f in features_csv.split(",") if f.strip() != ""]
                    for f in feats_here:
                        canonical_keys.add((barcode, "", f))

        # Iterate using tuple access to avoid pandas attribute-name sanitization
        if no_umi:
            selected_cols = [cell_col, feature_col, read_col]
        else:
            selected_cols = [cell_col, umi_col, feature_col, read_col]

        seen_rep: _Set[_Tuple[str, str, str]] = set()

        for rec in dfc[selected_cols].itertuples(index=False, name=None):
            if no_umi:
                barcode, feat, read = rec
                umi = ""
            else:
                barcode, umi, feat, read = rec

            # Normalize key components to strings for consistent matching
            barcode = str(barcode)
            feat = str(feat)
            umi = "" if no_umi else str(umi)

            key = (barcode, umi, feat)

            if key in canonical_keys:
                # Keep ONE representative per canonical triplet
                if key not in seen_rep:
                    if mark_only:
                        # Clear duplicate flag and optionally set a tag = 0
                        read.is_duplicate = False
                        if duplicate_tag is not None:
                            read.set_tag(duplicate_tag, 0, value_type="i", replace=True)
                        out.write(read)
                    else:
                        # Write only once
                        if duplicate_tag is not None:
                            read.set_tag(duplicate_tag, 0, value_type="i", replace=True)
                        out.write(read)
                    seen_rep.add(key)
                    kept += 1
                else:
                    # Additional reads within the same canonical EC are duplicates
                    if mark_only:
                        read.is_duplicate = True
                        if duplicate_tag is not None:
                            read.set_tag(duplicate_tag, 1, value_type="i", replace=True)
                        out.write(read)
                    dups += 1
            else:
                # Child EC (collapsed) or non-canonical: treat as duplicate
                if mark_only:
                    read.is_duplicate = True
                    if duplicate_tag is not None:
                        read.set_tag(duplicate_tag, 1, value_type="i", replace=True)
                    out.write(read)
                dups += 1

    out.close()
    src.close()
    # If you want an index, you can enable it where appropriate:
    # pysam.index(out_bam_fpn)

    return kept, dups



if __name__ == "__main__":
    # df contains: cell, feature, umi
    # e.g. df = pd.DataFrame({"cell":[...],"umi":[...],"feature":[...],"read_id":[...]} )
    from umiche.bam.Reader import ReaderChunk
    from umiche.util.Console import Console
    bam_fpn = "/mnt/d/Document/Programming/Python/umiche/umiche/data/r1/umicountr/Smartseq3.TTACCTGCCAGATTCG.bam"
    df_bam = ReaderChunk(
        bam_fpn=bam_fpn,
        bam_fields=['contig', 'pos', 'CIGAR', 'seq', 'read'],
        tag_whitelist=['BC', 'QU', 'UX', 'UB', 'GE'],
        verbose=True,
    ).todf(chunk_size=1_000_000)
    # print(df)
    spikecontig = "diySpike"
    spikename = "g_diySpike4"
    df_bam = df_bam[df_bam["UX"] != ""]
    df_bam = df_bam[df_bam["GE"] == spikename]
    df_bam = df_bam.reset_index(drop=True)
    df_bam = df_bam.drop("GE", axis=1)
    # print(df_bam)
    from umiche.deduplicate.spikein.SpikeUMI import SpikeUMI

    df_bam = SpikeUMI(verbose=True).extract_spike_dat(
        df=df_bam,
        match_seq_before_UMI="GAGCCTGGGGGAACAGGTAGG",
        match_seq_after_UMI="CTCGGAGGAGAAA",
    )
    Console(True).df_column_summary(df_bam)
    from umiche.deduplicate.method.trimer.Expand import Expand

    df_bam["htUMI"] = df_bam["UX"].apply(Expand().homotrimer_umi)
    mut_rate = 0.01
    from umiche.deduplicate.method.trimer.Error import Error

    df_bam["htUMI_" + str(mut_rate)] = df_bam["htUMI"].apply(lambda umi: Error().mutated(umi, mut_rate=mut_rate, mode="normal"))

    dedup = Irescue(
        cell_col="BC",
        umi_col="htUMI_" + str(mut_rate),
        feature_col="spikeUMI",
        read_id_col=None,
        max_hd=1,
        em_cycles=100,
        em_tol=1e-5,
        dump_ec=True,
        no_umi=False, # SMART-seq can be set as True
    )

    counts_long, ec_dump = dedup.fit_transform(df_bam)

    print(ec_dump)

    ec_dump.to_csv(
        '/mnt/d/Document/Programming/Python/umiche/umiche/data/r1/trumicount/10xn9k_10c/irescur.txt',
        header=True,
        sep="\t",
        index=True,
    )

    # # for matrix (cell × feature), please do pivot：
    # matrix = counts_long.pivot(index="feature", columns="cell", values="count").fillna(0.0)

    kept, dups = write_dedup_bam(
        df=df_bam,
        src_bam_fpn=bam_fpn,
        out_bam_fpn="/mnt/d/Document/Programming/Python/umiche/umiche/data/r1/umicountr/ccc.bam",
        cell_col="BC",
        umi_col="htUMI_" + str(mut_rate),
        feature_col="spikeUMI",
        read_col="read",
        read_id_col=None,
        max_hd=1,
        no_umi=False,
        mark_only=False,          # True -> keep all reads and mark duplicates instead of dropping
        duplicate_tag="PD",       # Optional tag to annotate duplicates (0/1)
        threads=4,
    )
    print("wrote dedup BAM:", kept, "kept,", dups, "duplicates")
