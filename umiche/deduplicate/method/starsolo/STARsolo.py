__version__ = "v1.0"
__copyright__ = "Copyright 2025"
__license__ = "GPL-3.0"
__developer__ = "Jianfeng Sun"
__maintainer__ = "Jianfeng Sun"
__email__="jianfeng.sunmt@gmail.com"

from typing import Dict, List, Iterable, Tuple, Optional, Any

import pandas as pd

from collections import defaultdict
from umiche.util.Console import Console


class STARsolo:

    def __init__(
            self,
            tie_break_lex=True,
            verbose=False,
    ):
        # When selecting "best UMI" inside a connected component:
        #  - primary key: highest count
        #  - tie-break: lexicographically smallest if tie_break_lex is True
        self.tie_break_lex = tie_break_lex
        self.console = Console()
        self.console.verbose = verbose

    def dedup_from_dataframe(
        self,
        df: pd.DataFrame,
        umi_col: str,
        group_cols: Optional[List[str]] = None,
        count_col: Optional[str] = None,
        dropna_umi: bool = True,
        return_per_group: bool = True,
    ) -> Dict[str, Any]:
        """
        Deduplicate UMIs from a DataFrame.

        Parameters
        ----------
        df : DataFrame with at least the UMI column. Each row is a read.
        umi_col : name of the UMI column (string of equal length; non-ACGT allowed but length must match within group)
        group_cols : list of column names to define grouping (e.g., ["cell", "gene"]). If None -> single global group.
        count_col : if provided, it is a per-row (read) weight; otherwise row count per UMI will be used as counts.
        dropna_umi : drop rows with NA/empty UMI.
        return_per_group : return per-group results under results["groups"].

        Returns
        -------
        dict with:
          overall, and if return_per_group=True then results["groups"][group_key] holds the same schema per group.
        """
        # basic input clean
        x = df
        if dropna_umi:
            x = x[x[umi_col].notna()]
        x = x.copy()

        if group_cols is None or len(group_cols) == 0:
            group_cols = ["__ALL__"]
            x["__ALL__"] = "ALL"

        results = {
            "overall": {
                "n_groups": 0,
                "n_before_sum": 0,
                "n_after_1mm_sum": 0,
                "n_after_directional_sum": 0,
            },
            "groups": {} if return_per_group else None,
        }

        for gkey, sub in x.groupby(group_cols, sort=False):
            # build counts per UMI
            if count_col is None:
                counts = sub.groupby(umi_col, sort=False).size().to_dict()
            else:
                counts = sub.groupby(umi_col, sort=False)[count_col].sum().to_dict()

            # adjacency by Hamming-1
            adjacency = self._build_hamming1_graph(list(counts.keys()))

            res = self.dedup_from_graph(adjacency=adjacency, counts=counts)

            # collect
            results["overall"]["n_groups"] += 1
            results["overall"]["n_before_sum"] += res["n_before"]
            results["overall"]["n_after_1mm_sum"] += res["n_after_1mm"]
            results["overall"]["n_after_directional_sum"] += res["n_after_directional"]

            if return_per_group:
                results["groups"][gkey] = res

        return results

    def dedup_from_graph(
            self,
            graph_adj: Dict[str, List[str]],
            counts: Dict[str, int],
            connected_components,
    ) -> Dict[str, Any]:
        """
        Deduplicate from an already-built 1MM graph_adj and UMI counts.

        Returns a dict (see class docstring).
        """
        # Normalize: ensure every node exists in graph_adj and counts
        nodes = set(counts.keys())
        adj = {u: list(set(v)) for u, v in graph_adj.items() if u in nodes}
        for u in nodes:
            if u not in adj:
                adj[u] = []

        # Connected components on 1MM graph (equivalent to final graph-coloring CCs)
        # components = self._connected_components(adj)
        # print(components)
        connected_components = [cc for cc in connected_components.values()]

        # Winner (best UMI) per component, and umi->winner mapping
        winner_per_component, map_to_comp_winner = self._pick_winners_and_map(connected_components, counts)

        # Directional collapse inside each component
        dir_clusters, map_dir, approvals, disapprovals = self._directional_clusters(adj, counts, connected_components)

        dir_clusters_ = {i: clus for i, clus in enumerate(dir_clusters)}

        return {
            "count": len(dir_clusters),
            "clusters": dir_clusters_,

            # "n_before": len(nodes),
            # "n_after_1mm": len(components),
            # "components": components,
            # "winner_per_component": winner_per_component,
            # "map_to_component_winner": map_to_comp_winner,
            # "map_directional": map_dir,
            # "approval_edges": approvals,
            # "disapproval_edges": disapprovals,
        }

    @staticmethod
    def _hamming1_buckets(umis: Iterable[str]) -> Dict[Tuple[int, str], List[str]]:
        """
        Build wildcard buckets to find Hamming-distance-1 pairs in O(N*L).
        Bucket key: (position, string with position replaced by '*')
        """
        buckets: Dict[Tuple[int, str], List[str]] = defaultdict(list)
        umis = list(umis)
        if not umis:
            return buckets
        L = len(umis[0])
        for u in umis:
            # skip different lengths (won't match any way)
            if len(u) != L:
                continue
            for i in range(L):
                key = (i, u[:i] + "*" + u[i+1:])
                buckets[key].append(u)
        return buckets

    def _build_hamming1_graph(self, umis: List[str]) -> Dict[str, List[str]]:
        """
        Build adjacency for all pairs with Hamming distance exactly 1.
        """
        adj: Dict[str, set] = {u: set() for u in umis}
        if not umis:
            return {u: [] for u in umis}
        L = len(umis[0])

        # group by wildcard buckets, and then connect pairs in each bucket that truly have Hamming distance=1
        buckets = self._hamming1_buckets(umis)
        for (_, wildcard), bucket_nodes in buckets.items():
            bn = bucket_nodes
            if len(bn) < 2:
                continue
            # connect all pairs that differ at exactly that wildcard position and are equal at others -> Hamming 1
            # since the bucket is per position, any two different UMIs in the same bucket have Hamming==1
            # but identical UMIs can appear; skip identical
            seen = set()
            for i, u in enumerate(bn):
                seen.add((i, u))
            # naive pair linking in this bucket
            for i in range(len(bn)):
                ui = bn[i]
                for j in range(i+1, len(bn)):
                    uj = bn[j]
                    if len(ui) != L or len(uj) != L:
                        continue
                    if ui != uj:
                        adj[ui].add(uj)
                        adj[uj].add(ui)

        return {u: sorted(v) for u, v in adj.items()}

    def _pick_winners_and_map(
        self,
        components: List[List[str]],
        counts: Dict[str, int],
    ) -> Tuple[List[str], Dict[str, str]]:
        winners: List[str] = []
        mapping: Dict[str, str] = {}

        for comp in components:
            # choose best by count, tie-break lexicographically if needed
            comp_sorted = sorted(comp, key=lambda u: (-counts.get(u, 0), u if self.tie_break_lex else ""))
            winner = comp_sorted[0]
            winners.append(winner)
            for u in comp:
                mapping[u] = winner
        return winners, mapping

    def _directional_clusters(
        self,
        adjacency: Dict[str, List[str]],
        counts: Dict[str, int],
        components: Optional[List[List[str]]] = None,
    ) -> Tuple[List[List[str]], Dict[str, str], Dict[str, List[Tuple[str, str]]], Dict[str, List[Tuple[str, str]]]]:
        """
        Directional collapse: recursively absorb neighbors n with
            count[root] >= 2*count[n] - 1
        starting from highest-count nodes within each CC.

        Returns clusters (list of nodes per root), loser->winner map, and approval/disapproval edges per root.
        """
        # if components is None:
        #     components = self._connected_components(adjacency)

        # order nodes within each CC by decreasing count
        clusters: List[List[str]] = []
        map_dir: Dict[str, str] = {}
        approvals: Dict[str, List[Tuple[str, str]]] = {}
        disapprovals: Dict[str, List[Tuple[str, str]]] = {}

        for comp in components:
            comp_counts = {u: counts.get(u, 0) for u in comp}
            # processing order (descending counts)
            order = sorted(comp, key=lambda u: (-comp_counts[u], u if self.tie_break_lex else ""))
            remaining = set(comp)

            for root in order:
                if root not in remaining:
                    continue
                # grow cluster from root using DFS with the 2*c-1 rule
                visited = set([root])
                approves: List[Tuple[str, str]] = []
                rejects: List[Tuple[str, str]] = []

                def dfs(u: str):
                    for v in adjacency[u]:
                        if v in visited:
                            continue
                        if v not in remaining:
                            continue
                        # apply directional rule exactly as in C++: if c_big >= 2*c_small -1
                        cu, cv = comp_counts.get(u, 0), comp_counts.get(v, 0)
                        # Try both directions, matching C++ symmetric check:
                        #  - if u not flagged AND cu > 2*cv -1 -> v duplicates (approve edge u->v)
                        #  - else if v not flagged AND cv > 2*cu -1 -> u duplicates (approve edge v->u)
                        # In graph-based DFS, we expand from the "winner side".
                        if cu >= 2 * cv - 1:
                            approves.append((u, v))
                            visited.add(v)
                            dfs(v)
                        elif cv >= 2 * cu - 1:
                            # v can absorb u — but since we expand from root, just mark rejection here
                            rejects.append((u, v))
                        else:
                            rejects.append((u, v))

                dfs(root)

                # finalize this directional cluster
                cluster_nodes = sorted(list(visited))
                clusters.append(cluster_nodes)

                # map losers to root (winner); keep winner->winner as well for completeness
                for u in cluster_nodes:
                    if u != root:
                        map_dir[u] = root
                    else:
                        map_dir[u] = root

                approvals[root] = approves
                disapprovals[root] = rejects

                # remove clustered nodes from remaining
                remaining -= visited

        return clusters, map_dir, approvals, disapprovals


if __name__ == "__main__":
    import pandas as pd
    from umiche.deduplicate.method.Cluster import Cluster as umiclust

    p = STARsolo()

    ### @@ data from UMI-tools
    # graph_adj = {
    #     'A': ['B', 'C', 'D'],
    #     'B': ['A', 'C'],
    #     'C': ['A', 'B'],
    #     'D': ['A', 'E', 'F'],
    #     'E': ['D'],
    #     'F': ['D'],
    # }
    # print("An adjacency list of a graph:\n{}".format(graph_adj))
    #
    # node_val_sorted = pd.Series({
    #     'A': 456,
    #     'E': 90,
    #     'D': 72,
    #     'B': 2,
    #     'C': 2,
    #     'F': 1,
    # })
    # print("Counts sorted:\n{}".format(node_val_sorted))

    ### @@ data from mine
    graph_adj = {
        'A': ['B', 'C', 'D'],
        'B': ['A', 'C'],
        'C': ['A', 'B'],
        'D': ['A', 'E', 'F'],
        'E': ['D', 'G'],
        'F': ['D', 'G'],
        'G': ['E', 'F'],
    }
    print("An adjacency list of a graph:\n{}".format(graph_adj))

    node_val_sorted = pd.Series({
        'A': 120,
        'D': 90,
        'E': 50,
        'G': 5,
        'B': 2,
        'C': 2,
        'F': 1,
    })

    ccs = umiclust().cc(graph_adj=graph_adj)
    print("Connected components:\n{}".format(ccs))

    res = p.dedup_from_graph(graph_adj=graph_adj, counts=node_val_sorted.to_dict(), connected_components=ccs)

    print(res)
