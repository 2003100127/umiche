__version__ = "v1.0"
__copyright__ = "Copyright 2025"
__license__ = "GPL-3.0"
__developer__ = "Jianfeng Sun"
__maintainer__ = "Jianfeng Sun"
__email__="jianfeng.sunmt@gmail.com"

from typing import Dict, Iterable, List, Tuple, Set, Optional

import itertools

from collections import defaultdict, Counter
from umiche.util.Console import Console

# ---------- Utilities ----------
def hamming(a: str, b: str) -> int:
    """Hamming distance for equal-length UMI strings."""
    if len(a) != len(b):
        raise ValueError("UMIs must have equal length for Hamming distance.")
    return sum(1 for x, y in zip(a, b) if x != y)

class DSU:
    """Disjoint Set Union (Union-Find) for connected components."""
    def __init__(self): self.p, self.r = {}, {}
    def find(self, x):
        if x not in self.p: self.p[x]=x; self.r[x]=0; return x
        if self.p[x]!=x: self.p[x]=self.find(self.p[x])
        return self.p[x]
    def union(self, x, y):
        rx, ry = self.find(x), self.find(y)
        if rx==ry: return
        if self.r[rx]<self.r[ry]: self.p[rx]=ry
        elif self.r[rx]>self.r[ry]: self.p[ry]=rx
        else: self.p[ry]=rx; self.r[rx]+=1

# ---------- Main class ----------
class Gencore:
    """
    Gencore-style UMI clustering: within the SAME coordinate group, connect UMIs if Hamming(u,v) ≤ t.
    Connected components are clusters.

    Inputs supported:
      1) DataFrame: auto-group by 'same coordinate', then build graph & cluster (includes coordinate judgment).
      2) Adjacency list: assumed to be a single 'same coordinate' subgraph; no coordinate judgment.

    Unified return format (pipeline style): {"count","clusters","apv","disapv","assigned"}
      - count: number of clusters
      - clusters: [{"labels":[...], "consensus_umi":str, "size":int}, ...]  # labels = node IDs / UMIs
      - apv: accepted undirected edges within clusters [("u","v"), ...]
      - disapv: rejected/cross-cluster edges (usually empty with this definition)
      - assigned: {label: cluster_id}
    """
    def __init__(
            self, t: int = 1, alphabet: str = "ACGTN",
            build_consensus: bool = True,
            verbose=False,
    ):
        self.t = t
        self.alphabet = alphabet
        self.build_consensus = build_consensus

        self.console = Console()
        self.console.verbose = verbose

    def fit_graph(
        self,
            graph_adj: Dict[str, Iterable[str]],   # {label: [neighbor_label,...]} already within the same coordinate
            connected_components,
            int_to_umi_dict: Dict[str, str],             # node label -> real UMI
            counts: Optional[Dict[str, int]] = None,  # {label: count} (optional; default=1)
            validate: bool = False,
            t: Optional[int] = None,
    ) -> Dict:
        """
        Cluster from a pre-built graph_adj list (assumed same coordinate group).
        Optionally validate each edge respects Hamming ≤ t on real UMIs.
        """
        adj = self._make_undirected(graph_adj)

        # Optional: verify edges respect the distance threshold on the REAL UMI sequences
        if validate:
            th = self.t if t is None else t
            for u, nbrs in adj.items():
                uu = int_to_umi_dict[u]
                for v in nbrs:
                    uv = int_to_umi_dict[v]
                    if len(uu) != len(uv): raise ValueError(f"UMI length mismatch: {u}={uu}, {v}={uv}")
                    if hamming(uu, uv) > th: raise ValueError(f"Edge violates distance ≤ {th}: {u}({uu})-{v}({uv})")

        # Connected components on labels
        _, assigned = self._connected_components(adj)

        clusters_sets = [set(cc) for cc in connected_components.values()]


        # Per-cluster consensus and size (weighted by label_counts if provided)
        if counts is None: label_counts = {lab: 1 for lab in adj.keys()}
        else: label_counts = {lab: int(counts.get(lab, 0)) for lab in adj.keys()}

        clusters = []
        for cid, members in enumerate(clusters_sets):
            members_sorted = sorted(members)
            umi_cnt = Counter()
            for lab in members:
                umi_cnt[int_to_umi_dict[lab]] += label_counts.get(lab, 0)
            cons = self._consensus_umi(umi_cnt) if self.build_consensus else None
            size = int(sum(label_counts.get(lab, 0) for lab in members)) or len(members)
            clusters.append({"labels": members_sorted, "consensus_umi": cons, "size": size})

        apv = self._edges_within_clusters(adj, assigned)
        disapv = []  # Usually empty for this clustering definition
        clusters = {i: clus['labels'] for i, clus in enumerate(clusters)}
        return {
            "count": len(clusters),
            "clusters": clusters,
            # "apv": apv,
            # "disapv": disapv,
            # "assigned": assigned,
        }

    # ---------------- Internals: graph build / components / consensus / edge collection ----------------
    def _build_graph_from_umis(self, umi_counts: Dict[str, int], t: int) -> Dict[str, Set[str]]:
        """
        Build a UMI graph within one coordinate group:
          - Nodes: UMIs
          - Edge(u,v): Hamming(u,v) ≤ t
        For t=1, use neighbor enumeration; for t≥2, use mask blocking to reduce candidates.
        """
        umis = list(umi_counts.keys())
        if not umis: return {}
        L = len(umis[0])
        if any(len(u)!=L for u in umis):
            raise ValueError("All UMIs in the same group must have equal length.")

        adj = {u:set() for u in umis}
        if t<=0: return adj

        # Fast path for t=1: enumerate single-base neighbors
        if t==1:
            s = set(umis)
            for u in umis:
                for i in range(L):
                    ori = u[i]
                    for b in self.alphabet:
                        if b==ori: continue
                        v = u[:i]+b+u[i+1:]
                        if v in s:
                            adj[u].add(v); adj[v].add(u)
            return adj

        # t>=2: mask blocking to limit candidate pairs, then verify distance
        positions = range(L)
        masks = list(itertools.combinations(positions, t))
        buckets = defaultdict(list)
        for u in umis:
            for mt in masks:
                arr = list(u)
                for p in mt: arr[p]="*"
                buckets["".join(arr)].append(u)
        seen = set()
        for bucket in buckets.values():
            if len(bucket)<2: continue
            for i in range(len(bucket)):
                ui = bucket[i]
                for j in range(i+1, len(bucket)):
                    uj = bucket[j]
                    a,b = (ui,uj) if ui<uj else (uj,ui)
                    if (a,b) in seen: continue
                    if hamming(ui, uj) <= t:
                        adj[ui].add(uj); adj[uj].add(ui); seen.add((a,b))
        return adj

    def _connected_components(self, adj: Dict[str, Set[str]]) -> Tuple[List[Set[str]], Dict[str,int]]:
        """Return list of components (as sets of labels) and a label->cluster_id mapping."""
        dsu = DSU()
        for u, nbrs in adj.items():
            dsu.find(u)
            for v in nbrs: dsu.union(u,v)
        comps = defaultdict(list)
        for u in adj.keys(): comps[dsu.find(u)].append(u)
        clusters_sets, assigned = [], {}
        for cid, (_, members) in enumerate(comps.items()):
            clusters_sets.append(set(members))
            for m in members: assigned[m]=cid
        return clusters_sets, assigned

    def _consensus_umi(self, umi_cnt: Counter) -> str:
        """
        Weighted per-position majority vote over A/C/G/T/N.
        Ties are broken by lexicographic order for determinism.
        """
        if not umi_cnt: return ""
        umis = list(umi_cnt.keys()); L = len(umis[0])
        cols = [Counter() for _ in range(L)]
        for u,w in umi_cnt.items():
            for i,ch in enumerate(u): cols[i][ch]+=w
        out=[]
        for c in cols:
            if not c: out.append("N"); continue
            m = max(c.values()); winners = sorted([b for b,v in c.items() if v==m])
            out.append(winners[0])
        return "".join(out)

    def _make_undirected(self, graph_adj: Dict[str, Iterable[str]]) -> Dict[str, Set[str]]:
        """Normalize to an undirected graph_adj (set-based) and ensure all nodes appear."""
        adj = defaultdict(set)
        for u, vs in graph_adj.items():
            for v in vs: adj[u].add(v); adj[v].add(u)
            if u not in adj: adj[u]=set()
        return adj

    def _edges_within_clusters(self, adj: Dict[str, Set[str]], assigned: Dict[str,int]) -> List[Tuple[str,str]]:
        """Collect unique undirected edges that fall within the same cluster."""
        apv_set=set()
        for u, nbrs in adj.items():
            for v in nbrs:
                if u==v: continue
                a,b = (u,v) if u<v else (v,u)
                if (a,b) in apv_set: continue
                if assigned[u]==assigned[v]: apv_set.add((a,b))
        return sorted(apv_set)


if __name__ == "__main__":
    from umiche.deduplicate.method.Cluster import Cluster as umiclust

    # Adjacency input example (assumed same-coordinate subgraph)
    graph_adj = {
        'A': ['B', 'C', 'D'],
        'B': ['A', 'C'],
        'C': ['A', 'B'],
        'D': ['A', 'E', 'F'],
        'E': ['D', 'G'],
        'F': ['D', 'G'],
        'G': ['E', 'F'],
    }
    label_counts = {'A':120, 'D':90, 'E':50, 'G':5, 'B':2, 'C':2, 'F':1}
    int_to_umi_dict = {
        'A': 'AGATCTCGCA','B': 'AGATCCCGCA','C': 'AGATCACGCA',
        'D': 'AGATCTGGCA','E': 'AGATCTGGGA','F': 'AGATCTGGCT','G': 'AGATCTGGGT',
    }

    ccs = umiclust().cc(graph_adj=graph_adj)
    print("Connected components:\n{}".format(ccs))

    deduper = Gencore(t=1, build_consensus=True)
    res_graph = deduper.fit_graph(
        graph_adj=graph_adj,
        int_to_umi_dict=int_to_umi_dict,
        connected_components=ccs,
        counts=label_counts,
        validate=True,
    )
    from pprint import pprint
    pprint(res_graph)

