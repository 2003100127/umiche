from pathlib import Path
from typing import Dict, Iterable
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt


def _read_dedup_cnt_one(path: str | Path, sep: str = "\t") -> pd.DataFrame:
    """
    Reads a single 'reads_with_tags_<er>_dedup_sum' file and returns two columns:
    ['transcript_id', 'dedup_cnt']. Explicit ID columns (transcript_id/target_id) are preferred;
    otherwise, the index or first column is used instead.
    """
    path = Path(path)

    # Read normally first
    df = pd.read_csv(path, sep=sep)

    # Find the ID column
    id_col = None
    for c in ["transcript_id", "target_id"]:
        if c in df.columns:
            id_col = c
            break

    # If dedup_cnt is not read, try to read the first column as the index again
    if "dedup_cnt" not in df.columns:
        df = pd.read_csv(path, sep=sep, index_col=0)
        # Try to get the id from the column or index
        if id_col is None:
            if df.index.name not in (None, ""):
                id_col = None # Use index
            elif len(df.columns) > 0 and df.columns[0] != "dedup_cnt":
                id_col = df.columns[0]

    # Unified structure output
    if "dedup_cnt" not in df.columns:
        # One step back: whitespace separation + first column index
        df = pd.read_csv(path, sep=r"\s+", engine="python", index_col=0)

    out = pd.DataFrame()
    out["dedup_cnt"] = pd.to_numeric(df["dedup_cnt"], errors="coerce")

    if id_col is not None and id_col in df.columns:
        out["transcript_id"] = df[id_col].astype(str)
    else:
        # Use index as ID (only makes sense if index is not 0..N)
        if not isinstance(df.index, pd.RangeIndex):
            out["transcript_id"] = df.index.astype(str)
        else:
            # Finally, find a non-dedup_cnt column as the ID
            candidates = [c for c in df.columns if c != "dedup_cnt"]
            if not candidates:
                raise ValueError(f"Cannot find transcript_id in: {path}")
            out["transcript_id"] = df[candidates[0]].astype(str)

    out = out.dropna(subset=["dedup_cnt"])
    return out[["transcript_id", "dedup_cnt"]]


def load_methods_errorrates(
    filemap: Dict[str, Dict[float, str]],
    sep: str = "\t",
) -> pd.DataFrame:
    """
    Read the mapping {method: {error_rate: filepath}} into a unified DataFrame.
    Return columns: ['method', 'error_rate', 'transcript_id', 'dedup_cnt'].
    """
    rows = []
    for method, er_dict in filemap.items():
        for er, fpn in er_dict.items():
            tmp = _read_dedup_cnt_one(fpn, sep=sep)
            tmp["method"] = method
            tmp["error_rate"] = float(er)
            rows.append(tmp)
    if not rows:
        raise RuntimeError("No files loaded. Check your filemap.")
    df = pd.concat(rows, ignore_index=True)
    return df[["method", "error_rate", "transcript_id", "dedup_cnt"]]


def load_methods_errorrates_with_relerr(
    filemap: dict, abundance_path: str | Path, sep: str = "\t", absolute: bool = False
) -> pd.DataFrame:
    """
    Wrapper: First use your existing load_methods_errorrates to read the full data,
    then apply relative error replacement. Do not change your original function
    signature/implementation, minimally invasive access.
    """
    df_all = load_methods_errorrates(filemap, sep=sep)
    print(df_all)
    ab = load_abundance(abundance_path)
    print(ab)
    df_all = apply_tpm_relative_error(df_all, ab, absolute=absolute)
    # print(df_all)
    return df_all


def load_abundance(abundance_path: str | Path) -> pd.DataFrame:
    """
    Reads the ground truth abundance.tsv and returns two columns: ['transcript_id', 'tpm'].
    Compatible with the table headers 'target_id' and 'tpm'/'TPM'; ignores 'cell'.
    """
    ab = pd.read_csv(abundance_path, sep="\t")
    # print(ab)
    ab = ab.rename(columns={
        "target_id": "transcript_id",
        "TPM": "tpm",
        "tpm": "tpm",
    })
    need = {"transcript_id", "tpm"}
    if not need.issubset(ab.columns):
        raise ValueError(f"abundance file must contain {need}, got {ab.columns.tolist()}")
    ab = ab[["transcript_id", "tpm"]].copy()
    ab["transcript_id"] = ab["transcript_id"].astype(str)
    ab["tpm"] = pd.to_numeric(ab["tpm"], errors="coerce")
    return ab


def apply_tpm_relative_error(
    df_all: pd.DataFrame,
    ab_df: pd.DataFrame,
    on: str = "transcript_id",
    col: str = "dedup_cnt",
    absolute: bool = False,
    keep_original: bool = True,
) -> pd.DataFrame:
    """
    Calculate a new dedup_cnt using the ground truth and replace the original column:
        new_dedup_cnt = (dedup_cnt - tpm) / dedup_cnt
        - dedup_cnt == 0 → Result is NaN (avoids division by zero)
        - absolute = True → Use |dedup_cnt - tpm| / dedup_cnt
        Other columns (method, error_rate, transcript_id, etc.) remain unchanged.
    """
    # print(1111)
    # print(ab_df)
    # print(df_all)
    m = df_all.merge(ab_df, on=on, how="left", validate="m:1")
    if keep_original and "dedup_cnt_orig" not in m.columns:
        m["dedup_cnt_orig"] = m[col]
    denom = pd.to_numeric(m[col], errors="coerce").replace(0, np.nan)
    diff = pd.to_numeric(m[col], errors="coerce") - m["tpm"]
    if absolute:
        diff = diff.abs()
    m[col] = diff / denom
    return m.drop(columns=["tpm"])


def make_distinct_colors(methods, cmap_name="tab20"):
    import matplotlib as mpl
    methods = list(methods)
    # base = plt.get_cmap(cmap_name).colors
    # print(base)
    import seaborn as sns
    cmap = sns.color_palette("cubehelix", as_cmap=True)  # colormap
    base = tuple([cmap(x)[:3] for x in np.linspace(0, 1, 20)])
    # print(base)
    order = [0,2,4,6,8,10,12,14,16,18,1,3,5,7,9,11,13,15,17,19]
    cols = [base[order[i % len(order)]] for i in range(len(methods))]
    to_hex = mpl.colors.to_hex
    return {m: to_hex(cols[i]) for i, m in enumerate(methods)}


def make_marker_map(
        methods: Iterable[str],
) -> Dict[str, str]:
    pool = ['o','s','^','D','v','P','X','*','h','H','>','<','d','p','8','x','+','1','2','3','4','|','_']
    methods = list(methods)
    return {m: pool[i % len(pool)] for i, m in enumerate(methods)}


def show_palette_strip(
        colors,
        figsize=(12, 1.1),
):
    """
    colors: dict(method -> hex)
    """
    from matplotlib.patches import Rectangle
    labels = list(colors.keys())
    fig, ax = plt.subplots(figsize=figsize, dpi=150)
    for i, lab in enumerate(labels):
        ax.add_patch(Rectangle((i, 0), 1, 1, color=colors[lab]))
    ax.set_xlim(0, len(labels)); ax.set_ylim(0, 1)
    ax.set_xticks(range(len(labels))); ax.set_xticklabels(labels, rotation=45, ha="right", fontsize=9)
    ax.set_yticks([]); ax.spines[:].set_visible(False)
    fig.tight_layout()
    return fig, ax


def plot_methods_vs_error_rates(
        df: pd.DataFrame,
        prefer: str = "min",
        agg: str = "median",
        interval: str = "iqr",
        methods_order: Iterable[str] | None = None,
        error_rates_order: Iterable[float] | None = None,
        xscale: str = "log",
        figsize=(11, 6),
        marker_size: float = 54,
        line_width: float = 2.2,
        alpha_band: float = 0.18,
        jitter_points: bool = False,

        annotate_lines: bool = False, # Turn offline endings (keep interfaces, but defaults to False)
        show_legend: bool = True, # use legend
        legend_outside: bool = False, # Whether to put it outside the picture
        legend_cols: int = 2,
        legend_fontsize: int = 10,
        # ★ Best way to mark stars; whether to add text to the stars
        show_best_star: bool = True,
        label_best_star: bool = False, # Set True to write the method name next to the star
        colors: dict | None = None,
        markers: dict | None = None,
        y_min: float | None = None,
        y_max: float | None = None,
        title: str | None = None,
        save_to: str | None = None,
        tick_exact: bool = True, # The x-axis is accurate using error_rates_order
):
    import seaborn as sns

    sns.set(font="Helvetica")
    sns.set_style("ticks")

    from matplotlib.ticker import ScalarFormatter
    need = {"method","error_rate","dedup_cnt"}
    if not need.issubset(df.columns):
        raise ValueError(f"df must contain {need}")

    df = df.copy()
    df["dedup_cnt"] = pd.to_numeric(df["dedup_cnt"], errors="coerce")
    df = df.dropna(subset=["dedup_cnt"])

    if methods_order is None:
        methods_order = sorted(df["method"].unique().tolist())
    if error_rates_order is None:
        error_rates_order = sorted(df["error_rate"].unique().tolist())

    if colors is None:
        colors = make_distinct_colors(methods_order)
    if markers is None:
        markers = make_marker_map(methods_order)

    g = df.groupby(["method","error_rate"])["dedup_cnt"]
    center = g.median() if agg=="median" else g.mean()
    if interval=="iqr":
        q25 = g.quantile(0.25); q75 = g.quantile(0.75)
    elif interval=="sem":
        m = g.mean(); s = g.std(ddof=1); n = g.size().clip(lower=1)
        sem = s/np.sqrt(n); q25, q75 = m-sem, m+sem
    else:
        raise ValueError("interval must be 'iqr' or 'sem'")

    fig, ax = plt.subplots(figsize=figsize, dpi=150)

    for method in methods_order:
        ers, y, ylo, yhi = [], [], [], []
        for er in error_rates_order:
            key = (method, er)
            if key in center.index:
                ers.append(er); y.append(center.loc[key])
                ylo.append(q25.loc[key]); yhi.append(q75.loc[key])
        if not ers:
            continue

        col = colors.get(method)
        mk  = markers.get(method, 'o')

        ax.fill_between(ers, ylo, yhi, alpha=alpha_band, facecolor=col, linewidth=0)
        ax.plot(ers, y, marker=mk, color=col, linewidth=line_width,
                markersize=marker_size**0.5, markerfacecolor='none', label=method)
        ax.spines['right'].set_color('none')
        ax.spines['top'].set_color('none')

        if jitter_points:
            rng = np.random.default_rng(0)
            for er in ers:
                ys = df.loc[(df["method"]==method)&(df["error_rate"]==er), "dedup_cnt"].values
                if ys.size:
                    xs = rng.uniform(er*0.98, er*1.02, size=ys.size) if xscale=="log" else \
                         rng.uniform(er-0.005, er+0.005, size=ys.size)
                    ax.scatter(xs, ys, s=10, alpha=0.25, edgecolors="none", color=col)

    # ★ for error rate best method
    if show_best_star:
        span = df["dedup_cnt"].max() - df["dedup_cnt"].min()
        dy = 0.012 * (span if span>0 else 1)
        for er in error_rates_order:
            items = [(m, center.loc[(m,er)]) for m in methods_order if (m,er) in center.index]
            if not items: continue
            best_m, best_v = (min if prefer=="min" else max)(items, key=lambda t:t[1])
            ax.scatter([er],[best_v], marker="*", s=160, zorder=5, color=colors.get(best_m))
            if label_best_star:
                ax.text(er, best_v+dy, best_m, ha="center", va="bottom", fontsize=8)

    # x ticks uses error rates
    ax.set_xlabel("Sequencing error rate")
    # ax.set_ylabel(f"dedup_cnt ({agg})")
    ax.set_ylabel(r'($\frac{N_e-N_t}{N_t}$)' + f" ({agg})")
    ax.set_xscale(xscale)
    if tick_exact:
        ax.set_xticks(error_rates_order)
        ax.get_xaxis().set_major_formatter(ScalarFormatter())
        ax.set_xticklabels([str(er) for er in error_rates_order])
    if xscale == "log":
        ax.set_xlim(min(error_rates_order)*0.95, max(error_rates_order)*1.15)
    if y_min is not None or y_max is not None:
        yl = ax.get_ylim()
        ax.set_ylim(y_min if y_min is not None else yl[0],
                    y_max if y_max is not None else yl[1])
    ax.grid(True, linestyle="--", alpha=0.35)

    # Legend
    if show_legend:
        if legend_outside:
            ax.legend(
                title="UMI deduplication method",
                frameon=False,
                ncol=legend_cols,
                bbox_to_anchor=(1.02, 1),
                loc="upper left",
                fontsize=legend_fontsize,
            )
            fig.tight_layout(rect=[0, 0, 0.82, 1])
        else:
            ax.legend(
                title="UMI deduplication method",
                frameon=False,
                ncol=legend_cols,
                fontsize=legend_fontsize,
            )
            fig.tight_layout()
    else:
        fig.tight_layout()
    if title: ax.set_title(title, fontsize=14)
    if save_to: fig.savefig(save_to, bbox_inches="tight")
    return fig, ax


def build_filemap_from_template(
        methods: Dict[str, str],
        error_rates: Iterable[float],
        template: str,
) -> Dict[str, Dict[float, str]]:
    """
    {method: {er: path}}
    template e.g. "/path/{method}/reads_with_tags_{er}_dedup_sum.txt"
    """
    fmap = {}
    for m_key, m_dir in methods.items():
        # print(m_dir)
        fmap[m_key] = {}
        for er in error_rates:
            # Keep the string format 0.001
            er_str = str(er)
            path = template.format(method=m_dir, er=er_str)
            fmap[m_key][er] = path
    # print(fmap)
    return fmap


if __name__ == "__main__":
    error_rates = [
        0,
        0.0005,
        0.001,
        0.01,
        0.05,
        0.1,
        # 0.2, 0.3,
    ]

    methods = {
        "UMI-tools: Adjacency": "adjacency",
        "UMI-tools: Cluster": "cluster",
        "UMI-tools: Directional": "directional",
        # "Set Cover1": "set_cover_spall",

        "Set Cover": "set_cover",
        "Majority Vote": "majority_vote",
        'UMICountR: adj': 'adj',
        'UMICountR: singleton': 'adj_singleton',
        'UMICountR: adj-direc': 'adj_direc',
        'mclUMI: MCL': 'mcl',
        'mclUMI: MCL-val': 'mcl_val',
        'mclUMI: MCL-ed': 'mcl_ed',
        'DBSCAN': 'dbscan',
        'Birch': 'birch',
        'Affinity Propagation': 'aprop',
        'STARsolo': 'starsolo',
        'Gencore': 'gencore',
        'IRescue': 'irescue',
        # 'UMIS': 'umis',
    }

    filemap = build_filemap_from_template(
        methods=methods,
        error_rates=error_rates,
        # template="D:/Document/Programming/Python/umiche/umiche/data/r1/simuread/tksm/ed1/{method}/reads_with_tags_{er}_dedup_sum.txt",
        # template="D:/Document/Programming/Python/umiche/umiche/data/r1/simuread/minnow/ed1/{method}/reads_with_tags_{er}_dedup_sum.txt",
        # template="D:/Document/Programming/Python/umiche/umiche/data/r1/simuread/asarusim/ed1/{method}/simulated_with_tags_{er}_dedup_sum.txt",
        template="D:/Document/Programming/Python/umiche/umiche/data/r1/simuread/screadsim/ed1/{method}/reads_with_tags_{er}_dedup_sum.txt",
        # template="D:/Document/Programming/Python/umiche/umiche/data/r1/simuread/tresor/ed1_ht_spl_pct10/{method}/tresor_{er}_dedup_sum.txt",

        # CMI
        # template="D:/Document/Programming/Python/umiche/umiche/data/r1/cmi/ed1/{method}/Aligned_GACTGCTACT_gene_sorted_primary_tagged_{er}_dedup_sum.txt",
    )

    df_all = load_methods_errorrates_with_relerr(
        filemap,
        # abundance_path="D:/Document/Programming/Python/umiche/umiche/data/r1/simuread/tksm/abundance.trimmed.tsv",
        # abundance_path="D:/Document/Programming/Python/umiche/umiche/data/r1/simuread/minnow/abundance.tsv",
        # abundance_path="D:/Document/Programming/Python/umiche/umiche/data/r1/simuread/asarusim/abundance.tsv",
        abundance_path="D:/Document/Programming/Python/umiche/umiche/data/r1/simuread/screadsim/abundance.tsv",
        # abundance_path="D:/Document/Programming/Python/umiche/umiche/data/r1/simuread/tresor/abundance.tsv",

        # CMI
        # abundance_path="D:/Document/Programming/Python/umiche/umiche/data/r1/cmi/abundance.tsv",
        sep="\t",
        absolute=False,
    )
    # print(df_all.columns)
    # print(df_all['transcript_id'])

    methods_order = sorted(df_all["method"].unique())
    error_rates_order = sorted(df_all["error_rate"].unique())

    # 2) Generate 17 color palette and view
    colors = make_distinct_colors(methods_order, cmap_name="tab20")
    # show_palette_strip(colors)
    # plt.show()

    markers = make_marker_map(methods_order)
    # print(df_all)
    fig, ax = plot_methods_vs_error_rates(
        df_all,
        prefer="min", # or 'max'
        agg="median", # or 'mean'
        interval="iqr", # or 'sem'
        methods_order=list(methods.keys()),
        error_rates_order=error_rates_order,
        xscale="log", # When the error rate span is large, log readings are clearer; you can also change to 'linear'
        figsize=(8, 6),
        marker_size=54,
        line_width=2.2,
        alpha_band=0.1,
        tick_exact=True,
        annotate_lines=False,
        show_legend=True, # Using legend
        legend_outside=False, # It is clearer if placed outside the image (can be changed to False)
        legend_cols=2, # 17 methods suggest 3-4 columns
        colors=colors, # <— Fixed color for each method
        markers=markers,
        jitter_points=False, # If you want to see the original points, press True (maybe very dense)
        show_best_star=False,
        y_min=None,
        y_max=None, # If you need to focus on the interval, you can set
        # title="TKSM (10X Long-read)",
        # title="Minnow (scRNA-seq)",
        # title="AsaruSim (Nanopore long-read)",
        title="scReadSim (ATAC-seq)",
        # title="Tresor",
        save_to=None, # e.g. "methods_vs_error_rates.png"
    )
    # plot_methods_vs_error_rates(
    #     df_all,
    #     prefer="min",  # or 'max'
    #     agg="median",
    #     interval="iqr",
    #     methods_order=methods_order,
    #     error_rates_order=error_rates_order,
    #     xscale="log",
    #     title="Dedup counts vs. error rate",
    # )
    plt.subplots_adjust(
        top=0.94,
        # bottom=0.06,
        # left=0.06,
        # right=0.98,
        # hspace=0.40,
        # wspace=0.15,
    )
    plt.show()