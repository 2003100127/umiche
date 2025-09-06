__version__ = "v1.0"
__copyright__ = "Copyright 2025"
__license__ = "GPL-3.0"
__developer__ = "Jianfeng Sun"
__maintainer__ = "Jianfeng Sun"
__email__="jianfeng.sunmt@gmail.com"

from typing import Dict, List, Sequence, Tuple, Optional

import os
import numpy as np
import pandas as pd
import seaborn as sns
import matplotlib.pyplot as plt
from math import ceil


def _try_read_summary(
        base_dir: str,
        method_key: str,
        seqtech: str,
        err_label: str,
        sep: str = "\t",
        ed_dir: str = "ed1",
) -> Optional[pd.DataFrame]:
    """Try multiple filename patterns; return DataFrame if found, else None."""
    candidates = [
        os.path.join(base_dir, seqtech, ed_dir, method_key, f"{seqtech}_{err_label}_dedup_sum.txt"),
        os.path.join(base_dir, method_key, f"{seqtech}_{err_label}_dedup_sum.txt"),
        os.path.join(base_dir, f"{seqtech}_{err_label}_dedup_sum.txt"),
        os.path.join(base_dir, method_key, f"{seqtech}-{err_label}-dedup-sum.txt"),
        os.path.join(base_dir, f"{seqtech}-{err_label}-dedup-sum.txt"),
    ]
    for fp in candidates:
        if os.path.isfile(fp):
            try:
                df = pd.read_csv(fp, sep=sep)
                if "dedup_cnt" not in df.columns:
                    print(f"[WARN] File found but no 'dedup_cnt' column: {fp}")
                    return None
                return df
            except Exception as e:
                print(f"[WARN] Failed to read {fp}: {e}")
                return None
    print(f"[MISS] No summary found for method={method_key}, seqtech={seqtech}, err={err_label}")
    return None


def load_aggregated_matrix(
        base_dir: str,
        seqtech: str,
        methods_map: Dict[str, str],
        error_rates: Sequence[str | float | int],
        agg: str = "median",
        ed_dir: str = "ed1",
        use_tpm_subtraction: bool = True,
        abun_filename: str = "abundance.tsv",
) -> Tuple[np.ndarray, np.ndarray]:
    method_disp = list(methods_map.keys())
    method_key = list(methods_map.values())
    err_labels = [str(e) for e in error_rates]
    vals = np.full((len(method_disp), len(err_labels)), np.nan, dtype=float)
    # Load abundance data (one-time read per sequencing technology)
    tpm_map = None
    if use_tpm_subtraction:
        tpm_series = _read_abundance(base_dir, seqtech, abun_filename=abun_filename, ed_dir=ed_dir)
        # Convert target_id to string before joining to handle mixed types.
        tpm_series.index = tpm_series.index.astype(str)
        tpm_map = tpm_series.to_dict()
    def _agg(series: pd.Series) -> float:
        if agg == "median":
            return float(series.median())
        elif agg == "mean":
            return float(series.mean())
        elif agg == "sum":
            return float(series.sum())
        else:
            raise ValueError(f"Unknown agg={agg}")
    for i, mkey in enumerate(method_key):
        for j, el in enumerate(err_labels):
            df = _try_read_summary(base_dir, mkey, seqtech, el, ed_dir=ed_dir)
            if df is None or "dedup_cnt" not in df.columns or len(df) == 0:
                continue
            # 1) If the result file contains a target_id/gene column, use it for alignment;
            # 2) Otherwise, use row numbers 0..N-1 as target_id to align with the 0..9 in abundance.tsv.
            tid_col = None
            for cand in ("target_id", "gene", "gid"):
                if cand in df.columns:
                    tid_col = cand
                    break
            if tid_col is not None:
                tids = df[tid_col].astype(str)
            else:
                tids = pd.Series(range(len(df)), name="target_id").astype(str)
            dedup = df["dedup_cnt"].astype(float).reset_index(drop=True)
            if use_tpm_subtraction and tpm_map:
                tpm_aligned = tids.map(tpm_map).astype(float).fillna(0.0).reset_index(drop=True)
                adj = (dedup - tpm_aligned) / tpm_aligned
            else:
                adj = dedup
            vals[i, j] = _agg(adj)
    mask = np.isnan(vals)
    return vals, mask


def _read_abundance(
        base_dir: str, seqtech: str,
        abun_filename: str = "abundance.tsv",
        ed_dir: str = "ed1",
):
    """
    Read abundance.tsv for a seqtech. Expected columns: target_id, tpm, cell
    Returns a Series mapping target_id -> tpm (float).
    """
    candidates = [
        # .../simu/{seqtech}/abundance.tsv
        os.path.join(base_dir, seqtech, abun_filename),
        # .../simu/{seqtech}/ed1/abundance.tsv (fallback)
        os.path.join(base_dir, seqtech, ed_dir, abun_filename),
    ]
    for fp in candidates:
        if os.path.isfile(fp):
            df = pd.read_csv(fp, sep="\s+", header=0)
            print(df)
            df.columns = [c.strip().lower() for c in df.columns]
            if not {"target_id", "tpm"}.issubset(set(df.columns)):
                raise ValueError(f"abundance file missing target_id/tpm columns: {fp}")
            s = df.set_index("target_id")["tpm"].astype(float).sort_index()
            return s
    raise FileNotFoundError(f"abundance.tsv not found for {seqtech} under {base_dir}")


class DedupSimu:

    def __init__(self, ):
        sns.set(font="Helvetica")
        sns.set_style("ticks")

    def plot_small_multiples_heatmaps(
            self,
            results: Dict[str, Tuple[np.ndarray, np.ndarray]],
            methods_display: List[str],
            error_rates: Sequence[str | float | int],
            suptitle: str = "UMI Dedup: median dedup_cnt per method × error rate",
            figsize_per_subplot: Tuple[float, float] = (4.0, 3.2),
            ncols: int = 5,
            cmap: str = "viridis",
            annotate: bool = False,
            sv_fpn: Optional[str] = None,
    ):
        """
        Draw one heatmap per seqtech, consistent color scale, shared colorbar.
        """
        seq_list = list(results.keys())
        n = len(seq_list)
        ncols = max(1, ncols)
        nrows = ceil(n / ncols)

        # Global color scale across all subplots (ignore NaNs)
        all_vals = []
        for seq in seq_list:
            m, _ = results[seq]
            if np.isfinite(m).any():
                all_vals.append(m[np.isfinite(m)])
        if all_vals:
            vmin = np.min([a.min() for a in all_vals])
            vmax = np.max([a.max() for a in all_vals])
        else:
            vmin, vmax = 0.0, 1.0

        # Figure & axes
        fig, axes = plt.subplots(
            nrows=nrows, ncols=ncols,
            figsize=(figsize_per_subplot[0] * ncols, figsize_per_subplot[1] * nrows),
            squeeze=False
        )

        err_labels = [str(e) for e in error_rates]
        for idx, seq in enumerate(seq_list):
            r, c = divmod(idx, ncols)
            ax = axes[r][c]
            matrix, mask = results[seq]
            im = ax.imshow(matrix, aspect="auto", cmap=cmap, vmin=vmin, vmax=vmax, interpolation="nearest")

            # Gridlines to improve readability
            ax.set_xticks(np.arange(len(err_labels)) + 0.5, minor=True)
            ax.set_yticks(np.arange(len(methods_display)) + 0.5, minor=True)
            ax.grid(which="minor", color="white", linewidth=0.6)
            ax.tick_params(which="minor", length=0)

            # Ticks & labels
            ax.set_xticks(np.arange(len(err_labels)))
            ax.set_yticks(np.arange(len(methods_display)))
            ax.set_xticklabels(err_labels, rotation=90, fontsize=8)
            ax.set_yticklabels(methods_display, fontsize=7)
            ax.set_xlabel("Error rate", fontsize=9)
            ax.set_title(seq, fontsize=10, pad=6)

            # Optional annotations
            if annotate and matrix.size <= 80:  # avoid clutter
                for i in range(matrix.shape[0]):
                    for j in range(matrix.shape[1]):
                        if np.isfinite(matrix[i, j]):
                            ax.text(j, i, f"{matrix[i, j]:.0f}", ha="center", va="center", fontsize=7)
                        else:
                            ax.text(j, i, "×", ha="center", va="center", fontsize=7, color="red")

        # Remove any empty axes
        for k in range(n, nrows * ncols):
            r, c = divmod(k, ncols)
            axes[r][c].axis("off")

        # Shared colorbar
        cax = fig.add_axes([0.92, 0.13, 0.015, 0.74])
        norm_mappable = plt.cm.ScalarMappable(cmap=cmap)
        norm_mappable.set_clim(vmin=vmin, vmax=vmax)
        cb = fig.colorbar(norm_mappable, cax=cax)
        cb.set_label("Median dedup_cnt", rotation=90)

        fig.suptitle(suptitle, fontsize=12)
        plt.subplots_adjust(right=0.90, wspace=0.30, hspace=0.35)
        if sv_fpn:
            plt.savefig(sv_fpn, dpi=300, bbox_inches="tight")
            print(f"[OK] Figure saved to: {sv_fpn}")

    def plot_small_multiples_bubbles(
            self,
            results: Dict[str, Tuple[np.ndarray, np.ndarray]],
            methods_display: List[str],
            error_rates: Sequence[str | float | int],
            suptitle: str = "UMI Dedup: median dedup_cnt (bubble size)",
            figsize_per_subplot: Tuple[float, float] = (4.0, 3.0),
            ncols: int = 5,
            s_range: Tuple[float, float] = (16, 520),
            use_log: bool = False,
            monochrome: bool = True,
            cmap: str = "viridis",
            sv_fpn: Optional[str] = None,
            annotate_topk: int = 0,
    ):
        """

        Parameters
        ----------
        results
        methods_display
        error_rates
        suptitle
        figsize_per_subplot
        ncols
        s_range
            Dot area range; can be scaled down to (14, 360) for a more compact layout.
        use_log
        monochrome
            Monochromatic strokes for cleaner look; set to False to enable colors.
        cmap
        sv_fpn
        annotate_topk

        Returns
        -------

        """
        import matplotlib.cm as cm

        seq_list = list(results.keys())
        n = len(seq_list)
        ncols = max(1, ncols)
        nrows = ceil(n / ncols)

        all_vals = []
        for seq in seq_list:
            m, _ = results[seq]
            if np.isfinite(m).any():
                all_vals.append(m[np.isfinite(m)])
        if all_vals:
            raw_min = float(np.min([a.min() for a in all_vals]))
            raw_max = float(np.max([a.max() for a in all_vals]))
        else:
            raw_min, raw_max = 0.0, 1.0

        def transform(x: np.ndarray) -> np.ndarray:
            if use_log:
                return np.log10(1.0 + np.clip(x, a_min=0, a_max=None))
            return x

        tmins, tmaxs = transform(np.array([raw_min])), transform(np.array([raw_max]))
        tmin, tmax = float(tmins), float(tmaxs)
        if tmax == tmin:
            tmax = tmin + 1e-6

        def map_size(v: float) -> float:
            tv = transform(np.array([v]))[0]
            frac = (tv - tmin) / (tmax - tmin)
            return s_range[0] + frac * (s_range[1] - s_range[0])

        norm = plt.Normalize(vmin=tmin, vmax=tmax)
        cmap_obj = cm.get_cmap(cmap)

        fig, axes = plt.subplots(
            nrows=nrows, ncols=ncols,
            figsize=(figsize_per_subplot[0] * ncols, figsize_per_subplot[1] * nrows),
            squeeze=False
        )
        err_labels = [str(e) for e in error_rates]
        n_err = len(err_labels)
        n_met = len(methods_display)

        legend_vals = np.linspace(raw_min, raw_max, num=3) if raw_max > raw_min else [raw_min]
        legend_handles, legend_labels = [], []
        for lv in legend_vals:
            h = plt.scatter([], [], s=map_size(lv),
                            facecolors="none" if monochrome else cmap_obj(norm(transform(np.array([lv]))[0])),
                            edgecolors="#444", linewidths=0.9)
            legend_handles.append(h)
            legend_labels.append(f"{lv:.0f}")

        for idx, seq in enumerate(seq_list):
            r, c = divmod(idx, ncols)
            ax = axes[r][c]
            matrix, mask = results[seq]

            for i in range(n_met):
                if i % 2 == 1:
                    ax.axhspan(i - 0.5, i + 0.5, color="#f4f4f4", zorder=0)

            for i in range(n_met):
                for j in range(n_err):
                    val = matrix[i, j]
                    if np.isfinite(val):
                        s = map_size(val)
                        if monochrome:
                            ax.scatter(j, i, s=s, facecolors="none", edgecolors="#2c2c2c",
                                       linewidths=0.9, zorder=3)
                        else:
                            ax.scatter(j, i, s=s, color=cmap_obj(norm(transform(np.array([val]))[0])),
                                       edgecolors="#2c2c2c", linewidths=0.4, zorder=3)
                    else:
                        ax.scatter(j, i, s=20, marker="x", color="red", linewidths=0.8, zorder=3)

            ax.set_xlim(-0.5, n_err - 0.5)
            ax.set_ylim(-0.5, n_met - 0.5)
            ax.set_xticks(range(n_err))
            ax.set_yticks(range(n_met))
            ax.set_xticklabels(err_labels, rotation=90, fontsize=8)
            ax.set_yticklabels(methods_display, fontsize=7)
            ax.grid(axis="x", linestyle=":", alpha=0.35)
            ax.set_title(seq, fontsize=10, pad=6)

            if annotate_topk > 0 and np.isfinite(matrix).any():
                flat = matrix.flatten()
                idxs = np.argpartition(-np.nan_to_num(flat, nan=-np.inf), kth=min(annotate_topk, flat.size - 1))[
                       :annotate_topk]
                for k in idxs:
                    i, j = divmod(int(k), n_err)
                    val = matrix[i, j]
                    if np.isfinite(val):
                        ax.text(j, i, f"{val:.0f}", ha="center", va="center",
                                fontsize=7, color="#1a1a1a", zorder=4)

        for k in range(n, nrows * ncols):
            r, c = divmod(k, ncols)
            axes[r][c].axis("off")

        leg = fig.legend(
            handles=legend_handles, labels=legend_labels, title="dedup_cnt",
            loc="center right", bbox_to_anchor=(0.98, 0.5), frameon=False
        )

        if not monochrome:
            cax = fig.add_axes([0.92, 0.13, 0.015, 0.74])
            cb = fig.colorbar(plt.cm.ScalarMappable(norm=norm, cmap=cmap_obj), cax=cax)
            cb.set_label("transformed value" if use_log else "value", rotation=90)

        fig.suptitle(suptitle, fontsize=12)
        plt.subplots_adjust(right=0.90, wspace=0.35, hspace=0.40)
        if sv_fpn:
            plt.savefig(sv_fpn, dpi=300, bbox_inches="tight")
            print(f"[OK] Figure saved to: {sv_fpn}")

    def plot_small_multiples_cellbars(
            self,
            results: Dict[str, Tuple[np.ndarray, np.ndarray]],
            methods_display: List[str],
            error_rates: Sequence[str | float | int],
            suptitle: str = "UMI dedup — median dedup_cnt (cell-bar length)",
            figsize_per_subplot: Tuple[float, float] = (3.6, 2.6),
            ncols: int = 5,
            use_log: bool = False,
            bar_face: str = "#3c3c3c",
            bar_alpha: float = 0.85,
            cell_bg: str = "#f5f6f8",
            sv_fpn: Optional[str] = None,
            seqtech_display_map=None,
            ylabel="Method",
            legend_label="dedup_cnt − TPM",
    ):
        """

        Parameters
        ----------
        results
        methods_display
        error_rates
        suptitle
        figsize_per_subplot
        ncols
        use_log
            If the panels vary significantly, change to True.
        bar_face
        bar_alpha
        cell_bg
        sv_fpn
        seqtech_display_map
        ylabel
        legend_label

        Returns
        -------

        """
        seq_list = list(results.keys())
        n = len(seq_list)
        ncols = max(1, ncols)
        nrows = ceil(n / ncols)

        # Global scale, ensuring a consistent length mapping across all panels
        all_vals = []
        for seq in seq_list:
            m, _ = results[seq]
            if np.isfinite(m).any():
                all_vals.append(m[np.isfinite(m)])
        if all_vals:
            vmin = float(np.min([a.min() for a in all_vals]))
            vmax = float(np.max([a.max() for a in all_vals]))
        else:
            vmin, vmax = 0.0, 1.0

        def transform(x: np.ndarray) -> np.ndarray:
            return np.log10(1.0 + np.clip(x, a_min=0, a_max=None)) if use_log else x

        tmin = float(transform(np.array([vmin])))
        tmax = float(transform(np.array([vmax])))
        if tmax == tmin:
            tmax = tmin + 1e-6

        def frac_len(v: float) -> float:
            tv = float(transform(np.array([v])))
            return max(0.0, min(1.0, (tv - tmin) / (tmax - tmin)))  # 0~1

        # canvas and axis
        fig, axes = plt.subplots(
            nrows=nrows,
            ncols=ncols,
            figsize=(figsize_per_subplot[0] * ncols, figsize_per_subplot[1] * nrows),
            sharey=True,
            # squeeze=False,
        )
        err_labels = [str(e) for e in error_rates]
        n_err = len(err_labels)
        n_met = len(methods_display)

        for idx, seq in enumerate(seq_list):
            r, c = divmod(idx, ncols)
            ax = axes[r][c]
            matrix, _ = results[seq]

            # Zebra striping to improve readability (subtle)
            for i in range(n_met):
                if i % 2 == 1:
                    ax.axhspan(i - 0.5, i + 0.5, color=cell_bg, zorder=0)

            # each cell as for a panel
            from matplotlib.patches import Rectangle
            cell_h = 0.66
            cell_w = 0.82
            for i in range(n_met):
                for j in range(n_err):
                    val = matrix[i, j]
                    cx_left = j - 0.41
                    cy_bottom = i - (cell_h / 2)

                    if np.isfinite(val):
                        w = cell_w * frac_len(val)
                        # bg frame
                        ax.add_patch(
                            Rectangle(
                                (cx_left, i - 0.5), 0.82, 1.0,
                                facecolor="none", edgecolor="#e9eaed",
                                linewidth=0.6, zorder=1
                            )
                        )
                        # actual
                        if w > 0:
                            ax.add_patch(
                                Rectangle(
                                    (cx_left, cy_bottom), w, cell_h,
                                    facecolor=bar_face, edgecolor="none",
                                    alpha=bar_alpha, zorder=3
                                )
                            )
                    else:
                        # red cross marker if missing
                        ax.scatter(j, i, s=18, marker="x", color="red", linewidths=0.9, zorder=4)
            # axis style
            ax.spines['right'].set_color('none')
            ax.spines['top'].set_color('none')
            ax.set_xlim(-0.5, n_err - 0.5)
            ax.set_ylim(-0.5, n_met - 0.5)
            ax.set_xticks(range(n_err))
            ax.set_yticks(range(n_met))
            ax.set_xticklabels(err_labels, rotation=30, fontsize=8)
            for lbl in ax.get_xticklabels():
                lbl.set_horizontalalignment('right')
            ax.set_yticklabels(methods_display, fontsize=8)
            if c == 0:
                ax.set_ylabel(ylabel, fontsize=11)
            ax.grid(axis="x", linestyle=":", alpha=0.3)
            title_txt = seqtech_display_map.get(seq, seq) if isinstance(seqtech_display_map, dict) else seq
            ax.set_title(title_txt, fontsize=12, pad=6)

        # clean empty axis
        for k in range(n, nrows * ncols):
            r, c = divmod(k, ncols)
            axes[r][c].axis("off")

        # min/median/max
        if vmax > vmin:
            scale_ax = fig.add_axes([0.905, 0.08, 0.085, 0.30])  # [left, bottom, width, height]
            scale_ax.axis("off")

            # title
            scale_ax.set_title(legend_label, fontsize=14, pad=8)

            nbars = 5
            vals = np.linspace(vmin, vmax, nbars)
            bar_h = 0.03  # The "thickness" (height) of the bar.
            top, bottom = 0.92, 0.10
            step = (top - bottom) / (nbars - 1)

            y = top
            for v in vals:
                L = 0.86 * frac_len(v)  # 0~0.86
                # Horizontal thin bars (solid color, cleaner look)
                scale_ax.add_patch(
                    Rectangle(
                        (0.05, y - bar_h / 2), L, bar_h,
                        facecolor=bar_face,
                        edgecolor="none",
                        alpha=bar_alpha
                    )
                )
                # Numerical labels slightly larger.
                scale_ax.text(0.05 + L + 0.025, y, f"{v:.0f}", va="center", fontsize=10)
                y -= step

        fig.suptitle(suptitle, fontsize=12)
        # plt.subplots_adjust(right=0.88, wspace=0.35, hspace=0.40)

        try:
            # plt ≥ 3.4
            fig.supxlabel("Sequencing error rate", fontsize=14, y=0.05)
        except AttributeError:
            # fallback for older Matplotlib
            fig.text(
                x=0.5,
                y=0.001,
                s="Sequencing error rate",
                ha="center",
                va="center",
                fontsize=14,
            )
        if sv_fpn:
            plt.savefig(sv_fpn, dpi=300, bbox_inches="tight")
            print(f"[OK] Figure saved to: {sv_fpn}")


if __name__ == "__main__":
    p = DedupSimu()
    seqtechs = {
        'vasaseq': 'VASA-seq',
        'pipseq': 'PIP-seq',
        'scifiseq': 'scifi-seq',
        'scrbseq': 'SCRB-seq/mcSCRB-seq',
        'shareseq': 'SHARE-seq',
        'snareseq': 'SNARE-seq',
        'splitseq': 'SPLiT-seq',
        'strtseqc1': 'STRT-seq-C1',
        'strtseq2i': 'STRT-seq-2i',
        'quartzseq2': 'Quartz-seq2',
        'petipseq': 'PETIP-seq',
        'marsseq2': 'MARS-seq2',
        'pairedseq': 'Paired-seq',
        'issaacseq': 'ISSAAC-seq',
        'indrop': 'inDrop',
        'dropseq': 'Drop-seq',
        '10xv2': '10x Chromium V2',
        '10xv3': '10x Chromium V3',
        'celseq2': 'CEL-seq2',
        'flashseq': 'FLASH-seq-UMI',
    }

    seq_keys = list(seqtechs.keys())

    methods = {
        "UMI-tools: Adjacency": "adjacency",
        "UMI-tools: Cluster": "cluster",
        "UMI-tools: Directional": "directional",
        "Set Cover": "set_cover",
        "Majority Vote": "majority_vote",
        "UMICountR: adj": "adj",
        "UMICountR: singleton": "adj_singleton",
        "UMICountR: adj-direc": "adj_direc",
        "mclUMI: MCL": "mcl",
        "mclUMI: MCL-val": "mcl_val",
        "mclUMI: MCL-ed": "mcl_ed",
        "DBSCAN": "dbscan",
        "Birch": "birch",
        "Affinity Propagation": "aprop",
        "STARsolo": "starsolo",
        "Gencore": "gencore",
        'IRescue': 'irescue',
        'UMIS': 'umis',
    }

    # Set your base directory containing results.
    # Example layout supported:
    #   base_dir/method_key/seqtech_err_dedup_sum.txt
    #   base_dir/seqtech_err_dedup_sum.txt
    base_dir = "D:/Document/Programming/Python/umiche/umiche/data/r1/seqtech/simu/"
    # base_dir = "D:/Document/Programming/Python/umiche/umiche/data/r1/seqtech/simu/pipseq/ed1/directional/pipseq_0_dedup_sum.txt"
    # ed_dir = "ed1"
    # ed_dir = "ed1_ht"
    ed_dir = "ed1_ht_spl_pct10"

    # Provide the list of error rates that appear in filenames (as strings or numbers)
    # e.g., ["0", "0.01", "0.05", "0.1"]
    # error_rates = ["0"]  # <-- change to your actual list, e.g. ["0","0.01","0.05","0.1"]
    error_rates = ["0.05", "0.1"]

    results = {}
    for st in seqtechs:
        mat, mask = load_aggregated_matrix(
            base_dir=base_dir,
            seqtech=st,
            methods_map=methods,
            error_rates=error_rates,
            agg="median",  # or "mean" / "sum"
            use_tpm_subtraction=True,  # (dedup_cnt - tpm)
            ed_dir=ed_dir,
            abun_filename="abundance.tsv",
        )
        results[st] = (mat, mask)

    # p.plot_small_multiples_heatmaps(
    #     results=results,
    #     methods_display=list(methods.keys()),
    #     error_rates=error_rates,
    #     suptitle="UMI deduplication — median dedup_cnt by method × error rate",
    #     figsize_per_subplot=(4.2, 3.0),
    #     ncols=5, # 5 columns × ceil(len(seqtechs)/5) rows
    #     cmap="magma", # or "magma", "plasma", "cividis"
    #     annotate=False,
    #     sv_fpn=None,
    # )

    # p.plot_small_multiples_bubbles(
    #     results=results,
    #     methods_display=list(methods.keys()),
    #     error_rates=error_rates,
    #     suptitle="UMI dedup — median dedup_cnt (bubble size)",
    #     figsize_per_subplot=(3.8, 2.8),
    #     ncols=5,
    #     s_range=(20, 520),
    #     use_log=False,
    #     monochrome=True,
    #     cmap="viridis",
    #     sv_fpn=None,
    #     annotate_topk=0,
    # )

    p.plot_small_multiples_cellbars(
        results=results,
        methods_display=list(methods.keys()),
        error_rates=error_rates,
        suptitle="",
        figsize_per_subplot=(3.4, 2.9),
        ncols=5,
        use_log=False,
        bar_face="#30353b",  # "#5b5f66"
        bar_alpha=0.9,
        sv_fpn=None, # "umi_dedup_cellbars.png"
        seqtech_display_map=seqtechs,
        ylabel="UMI deduplication method",
        legend_label=r'($\frac{N_e-N_t}{N_t}$)' + "(median)",
    )

    plt.subplots_adjust(
        top=0.94,
        # bottom=0.06,
        left=0.10,
        # right=0.98,
        hspace=0.30,
        wspace=0.15,
    )
    plt.show()
