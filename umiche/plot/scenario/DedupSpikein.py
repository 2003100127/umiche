__version__ = "v1.0"
__copyright__ = "Copyright 2025"
__license__ = "GPL-3.0"
__developer__ = "Jianfeng Sun"
__maintainer__ = "Jianfeng Sun"
__email__="jianfeng.sunmt@gmail.com"


from typing import Dict

import math
import numpy as np
import pandas as pd
import seaborn as sns
import matplotlib.pyplot as plt
from umiche.simu.Parameter import Parameter as params
from umiche.deduplicate.io.Stat import Spikein as dedupstat

from umiche.util.Reader import Reader as freader
from umiche.util.Console import Console


class DedupSpikein:

    def __init__(
            self,
            scenarios: Dict,
            methods: Dict,
            token: str,
            umi_gt_cnt: int = 50,
            param_fpn : str = None,
            verbose : bool = False,
    ):
        self.scenarios = scenarios
        self.methods = methods
        self.umi_gt_cnt = umi_gt_cnt
        self.token = token
        self.param_fpn = param_fpn
        self.freader = freader()
        self.params = params(param_fpn=self.param_fpn)

        self.dedupstat = dedupstat(
            scenarios=self.scenarios,
            methods=self.methods,
            token=self.token,
            param_fpn=self.param_fpn,
        )

        self.df_dedup = self.dedupstat.df_dedup

        sns.set(font="Helvetica")
        sns.set_style("ticks")

        self.console = Console()
        self.console.verbose = verbose

    def plot_dedup_by_method_per_scenario(
            self,
            value_col: str = "dedup_cnt",
            method_col: str = "method",
            scenario_col: str = "scenario",
            method_order: list | None = None,
            scenario_order: list | None = None,
            ncols: int = 3,  # Number of subplots per row
            figsize=(21, 9),
            sharey: bool = True,  # Whether all subplots share the y-axis for easier comparison
            jitter_width: float = 0.18,  # Half-width for horizontal jitter of scatter points around the box center
            point_size: float = 14,
            point_alpha: float = 0.6,
            seed: int = 0,
            showfliers: bool = False,  # Whether to show outliers in the boxplot; usually disabled when plotting custom scatter points
            title: str | None = None,
            save_to: str | None = None,
    ):
        """
        Given a DataFrame containing [dedup_cnt, method, scenario],
        create one subplot per scenario; in each subplot compare dedup_cnt for each method
        using a boxplot + uniformly jittered scatter points.
        """
        # Basic validation
        for col in (value_col, method_col, scenario_col):
            if col not in self.df_dedup.columns:
                raise ValueError(f"Column '{col}' not in df. Available: {list(self.df_dedupcolumns)}")

        # Ensure value column is numeric
        df = self.df_dedup.copy()
        df[value_col] = pd.to_numeric(df[value_col], errors="coerce")

        # Sorting: if not specified, use alphabetical order of unique values
        if scenario_order is None:
            scenario_order = sorted(df[scenario_col].dropna().unique().tolist())
        if method_order is None:
            method_order = sorted(df[method_col].dropna().unique().tolist())

        # Subplot grid
        n_scen = len(scenario_order)
        ncols = max(1, int(ncols))
        nrows = math.ceil(n_scen / ncols)

        fig, axes = plt.subplots(nrows, ncols, figsize=figsize, sharey=sharey)
        if isinstance(axes, np.ndarray):
            axes = axes.ravel()
        else:
            axes = np.array([axes])

        # Set unified y-axis limits for comparability
        if sharey:
            y_min = np.nanmin(df[value_col].values)
            y_max = np.nanmax(df[value_col].values)
            # Add some padding to the y-axis range
            pad = 0.02 * (y_max - y_min if np.isfinite(y_max - y_min) else 1.0)
            y_lim = (y_min - pad, y_max + pad)
        else:
            y_lim = None

        rng = np.random.default_rng(seed)

        # Plot each scenario
        for i, scen in enumerate(scenario_order):
            ax = axes[i]
            sub = df.loc[df[scenario_col] == scen].copy()

            # Keep only methods that appear in this scenario, preserving method_order
            local_methods = [m for m in method_order if (sub[method_col] == m).any()]
            if len(local_methods) == 0:
                ax.set_axis_off()
                ax.set_title(f"{scen} (no data)")
                continue

            data = [sub.loc[sub[method_col] == m, value_col].dropna().to_numpy() for m in local_methods]
            positions = np.arange(1, len(local_methods) + 1)

            # Boxplot (no fill color, simple lines; no fliers because we plot scatter points ourselves)
            bp = ax.boxplot(
                data,
                positions=positions,
                widths=0.55,
                showfliers=showfliers,
                patch_artist=False,
                medianprops=dict(linewidth=1.5),
                whiskerprops=dict(linewidth=1.2),
                capprops=dict(linewidth=1.2),
                boxprops=dict(linewidth=1.2),
            )
            ax.spines['right'].set_color('none')
            ax.spines['top'].set_color('none')
            # Uniform jitter scatter points: spread points around each boxplot position
            for x0, ys in zip(positions, data):
                if len(ys) == 0:
                    continue
                xs = rng.uniform(x0 - jitter_width, x0 + jitter_width, size=len(ys))
                ax.scatter(xs, ys, s=point_size, alpha=point_alpha, edgecolors="none", linewidths=0)

            # Axis and labels
            ax.set_xticks(positions)
            ax.set_xticklabels(local_methods, rotation=30, ha="right")

            if (i % ncols) == 0:  # Add y-axis label to the first subplot of each row
                ax.set_ylabel(value_col)

            ax.set_title(str(scen))
            ax.grid(axis="y", linestyle="--", alpha=0.35)
            ax.set_xlim(0.5, len(local_methods) + 0.5)
            if y_lim is not None:
                ax.set_ylim(*y_lim)

        # Turn off extra empty axes if the subplot count is not a perfect multiple of nrows*ncols
        for j in range(n_scen, len(axes)):
            axes[j].set_axis_off()

        if title:
            fig.suptitle(title, y=0.995)

        fig.tight_layout()
        if title:
            # suptitle and tight_layout may conflict
            plt.subplots_adjust(top=0.92)

        if save_to:
            fig.savefig(save_to, bbox_inches="tight", dpi=150)

        return fig, axes

    def plot_dedup_by_method_per_scenario_control(
            self,
            # df: pd.DataFrame,
            value_col: str = "dedup_cnt",
            method_col: str = "method",
            scenario_col: str = "scenario",
            method_order: list | None = None,
            scenario_order: list | None = None,
            ncols: int = 3,
            figsize=(16, 9),
            sharey: bool = False,
            sharex: bool = True,
            jitter_width: float = 0.18,
            point_size: float = 14,
            point_alpha: float = 0.6,
            seed: int = 0,
            showfliers: bool = False,
            title: str | None = None,
            save_to: str | None = None,

            y_min: float | None = None,
            y_max: float | None = None,
            clip_scatter_above: float | None = None,
            clip_scatter_below: float | None = None,
            # The box plot statistics still use the full data.
            truncate_for_boxplot: bool = False,
    ):
        """
        One subplot per scenario; in each subplot, plot dedup_cnt
        by method (box plot + uniform scatter points).

        New:
          - y_min / y_max: y-axis
          - clip_scatter_above / clip_scatter_below
          - truncate_for_boxplot：True, y_min/y_max
        """
        for col in (value_col, method_col, scenario_col):
            if col not in self.df_dedup.columns:
                raise ValueError(f"Column '{col}' not in df. Available: {list(df.columns)}")

        df = self.df_dedup.copy()
        df[value_col] = pd.to_numeric(df[value_col], errors="coerce")
        df[value_col] = df[value_col].apply(lambda x: x - 1)
        # print(df[value_col])

        if scenario_order is None:
            scenario_order = sorted(df[scenario_col].dropna().unique().tolist())
        if method_order is None:
            method_order = df[method_col].dropna().unique().tolist()
            # method_order = sorted(df[method_col].dropna().unique().tolist())
        # palette_name = "Set3"  # 可改为 "Set2" / "tab20" / "Set3" husl 等
        # pal = sns.color_palette(palette_name, n_colors=len(method_order))
        pal = ["grey"] * len(method_order)
        method_color = dict(zip(method_order, pal))

        n_scen = len(scenario_order)
        ncols = max(1, int(ncols))
        nrows = math.ceil(n_scen / ncols)

        fig, axes = plt.subplots(nrows, ncols, figsize=figsize, sharex=sharex, sharey=sharey)
        if isinstance(axes, np.ndarray):
            axes = axes.ravel()
        else:
            axes = np.array([axes])

        # y-axis y_min/y_max
        if sharey:
            if y_min is None:
                y_min_auto = np.nanmin(df[value_col].values)
            else:
                y_min_auto = y_min
            if y_max is None:
                y_max_auto = np.nanmax(df[value_col].values)
            else:
                y_max_auto = y_max
            pad = 0.02 * (y_max_auto - y_min_auto if np.isfinite(y_max_auto - y_min_auto) else 1.0)
            y_lim_global = (
                (y_min if y_min is not None else y_min_auto) - (0 if (y_min is not None) else pad),
                (y_max if y_max is not None else y_max_auto) + (0 if (y_max is not None) else pad),
            )
        else:
            y_lim_global = None

        rng = np.random.default_rng(seed)

        for i, scen in enumerate(scenario_order):
            ax = axes[i]
            sub = df.loc[df[scenario_col] == scen].copy()
            box_data = []
            scatter_data = []
            local_methods = []
            for m in method_order:
                ys_full = sub.loc[sub[method_col] == m, value_col].dropna().to_numpy()
                if ys_full.size == 0:
                    continue

                if truncate_for_boxplot and (y_min is not None or y_max is not None):
                    mask_box = np.ones_like(ys_full, dtype=bool)
                    if y_min is not None:
                        mask_box &= ys_full >= y_min
                    if y_max is not None:
                        mask_box &= ys_full <= y_max
                    ys_box = ys_full[mask_box]
                else:
                    ys_box = ys_full

                if ys_box.size == 0:
                    continue

                ys_scat = ys_full.copy()
                if clip_scatter_below is not None:
                    ys_scat = ys_scat[ys_scat >= clip_scatter_below]
                elif y_min is not None:
                    # clip_scatter_below, y_min by default
                    ys_scat = ys_scat[ys_scat >= y_min]

                if clip_scatter_above is not None:
                    ys_scat = ys_scat[ys_scat <= clip_scatter_above]
                elif y_max is not None:
                    # clip_scatter_above, y_max by default
                    ys_scat = ys_scat[ys_scat <= y_max]
                local_methods.append(m)
                box_data.append(ys_box)
                scatter_data.append(ys_scat)
            if len(local_methods) == 0:
                ax.set_axis_off()
                ax.set_title(f"{scen} (no data)")
                continue
            positions = np.arange(1, len(local_methods) + 1)
            # boxplot
            ax.boxplot(
                box_data,
                positions=positions,
                widths=0.55,
                showfliers=showfliers,
                patch_artist=False,
                medianprops=dict(linewidth=1.5),
                whiskerprops=dict(linewidth=1.2),
                capprops=dict(linewidth=1.2),
                boxprops=dict(linewidth=1.2),
            )
            # scatter jitter evenly
            for x0, ys, m in zip(positions, scatter_data, local_methods):
                if len(ys) == 0:
                    continue
                xs = rng.uniform(x0 - jitter_width, x0 + jitter_width, size=len(ys))
                ax.scatter(
                    xs, ys,
                    s=point_size,
                    alpha=point_alpha,
                    edgecolors="none",
                    linewidths=0,
                    color=method_color[m],
                )
            # style
            ax.set_xticks(positions)
            ax.set_xticklabels(local_methods, rotation=20, ha="right", fontsize=14)
            if (i % ncols) == 0:
                # ax.set_ylabel(value_col, fontsize=14)
                ax.set_ylabel(r'$\frac{N_e-N_t}{N_t}$', fontsize=18)
            ax.set_title(str(scen), fontsize=16)
            ax.grid(axis="y", linestyle="--", alpha=0.35)
            ax.set_xlim(0.5, len(local_methods) + 0.5)
            ax.spines['right'].set_color('none')
            ax.spines['top'].set_color('none')
            # y-axis
            if sharey:
                ax.set_ylim(*y_lim_global)
            else:
                # y_min/y_max
                if y_min is not None or y_max is not None:
                    yl = ax.get_ylim()
                    ax.set_ylim(
                        y_min if y_min is not None else yl[0],
                        y_max if y_max is not None else yl[1],
                    )
        # off empty axis
        for j in range(n_scen, len(axes)):
            axes[j].set_axis_off()
        if title:
            fig.suptitle(title, y=0.995, fontsize=14)
        fig.tight_layout()
        if save_to:
            fig.savefig(save_to, bbox_inches="tight", dpi=150)
        return fig, axes


if __name__ == "__main__":
    from umiche.path import to

    p = DedupSpikein(
        scenarios={
            # '0.0005': 'Error 0.0005',
            # '0.001': 'Error 0.001',
            '0.01': 'Sequencing error rate 0.01',
            # '0.05': 'Error 0.05',
            # '0.1': 'Error 0.1',
            # '0.2': 'Error 0.2',
            # '0.3': 'Error 0.3',
            # '0.4': 'Error 0.4',
            # '0.5': 'Error 0.5',
            '0': 'Oringinal',
        },
        methods={
            'unique': 'Unique',
            'cluster': 'UMI-tools:Cluster',
            'adjacency': 'UMI-tools:Adjacency',
            'directional': 'UMI-tools:Directional',
            'adj': 'UMICountR:adj',
            'adj_singleton': 'UMICountR:singleton',
            'adj_direc': 'UMICountR:adj-direc',
            'mcl': 'mclUMI:MCL',
            'mcl_val': 'mclUMI:MCL-val',
            'mcl_ed': 'mclUMI:MCL-ed',
            'dbscan': 'DBSCAN',
            'birch': 'Birch',
            'aprop': 'Affinity Propagation',
            'set_cover': 'Set Ccover',
            'majority_vote': 'Majority Vote',
            'dropest': 'dropEst',
            'starsolo': 'STARsolo',
            'irescue': 'IRescue',
            'gencore': 'Gencore',
            # 'umis': 'UMIS',
        },
        # token='Aligned_GACTGCTACT_gene_sorted_primary_tagged',
        token='Smartseq3.TTACCTGCCAGATTCG',
        param_fpn=to('data/params_spikein.yml'),
    )

    # fig, axes = p.plot_dedup_by_method_per_scenario(
    #     value_col="dedup_cnt",
    #     method_col="method",
    #     scenario_col="scenario",
    #     # Optional: fix the display order of methods and scenarios (if not provided, alphabetical order is used)
    #     # method_order=["MethodA", "MethodB", ..., "MethodL"],
    #     # scenario_order=["Error 0.05", "Error 0.1", "Error 0.2", "Error 0.5", "Error 1.0"],
    #     ncols=1,
    #     figsize=(21, 10),
    #     jitter_width=0.20,
    #     point_size=16,
    #     point_alpha=0.6,
    #     sharey=True,
    #     showfliers=False,
    #     title="",
    #     save_to=None, # e.g. "dedup_cnt_box_scatter_by_scenario.png"
    # )
    # plt.show()

    p.plot_dedup_by_method_per_scenario_control(
        # y_max=140, # max to 200
        # y_max=20, # max to 200
        ncols=1,
        truncate_for_boxplot=False,
        figsize=(16, 7),
        title=None,
    )
    plt.subplots_adjust(
        top=0.94,
        # bottom=0.06,
        left=0.06,
        # right=0.98,
        # hspace=0.40,
        # wspace=0.15,
    )
    plt.show()
