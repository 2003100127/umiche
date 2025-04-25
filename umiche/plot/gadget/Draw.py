__author__ = "Jianfeng Sun"
__version__ = "v1.0"
__copyright__ = "Copyright 2025"
__license__ = "GPL v3.0"
__email__ = "jianfeng.sunmt@gmail.com"
__maintainer__ = "Jianfeng Sun"

from functools import wraps
from matplotlib.ticker import FuncFormatter


class Draw:

    def __init__(self, *args, **kwargs):
        self.args = args
        self.kwargs = kwargs
        self.type = self.kwargs['type']

    def __call__(self, deal):
        if self.type == 'bar_individual':
            draw = self.bar_individual
        # if self.type == 'bar_grouped':
        #     draw = self.bar_grouped
        if self.type == 'hist_grouped':
            draw = self.hist_grouped
        if self.type == 'line':
            draw = self.line
        if self.type == 'line_scatter':
            draw = self.line_scatter

        @wraps(deal)
        def config(self, *args, **kwargs):
            params_all = deal(self, *args, **kwargs)
            # print(kwargs)
            # print(params_all)
            kwarg_keys = kwargs.keys()
            if "color" not in kwarg_keys:
                color = None
            else:
                if kwargs['color'] == "by default":
                    color = 'grey'
                else:
                    color = kwargs['color']

            if "edgecolor" not in kwarg_keys:
                edgecolor = None
            else:
                if kwargs['edgecolor'] == "by default":
                    edgecolor = 'grey'
                else:
                    edgecolor = kwargs['edgecolor']

            if "facecolor" not in kwarg_keys:
                facecolor = None
            else:
                if kwargs['facecolor'] == "by default":
                    facecolor = 'grey'
                else:
                    facecolor = kwargs['facecolor']

            if "width" not in kwarg_keys:
                width = 0.8
            else:
                if kwargs['width'] == "by default":
                    width = 0.4
                else:
                    width = kwargs['width']

            if "linewidth" not in kwarg_keys:
                linewidth = None
            else:
                if kwargs['linewidth'] == "by default":
                    linewidth = 1
                else:
                    linewidth = kwargs['linewidth']

            if "alpha" not in kwarg_keys:
                alpha = None
            else:
                if kwargs['alpha'] == "by default":
                    alpha = 1
                else:
                    alpha = kwargs['alpha']

            if "label" not in kwarg_keys:
                label = None
            else:
                if kwargs['label'] == "by default":
                    label = ""
                else:
                    label = kwargs['label']

            if "x_label" not in kwarg_keys:
                x_label = None
            else:
                if kwargs['x_label'] == "by default":
                    x_label = ""
                else:
                    x_label = kwargs['x_label']

            if "y_label" not in kwarg_keys:
                y_label = None
            else:
                if kwargs['y_label'] == "by default":
                    y_label = ""
                else:
                    y_label = kwargs['y_label']

            if "title" not in kwarg_keys:
                title = None
            else:
                if kwargs['title'] == "by default":
                    title = ""
                else:
                    title = kwargs['title']

            if "x_ticks" not in kwarg_keys:
                x_ticks = None
            else:
                if kwargs['x_ticks'] == "by default":
                    x_ticks = None
                else:
                    x_ticks = kwargs['x_ticks']

            if "y_ticks" not in kwarg_keys:
                y_ticks = None
            else:
                if kwargs['y_ticks'] == "by default":
                    y_ticks = None
                else:
                    y_ticks = kwargs['y_ticks']

            if "x_ticklabels" not in kwarg_keys:
                x_ticklabels = None
            else:
                if kwargs['x_ticklabels'] == "by default":
                    x_ticklabels = None
                else:
                    x_ticklabels = kwargs['x_ticklabels']

            if "x_label_rotation" not in kwarg_keys:
                x_label_rotation = 0
            else:
                if kwargs['x_label_rotation'] == "by default":
                    x_label_rotation = 0
                else:
                    x_label_rotation = kwargs['x_label_rotation']

            if "x_label_rotation_align" not in kwarg_keys:
                x_label_rotation_align = 'center'
            else:
                if kwargs['x_label_rotation_align'] == "by default":
                    x_label_rotation_align = 'center'
                else:
                    x_label_rotation_align = kwargs['x_label_rotation_align']

            if "y_ticklabels" not in kwarg_keys:
                y_ticklabels = None
            else:
                if kwargs['y_ticklabels'] == "by default":
                    y_ticklabels = None
                else:
                    y_ticklabels = kwargs['y_ticklabels']

            if "x_ticklabel_fs" not in kwarg_keys:
                x_ticklabel_fs = None
            else:
                if kwargs['x_ticklabel_fs'] == "by default":
                    x_ticklabel_fs = None
                else:
                    x_ticklabel_fs = kwargs['x_ticklabel_fs']

            if "y_ticklabel_fs" not in kwarg_keys:
                y_ticklabel_fs = None
            else:
                if kwargs['y_ticklabel_fs'] == "by default":
                    y_ticklabel_fs = None
                else:
                    y_ticklabel_fs = kwargs['y_ticklabel_fs']

            if "x_label_fs" not in kwarg_keys:
                x_label_fs = None
            else:
                if kwargs['x_label_fs'] == "by default":
                    x_label_fs = None
                else:
                    x_label_fs = kwargs['x_label_fs']

            if "y_label_fs" not in kwarg_keys:
                y_label_fs = None
            else:
                if kwargs['y_label_fs'] == "by default":
                    y_label_fs = None
                else:
                    y_label_fs = kwargs['y_label_fs']

            if "title_fs" not in kwarg_keys:
                title_fs = None
            else:
                if kwargs['title_fs'] == "by default":
                    title_fs = None
                else:
                    title_fs = kwargs['title_fs']

            if "x_format" not in kwarg_keys:
                x_format = None
            else:
                if kwargs['x_format'] == "by default":
                    x_format = None
                else:
                    x_format = kwargs['x_format']

            if "y_format" not in kwarg_keys:
                y_format = None
            else:
                if kwargs['y_format'] == "by default":
                    y_format = None
                else:
                    y_format = kwargs['y_format']

            if "decoration_mark" not in kwarg_keys:
                decoration_mark = False
            else:
                if kwargs['decoration_mark'] == "by default":
                    decoration_mark = False
                else:
                    decoration_mark = kwargs['decoration_mark']

            if "marker_list" not in kwarg_keys:
                marker_list = None
            else:
                if kwargs['marker_list'] == "by default":
                    marker_list = None
                else:
                    marker_list = kwargs['marker_list']

            if "marker_label_list" not in kwarg_keys:
                marker_label_list = None
            else:
                if kwargs['marker_label_list'] == "by default":
                    marker_label_list = None
                else:
                    marker_label_list = kwargs['marker_label_list']

            if "marker" not in kwarg_keys:
                marker = None
            else:
                if kwargs['marker'] == "by default":
                    marker = None
                else:
                    marker = kwargs['marker']

            if "marker_size" not in kwarg_keys:
                marker_size = None
            else:
                if kwargs['marker_size'] == "by default":
                    marker_size = None
                else:
                    marker_size = kwargs['marker_size']

            if "xl_mark" not in kwarg_keys:
                xl_mark = False
            else:
                if kwargs['xl_mark'] == "by default":
                    xl_mark = False
                else:
                    xl_mark = kwargs['xl_mark']

            if "yl_mark" not in kwarg_keys:
                yl_mark = False
            else:
                if kwargs['yl_mark'] == "by default":
                    yl_mark = False
                else:
                    yl_mark = kwargs['yl_mark']

            if "line_color" not in kwarg_keys:
                line_color = None
            else:
                if kwargs['line_color'] == "by default":
                    line_color = None
                else:
                    line_color = kwargs['line_color']

            if "marker_color" not in kwarg_keys:
                marker_color = None
            else:
                if kwargs['marker_color'] == "by default":
                    marker_color = None
                else:
                    marker_color = kwargs['marker_color']

            if "marker_edge_width" not in kwarg_keys:
                marker_edge_width = None
            else:
                if kwargs['marker_edge_width'] == "by default":
                    marker_edge_width = None
                else:
                    marker_edge_width = kwargs['marker_edge_width']

            if "marker_face_color" not in kwarg_keys:
                marker_face_color = 'none'
            else:
                if kwargs['marker_face_color'] == "by default":
                    marker_face_color = 'none'
                else:
                    marker_face_color = kwargs['marker_face_color']

            if "legend_fs" not in kwarg_keys:
                legend_fs = 'none'
            else:
                if kwargs['legend_fs'] == "by default":
                    legend_fs = 'none'
                else:
                    legend_fs = kwargs['legend_fs']

            if "legend_loc" not in kwarg_keys:
                legend_loc = None
            else:
                if kwargs['legend_loc'] == "by default":
                    legend_loc = None
                else:
                    legend_loc = kwargs['legend_loc']

            canvas = draw(
                ax=kwargs['ax'],
                x=kwargs['x'],
                y=kwargs['y'],
                xl_mark=xl_mark,
                yl_mark=yl_mark,
                title=title,
                marker=marker,
                line_color=line_color,
                marker_color=marker_color,
                marker_size=marker_size,
                decoration_mark=decoration_mark,
                marker_list=marker_list,
                marker_label_list=marker_label_list,
                marker_edge_width=marker_edge_width,
                marker_face_color=marker_face_color,
                legend_loc=legend_loc,
                x_label=x_label,
                y_label=y_label,
                color=color,
                width=width,
                edgecolor=edgecolor,
                facecolor=facecolor,
                linewidth=linewidth,
                label=label,
                x_format=x_format,
                y_format=y_format,
                x_ticks=x_ticks,
                y_ticks=y_ticks,
                x_ticklabels=x_ticklabels,
                y_ticklabels=y_ticklabels,

                x_ticklabel_fs=x_ticklabel_fs,
                x_label_rotation=x_label_rotation,
                x_label_rotation_align=x_label_rotation_align,
                x_label_fs=x_label_fs,
                y_ticklabel_fs=y_ticklabel_fs,
                y_label_fs=y_label_fs,
                title_fs=title_fs,
                legend_fs=legend_fs,

                alpha=alpha,
            )
            params_all["canvas"] = canvas
            return params_all
        return config

    def bar_individual(self, **kwargs):
        canvas = kwargs["ax"].bar(
            x=kwargs["x"],
            height=kwargs["y"],
            width=kwargs["width"],
            color=kwargs["color"],
            label=kwargs["label"],
            edgecolor=kwargs["edgecolor"],
            facecolor=kwargs["facecolor"],
            linewidth=kwargs["linewidth"],
            alpha=kwargs["alpha"],
        )
        if kwargs["y_format"] == 'int':
            kwargs["ax"].yaxis.set_major_formatter(FuncFormatter(lambda y, _: '{:d}'.format(int(y))))
        elif kwargs["y_format"] == 'pct':
            kwargs["ax"].yaxis.set_major_formatter(FuncFormatter(lambda y, _: '{:.0%}'.format(y)))
        elif kwargs["y_format"] == 'exp':
            kwargs["ax"].yaxis.set_major_formatter(FuncFormatter(lambda y, _: '{:.2e}'.format(y)))
        else:
            kwargs["ax"].yaxis.set_major_formatter(FuncFormatter(lambda y, _: '{}'.format(y)))
        kwargs["ax"].spines['right'].set_color('none')
        kwargs["ax"].spines['top'].set_color('none')
        if kwargs["x_ticks"]:
            kwargs["ax"].set_xticks(kwargs["x"])
        if kwargs["y_ticks"]:
            kwargs["ax"].set_yticks(kwargs["y_ticks"])
        if kwargs["x_ticklabels"] is not None or kwargs["x_ticklabels"] == []:
            kwargs["ax"].set_xticklabels(kwargs["x_ticklabels"], rotation=15, ha='right', fontsize=kwargs["x_ticklabel_fs"])
        if kwargs["y_ticklabels"] is not None or kwargs["y_ticklabels"] == []:
            kwargs["ax"].set_yticklabels(kwargs["y_ticklabels"], fontsize=kwargs["y_ticklabel_fs"])
        if kwargs["xl_mark"]:
            kwargs["ax"].set_xlabel(kwargs["x_label"], fontsize=kwargs["x_label_fs"])
            kwargs['ax'].set_xticklabels(
                kwargs['x'],
                rotation=kwargs['x_label_rotation'],
                ha=kwargs['x_label_rotation_align'],
                fontsize=kwargs['x_ticklabel_fs'],
            )
        if kwargs["yl_mark"]:
            kwargs["ax"].set_ylabel(kwargs["y_label"], fontsize=kwargs["y_label_fs"])
        kwargs["ax"].set_title(kwargs["title"], fontsize=kwargs["title_fs"])
        return canvas

    def hist_grouped(self, **kwargs):
        # tc = plt.cm.tab20b(np.linspace(0, 1, len(y.keys()))).tolist()
        import seaborn as sns
        tc = sns.color_palette("Set2")
        for i, (label, y1) in enumerate(kwargs['y'].items()):
            kwargs['ax'].plot(
                kwargs['x'],
                y1,
                label=label,
                c=tc[i],
                alpha=kwargs['alpha'],
            )
        kwargs['ax'].spines['right'].set_color('none')
        kwargs['ax'].spines['top'].set_color('none')
        # if kwargs['x_ticks']:
        #     kwargs['ax'].set_xticks(kwargs['x_ticks'])
        # if kwargs['y_ticks']:
        #     kwargs['ax'].set_yticks(kwargs['y_ticks'])
        # if kwargs['x_ticklabels'] is not None or kwargs['x_ticklabels'] == []:
        #     kwargs['ax'].set_xticklabels(kwargs['x_ticklabels'], rotation=15, ha='right', fontsize=kwargs['x_ticklabel_fs'])
        # if kwargs['y_ticklabels'] is not None or kwargs['y_ticklabels'] == []:
        #     kwargs['ax'].set_yticklabels(kwargs['y_ticklabels'], fontsize=kwargs['y_ticklabel_fs'])
        if kwargs['xl_mark']:
            print(True)
            kwargs['ax'].set_xlabel(kwargs['x_label'], fontsize=kwargs['x_label_fs'])
            kwargs['ax'].set_xticks(kwargs['x_ticks'])
            kwargs['ax'].set_xticklabels(
                kwargs['x_ticklabels'],
                rotation=kwargs['x_label_rotation'],
                ha=kwargs['x_label_rotation_align'],
                fontsize=kwargs['x_ticklabel_fs'],
            )
        if kwargs['yl_mark']:
            kwargs['ax'].set_ylabel(kwargs['y_label'], fontsize=kwargs['y_label_fs'])
        kwargs['ax'].set_title(kwargs['title'], fontsize=kwargs['title_fs'])
        kwargs['ax'].legend(loc=kwargs['legend_loc'], ncol=1, fontsize=kwargs['legend_fs'])
        return

    def line(self, **kwargs):
        kwargs['ax'].plot(
            kwargs['x'],
            kwargs['y'],
            color=kwargs['line_color'],
            label=kwargs['label'],  # " ".join(ds_key.split("_"))
            linewidth=kwargs['linewidth'],
            marker=kwargs['marker'],
            markersize=kwargs['marker_size'],
            markeredgewidth=kwargs['marker_edge_width'],
            markerfacecolor=kwargs['marker_face_color'],
            # linestyle=line_styles[configs[ds_key]['position']],
        )
        kwargs['ax'].spines['right'].set_color('none')
        kwargs['ax'].spines['top'].set_color('none')
        if kwargs['decoration_mark']:
            kwargs['ax'].set_xlabel(kwargs['x_label'], fontsize=kwargs['x_label_fs'])
            kwargs['ax'].set_ylabel(kwargs['y_label'], fontsize=kwargs['y_label_fs'])
            if kwargs['title']:
                kwargs['ax'].set_title(kwargs['title'], fontsize=kwargs['title_fs'])
            kwargs['ax'].legend(loc=kwargs['legend_loc'], ncol=1, fontsize=kwargs['legend_fs'])
            kwargs['ax'].tick_params(axis='both', which='major', labelsize=14)
            kwargs['ax'].tick_params(axis='both', which='minor', labelsize=14)

        if kwargs["x_ticks"] and ["xl_mark"] is True:
            kwargs["ax"].set_xticks(kwargs["x_ticks"])
        if kwargs["y_ticks"] and ["xl_mark"] is True:
            kwargs["ax"].set_yticks(kwargs["y_ticks"])
        if kwargs["x_ticklabels"] is not None or kwargs["x_ticklabels"] == []:
            kwargs["ax"].set_xticklabels([], rotation=15, ha='right', fontsize=kwargs["x_ticklabel_fs"])
        if kwargs["y_ticklabels"] is not None or  kwargs["y_ticklabels"] == []:
            kwargs["ax"].set_yticklabels(kwargs["y_ticklabels"], fontsize=kwargs["y_ticklabel_fs"])
        if kwargs["xl_mark"]:
            kwargs["ax"].set_xlabel(kwargs["x_label"], fontsize=kwargs["x_label_fs"])
            kwargs["ax"].set_xticks(kwargs["x_ticks"])
            kwargs['ax'].set_xticklabels(
                kwargs['x_ticklabels'],
                rotation=kwargs['x_label_rotation'],
                ha=kwargs['x_label_rotation_align'],
                fontsize=kwargs['x_ticklabel_fs'],
            )
        if kwargs["yl_mark"]:
            kwargs["ax"].set_ylabel(kwargs["y_label"], fontsize=kwargs["y_label_fs"])
        kwargs["ax"].set_title(kwargs["title"], fontsize=kwargs["title_fs"])
        kwargs['ax'].legend(
            loc=kwargs['legend_loc'],
            ncol=1,
            fontsize=kwargs['legend_fs'],
            frameon=False,
        )

    def line_scatter(self, **kwargs):
        kwargs['ax'].plot(
            kwargs['x'],
            kwargs['y'],
            color=kwargs['line_color'],
            label=kwargs['label'],
            linewidth=kwargs['linewidth'],
            alpha=kwargs['alpha'],
            # marker="o" if configs[key]['position'] == "Cells" else "*",
            # markersize=5,
            # linestyle=line_styles[configs[key]['position']],
        )
        kwargs['ax'].scatter(
            kwargs['x'],
            kwargs['y'],
            color=kwargs['marker_color'],
            marker=kwargs['marker'],
            s=kwargs['marker_size'],
            facecolors='none',
            alpha=kwargs['alpha'],
        )
        if kwargs['decoration_mark']:
            print(True)
            kwargs['ax'].spines['right'].set_color('none')
            kwargs['ax'].spines['top'].set_color('none')
            kwargs["ax"].set_xticks(kwargs["x"])
            kwargs['ax'].set_xticklabels(
                kwargs['x'],
                rotation=kwargs['x_label_rotation'],
                ha=kwargs['x_label_rotation_align'],
                fontsize=kwargs['x_ticklabel_fs'],
            )

            kwargs['ax'].set_title(kwargs['title'], fontsize=kwargs['title_fs'])
            kwargs['ax'].set_xlabel(kwargs['x_label'], fontsize=kwargs['x_label_fs'])
            kwargs['ax'].set_ylabel(kwargs['y_label'], fontsize=kwargs['y_label_fs'])
            legend_main = kwargs['ax'].legend(
                loc='center left',
                ncol=1,
                bbox_to_anchor=(1, 0.5),
                frameon=False,
            )
            kwargs['ax'].add_artist(legend_main)
            import matplotlib.lines as mlines
            handles1 = [
                mlines.Line2D(
                    [],
                    [],
                    marker=marker,
                    markeredgecolor='k',
                    mfc='w',
                    ls='',
                ) for marker in kwargs['marker_list']
            ]
            legend_sub1 = kwargs['ax'].legend(
                handles1,
                kwargs['marker_label_list'],
                ncols=2,
                loc=kwargs['legend_loc'],
                fontsize=kwargs['legend_fs'],
            )
            kwargs['ax'].add_artist(legend_sub1)
            # handles2 = [mlines.Line2D([], [], linestyle=line_style, color='grey', mfc='w',) for line_style in [*line_styles.values()]]
            # ax[i].legend(handles2, [*line_styles.keys()], loc='lower right', ncols=3, fontsize=10)

            # handles = [mlines.Line2D([], [], marker=marker, mec='k', mfc='w', ls='') for marker in ['o', '^']]
            # ax[i].legend(handles, ['Radial', 'Transit'], loc='upper left', title="Detection")
        return