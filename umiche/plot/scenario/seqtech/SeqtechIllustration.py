__version__ = "v1.0"
__copyright__ = "Copyright 2025"
__license__ = "GPL-3.0"
__developer__ = "Jianfeng Sun"
__maintainer__ = "Jianfeng Sun"
__email__="jianfeng.sunmt@gmail.com"

import colorsys
import matplotlib.colors as mcolors
import matplotlib.pyplot as plt
import seaborn as sns
from matplotlib.patches import Rectangle, FancyArrowPatch

sns.set(font="Helvetica")
sns.set_style("ticks")

# base = sns.color_palette("Set3", 8)
base = sns.color_palette("Set3", 11)[3:]


def to_hex(c):
    return mcolors.to_hex(c)

def adjust_lightness(color, amount=1.0):
    r,g,b = mcolors.to_rgb(color)
    h,l,s = colorsys.rgb_to_hls(r,g,b)
    l = max(0, min(1, l * amount))
    r,g,b = colorsys.hls_to_rgb(h,l,s)
    return (r,g,b)

PALETTE = {
    'umi'     : to_hex(base[2]),
    'barcode' : to_hex(base[4]),
    'barcode2': to_hex(adjust_lightness(base[4], 0.85)),
    'barcode3': to_hex(adjust_lightness(base[4], 0.70)),
    'polyT'   : to_hex(base[5]),
    'spacer'  : to_hex(base[1]),
    'primer'  : to_hex(base[0]),
    'tail'    : to_hex(base[6]),
    'rbG'     : to_hex(base[7]),
    'constant': to_hex(base[7]),
    'insert'  : '#9E9E9E',
    'other'   : to_hex(base[3]),
}


def family_and_color(seg_name: str, idx_in_family: int = 1):
    """barcode1/2/3, dT(s), spacer, primer, tail, rb(G)s"""
    n = seg_name.lower()
    if n.startswith('barcode'):
        if n.startswith('barcode1'):
            return 'barcode', PALETTE['barcode']
        if n.startswith('barcode2'):
            return 'barcode', PALETTE['barcode2']
        if n.startswith('barcode3'):
            return 'barcode', PALETTE['barcode3']
        return 'barcode', PALETTE['barcode']
    if 'umi' in n:
        return 'umi', PALETTE['umi']
    if 'dt(' in n or 'polyt' in n:
        return 'polyT', PALETTE['polyT']
    if 'spacer' in n:
        return 'spacer', PALETTE['spacer']
    if 'primer' in n:
        return 'primer', PALETTE['primer']
    if 'tail' in n:
        return 'tail', PALETTE['tail']
    if 'rb(' in n or 'ggg' in n:
        return 'rbG', PALETTE['rbG']
    if 'const' in n or 'linker' in n:
        return 'constant', PALETTE['constant']
    return 'other', PALETTE['other']


def draw_one(ax, tech_name: str, scheme: dict):
    BAR_X, BAR_Y, BAR_H, BAR_W = 0.08, 0.20, 0.40, 0.84
    MIN_W_FOR_INSIDE = 0.10
    NAME_FS = 10
    LEN_FS  = 8
    ABOVE_DY = 0.05
    BELOW_DY = 0.06
    TITLE_Y = 0.96
    ARROW_PAD = 0.06   # dist from arrow to compo

    total = sum(int(v.get('len', 0)) for v in scheme.values()) or 1

    x, y, h, W = BAR_X, BAR_Y, BAR_H, BAR_W
    bar_start = x

    for seg_name, seg in scheme.items():
        L = int(seg.get('len', 0)) or 1
        w = W * L / total
        fam, color = family_and_color(seg_name)

        ax.add_patch(Rectangle((x, y), w, h, facecolor=color, edgecolor='black', linewidth=1.0))

        if w >= MIN_W_FOR_INSIDE:
            ax.text(x + w/2, y + h/2, seg_name, ha='center', va='center',
                    fontsize=NAME_FS, color='black')
        else:
            ax.text(x + w/2, y + h + ABOVE_DY, seg_name, ha='center', va='bottom',
                    fontsize=NAME_FS-1, color='black')
        # bp put below the bar
        ax.text(x + w/2, max(0.02, y - BELOW_DY), f"{L} bp",
                ha='center', va='top', fontsize=LEN_FS, color='black')
        x += w

    bar_end = x
    mid_y = y + h/2

    # title stays far above
    ax.text(0.02, TITLE_Y, tech_name, fontsize=12, ha='left', va='center', fontweight='bold')

    # 5'… ——> start
    left_text_x = max(0.02, bar_start - ARROW_PAD)
    ax.add_patch(FancyArrowPatch(
        (bar_start - ARROW_PAD, mid_y), (bar_start, mid_y),
        arrowstyle='->', mutation_scale=14, lw=1.4, color='grey',
        zorder=10, clip_on=False
    ))
    ax.text(left_text_x - 0.005, mid_y, "5′ …", ha="right", va="center", fontsize=10)

    # end ——> …3′
    right_text_x = min(0.98, bar_end + ARROW_PAD)
    ax.add_patch(FancyArrowPatch(
        (bar_end, mid_y), (bar_end + ARROW_PAD, mid_y),
        arrowstyle='->', mutation_scale=14, lw=1.4, color='grey',
        zorder=10, clip_on=False
    ))
    ax.text(right_text_x + 0.005, mid_y, "… 3′", ha="left", va="center", fontsize=10)

    ax.set_xlim(0, 1); ax.set_ylim(0, 1); ax.axis('off')


def pl(params, ):
    techs = list(params.keys())
    # ncols, nrows = 2, 10
    # fig, axes = plt.subplots(nrows, ncols, figsize=(10, 8))
    ncols, nrows = 4, 5
    fig, axes = plt.subplots(nrows, ncols, figsize=(19, 5))
    plt.subplots_adjust(
        left=0.04,
        right=0.98,
        top=0.98,
        bottom=0.04,
        wspace=0.12,
        hspace=0.35,
    )
    for i, tech in enumerate(techs):
        r = i % nrows
        # first left and then right
        c = i // nrows
        draw_one(axes[r, c], tech, params[tech])
    # legend_fig, legend_ax = plt.subplots(figsize=(10, 1.2))
    # legend_ax.axis('off')
    # items = [
    #     ('UMI', PALETTE['umi']),
    #     ('Barcode', PALETTE['barcode']),
    #     ('Barcode2', PALETTE['barcode2']),
    #     ('Barcode3', PALETTE['barcode3']),
    #     ('polyT', PALETTE['polyT']),
    #     ('Spacer', PALETTE['spacer']),
    #     ('Primer', PALETTE['primer']),
    #     ('Tail', PALETTE['tail']),
    #     ('rb(G)s', PALETTE['rbG']),
    # ]
    # x0 = 0.02
    # for label, col in items:
    #     legend_ax.add_patch(Rectangle((x0, 0.2), 0.08, 0.6, facecolor=col, edgecolor='black'))
    #     legend_ax.text(x0 + 0.04, 0.1, label, ha='center', va='bottom', fontsize=10)
    #     x0 += 0.10
    # legend_fig.savefig("read1_legend.png", dpi=200, bbox_inches="tight")
    plt.show()


if __name__ == "__main__":
    params = {
        'VASA-seq': {
            'UMI': {'len': 6},
            'barcode': {'len': 6},
            'dT(s)': {'len': 24, 'seq': 'TTTTTTTTTTTTTTTTTTTTTTTT'},
        },
        'PIP-seq': {
            'BC1': {'len': 8},
            'spacer1': {'len': 7, 'seq': 'ATGCATC'},
            'BC2': {'len': 8},
            'spacer2': {'len': 7, 'seq': 'CCTCGAG'},
            'BC3': {'len': 8},
            'UMI': {'len': 12},
            'dT(s)': {'len': 19, 'seq': 'TTTTTTTTTTTTTTTTTTT'},
        },
        'SCIFI-seq': {
            'UMI': {'len': 8},
            'barcode': {'len': 13},
            'dT(s)': {'len': 30, 'seq': 'TTTTTTTTTTTTTTTTTTTTTTTTTTTTTT'},
        },
        'SCRB-seq': {
            'BC': {'len': 6},
            'UMI': {'len': 10},
            'dT(s)': {'len': 30, 'seq': 'TTTTTTTTTTTTTTTTTTTTTTTTTTTTTT'},
        },
        'SHARE-seq': {
            'UMI': {'len': 10},
            'BiodT/(dT)': {'len': 10},
        },
        'SNARE-seq': {
            'barcode': {'len': 12},
            'UMI': {'len': 8},
            'dT(s)': {'len': 30, 'seq': 'TTTTTTTTTTTTTTTTTTTTTTTTTTTTTT'},
        },
        'SPLiT-seq': {
            'UMI': {'len': 10},
            'barcode': {'len': 8},
            'dT(s)': {'len': 15, 'seq': 'GTGGCCGATGTTTCG'},
        },
        'STRT-seq-C1': {
            'UMI': {'len': 5},
            'GGG': {'len': 3},
        },
        'STRT-seq-2i': {
            'UMI': {'len': 6},
            'GGG': {'len': 3},
        },
        'Quartz-seq2': {
            'barcode': {'len': 15},
            'UMI': {'len': 8},
            'dT(s)': {'len': 24, 'seq': 'TTTTTTTTTTTTTTTTTTTTTTTT'},
        },
        'PETRI-seq': {
            'UMI': {'len': 7},
            'barcode': {'len': 7},
            'tail': {'len': 8, 'seq': 'GGTCCTTG'},
        },
        'MARS-seq2': {
            'barcode': {'len': 6},
            'UMI': {'len': 4},
            'dT(s)': {'len': 20, 'seq': 'TTTTTTTTTTTTTTTTTTTT'},
        },
        'Paired-seq': {
            'primer': {'len': 22, 'seq': 'TCTAGCCTTCTCGTGTGCAGAC'},
            'UMI': {'len': 10},
            'barcode': {'len': 7},
        },
        'ISSAAC-seq': {
            'UMI': {'len': 10},
            'dT(s)': {'len': 30, 'seq': 'TTTTTTTTTTTTTTTTTTTTTTTTTTTTTT'},
        },
        'inDrop': {
            'barcode': {'len': 8},
            'UMI': {'len': 6},
            'dT(s)': {'len': 18, 'seq': 'TTTTTTTTTTTTTTTTTT'},
        },
        'Drop-seq': {
            'barcode': {'len': 12},
            'UMI': {'len': 8},
            'dT(s)': {'len': 30},
        },
        '10X Chromium V2': {
            'barcode': {'len': 16},
            'UMI': {'len': 10},
            'dT(s)': {'len': 30, 'seq': 'TTTTTTTTTTTTTTTTTTTTTTTTTTTTTT'},
        },
        '10X Chromium V3': {
            'barcode': {'len': 16},
            'UMI': {'len': 12},
            'dT(s)': {'len': 30, 'seq': 'TTTTTTTTTTTTTTTTTTTTTTTTTTTTTT'},
        },
        'CEL-seq2': {
            'UMI': {'len': 6},
            'barcode': {'len': 6},
            'dT(s)': {'len': 24, 'seq': 'TTTTTTTTTTTTTTTTTTTTTTTT'},
        },
        'FLASH-seq-UMI': {
            'UMI': {'len': 8},
            'spacer': {'len': 5, 'seq': 'CTAAC'},
            'rb(G)s': {'len': 3, 'seq': 'GGG'},
        },
    }
    pl(params)
