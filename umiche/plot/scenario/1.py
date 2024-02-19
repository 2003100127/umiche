import numpy as np
import matplotlib.pyplot as plt
import scipy.stats as st

from palmerpenguins import load_penguins

penguins = load_penguins()
print(penguins)
penguins = penguins.dropna()
print(penguins)


# Get the species, sorted alphabetically
species = sorted(penguins["species"].unique())

# y_data is a list of length 3 containing the bill_length_mm values for each specie
y_data = [penguins[penguins["species"] == specie]["bill_length_mm"].values for specie in species]

# Create jittered version of "x" (which is only 0, 1, and 2)
# More about this in the bonus track!
x_data = [np.array([i] * len(d)) for i, d in enumerate(y_data)]
print(x_data)
x_jittered = [x + st.t(df=6, scale=0.04).rvs(len(x)) for x in x_data]
print(x_jittered)
# Colors
BG_WHITE = "#fbf9f4"
GREY_LIGHT = "#b4aea9"
GREY50 = "#7F7F7F"
BLUE_DARK = "#1B2838"
BLUE = "#2a475e"
BLACK = "#282724"
GREY_DARK = "#747473"
RED_DARK = "#850e00"

# Colors taken from Dark2 palette in RColorBrewer R library
COLOR_SCALE = ["#1B9E77", "#D95F02", "#7570B3"]

# Horizontal positions for the violins.
# They are arbitrary numbers. They could have been [-1, 0, 1] for example.
POSITIONS = [0, 1, 2]

# Horizontal lines
HLINES = [40, 50, 60]

fig, ax = plt.subplots(figsize=(14, 10))

# Some layout stuff ----------------------------------------------
# Background color
# fig.patch.set_facecolor(BG_WHITE)
# ax.set_facecolor(BG_WHITE)

# Horizontal lines that are used as scale reference
# for h in HLINES:
#     ax.axhline(h, color=GREY50, ls=(0, (5, 5)), alpha=0.8, zorder=0)

# Add violins ----------------------------------------------------
# bw_method="silverman" means the bandwidth of the kernel density
# estimator is computed via Silverman's rule of thumb.
# More on this in the bonus track ;)

# The output is stored in 'violins', used to customize their appearence
print(y_data)
violins = ax.violinplot(
    y_data,
    positions=POSITIONS,
    widths=0.45,
    bw_method="silverman",
    showmeans=False,
    showmedians=False,
    showextrema=False
)

# Customize violins (remove fill, customize line, etc.)
for pc in violins["bodies"]:
    pc.set_facecolor("none")
    pc.set_edgecolor(BLACK)
    pc.set_linewidth(1.4)
    pc.set_alpha(1)

# Add boxplots ---------------------------------------------------
# Note that properties about the median and the box are passed
# as dictionaries.

medianprops = dict(
    linewidth=4,
    color=GREY_DARK,
    solid_capstyle="butt"
)
boxprops = dict(
    linewidth=2,
    color=GREY_DARK
)

ax.boxplot(
    y_data,
    positions=POSITIONS,
    showfliers=False,  # Do not show the outliers beyond the caps.
    showcaps=False,  # Do not show the caps
    medianprops=medianprops,
    whiskerprops=boxprops,
    boxprops=boxprops
)

# Add jittered dots ----------------------------------------------
for x, y, color in zip(x_jittered, y_data, COLOR_SCALE):
    ax.scatter(x, y, s=100, color=color, alpha=0.4)

plt.show()