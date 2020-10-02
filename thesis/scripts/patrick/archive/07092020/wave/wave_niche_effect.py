import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import seaborn as sns
import random
from my_debug import message
from scipy import spatial

#load the df

no_of_runs = 10
base_path = "/extra/brunner/thesis/kinetic/q_fraction_medium_g_scan_multi/"
cell_cell_distance = 20
gammas = [0.1, 0.5, 2, 10]

avg_hist_df = pd.read_hdf(base_path + str(no_of_runs) + '_runs_avg_hist_df.h5', mode="r")
message("loaded avg_hist_df from " + base_path)


#############################################################################################################################################

# load average secretory distance and norm to cell-cell-distance
mean_sec_dist = float(np.loadtxt(base_path + str(no_of_runs) + "_runs_avg_sec_dist.txt"))/cell_cell_distance

niche_size_at_gamma = []

# get first minimum as niche size
for i, gamma in enumerate(gammas):
    for i in range(len(avg_hist_df[gamma][39])):
        try:
            if avg_hist_df[gamma][39][i] - avg_hist_df[gamma][39][i+1] < 0:
                niche_size_at_gamma.append(i)
                break
        except IndexError:
            niche_size_at_gamma.append(i)

#take the last timepoints user defined niche and get the average concentration within that niche

niche_effect_at_gamma = []
for i,gamma in enumerate(gammas):
    hist = avg_hist_df[gamma][39]
    avg_c = []
    for k in range(niche_size_at_gamma[i]):
        avg_c.append(hist[k])
    niche_effect_at_gamma.append(np.mean(avg_c))

# ax13 = sns.lineplot(points, niche, label="hist_niche", legend=False)
# ax13 = sns.lineplot(points[:], AC_niche[:], label="AC_niche")
# ax13.set(xlabel="time", ylabel="cell-cell-distance", yscale="linear", xscale="linear", title="Niche size at threshold " + str(threshold))
# plt.show()

#############################################################################################################################################

# Plotting
#

# sns.set(rc={'figure.figsize':(12,10)})
# sns.set(palette="Greens")
# sns.set_style("ticks")
# sns.set_context("talk", font_scale=1, rc={"lines.linewidth": 3})


# sns.set_context("talk", font_scale=1.7, rc={"lines.linewidth": 5, 'lines.markersize': 10})

# fig.savefig("plots/" + "histogram" + ".png", bbox_inches='tight')


palette = sns.color_palette("bwr", 4)
#
x = niche_effect_at_gamma
y = niche_size_at_gamma
 # fig.add_subplot(a_x,a_y,3)
# ax13 = sns.lineplot(np.array(points)*10, niche, label="hist_niche", legend=False)
# ax13 = sns.lineplot(np.array(points[:])*10, AC_niche[:], label="AC_niche")
# ax13.set(xlabel="time in h", ylabel="cell-cell-distance", yscale=yscale, xscale=xscale, title="Niche size at threshold " + str(round(threshold,1)))
#
# fig.add_subplot(a_x,a_y,4)
# plt.plot(np.array(points)*10, niche_avg_c, label="niche avg. c")
# fig.savefig("/home/brunner/Documents/Current work/28082020/" + "niche_" + "f=" + str(cell_df["IL-2_fraction"].unique()[0]) + "_t=" + str(threshold) + ".pdf", bbox_inches='tight')
# fig.savefig("/home/brunner/Documents/Current work/28082020/" + "niche_" + "f=" + str(cell_df["IL-2_fraction"].unique()[0]) + "_t=" + str(threshold) + ".png", bbox_inches='tight')

for i in range(len(x)):
    plt.scatter(x[i], y[i], color=palette[3-i])
plt.xlabel("niche effect")
plt.ylabel("niche size")
import matplotlib.lines as mlines

zero_line = mlines.Line2D([], [], color=palette[0], label='lightcoral line')
one_line = mlines.Line2D([], [], color=palette[1], label='red line')
two_line = mlines.Line2D([], [], color=palette[2], label='lightsteelblue line')
three_line = mlines.Line2D([], [], color=palette[3], label='blue line')

grey_line = mlines.Line2D([], [], color="grey", label='grey line')
grey_dashed_line = mlines.Line2D([], [], color="grey", label='grey dashed line', linestyle="dashed")
white_line = mlines.Line2D([], [], color="white", label='white line')


new_handles = [three_line, zero_line]
new_labels = ["neg. feedback", "pos. feedback"]
# new_handles = [three_line, two_line, black_line, one_line, zero_line]
# new_labels = ["strong negative feedback", "negative feedback", "no feedback", "positive feedback", "strong positive feedback"]
plt.legend(handles=new_handles, labels=new_labels, loc='center right', fancybox=True)

plt.show()
