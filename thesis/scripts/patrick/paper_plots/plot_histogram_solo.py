import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import seaborn as sns
import random
from my_debug import message
from scipy import spatial

no_of_runs = 10
timepoints = [0,2,5,10,20]
gammas = [0.1, 0.5, 2, 10]
cell_cell_distance = 20
get_dataframes = []

#create multiindex df
array1 = np.array([[x]*len(gammas) for x in range(no_of_runs)]).flatten()
array2 = np.array([[g for g in gammas]*no_of_runs]).flatten()
df_columns = [array1,array2]

df_indices = timepoints.copy()
df_indices.append("sec_dists")

runs_hist_df = pd.DataFrame(columns=df_columns, index=df_indices)
base_path = "/extra/brunner/thesis/kinetic/q_fraction_medium_g_scan_multi/"

try:
    runs_hist_df = pd.read_hdf(base_path + str(no_of_runs) + "_runs_hist_df.h5", mode="r")
    message("loaded runs_hist_df from " + base_path)
except FileNotFoundError:
    # message("runs_hist_df not found!")
    assert runs_hist_df[0][gammas[0]].isnull().values[0] == False, "runs_hist_df not found/loaded!"

#############################################################################################################################################

# Avergaging runs. columns = gamma values, indices = timepoints
# Transpose df for easier handling
transposed_runs_hist_df = runs_hist_df.transpose().copy()
avg_hist_df = pd.DataFrame(columns=gammas, index=timepoints)

for timepoint in timepoints:
    for gamma in gammas:
        avg_hist = np.zeros(len(transposed_runs_hist_df[timepoint][0][gamma])) # taking the length of the zero run as reference
        temp_runs_list = []
        for run in range(transposed_runs_hist_df.index[-1][0]+1):
            temp_runs_list.append(transposed_runs_hist_df[timepoint][run][gamma])
        temp_runs_list = np.transpose(temp_runs_list)
        for j,entry in enumerate(temp_runs_list):
            avg_hist[j] = np.mean(entry)

            avg_hist_df[gamma][timepoint] = avg_hist

# Averaging secretory cells distances
temp_sec_dists = []
for i in range(transposed_runs_hist_df["sec_dists"].index[-1][0]+1):
    temp_sec_dists.append(transposed_runs_hist_df["sec_dists"][i][0.1])
temp_sec_dists = [item for sublist in temp_sec_dists for item in sublist]
sec_dists = np.mean(temp_sec_dists)
# Write averaged dfs to file
if save_avg_hist_df == True:
    np.savetxt(base_path + str(no_of_runs) + "_runs_avg_sec_dist.txt", [sec_dists])
    avg_hist_df.to_hdf(base_path + str(no_of_runs) + '_runs_avg_hist_df.h5', key="data", mode="w")



#############################################################################################################################################

# Plotting

# sns.set(rc={'figure.figsize':(12,10)})
# sns.set(palette="Greens")
# sns.set_style("ticks")
# sns.set_context("talk", font_scale=1, rc={"lines.linewidth": 3})

fig = plt.figure(figsize=(8,8))
# plt.figure(figsize=(12,10))
# sns.set_context("talk", font_scale=1.1, rc={"lines.linewidth": 2.5})
# sns.set(rc={'figure.figsize':(8,8)})
sns.set_style("ticks")
sns.set_context("talk", font_scale=1.5, rc={"lines.linewidth": 5})
fig,ax = plt.subplots(figsize=(6,5))
# plt.subplots_adjust(wspace=.4, hspace=0.6)
a_x = 2
a_y = 2

xscale = "linear"
yscale = "linear"

mean_sec_dist = np.mean(sec_dists)/cell_cell_distance

#
# fig = plt.figure(figsize=(9,8))
# sns.set_context("talk", font_scale=1.7, rc={"lines.linewidth": 5})
#
# # plt.figure(figsize=(12,10))
# # sns.set_context("talk", font_scale=1.1, rc={"lines.linewidth": 2.5})
# # sns.set(rc={'figure.figsize':(8,8)})
# sns.set_style("ticks")
# sns.set_context("talk", font_scale=1.7, rc={"lines.linewidth": 5, 'lines.markersize': 10})


time = 20
palette = sns.color_palette("bwr", 4)
gammas = [10]

from scipy.optimize import curve_fit

def func(x, a,b,c):
    return a*np.exp(-x/b) + c

mean_sec_dist = float(np.loadtxt(base_path + str(no_of_runs) + "_runs_avg_sec_dist.txt"))/cell_cell_distance


sns.set_context(rc={"lines.linewidth": 3, 'lines.markersize': 6.5})
for i,gamma in enumerate(gammas):
    for run in [3,2,0,1]: #range(10):
        sns.lineplot(np.arange(0,len(runs_hist_df.loc[20,(run,gamma)][:-1])), runs_hist_df.loc[20,(run,gamma)][:-1], label="strong neg.", color=palette[0], marker="o", legend=None)
        # plt.axvline(1,ymin = 0,ymax = (runs_hist_df.loc[20,(run,gamma)][:-1][1]-7.5)/13.3, label="sec. avg. dist.", color="black")
    # ax.fill_between([0,1],[runs_hist_df.loc[20,(run,gamma)][:-1][0],runs_hist_df.loc[20,(run,gamma)][:-1][1]], color="silver")
    # plt.axvline(mean_sec_dist,ymax = 0.1, label="sec. avg. dist.", color="black")
    # sns.set_context(rc={"lines.linewidth": 5, 'lines.markersize': 8})
    # sns.lineplot(np.arange(0,len(avg_hist_df[gamma][time][:-1])), avg_hist_df[gamma][time][:-1], label="strong neg.", color=palette[i-1], marker="o", legend=None)
# popt = curve_fit(func, np.arange(0,len(avg_hist_df[0.1][time])), avg_hist_df[0.1][time], maxfev=10000)[0]
# sns.lineplot(np.arange(0,len(avg_hist_df[0.1][time])), func(np.arange(0,len(avg_hist_df[0.1][time])), *popt), color=palette[3], linewidth=2)
for i, gamma in enumerate([0.1]):
    for run in [1,6,3]:  # range(10):
        sns.lineplot(np.arange(0, len(runs_hist_df.loc[20, (run, gamma)][:-1])),
                     runs_hist_df.loc[20, (run, gamma)][:-1], label="strong neg.", color=palette[-1], marker="o",
                     legend=None)
for i, gamma in enumerate([0.1]):
    for run in [0]:  # range(10):
        sns.lineplot(np.arange(0, len(runs_hist_df.loc[0, (run, gamma)][:-1])),
                     runs_hist_df.loc[0, (run, gamma)][:-1], label="strong neg.", color="black", marker="o",
                     legend=None)

palette = sns.color_palette("bwr", 4)
# ax.set(ylim=(7.5,28))

# include exponential fit with every curve
#
# sns.lineplot(np.arange(0,len(avg_hist_df[0.5][time])), avg_hist_df[0.5][time], label="neg.", color=palette[2], marker="o")
# popt = curve_fit(func, np.arange(0,len(avg_hist_df[0.5][time])), avg_hist_df[0.5][time], maxfev=10000)[0]
# sns.lineplot(np.arange(0,len(avg_hist_df[0.5][time])), func(np.arange(0,len(avg_hist_df[0.5][time])), *popt), color=palette[2], linewidth=2)
#
#
# sns.set_context("talk", font_scale=1.7, rc={"lines.linewidth": 3, 'lines.markersize': 5})
# sns.lineplot(np.arange(0,len(avg_hist_df[0.5][0])), avg_hist_df[0.5][0], label="no", color="black", marker="o")
# sns.set_context("talk", font_scale=1.7, rc={"lines.linewidth": 5, 'lines.markersize': 10})
#
#
# sns.lineplot(np.arange(0,len(avg_hist_df[2][time])), avg_hist_df[2][time], label="pos.", color=palette[1], marker="o")
# popt = curve_fit(func, np.arange(0,len(avg_hist_df[2][time])), avg_hist_df[2][time], maxfev=10000)[0]
# sns.lineplot(np.arange(0,len(avg_hist_df[2][time])), func(np.arange(0,len(avg_hist_df[2][time])), *popt), color=palette[1], linewidth=2)
#
# sns.lineplot(np.arange(0,len(avg_hist_df[10][time])), avg_hist_df[10][time], label="strong pos.", color=palette[0], marker="o")
# popt = curve_fit(func, np.arange(0,len(avg_hist_df[10][time])), avg_hist_df[10][time], maxfev=10000)[0]
# sns.lineplot(np.arange(0,len(avg_hist_df[10][time])), func(np.arange(0,len(avg_hist_df[10][time])), *popt), color=palette[0], linewidth=2)
#

# customise plot

plt.xlabel("cell-cell-distance")
plt.ylabel("Avg. c. (pM)")
plt.xticks([0,2,4,6,8,10])

# create custom legend
import matplotlib.lines as mlines
# zero-three refers to palette entries
zero_line = mlines.Line2D([], [], color=palette[0], label='lightcoral line')
one_line = mlines.Line2D([], [], color=palette[1], label='red line')
black_line = mlines.Line2D([], [], color='black', label='Black line')
two_line = mlines.Line2D([], [], color=palette[2], label='lightsteelblue line')
three_line = mlines.Line2D([], [], color=palette[3], label='blue line')

grey_line = mlines.Line2D([], [], color="grey", label='grey line')
grey_dashed_line = mlines.Line2D([], [], color="grey", label='grey dashed line', linestyle="dashed")
white_line = mlines.Line2D([], [], color="white", label='white line')


# new_handles = [three_line, black_line, zero_line]
# new_labels = ["neg. feedback", "no feedback", "pos. feedback"]
new_handles = [three_line]
new_labels = ["neg. feedback"]
# new_handles = [three_line, two_line, black_line, one_line, zero_line]
# new_labels = ["strong negative feedback", "negative feedback", "no feedback", "positive feedback", "strong positive feedback"]
# plt.legend(handles=new_handles, labels=new_labels,loc='upper right', fancybox=True)
fig.savefig("plots/" + "histogram_triple" + ".svg", bbox_inches='tight')
plt.show()