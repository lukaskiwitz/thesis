import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import seaborn as sns
import random
from my_debug import message
from scipy import spatial

load_runs_hist_df = True
save_runs_hist_df = True

save_avg_hist_df = True

verbose = False

no_of_runs = 10
timepoints = [0,2,5,10,20]
gammas = [1/9]
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
if load_runs_hist_df == True:
    try:
        runs_hist_df = pd.read_hdf(base_path + str(no_of_runs) + "_runs_hist_df.h5", mode="r")
        print("CHANGED RUNS HIST DF TO NEG")
        runs_hist_df = pd.read_hdf(
            "/extra/brunner/thesis/kinetic/q_fraction_medium_neg_g_scan_multi/" + "6_runs_hist_df.h5", mode="r")
        message("loaded runs_hist_df from " + base_path)
        df_loaded = True
    except FileNotFoundError:
        message("runs_hist_df not found, commencing calculation")
        df_loaded = False

if load_runs_hist_df == False or df_loaded == False:
    for run in range(no_of_runs):
        message("calculating run " + str(run))
        for gamma in gammas:
            if verbose == True:
                message("calculating gamma = " + str(gamma))
            get_dataframes.append([])

            path = base_path + "run" + str(run) + "/"

            global_df = pd.read_hdf(path + 'global_df.h5', mode="r")
            cell_df = pd.read_hdf(path + 'cell_df.h5', mode="r")


            get_dataframes[run] = [global_df, cell_df]

            variable = "time"
            cell_df = cell_df.loc[cell_df["IL-2_fraction"] == 0.25]
            cell_df = cell_df.loc[cell_df["IL-2_gamma"] == gamma]

            #for which timepoints to calculate for
            points = [cell_df[variable].unique()[i] for i in timepoints]

            #create multiindex dfs
            hist_df = pd.DataFrame(columns=[str(x) for x in points], index=cell_df.loc[(cell_df[variable] == points[0]) , "id"].values)

            # prepare df's to build kd_tree
            cell_df_split = cell_df.loc[cell_df[variable] == points[0]].copy()
            cell_points = np.array([cell_df_split["x"].values, cell_df_split["y"].values, cell_df_split["z"].values]).T
            kd_tree = spatial.KDTree(cell_points)

            max_distance = 0

            #get the distances from the KDTree at the given time points
            for entry in points: #cell_df["time"].unique():
                if verbose == True:
                    message("calculating KDTree at " + variable + " = " + str(entry))
                # create containers
                cell_df_split = cell_df.loc[cell_df[variable] == entry].copy()
                # iterate over cells
                for il2cells in cell_df.loc[(cell_df[variable] == points[0]), "id"]:

                    kd_tree_data = kd_tree.query(cell_points[int(il2cells)], len(cell_points))
                    hist_df[str(entry)][str(il2cells)] = [np.zeros(len(cell_points)), np.zeros(len(cell_points))]

                    #write concentrations
                    hist_df[str(entry)][str(il2cells)][0][kd_tree_data[1]] = cell_df_split.iloc[kd_tree_data[1]]["IL-2_surf_c"].mul(1e3).values
                    #write distances
                    hist_df[str(entry)][str(il2cells)][1][kd_tree_data[1]] = kd_tree_data[0]
                    #get the max distance between two cells while we are at it
                    if max_distance < kd_tree_data[0].max():
                        max_distance = kd_tree_data[0].max()

            #############################################################################################################################################

            # Histogram building

            # build various containers
            hist_timepoints = [str(x) for x in points]
            hist_data_array = []
            c_data_array = []
            conv_data_array = []
            auto_corr_data_array = []

            sec_cell_corr_array = []

            niche = []

            #iterate over timepoints:
            for hist_timepoint in hist_timepoints:
                if verbose == True:
                    message("creating hist at time = " + hist_timepoint)

                # calculate necessary bins
                bins = np.arange(0,(max_distance//cell_cell_distance + 1)*cell_cell_distance, cell_cell_distance)
                center = bins[:-1]/cell_cell_distance

                # create container for correlation computation
                c_bin_array = [[] for i in range((len(bins)-1))]
                multi_cell_corr_array = [[] for i in range(len(hist_df[hist_timepoint][cell_df.loc[(cell_df["type_name"] == "changed") & (cell_df["time"] == float(hist_timepoint)), "id"]]))]


                # iterate over all secreting cells in hist_df to average over them
                # grab the secreting cells average distance from each other while we are at it
                sec_dists = []
                cell_tracker = 0
                no_of_cells = len(hist_df[hist_timepoint][cell_df.loc[(cell_df["type_name"] == "changed") & (cell_df["time"] == float(hist_timepoint)), "id"]])
                for k,cell in enumerate(hist_df[hist_timepoint][cell_df.loc[(cell_df["type_name"] == "changed") & (cell_df["time"] == float(hist_timepoint)), "id"]]):  #.loc[hist_df[plot_time_point].index == str(2)].dropna():
                    #sort the arrays and corresponding concentrations by distance via zipping und unzipping
                    cell_tracker += 1

                    zipped = zip(cell[0], cell[1])
                    sorted_zipped = sorted(list(zipped), key=lambda x: x[1])

                    c, distances = zip(*sorted_zipped)
                    c = list(c)
                    distances = list(distances)

                    # since "c" and "distances" are sorted, we can use np.where's first and last argument as a range for the entries of the bin
                    # this way we only iterate over the bins, the iteration over the distances is done by numpy.where
                    for i,bin in enumerate(bins[1:]):
                        if i == 0:
                            current_range = np.where(distances < bins[1:][i])[0]
                        else:
                            current_range = np.where((distances<bins[1:][i]) & (distances>bins[1:][i-1]))[0]

                        if len(current_range) == 0:
                            pass
                        elif len(current_range) == 1:
                            c_bin_array[i].extend([c[current_range[0]]])
                        else:
                            c_bin_array[i].extend(c[current_range[0]:current_range[-1]])
                            # print(0 + "test")

                    #get distances between secs
                    ids = [int(x) for x in cell_df.loc[(cell_df["type_name"] == "changed") & (cell_df["time_index"] == 0), "id"].values]
                    for i in cell[1][ids]:
                        if i != 0.0:
                            sec_dists.append(i)

                # calc mean and std of bins
                hist = [np.mean(k) for k in c_bin_array]
                # write histogram into df
                runs_hist_df.loc[:, (run, gamma)][int(float(hist_timepoint))] = hist
                runs_hist_df.loc[:, (run, gamma)]["sec_dists"] = sec_dists

                #############################################################################################################################################

                # Autocorrelation with space instead of time:

                # auto_corr = []
                # for r_dash in range(len(hist)-1):
                #     R = 0
                #     for r in range(len(hist)):
                #         try:
                #             R += hist[r + r_dash]*hist[r]
                #         except IndexError:
                #             pass
                #     auto_corr.append(R)
                # auto_corr_data_array.append(auto_corr)

    # write runs_hist_df to file
    if save_runs_hist_df == True:
        runs_hist_df.to_hdf(base_path + str(no_of_runs) + '_runs_hist_df.h5', key="data", mode="w")
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
    temp_sec_dists.append(transposed_runs_hist_df["sec_dists"][i][1/9])
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

from scipy.optimize import curve_fit

def func(x, a,b,c):
    return a*np.exp(-x/b) + c

niche_at_0_1 = [3, 3, 4, 4, 3, 3,3]
niche_effect_at_0_1 = np.array([1.248354058933971,
1.374301007199112,
  1.4110233480906524,
  1.3712516856247705,
  1.253301502010431,
  1.3139483268186123,
  1.2679607881211836,
  1.377898802550747,
  1.2393480407747517,
  1.395052223882936])

# palette = palette
for i,gamma in enumerate([1/9]):
    sns.set_context(rc={"lines.linewidth": 1, 'lines.markersize': 3})
    for run in [3]: #range(10):
        sns.set_context(rc={"lines.linewidth": 5, 'lines.markersize': 8})
        sns.lineplot(np.arange(0,len(runs_hist_df.loc[20,(run,gamma)][:-1])), runs_hist_df.loc[20,(run,gamma)][:-1], label="strong neg.", color=palette[i-1], marker="o", legend=None)
        sns.set_context(rc={"lines.linewidth": 1, 'lines.markersize': 7})
        # plt.scatter(x=niche_at_0_1[run], y=runs_hist_df.loc[20,(run,gamma)][:-1][niche_at_0_1[run]], color="black")

        # plt.axhline(xmin=0, xmax=niche_at_0_1[run]/10, y=niche_effect_at_0_1[run]*np.mean(runs_hist_df.loc[20, (run, gamma)][:-1]), color="black")
    # plt.axvline(mean_sec_dist, label="sec. avg. dist.", color="black")
    sns.set_context(rc={"lines.linewidth": 5, 'lines.markersize': 8})
    # sns.lineplot(np.arange(0,len(avg_hist_df[gamma][time][:-1])), avg_hist_df[gamma][time][:-1], label="strong neg.", color=palette[i-1], marker="o", legend=None)
# popt = curve_fit(func, np.arange(0,len(avg_hist_df[0.1][time])), avg_hist_df[0.1][time], maxfev=10000)[0]
# sns.lineplot(np.arange(0,len(avg_hist_df[0.1][time])), func(np.arange(0,len(avg_hist_df[0.1][time])), *popt), color=palette[3], linewidth=2)
palette = sns.color_palette("bwr", 4)
# ax.set(ylim=(62,132.2))

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


plt.xlabel("cell-cell-distance")
plt.ylabel("Avg. c. (pM)")
plt.xticks([0,2,4,6,8,10])

# ax2 = sns.lineplot(center[:len(hist_data_array[1])], hist_data_array[1], label=str(float(time_list[1])*10)+"h", legend=False)
# ax3 = sns.lineplot(center[:len(hist_data_array[2])], hist_data_array[2], label=str(float(time_list[2])*10)+"h", legend=False)
import matplotlib.lines as mlines

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
fig.savefig("plots/" + "histogram_0.1_aufgefaechert" + ".svg", bbox_inches='tight')
plt.show()
# plt.plot(niche, label="c_niche")
# plt.plot(AC_niche, label="AC_niche")
# plt.legend()
# plt.show()

# # fig.add_subplot(a_x, a_y, 2)
# # for i in range(len(c_data_array)):
# #     ax4 = sns.lineplot(center[:len(c_data_array[i])], c_data_array[i], label=variable + "=" + str(float(points[i])), legend=False)
# # # ax5 = sns.lineplot(center[:len(c_data_array[1])], c_data_array[1], label=str(float(time_list[1])*10)+"h", legend=False)
# # # ax6 = sns.lineplot(center[:len(c_data_array[2])], c_data_array[2], label=str(float(time_list[2])*10)+"h", legend=False)
# # plt.axvline(mean_sec_dist, label="sec. avg. dist.", color="black")
# # ax4.set(xlabel="cell-cell-distance", yscale=yscale, xscale=xscale, title="c_0*c_r", xticks=[0,5,10])
#
# # fig.add_subplot(a_x,a_y,3)
# # for i in range(len(conv_data_array)):
# #     ax7 = sns.lineplot(center[:len(conv_data_array[i])], conv_data_array[i], label=variable + "=" + str(float(points[i])), legend=False)
# # # ax8 = sns.lineplot(center[:len(conv_data_array[1])], conv_data_array[1], label=str(float(time_list[1])*10)+"h", legend=False)
# # # ax9 = sns.lineplot(center[:len(conv_data_array[2])], conv_data_array[2], label=str(float(time_list[2])*10)+"h", legend=False)
# # plt.axvline(mean_sec_dist, label="sec. avg. dist.", color="black")
# # ax7.set(xlabel="cell-cell-distance", ylabel="Convolution", yscale=yscale, xscale=xscale, title="conv(c_r * c_r)", xticks=[0,5,10])
# #
# fig.add_subplot(a_x,a_y,2)
# for i in range(len(auto_corr_data_array)):
#     ax10 = sns.lineplot(center[:len(auto_corr_data_array[i])], auto_corr_data_array[i], label=variable + "=" + str(float(points[i])*10) + "h", legend=False)
# # ax11 = sns.lineplot(center[:len(auto_corr_data_array[1])], auto_corr_data_array[1], label=str(float(time_list[1])*10)+"h", legend=False)
# # ax12 = sns.lineplot(center[:len(auto_corr_data_array[-1])], auto_corr_data_array[-1], label=str(float(time_list[-1])*10)+"h", legend=False)
#
# ax10.set(xlabel="cell-cell-distance", ylabel="G(r)", yscale=yscale, xscale=xscale, title="AC with r")
# # plt.axvline(mean_sec_dist, label="sec. avg. dist.", color="black")
# #
# # fig.add_subplot(a_x, a_y, 3)
# # for i in range(len(sec_cell_corr_array)):
# #     ax14 = sns.lineplot(center[:len(sec_cell_corr_array[i])], sec_cell_corr_array[i], label=variable + "=" + str(float(points[i])), legend=False)
# # # ax14.set(xlabel="cell-cell-distance", ylabel="G(r)", yscale=yscale, xscale=xscale, title="iterating AC", xticks=[0,5,10])
# #
# # plt.legend(bbox_to_anchor=(1.04,0.5), loc="center left", borderaxespad=0)
# #
# #
# fig.add_subplot(a_x,a_y,3)
# ax13 = sns.lineplot(np.array(points)*10, niche, label="hist_niche", legend=False)
# ax13 = sns.lineplot(np.array(points[:])*10, AC_niche[:], label="AC_niche")
# ax13.set(xlabel="time in h", ylabel="cell-cell-distance", yscale=yscale, xscale=xscale, title="Niche size at threshold " + str(threshold))
#
#
# # fig.savefig("/home/brunner/Documents/Current work/28082020/" + "niche_" + "f=" + str(cell_df["IL-2_fraction"].unique()[0]) + "_t=" + str(threshold) + ".pdf", bbox_inches='tight')
# # fig.savefig("/home/brunner/Documents/Current work/28082020/" + "niche_" + "f=" + str(cell_df["IL-2_fraction"].unique()[0]) + "_t=" + str(threshold) + ".png", bbox_inches='tight')
#
# plt.show()
