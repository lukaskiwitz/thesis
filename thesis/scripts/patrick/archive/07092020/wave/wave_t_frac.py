import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import seaborn as sns
import random
from my_debug import message
from scipy import spatial

#load the cell_df
# from parameters_q_fraction import path_kinetic
# path = path_kinetic
# path = "/extra/brunner/thesis/kinetic/q_fraction/"
# path = "/extra/brunner/thesis/kinetic/q_fraction_wave_gamma_scan/"
# path = "/extra/brunner/thesis/kinetic/q_fraction_wave_multi/q_fraction_wave_multi_0/"
# path = "/extra/brunner/thesis/kinetic/q_fraction_large_gamma_scan_new_paras/"
# path = "/extra/brunner/thesis/kinetic/q_fraction_test/"
# path = "/extra/brunner/thesis/kinetic/q_fraction_large_gamma_scan_g_0.01/"
path = "/extra/brunner/thesis/kinetic/q_fraction_medium_g_scan_multi/run0/"

cell_df = pd.read_hdf(path + 'cell_df.h5', mode="r")
global_df = pd.read_hdf(path + 'global_df.h5', mode="r")

print("cell types: ", cell_df["type_name"].unique())
print("gammas:", global_df["IL-2_gamma"].unique())
print("time:" , cell_df["time"].unique())
try:
    print("fractions: ", global_df["IL-2_fraction"].unique())
except:
    print("no fraction")

# print(global_df["Concentration"])
# cell_df = cell_df.loc[cell_df["scan_index"] <15]

# variable = "IL-2_fraction"
variable = "time"
cell_df = cell_df.loc[cell_df["IL-2_fraction"] == 0.25]
cell_df = cell_df.loc[cell_df["IL-2_gamma"] == 2]
# cell_df = cell_df.loc[cell_df["x"] == 20]
# cell_df = cell_df.reset_index()
# for time in cell_df["time_index"].unique():
#     for i in range(len(cell_df.loc[cell_df["time_index"] == 0])):
#         cell_df.loc[cell_df["time_index"] == time, "id"] = pd.Series(range(121))


#for which timepoints to calculate for
points = [cell_df[variable].unique()[i] for i in [0,5,10,15]]#[0., 10., 19.]


#create multiindex dfs
hist_df = pd.DataFrame(columns=[str(x) for x in points], index=cell_df.loc[(cell_df[variable] == points[0]) , "id"].values)
# hist_df = pd.DataFrame(columns=[str(x) for x in time_points], index=cell_df.loc[(cell_df["type_name"] == "changed") & (cell_df["time"] == 0) , "id"].values)

# normalising_df = pd.DataFrame(columns=[str(x) for x in time_points], index=cell_df.loc[(cell_df["type_name"] == "changed") & (cell_df["time"] == 0), "id"].values)
# exit()
# small_cell_sample = random.sample(cell_df.loc[(cell_df["type_name"] == "changed") & (cell_df["time"] == time_points[0]), "id"].values.tolist(), 1)
# for desired time_points
cell_df_split = cell_df.loc[cell_df[variable] == points[0]].copy()
cell_points = np.array([cell_df_split["x"].values, cell_df_split["y"].values, cell_df_split["z"].values]).T
kd_tree = spatial.KDTree(cell_points)
max_distance = 0

for entry in points: #cell_df["time"].unique():
    message(variable + "= " + str(entry))
    # create containers
    cell_df_split = cell_df.loc[cell_df[variable] == entry].copy()

    # for il2cells in cell_df.loc[(cell_df["type_name"] == "changed") & (cell_df["time_index"] == 0), "id"]: #small_cell_sample: #
    for il2cells in cell_df.loc[(cell_df[variable] == points[0]), "id"]:

        # kd_tree_data = kd_tree.query(cell_points[int(il2cells)], len(cell_points))
        kd_tree_data = kd_tree.query(cell_points[int(il2cells)], len(cell_points))
        hist_df[str(entry)][str(il2cells)] = [np.zeros(len(cell_points)), np.zeros(len(cell_points))]

        #write concentrations
        hist_df[str(entry)][str(il2cells)][0][kd_tree_data[1]] = cell_df_split.iloc[kd_tree_data[1]]["IL-2_surf_c"].mul(1e3).values
        #write distances
        hist_df[str(entry)][str(il2cells)][1][kd_tree_data[1]] = kd_tree_data[0]
        if max_distance < kd_tree_data[0].max():
            max_distance = kd_tree_data[0].max()



# sns.set(rc={'figure.figsize':(11,8.27)})
# fig = plt.figure()
# plt.subplots_adjust(wspace=.3)
# a_x = 1
# a_y = 2

plot_list = [str(x) for x in points]
hist_data_array = []
c_data_array = []
conv_data_array = []
auto_corr_data_array = []

sec_cell_corr_array = []

niche = []

#iterate over time points:

# import time
for plot_point in plot_list:
    # calculate necessary bins
    message(plot_point)
    cell_cell_distance = 20
    # max_distance = np.array([x for x in hist_df[plot_time_point]]).flatten().max()
    # max_distance = np.array([x for x in normalising_df[plot_time_point][10]]).max()
    bins = np.arange(0,(max_distance//cell_cell_distance + 1)*cell_cell_distance, cell_cell_distance)
    center = bins[:-1]/cell_cell_distance

    # create container for histogram computation
    c_bin_array = [[] for i in range((len(bins)-1))]
    multi_cell_corr_array = [[] for i in range(len(hist_df[plot_point][cell_df.loc[(cell_df["type_name"] == "changed") & (cell_df["time"] == float(plot_point)), "id"]]))]
    # multi_cell_corr_array = np.empty((len(hist_df[plot_point][cell_df.loc[(cell_df["type_name"] == "changed") & (cell_df["time_index"] == 0), "id"]]), 100))
    # multi_cell_corr_array[:] = np.NaN
    # multi_cell_corr_array = {}

    # iterate over all secreting cells in hist_df to average over them
    # grab the secreting cells average distance from each other while we are at it
    sec_dists = []
    # hist_df[plot_time_point] = hist_df[plot_time_point].loc[hist_df[plot_time_point].index == str(56)]
    cell_tracker = 0
    no_of_cells = len(hist_df[plot_point][cell_df.loc[(cell_df["type_name"] == "changed") & (cell_df["time"] == float(plot_point)), "id"]])
    for k,cell in enumerate(hist_df[plot_point][cell_df.loc[(cell_df["type_name"] == "changed") & (cell_df["time"] == float(plot_point)), "id"]]):  #.loc[hist_df[plot_time_point].index == str(2)].dropna():
        #sort the arrays and corresponding concentrations by distance via zipping und unzipping

        # print(str(cell_tracker) + "/" + str(no_of_cells))
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

    #     single_cell_correlation = []
    #     for j in range(len(c_bin_array)):
    #         single_cell_correlation.append(np.mean([c[0]*i for i in c_bin_array[j]])/len(c_bin_array[j]))
    #     multi_cell_corr_array[k] = single_cell_correlation
    #
    #
    # if len(multi_cell_corr_array) == 1:
    #     sec_cell_corr_array.append([multi_cell_corr_array[0]])
    # else:
    #     sec_cell_corr_array.append([np.nanmean(k) for k in np.array(multi_cell_corr_array).T])


    # calc mean and std of bins
    hist = [np.mean(k) for k in c_bin_array]   #/global_df.loc[(global_df["time"] == float(plot_point)) & (global_df["IL-2_fraction"] == cell_df["IL-2_fraction"].values[0]), "surf_c"].values[0]
    # plt.plot(hist)
    # plt.show()

    hist_data_array.append(hist)
    c_data_array.append([hist[0]*hist[i] for i in range(len(hist))])
    # error = [np.std(k, dtype=np.float64) for k in c_bin_array]
    # lolims = [0 for k in error]
    g0_hist = [k for k in hist if np.isnan(k) == False]

    #############################################################################################################################################

    # # for convolution purposes get the hist for every secreting cell, convolute with the sum (g0_hist) and average over it?
    # # DISCUSS if this is the way to go
    # sec_cells_hist_array = [[] for i in range(len(ids))]
    # conv_array = [[] for i in range(len(ids))]
    # # for cell in hist_df[plot_time_point][cell_df.loc[(cell_df["type_name"] == "default") & (cell_df["time_index"] == 0), "id"]]:
    # for j, cell in enumerate(hist_df[plot_point][ids]): #.loc[hist_df[plot_time_point].index == str(11)].dropna():
    #     #sort the arrays and corresponding concentrations by distance via zipping und unzipping
    #     zipped = zip(cell[0], cell[1])
    #     sorted_zipped = sorted(list(zipped), key=lambda x: x[1])
    #
    #     c, distances = zip(*sorted_zipped)
    #     c = list(c)
    #     distances = list(distances)
    #
    #     # since "c" and "distances" are sorted, we can use np.where's first and last argument as a range for the entries of the bin
    #     # this way we only iterate over the bins, the iteration over the distances is done by numpy
    #     for i,bin in enumerate(bins[1:]):
    #         if i == 0:
    #             current_range = np.where(distances < bins[1:][i])[0]
    #         else:
    #             current_range = np.where((distances<bins[1:][i]) & (distances>bins[1:][i-1]))[0]
    #
    #         if len(current_range) == 0:
    #             pass
    #         elif len(current_range) == 1:
    #             c_bin_array[i].extend([c[current_range[0]]])
    #         else:
    #             c_bin_array[i].extend(c[current_range[0]:current_range[-1]])
    #     sec_cells_hist_array[j] = [np.mean(k) for k in c_bin_array]
    #     g1_hist = [np.mean(k) for k in c_bin_array if len(k) != 0]
    #     conv_array[j] = np.convolve(g0_hist, g0_hist, "same")
    #
    # # g1_hist = [np.mean(k) for k in np.transpose(sec_cells_hist_array)]
    # # g1_hist = [np.mean(k) for k in c_bin_array]
    # # plt.semilogy(np.convolve(g1_hist[:-1], g0_hist[:-1], "same"), label="Convolution")
    # # plt.plot(np.convolve(g1_hist[0], np.array(g0_hist)), "-o", label=str(float(plot_time_point)*10)+"h")
    # conv_data_array.append([np.mean(k) for k in np.transpose(conv_array)])
    # # fig.add_subplot(a_x, a_y, 2)
    # ax2 = plt.plot([np.mean(k) for k in np.transpose(conv_array)], "-o", label=str(float(plot_point) * 10) + "h")
    # plt.show()
    #############################################################################################################################################

    # Autocorrelation with space instead of time:

    auto_corr = []
    for r_dash in range(len(g0_hist)-1):
        R = 0
        for r in range(len(g0_hist)):
            try:
                R += g0_hist[r + r_dash]*g0_hist[r]
            except IndexError:
                pass
        auto_corr.append(R)
    auto_corr_data_array.append(auto_corr)

#############################################################################################################################################

# Niche
# Fitting

AC_niche = []

def func(x, a, b, c):
    return a * b * np.exp(-x/b) + c
from scipy.optimize import curve_fit
for i in range(len(auto_corr_data_array)):
    popt, pcov = curve_fit(func, np.arange(len(auto_corr_data_array[i])), auto_corr_data_array[i], maxfev=10000)
    AC_niche.append(popt[1])


# Thresholding

niche = []
niche_avg_c = []

for j, hist in enumerate(hist_data_array):
    if j == 0:
        threshold = hist[1] #200
    for i,c in enumerate(hist):
        if c < threshold:
            niche.append(i-1)
            avg_c = []
            for k in range(i):
                avg_c.append(hist[k])
            niche_avg_c.append(np.mean(avg_c))
            break
        if i == len(hist)-1:
            niche.append(i-1)
            avg_c = []
            for k in range(i):
                avg_c.append(hist[k])
            niche_avg_c.append(np.mean(avg_c))

#take the last timepoints user defined niche and get the average concentration within that niche
user_defined_niche = 3
# user_defined_niche_avg_c = []
hist = hist_data_array[-1]
avg_c = []
for k in range(i):
    avg_c.append(hist[k])
user_defined_niche_avg_c = np.mean(avg_c)

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

fig = plt.figure(figsize=(8,8))
# plt.figure(figsize=(12,10))
# sns.set_context("talk", font_scale=1.1, rc={"lines.linewidth": 2.5})
# sns.set(rc={'figure.figsize':(8,8)})
sns.set_style("ticks")
sns.set_context("talk", font_scale=1.7, rc={"lines.linewidth": 5, 'lines.markersize': 10})
plt.subplots_adjust(wspace=.4, hspace=0.4)
a_x = 2
a_y = 2

xscale = "linear"
yscale = "linear"

mean_sec_dist = np.mean(sec_dists)/cell_cell_distance
# if np.isnan(mean_sec_dist) == True:
#     mean_sec_dist = 0

variable = "t"

# fig.add_subplot(a_x, a_y, 1)
palette = sns.color_palette("binary", len(hist_data_array))
for i in range(len(hist_data_array)):
    ax1 = sns.lineplot(center[:len(hist_data_array[i])], hist_data_array[i], label=variable + "=" + str(float(points[i])*10) + "h", legend="brief", color=palette[i], marker="o")
# ax2 = sns.lineplot(center[:len(hist_data_array[1])], hist_data_array[1], label=str(float(time_list[1])*10)+"h", legend=False)
# ax3 = sns.lineplot(center[:len(hist_data_array[2])], hist_data_array[2], label=str(float(time_list[2])*10)+"h", legend=False)
plt.axvline(mean_sec_dist, color="black", linestyle="dashed")
ax1.set(xlabel="cell-cell-distance", ylabel="Avg. c.", yscale=yscale, xscale=xscale, xticks=[0,5,10])
plt.xticks([0,2,4,6,8,10,12])
plt.show()
# fig.savefig("plots/" + "histogram" + ".png", bbox_inches='tight')



# for i in range(len(auto_corr_data_array)):
#     ax10 = sns.lineplot(center[:len(auto_corr_data_array[i])], auto_corr_data_array[i], label=variable + "=" + str(float(points[i])*10) + "h", legend=False)
# # ax11 = sns.lineplot(center[:len(auto_corr_data_array[1])], auto_corr_data_array[1], label=str(float(time_list[1])*10)+"h", legend=False)
# # ax12 = sns.lineplot(center[:len(auto_corr_data_array[-1])], auto_corr_data_array[-1], label=str(float(time_list[-1])*10)+"h", legend=False)
#
# ax10.set(xlabel="cell-cell-distance", ylabel="G(r)", yscale=yscale, xscale=xscale, title="AC with r")
# # plt.axvline(mean_sec_dist, label="sec. avg. dist.", color="black")





# fig.add_subplot(a_x,a_y,3)
# ax13 = sns.lineplot(np.array(points)*10, niche, label="hist_niche", legend=False)
# ax13 = sns.lineplot(np.array(points[:])*10, AC_niche[:], label="AC_niche")
# ax13.set(xlabel="time in h", ylabel="cell-cell-distance", yscale=yscale, xscale=xscale, title="Niche size at threshold " + str(round(threshold,1)))
#
# fig.add_subplot(a_x,a_y,4)
# plt.plot(np.array(points)*10, niche_avg_c, label="niche avg. c")
# # fig.savefig("/home/brunner/Documents/Current work/28082020/" + "niche_" + "f=" + str(cell_df["IL-2_fraction"].unique()[0]) + "_t=" + str(threshold) + ".pdf", bbox_inches='tight')
# # fig.savefig("/home/brunner/Documents/Current work/28082020/" + "niche_" + "f=" + str(cell_df["IL-2_fraction"].unique()[0]) + "_t=" + str(threshold) + ".png", bbox_inches='tight')
#
# plt.show()
