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
path = "/extra/brunner/thesis/kinetic/q_fraction_test/"

cell_df = pd.read_hdf(path + 'cell_df.h5', mode="r")
global_df = pd.read_hdf(path + 'global_df.h5', mode="r")

cell_df = cell_df.loc[cell_df["IL-2_fraction"] == 0.25]
# cell_df = cell_df.loc[cell_df["IL-2_gamma"] == 0.5]
# cell_df = cell_df.loc[cell_df["z"] == 60]
# cell_df = cell_df.reset_index()
# for time in cell_df["time_index"].unique():
#     for i in range(len(cell_df.loc[cell_df["time_index"] == 0])):
#         cell_df.loc[cell_df["time_index"] == time, "id"] = pd.Series(range(121))


#for which timepoints to calculate for
time_points = [0., 1., 2.] #cell_df["time"].unique()#[0., 10., 19.]


#create multiindex dfs
hist_df = pd.DataFrame(columns=[str(x) for x in time_points], index=cell_df.loc[(cell_df["time"] == 0) , "id"].values)
# hist_df = pd.DataFrame(columns=[str(x) for x in time_points], index=cell_df.loc[(cell_df["type_name"] == "changed") & (cell_df["time"] == 0) , "id"].values)

# normalising_df = pd.DataFrame(columns=[str(x) for x in time_points], index=cell_df.loc[(cell_df["type_name"] == "changed") & (cell_df["time"] == 0), "id"].values)
# exit()
# small_cell_sample = random.sample(cell_df.loc[(cell_df["type_name"] == "changed") & (cell_df["time"] == time_points[0]), "id"].values.tolist(), 1)
# for desired time_points
cell_df_timed = cell_df.loc[cell_df["time_index"] == 0].copy()
cell_points = np.array([cell_df_timed["x"].values, cell_df_timed["y"].values, cell_df_timed["z"].values]).T
kd_tree = spatial.KDTree(cell_points)
max_distance = 0

for time in time_points: #cell_df["time"].unique():
    message("time = " + str(time))
    # create containers
    cell_df_timed = cell_df.loc[cell_df["time"] == time].copy()

    # for il2cells in cell_df.loc[(cell_df["type_name"] == "changed") & (cell_df["time_index"] == 0), "id"]: #small_cell_sample: #
    for il2cells in cell_df.loc[(cell_df["time_index"] == 0), "id"]:

        kd_tree_data = kd_tree.query(cell_points[int(il2cells)], len(cell_points))
        hist_df[str(time)][str(il2cells)] = [np.zeros(len(cell_points)), np.zeros(len(cell_points))]

        #write concentrations
        hist_df[str(time)][str(il2cells)][0][kd_tree_data[1]] = cell_df_timed.iloc[kd_tree_data[1]]["IL-2_surf_c"].mul(1e3).values
        #write distances
        hist_df[str(time)][str(il2cells)][1][kd_tree_data[1]] = kd_tree_data[0]
        if max_distance < kd_tree_data[0].max():
            max_distance = kd_tree_data[0].max()


sns.set(palette="Greens")
sns.set_style("ticks")
sns.set_context("talk", font_scale=1, rc={"lines.linewidth": 3})
# sns.set(rc={'figure.figsize':(11,8.27)})
# fig = plt.figure()
# plt.subplots_adjust(wspace=.3)
# a_x = 1
# a_y = 2

time_list = [str(x) for x in time_points]
hist_data_array = []
c_data_array = []
conv_data_array = []
auto_corr_data_array = []

sec_cell_corr_array = []

niche = []

#iterate over time points:
for plot_time_point in time_list:
    # calculate necessary bins
    message(plot_time_point)
    cell_cell_distance = 20
    # max_distance = np.array([x for x in hist_df[plot_time_point]]).flatten().max()
    # max_distance = np.array([x for x in normalising_df[plot_time_point][10]]).max()
    bins = np.arange(0,(max_distance//cell_cell_distance + 1)*cell_cell_distance, cell_cell_distance)
    center = bins[:-1]/cell_cell_distance

    # create container for histogram computation
    c_bin_array = [[] for i in range((len(bins)-1))]

    # iterate over all secreting cells in hist_df to average over them
    # grab the secreting cells average distance from each other while we are at it
    sec_dists = []
    # hist_df[plot_time_point] = hist_df[plot_time_point].loc[hist_df[plot_time_point].index == str(56)]
    for cell in hist_df[plot_time_point][cell_df.loc[(cell_df["type_name"] == "changed") & (cell_df["time_index"] == 0), "id"]]:  #.loc[hist_df[plot_time_point].index == str(2)].dropna():
        #sort the arrays and corresponding concentrations by distance via zipping und unzipping
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

        #get distances between secs
        ids = [int(x) for x in cell_df.loc[(cell_df["type_name"] == "changed") & (cell_df["time_index"] == 0), "id"].values]
        for i in cell[1][ids]:
            if i != 0.0:
                sec_dists.append(i)
        single_cell_correlation = []
        for j in range(len(c_bin_array)):
            single_cell_correlation.append(np.mean([c[0]*i for i in c_bin_array[j]])/len(c_bin_array[j]))
        sec_cell_corr_array.append(single_cell_correlation)

    if len(sec_cell_corr_array) == 1:
        test = sec_cell_corr_array[0]
    else:
        test = [np.nanmean(k) for k in np.array(sec_cell_corr_array).T]
    plt.plot(test)
    plt.show()

    # calc mean and std of bins
    hist = [np.mean(k) for k in c_bin_array]
    hist_data_array.append(hist)
    c_data_array.append([hist[0]*hist[i] for i in range(len(hist))])
    # error = [np.std(k, dtype=np.float64) for k in c_bin_array]
    # lolims = [0 for k in error]
    g0_hist = [k for k in hist if np.isnan(k) == False]

    threshold = 80
    for i,c in enumerate(g0_hist):
        if c < threshold:
            niche.append(i)
            break
        if i == len(g0_hist)-1:
            niche.append(i)


    #############################################################################################################################################

    # for convolution purposes get the hist for every secreting cell, convolute with the sum (g0_hist) and average over it?
    # DISCUSS if this is the way to go
    sec_cells_hist_array = [[] for i in range(len(ids))]
    conv_array = [[] for i in range(len(ids))]
    # for cell in hist_df[plot_time_point][cell_df.loc[(cell_df["type_name"] == "default") & (cell_df["time_index"] == 0), "id"]]:
    for j, cell in enumerate(hist_df[plot_time_point][ids]): #.loc[hist_df[plot_time_point].index == str(11)].dropna():
        #sort the arrays and corresponding concentrations by distance via zipping und unzipping
        zipped = zip(cell[0], cell[1])
        sorted_zipped = sorted(list(zipped), key=lambda x: x[1])

        c, distances = zip(*sorted_zipped)
        c = list(c)
        distances = list(distances)

        # since "c" and "distances" are sorted, we can use np.where's first and last argument as a range for the entries of the bin
        # this way we only iterate over the bins, the iteration over the distances is done by numpy
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
        sec_cells_hist_array[j] = [np.mean(k) for k in c_bin_array]
        g1_hist = [np.mean(k) for k in c_bin_array if len(k) != 0]
        conv_array[j] = np.convolve(g0_hist, g0_hist, "same")

    # g1_hist = [np.mean(k) for k in np.transpose(sec_cells_hist_array)]
    # g1_hist = [np.mean(k) for k in c_bin_array]
    # plt.semilogy(np.convolve(g1_hist[:-1], g0_hist[:-1], "same"), label="Convolution")
    # plt.plot(np.convolve(g1_hist[0], np.array(g0_hist)), "-o", label=str(float(plot_time_point)*10)+"h")
    conv_data_array.append([np.mean(k) for k in np.transpose(conv_array)])
    # fig.add_subplot(a_x, a_y, 2)
    # ax2 = plt.plot([np.mean(k) for k in np.transpose(conv_array)], "-o", label=str(float(plot_time_point) * 10) + "h")

    #############################################################################################################################################

    # Autocorrelation with space instead of time:

    # auto_corr = []
    # for r_dash in range(len(g0_hist)-1):
    #     R = 0
    #     for r in range(len(g0_hist)):
    #         try:
    #             R += g0_hist[r + r_dash]*g0_hist[r]
    #         except IndexError:
    #             pass
    #     auto_corr.append(R)
    # auto_corr_data_array.append(auto_corr)

    #############################################################################################################################################

    # Build matrix from cell_df
    R = []
    matrix = np.zeros((len(cell_df["y"].unique()),len(cell_df["x"].unique())))

    cell_df_timed = cell_df.loc[cell_df["time"] == float(plot_time_point)].copy()
    for row in cell_df_timed.iterrows():
        matrix[int(row[1]["x"]/20)-1][int(row[1]["y"]/20)-1] = row[1]["IL-2_surf_c"]
    orig = matrix
    origFT = np.fft.fft(matrix, axis=1)

    # matrix = np.array([[1,2,3],[4,5,6],[7,8,9]])
    shifted = matrix
    for i in range(len(matrix)):
        # fade right and top
        deleted = np.delete(np.delete(shifted,0,0),-1,1)
        shifted = np.insert(np.insert(deleted, len(deleted)-1,0, axis=0), 0,0,axis=1)

        # halfed_hill = shifted[1:,:-1]
        # first = np.append(halfed_hill, [shifted[1:,-1]], 0)
        # second = np.append([[shifted[0,i]] for i in range(len(shifted[0]))], first, 1)
        # shifted = second

        # roll through
        # shifted = np.roll(np.roll(shifted, -1, axis=0), 1, axis=1)

        # fade right and bot
        # deleted = np.delete(np.delete(shifted,-1,0),-1,1)
        # shifted = np.insert(np.insert(deleted, 0,0, axis=0), 0,0,axis=1)

        # fade only right
        # deleted = np.delete(shifted,-1,1)
        # shifted = np.insert(deleted, 0,0, axis=1)

        # plt.matshow(matrix)
        # plt.show()
        # if plot_time_point == '15.0':
        #     plt.matshow(shifted)
        #     plt.show()
        # shifted = np.delete(np.delete(shifted,0,0),-1,1)

        shiftedFT = np.fft.fft(shifted, axis=1)
        dataAC = np.fft.ifft(origFT * np.conjugate(shiftedFT), axis=1).real
        R.append(np.sum(dataAC))
    auto_corr_data_array.append(R)


#############################################################################################################################################

# Fitting

AC_niche = []

def func(x, a, b, c):
    return a * b * np.exp(-x/b) + c
from scipy.optimize import curve_fit
for i in range(len(auto_corr_data_array)):
    popt, pcov = curve_fit(func, np.arange(len(auto_corr_data_array[i])), auto_corr_data_array[i], maxfev=10000)
    AC_niche.append(popt[1])
# plt.plot(np.arange(len(auto_corr_data_array[i])), func(np.arange(len(auto_corr_data_array[i])), *popt), 'r-', label='fit: a=%5.3f, b=%5.3f, c=%5.3f' % tuple(popt))
# plt.plot(auto_corr_data_array[i], label="R")
# plt.legend()
# plt.show()
# [k*10 for k in time_points]

#############################################################################################################################################

# Plotting

sns.set(rc={'figure.figsize':(8,8)})
fig = plt.figure()
plt.subplots_adjust(wspace=.3)
a_x = 3
a_y = 2

xscale = "linear"
yscale = "linear"

mean_sec_dist = np.mean(sec_dists)/cell_cell_distance
# if np.isnan(mean_sec_dist) == True:
#     mean_sec_dist = 0
fig.add_subplot(a_x, a_y, 1)
for i in range(len(hist_data_array)):
    ax1 = sns.lineplot(center[:len(hist_data_array[i])], hist_data_array[i], label=str(float(time_list[i])*10)+"h", legend=False)
# ax2 = sns.lineplot(center[:len(hist_data_array[1])], hist_data_array[1], label=str(float(time_list[1])*10)+"h", legend=False)
# ax3 = sns.lineplot(center[:len(hist_data_array[2])], hist_data_array[2], label=str(float(time_list[2])*10)+"h", legend=False)
plt.axvline(mean_sec_dist, color="black")
ax1.set(xlabel="cell-cell-distance", ylabel="Avg. c.", yscale=yscale, xscale=xscale, title="c_r", xticks=[0,5,10])

fig.add_subplot(a_x, a_y, 2)
for i in range(len(c_data_array)):
    ax4 = sns.lineplot(center[:len(c_data_array[i])], c_data_array[i], label=str(float(time_list[i])*10)+"h", legend=False)
# ax5 = sns.lineplot(center[:len(c_data_array[1])], c_data_array[1], label=str(float(time_list[1])*10)+"h", legend=False)
# ax6 = sns.lineplot(center[:len(c_data_array[2])], c_data_array[2], label=str(float(time_list[2])*10)+"h", legend=False)
plt.axvline(mean_sec_dist, label="sec. avg. dist.", color="black")
ax4.set(xlabel="cell-cell-distance", yscale=yscale, xscale=xscale, title="c_0*c_r", xticks=[0,5,10])

fig.add_subplot(a_x,a_y,3)
for i in range(len(conv_data_array)):
    ax7 = sns.lineplot(center[:len(conv_data_array[i])], conv_data_array[i], label=str(float(time_list[i])*10)+"h", legend=False)
# ax8 = sns.lineplot(center[:len(conv_data_array[1])], conv_data_array[1], label=str(float(time_list[1])*10)+"h", legend=False)
# ax9 = sns.lineplot(center[:len(conv_data_array[2])], conv_data_array[2], label=str(float(time_list[2])*10)+"h", legend=False)
plt.axvline(mean_sec_dist, label="sec. avg. dist.", color="black")
ax7.set(xlabel="cell-cell-distance", ylabel="Convolution", yscale=yscale, xscale=xscale, title="conv(c_r * c_r)", xticks=[0,5,10])

fig.add_subplot(a_x,a_y,4)
for i in range(len(auto_corr_data_array)):
    ax10 = sns.lineplot(center[:len(auto_corr_data_array[i])], auto_corr_data_array[i], label=str(float(time_list[i])*10)+"h", legend=False)
# ax11 = sns.lineplot(center[:len(auto_corr_data_array[1])], auto_corr_data_array[1], label=str(float(time_list[1])*10)+"h", legend=False)
# ax12 = sns.lineplot(center[:len(auto_corr_data_array[-1])], auto_corr_data_array[-1], label=str(float(time_list[-1])*10)+"h", legend=False)
ax10.set(xlabel="cell-cell-distance", ylabel="Auto-corr", yscale=yscale, xscale=xscale, title="Auto-correlation", xticks=[0,5,10])
plt.axvline(mean_sec_dist, label="sec. avg. dist.", color="black")
plt.legend(bbox_to_anchor=(1.04,0.5), loc="center left", borderaxespad=0)

fig.add_subplot(a_x,a_y,5)
ax13 = sns.lineplot([k*10 for k in time_points], niche, label="hist_niche", legend=False)
ax13 = sns.lineplot([k*10 for k in time_points], AC_niche, label="AC_niche")
ax13.set(xlabel="time in h", ylabel="cell-cell-distance", yscale=yscale, xscale=xscale, title="Niche size at threshold " + str(threshold), xticks=[0,5,10])


# fig.savefig("/home/brunner/Documents/Current work/31072020/" + "niche_" + str(cell_df["IL-2_fraction"].unique()[0]) + ".png", bbox_inches='tight')
plt.show()



#     # Plot the first histogram, just the concentration bins
#
#     width = 0.9 * (bins[1] - bins[0])
#     center = (bins[:-1] + bins[1:]) / 2
#
#     # plt.errorbar(center[:13], hist[:13], yerr=error[:13], label=str(float(plot_time_point)*10)+"h", fmt='-o', capsize=8)
#     plt.plot(center, hist,'-o', label=str(float(plot_time_point)*10)+"h")
#     plt.legend(loc=1)
# #     # plt.title("time = " + str(int((float(plot_time_point)+0)*10)) + "h")
#     plt.title("gamma = " +  str(cell_df["IL-2_gamma"].unique()[0]))
#     plt.xlabel("Cell-cell-distances to secreting cell")
#     plt.ylabel("Mean c. in pM")
#     # plt.ylim((0,380))
#     plt.xticks(bins[:]-10, [int(x/20) for x in bins][:])
#     plt.xlim((0,240))
# fig.savefig("/home/brunner/Documents/Current work/19062020/" + "wave_c01" + ".png", bbox_inches='tight')
# plt.show()
