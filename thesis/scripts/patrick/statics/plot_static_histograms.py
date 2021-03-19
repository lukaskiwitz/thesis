import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import seaborn as sns
import random
from my_debug import message
from scipy import spatial
from scipy.optimize import curve_fit
from scipy.constants import N_A

try:
    from new_AC import AC_1_1, AC_1_2, AC_2, AC_2_1, AC_2_2
except ModuleNotFoundError:
    exec(open('/home/brunner/thesis/thesis/scripts/patrick/niches/new_AC.py').read())


def histogram_plotting(avg_hist_df, gamma, time, color, label, exp_fit=False):
    plt.plot(avg_hist_df[gamma].loc[np.abs(avg_hist_df.index - float(time)) < 0.001].iloc[0], label=label, color=color, marker="o")
    if exp_fit == True:
        data = avg_hist_df[gamma].loc[np.abs(avg_hist_df.index - float(time)) < 0.001].iloc[0]
        data = data[~np.isnan(data)]
        def func(x,a,b,c):
            return a/x * np.exp(-x / b) + c
        myRange = np.arange(1, len(data)+1, 1)
        popt, pcov = curve_fit(func, myRange, data, maxfev=10000, p0=[15,10, 190])
        plotting_range = np.arange(1, len(data)+1, 1)
        plt.plot(plotting_range-1, func(plotting_range, *popt), color="black", lw=1)
    return pcov
        #
        #
        # def func(x,a,b,c):
        #     return a/x +0*b*c
        # plt.plot(np.arange(1, len(data), 0.01), func(np.arange(1, len(data), 0.01), 200, 1.3e100, 2.6e2), color="black", lw=1)
        # plt.plot(avg_hist_df[gamma].loc[np.abs(avg_hist_df.index - float(time)) < 0.001].iloc[0],
        #          color="red", marker="o")
        # plt.show()


#load the cell_df

base_path = "/extra/brunner/thesis/static/static_exp_decay/"
global_df = pd.read_hdf(base_path + "global_df.h5", mode="r")
cell_df = pd.read_hdf(base_path + "cell_df.h5", mode="r")

verbose = False
save_avg_hist_df = True
no_of_runs = 1
dt = 1
cell_cell_distance = 20
max_distance = 13
pSTAT5 = True

scan_variables = np.arange(0.01,0.25,0.02)[:]
df_string = "IL-2_Tsec_fraction"
saving_string = "_Tsec_fraction"

AC_niches_over_sv = pd.DataFrame(columns=np.round(scan_variables,2), index=["1_2"])
EC50_niches_sv = [[] for f in scan_variables]
EC50_sv = [[] for f in scan_variables]

unavg_EC50_niche_over_sv = [[] for f in scan_variables]

AC_1_2_autocorr_sv = [[] for x in scan_variables]

for sv,scan_v in enumerate(scan_variables):
    if len(str(scan_v)) > 5:
        scan_v = round(scan_v,2)
    # base_path = "/extra/brunner/thesis/kinetic/q_fraction_medium_pos_g_scan_multi/"

    runs_hist_df = pd.read_hdf(base_path + str(no_of_runs) + '_runs_' + str(scan_v) + saving_string + '_hist_df.h5', mode="r")
    avg_hist_df = pd.read_hdf(base_path + str(no_of_runs) + '_runs_' + str(scan_v) + saving_string + '_avg_hist_df.h5', mode="r")
    mean_sec_dist = float(np.loadtxt(base_path + str(no_of_runs) + '_runs_' + str(scan_v) + saving_string + '_avg_sec_dists.txt'))/cell_cell_distance

    fractioned_cell_df = cell_df.loc[(np.abs(cell_df[df_string] - scan_v) < 0.0001)]

    if pSTAT5 == True:
        pSTAT5_sec_histograms = pd.read_hdf(base_path + "pSTAT5_hist_dfs/sec_pSTAT5" + "_" + str(scan_v) + saving_string +  "_histograms.h5", mode="r")
    sec_histograms = pd.read_hdf(base_path + "hist_dfs/sec" + "_" + str(scan_v) + saving_string + "_histograms.h5", mode="r")

    message("loaded runs_hist_df from " + base_path)

    timepoints = list(runs_hist_df.index[:-1])
    # timepoints = [timepoints[x] for x in [3,4,5,7,11,15,19]]
    timepoints = [str(2.0)]
    # gammas = [runs_hist_df.columns[x][1] for x in range(len(runs_hist_df.columns))]
    gammas = [1.]
    gamma = 1.

    #############################################################################################################################################

    # Thresholding with EC50
    for t, time in enumerate(timepoints):
        # for niche_timepoint in [runs_hist_df.index[-2]]:
        niche = 0
        # for i in range(len(runs_hist_df.loc[niche_timepoint, (run,gamma)][:cut_off])):
        mean_R = fractioned_cell_df.loc[
            (np.abs(fractioned_cell_df["time"] - float(time)) < 0.0001) & (
                    fractioned_cell_df["IL-2_gamma"] == gamma) & (
                    fractioned_cell_df["type_name"] == "Treg"), "IL-2_R"].mean()
        mean_R = mean_R * N_A ** -1 * 1e9  # into M
        R_k_max = 4e-11  # max pSTAT5 EC50 in M
        R_k_min = 2e-12  # min pSTAT55 EC50 in M
        R_k = 1e4 * N_A ** -1 * 1e9  # Receptors at half pSTAT5 maximum
        R_N = 1.2  # hill coefficient
        EC50 = (R_k_max * R_k ** R_N + R_k_min * mean_R ** R_N) / (R_k ** R_N + mean_R ** R_N)
        EC50 *= 1e12  # into pM
        EC50_sv[sv].append(EC50)
        for i in range(len(avg_hist_df[gammas[0]][float(time)]) + 1):
            try:
                # if runs_hist_df.loc[niche_timepoint, (run, gamma)][i] - runs_hist_df.loc[niche_timepoint, (run, gamma)][i + 1] < slope_threshold * runs_hist_df.loc[niche_timepoint, (run, gamma)][i + 1]:
                if avg_hist_df[gamma][float(time)][i] < EC50:
                    EC50_niches_sv[sv].append(i)
                    niche = i
                    break
                # if runs_hist_df.loc[niche_timepoint, (run,gamma)][i+1] - runs_hist_df.loc[niche_timepoint, (run,gamma)][i] < slope_threshold* runs_hist_df.loc[niche_timepoint, (run,gamma)][i]:
                #     niche_size_at_gamma[g].append(i+1)
                #     niche = i
                #     break
            except IndexError:
                EC50_niches_sv[sv].append(i)
                niche = i



    def func(x, a, b, c):
        return a * np.exp(-x / b) + c
    def func2(x, a, b):
        return a - x / b

    # calculate AC niches:
    if pSTAT5 == True:
        # AC_over_gamma[gamma]["1_1"] = AC_1_1(pSTAT5_sec_histograms, gammas, timepoints, cell_cell_distance, max_distance, func)[0][0]
        AC_1_2_niche, AC_1_2_autocorr = AC_1_2(pSTAT5_sec_histograms, gammas, timepoints, cell_cell_distance, max_distance, func2, plot=False, linear=True, p0=None)
        AC_niches_over_sv[scan_v]["1_2"] = AC_1_2_niche[0][0]
        AC_1_2_autocorr_sv[sv].append(AC_1_2_autocorr[0][0])

    
    for t, time in enumerate(timepoints):
        sec_niches = []
        for i in range(len(sec_histograms[gamma][str(time)])):
            for j in range(len(sec_histograms[gamma][str(time)].iloc[i]) + 1):
                try:
                    # if runs_hist_df.loc[niche_timepoint, (run, gamma)][i] - runs_hist_df.loc[niche_timepoint, (run, gamma)][i + 1] < slope_threshold * runs_hist_df.loc[niche_timepoint, (run, gamma)][i + 1]:
                    if sec_histograms[gamma][str(time)].iloc[i][j] < EC50_sv[sv][t]:
                        sec_niches.append(j)
                        break
                except IndexError:
                    sec_niches.append(j)
        unavg_EC50_niche_over_sv[sv].append(np.mean(sec_niches))
#############################################################################################################################################

# AC_niches_over_sv[sv] = AC_1_2(pSTAT5_sec_histograms, gammas, timepoints, cell_cell_distance, max_distance, func, plot=True)[1][0][0]
#
# colors = ["blue","orange","green", "red"]
#
# for e,entry in enumerate(AC_1_2_autocorr_sv):
#     popt, pcov = curve_fit(func, np.arange(len(entry)), entry, maxfev=10000)
#     plt.semilogy(np.arange(len(entry)), func(np.arange(len(entry)), *popt), "--", color= colors[e], lw=1)
#     plt.semilogy(entry, "-o", label=scan_variables[e], color=colors[e])
#     plt.xlabel("c-c-d")
#     plt.ylabel("AC")
# plt.legend()
# # plt.savefig("/home/brunner/Documents/Current work/16102020/static_AC_fit_log.png", bbox_inches='tight')
#
# plt.show()


# Plotting
AC_niches_over_sv_T = AC_niches_over_sv.T


fig = plt.figure(figsize=(20,6))
sns.set_style("ticks")
sns.set_context("talk", font_scale=1, rc={"lines.linewidth": 3, 'lines.markersize': 5})
plt.subplots_adjust(wspace=.4, hspace=0.6)
a_x = 1
a_y = 2

xscale = "linear"
yscale = "linear"

if np.isnan(mean_sec_dist) == True:
    mean_sec_dist = 0

variable = "t"

fig.add_subplot(a_x, a_y, 2)
for AC_method in ["1_2"]:
    if len(AC_niches_over_sv_T.index.values) == 1:
        plt.scatter(AC_niches_over_sv_T.index.values, AC_niches_over_sv_T[AC_method].values, label=AC_method)
    else:
        plt.plot(AC_niches_over_sv_T.index.values, AC_niches_over_sv_T[AC_method].values, label="AC")
# plt.plot(scan_variables, EC50_niches_sv, label="EC50 Threshold", color="black")
plt.plot(scan_variables, unavg_EC50_niche_over_sv, label="Niche size by EC50", color="green")
plt.ylabel("niche size (c-c-d)")
plt.xlabel(df_string)
# plt.ylim((-0.1,8.5))
plt.title("pSTAT5 related")
plt.legend()


fig.add_subplot(a_x, a_y, 1)
scan_variables = [scan_variables[x] for x in [0,1,5,10]]
palette = sns.color_palette("bwr", len(scan_variables)*2)
palette = palette[len(scan_variables):]
# palette = ["red"]
perr_list = []
for sv,scan_v in enumerate(scan_variables):
    if len(str(scan_v)) > 5:
        scan_v = round(scan_v,2)
    avg_hist_df = pd.read_hdf(base_path + str(no_of_runs) + '_runs_' + str(scan_v) + saving_string + '_avg_hist_df.h5', mode="r")
    pcov = histogram_plotting(avg_hist_df,gammas[0],float(timepoints[0]),palette[sv], label="f = " + str(scan_v), exp_fit=True)
    perr = np.sqrt(np.diag(pcov)[1])
    perr_list.append(perr)

# ax2 = sns.lineplot(center[:len(hist_data_array[1])], hist_data_array[1], label=str(float(time_list[1])*10)+"h", legend=False)
# ax3 = sns.lineplot(center[:len(hist_data_array[2])], hist_data_array[2], label=str(float(time_list[2])*10)+"h", legend=False)
# plt.axvline(mean_sec_dist, color="black", linestyle="dashed")
plt.title("Histograms")
# ax1.set(xlabel="cell-cell-distance", ylabel="Avg. c. in pM", yscale=yscale, xscale=xscale, xticks=[0,5,10], title="g = " + str(gamma)+" , "+gamma_str)
plt.xticks([0,2,4,6,8,10,12])
plt.xlabel("c-c-d")
plt.ylabel("Concentration (pM)")
plt.legend()
# plt.savefig("/home/brunner/Documents/Current work/30102020/static_exp_fit.png", bbox_inches='tight')
plt.show()