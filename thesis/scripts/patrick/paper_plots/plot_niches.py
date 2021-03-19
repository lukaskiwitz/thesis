import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import seaborn as sns
import random
from thesis.main.my_debug import message
from scipy import spatial
from scipy.optimize import curve_fit
from scipy.constants import N_A

#load the cell_df
verbose = False
save_avg_hist_df = True

no_of_runs = 1
dt = 1

cell_cell_distance = 20

# base_path = "/extra/brunner/thesis/kinetic/q_fraction_medium_pos_g_scan_multi/"
base_path = "/extra/brunner/thesis/kinetic/kin_large_new_hill_4/"


runs_hist_df = pd.read_hdf(base_path + str(no_of_runs) + '_runs_0.25_f_hist_df.h5', mode="r")
avg_hist_df = pd.read_hdf(base_path + str(no_of_runs) + '_runs_0.25_f_avg_hist_df.h5', mode="r")
mean_sec_dist = float(np.loadtxt(base_path + str(no_of_runs) + "_runs_0.25_f_avg_sec_dists.txt"))/cell_cell_distance

global_df = pd.read_hdf(base_path + "global_df.h5", mode="r")
cell_df = pd.read_hdf(base_path + "cell_df.h5", mode="r")
fractioned_cell_df = cell_df.loc[(cell_df["IL-2_Tsec_fraction"] == 0.25)]
message("loaded dataframes from " + base_path)

timepoints = list(runs_hist_df.index[:-1])[:]
timepoints = [timepoints[x] for x in [5,20,37]]
# timepoints = np.array([0])*dt
# gammas = [runs_hist_df.columns[x][1] for x in range(len(runs_hist_df.columns))]
# gammas = [100,50,20,10]
gammas = [10.0]
thresholds = [10]

#############################################################################################################################################

# Niche
# Fitting

# AC_niche = []
#
# def func(x, a, b, c):
#     return a * b * np.exp(-x/b) + c
# from scipy.optimize import curve_fit
# for i in range(len(auto_corr_data_array)):
#     popt, pcov = curve_fit(func, np.arange(len(auto_corr_data_array[i])), auto_corr_data_array[i], maxfev=10000)
#     AC_niche.append(popt[1])


# Thresholding

# niche = []
#
# for j, hist in enumerate(hist_data_array):
#     if j == 0:
#         threshold = 200
#     for i,c in enumerate(hist):
#         if c < threshold:
#             niche.append(i)
#             break
#         if i == len(hist)-1:
#             niche.append(i)

# ax13 = sns.lineplot(points, niche, label="hist_niche", legend=False)
# ax13 = sns.lineplot(points[:], AC_niche[:], label="AC_niche")
# ax13.set(xlabel="time", ylabel="cell-cell-distance", yscale="linear", xscale="linear", title="Niche size at threshold " + str(threshold))
# plt.show()

#############################################################################################################################################

# Plotting of the niche selection process

# sns.set(rc={'figure.figsize':(12,10)})
# sns.set(palette="Greens")
# sns.set_style("ticks")
# sns.set_context("talk", font_scale=1, rc={"lines.linewidth": 3})

fig = plt.figure(figsize=(8,8))
# plt.figure(figsize=(12,10))
# sns.set_context("talk", font_scale=1.1, rc={"lines.linewidth": 2.5})
# sns.set(rc={'figure.figsize':(8,8)})
sns.set_style("ticks")
sns.set_context("talk", font_scale=1, rc={"lines.linewidth": 3, 'lines.markersize': 5})
plt.subplots_adjust(wspace=.4, hspace=0.6)
from math import ceil
if len(gammas) > 1:
    # a_x = ceil(len(slope_thresholds)/3)
    # a_y = ceil(len(slope_thresholds)/a_x)
    a_x = ceil(np.sqrt(len(gammas)))
    a_y = ceil(len(gammas)/a_x)
else:
    a_x = 1
    a_y = 1

xscale = "linear"
yscale = "linear"


if np.isnan(mean_sec_dist) == True:
    mean_sec_dist = 0

variable = "t"
# avg_hist_df = avg_hist_df.drop(20)

def func(x, a, b, c):
    return a * b * np.exp(-x / b) + c

AC_niche_over_gamma = [[] for gamma in gammas]
for g, gamma in enumerate(gammas):
    for t,time in enumerate(timepoints):
        hist = avg_hist_df[gamma][time]
        auto_corr = []
        for r_dash in range(len(hist) - 1):
            R = 0
            for r in range(len(hist)):
                try:
                    R += hist[r + r_dash] * hist[r]
                except IndexError:
                    pass
            auto_corr.append(R)

        popt, pcov = curve_fit(func, np.arange(len(auto_corr)), auto_corr, maxfev=10000)
        AC_niche_over_gamma[g].append(popt[1]/cell_cell_distance)
        plt.plot(auto_corr, label=str(gamma) + " " + str(round(time/3600, 1)) + "h")
        plt.plot(np.arange(len(auto_corr)), func(np.arange(len(auto_corr)), *popt), color="red", lw = 1)
        plt.xlabel("lag (c-c-d)")
        plt.ylabel("AC")

plt.legend()
# plt.show()

exponential_fit_over_gamma = [[] for gamma in gammas]
for g, gamma in enumerate(gammas):
    exp_niche = []
    for t,time in enumerate(timepoints):
        hist = avg_hist_df[gamma][time]
        popt, pcov = curve_fit(func, np.arange(len(hist)), hist, maxfev=10000)
        exponential_fit_over_gamma[g].append(popt[1])



niche_timepoint = runs_hist_df.index[-2]
cut_off = None

slope_thresholds = [0.03] #np.arange(0.01,0.08,0.01)

ylim = (None,None)

palette = sns.color_palette("bwr", len(slope_thresholds)*2)
# palette = palette[len(slope_thresholds):]

for s,slope_threshold in enumerate(slope_thresholds):
    slope_threshold = round(slope_threshold,2)

    niche_size_at_gamma = [[] for gamma in gammas]
    niche_effect_at_gamma = [[] for gamma in gammas]
    for g,gamma in enumerate(gammas):
        # niche_size_at_gamma.append([])
        # for run in range(runs_hist_df.columns[-1][0]+1):
        for t, time in enumerate(timepoints):
        # for niche_timepoint in [runs_hist_df.index[-2]]:
            niche = 0
            # for i in range(len(runs_hist_df.loc[niche_timepoint, (run,gamma)][:cut_off])):
            for i in range(len(avg_hist_df[gamma][time])):
                try:
                    # if runs_hist_df.loc[niche_timepoint, (run, gamma)][i] - runs_hist_df.loc[niche_timepoint, (run, gamma)][i + 1] < slope_threshold * runs_hist_df.loc[niche_timepoint, (run, gamma)][i + 1]:
                    if avg_hist_df[gamma][time][i] - avg_hist_df[gamma][time][i + 1] < slope_threshold * avg_hist_df[gamma][time][i + 1]:
                        niche_size_at_gamma[g].append(i)
                        niche = i
                        break
                    # if runs_hist_df.loc[niche_timepoint, (run,gamma)][i+1] - runs_hist_df.loc[niche_timepoint, (run,gamma)][i] < slope_threshold* runs_hist_df.loc[niche_timepoint, (run,gamma)][i]:
                    #     niche_size_at_gamma[g].append(i+1)
                    #     niche = i
                    #     break
                except IndexError:
                    niche_size_at_gamma[g].append(i)
                    niche = i

    # plt.plot(avg_hist_df[gamma].index.values, AC_niche_over_gamma[g], label=gammas[g], color="blue")
#     plt.plot(np.array(timepoints), niche_size_at_gamma[g], label=slope_threshold, color=palette[s*2])
#     # if gamma == gammas[-1] or gamma == gammas[-2]:
#     plt.xlabel("time in s")
#     plt.ylabel("Niche in c-c-d")
#     plt.title(gamma)
#     plt.legend()
# plt.show()



# gammas =[1/10]
palette = sns.color_palette("bwr", len(thresholds)*2)
for th,threshold in enumerate(thresholds):

    threshold_niche_size_at_gamma = [[] for gamma in gammas]
    threshold_niche_effect_at_gamma = [[] for gamma in gammas]
    for g,gamma in enumerate(gammas):
        # niche_size_at_gamma.append([])
        # for run in range(runs_hist_df.columns[-1][0]+1):
        for t, time in enumerate(timepoints):
        # for niche_timepoint in [runs_hist_df.index[-2]]:
            niche = 0
            # for i in range(len(runs_hist_df.loc[niche_timepoint, (run,gamma)][:cut_off])):
            for i in range(len(avg_hist_df[gamma][time])+1):
                try:
                    # if runs_hist_df.loc[niche_timepoint, (run, gamma)][i] - runs_hist_df.loc[niche_timepoint, (run, gamma)][i + 1] < slope_threshold * runs_hist_df.loc[niche_timepoint, (run, gamma)][i + 1]:
                    if avg_hist_df[gamma][time][i] < threshold:
                        threshold_niche_size_at_gamma[g].append(i-1)
                        niche = i-1
                        break
                    # if runs_hist_df.loc[niche_timepoint, (run,gamma)][i+1] - runs_hist_df.loc[niche_timepoint, (run,gamma)][i] < slope_threshold* runs_hist_df.loc[niche_timepoint, (run,gamma)][i]:
                    #     niche_size_at_gamma[g].append(i+1)
                    #     niche = i
                    #     break
                except IndexError:
                    threshold_niche_size_at_gamma[g].append(i)
                    niche = i-1

        # plt.plot(avg_hist_df[gamma].index.values, AC_niche_over_gamma[g], label=gammas[g], color="blue")
        # plt.plot(threshold_niche_size_at_gamma[g], label=threshold, color=palette[th*2])
        # # if gamma == gammas[-1] or gamma == gammas[-2]:
        # plt.xlabel("time in s")
        # plt.ylabel("Niche in c-c-d")
        # plt.title(gamma)
        # plt.legend()
        # plt.show()

# exit()



EC50_niche_size_at_gamma = [[] for gamma in gammas]
EC50_list = [[] for gamma in gammas]

for g,gamma in enumerate(gammas):
    # niche_size_at_gamma.append([])
    # for run in range(runs_hist_df.columns[-1][0]+1):
    for t, time in enumerate(timepoints):
    # for niche_timepoint in [runs_hist_df.index[-2]]:
        niche = 0
        # for i in range(len(runs_hist_df.loc[niche_timepoint, (run,gamma)][:cut_off])):
        mean_R = fractioned_cell_df.loc[(np.abs(fractioned_cell_df["time"] - time) < 0.001) & (fractioned_cell_df["IL-2_gamma"] == gamma) &  (fractioned_cell_df["type_name"] == "Tsec"), "IL-2_R"].mean()
        mean_R = mean_R * N_A ** -1 * 1e9 #into M
        R_k_max = 4e-11  # max pSTAT5 EC50 in M
        R_k_min = 2e-12  # min pSTAT55 EC50 in M
        R_k = 1e4 * N_A ** -1 * 1e9  # Receptors at half pSTAT5 maximum
        R_N = 1.2  # hill coefficient
        EC50 = (R_k_max * R_k ** R_N + R_k_min * mean_R ** R_N) / (R_k ** R_N + mean_R ** R_N)
        EC50 *= 1e12 #into pM
        EC50_list[g].append(EC50)
        print(EC50)
        for i in range(len(avg_hist_df[gamma][time])+1):
            try:
                # if runs_hist_df.loc[niche_timepoint, (run, gamma)][i] - runs_hist_df.loc[niche_timepoint, (run, gamma)][i + 1] < slope_threshold * runs_hist_df.loc[niche_timepoint, (run, gamma)][i + 1]:
                if avg_hist_df[gamma][time][i] < EC50:
                    EC50_niche_size_at_gamma[g].append(i)
                    niche = i
                    break
                # if runs_hist_df.loc[niche_timepoint, (run,gamma)][i+1] - runs_hist_df.loc[niche_timepoint, (run,gamma)][i] < slope_threshold* runs_hist_df.loc[niche_timepoint, (run,gamma)][i]:
                #     niche_size_at_gamma[g].append(i+1)
                #     niche = i
                #     break
            except IndexError:
                EC50_niche_size_at_gamma[g].append(i)
                niche = i

    #
    # plt.plot(timepoints, EC50_list[g], label=gamma, color=palette[th*2])
    # # if gamma == gammas[-1] or gamma == gammas[-2]:
    # plt.xlabel("time in s")
    # plt.ylabel("Niche in c-c-d")
    # plt.title(gamma)
    # plt.legend()
    # plt.show()



############ Plotting #############
fig = plt.figure(figsize=(8,8))
sns.set_style("ticks")
sns.set_context("talk", font_scale=0.7, rc={"lines.linewidth": 3, 'lines.markersize': 5})
plt.subplots_adjust(wspace=.4, hspace=0.6)
if len(gammas) > 1:
    # a_x = ceil(len(slope_thresholds)/3)
    # a_y = ceil(len(slope_thresholds)/a_x)
    a_x = ceil(np.sqrt(len(gammas)))
    a_y = ceil(len(gammas)/a_x)
else:
    a_x = 1
    a_y = 1

xscale = "linear"
yscale = "linear"


for g,gamma in enumerate(gammas):
    fig.add_subplot(a_x, a_y, g + 1)
    plt.plot(np.array(timepoints)/3600,AC_niche_over_gamma[g], label="AC", color="blue")
    plt.plot(np.array(timepoints)/3600, niche_size_at_gamma[g], label="slope", color="red")
    plt.plot(np.array(timepoints) / 3600, threshold_niche_size_at_gamma[g], label="threshold", color="green")
    plt.plot(np.array(timepoints) / 3600, EC50_niche_size_at_gamma[g], label="EC50", color="orange")
    plt.plot(np.array(timepoints) / 3600, exponential_fit_over_gamma[g], label="exp fit", color="black")
    # if gamma == gammas[-1] or gamma == gammas[-2]:
    #     plt.xlabel("time in s")
    # plt.ylabel("Niche in c-c-d")
    if gamma < 1:
        plt.title("g = 1/" + str(int(1 / gamma)))
    else:
        plt.title("g = " + str(gamma))
    if gamma == gammas[-1]:
        plt.legend()
# plt.show()

print("AC", AC_niche_over_gamma[0])
print("slope", niche_size_at_gamma[0])
print("threshold", threshold_niche_size_at_gamma[0])
print("EC50", EC50_niche_size_at_gamma[0])
print("exp", exponential_fit_over_gamma[0])

