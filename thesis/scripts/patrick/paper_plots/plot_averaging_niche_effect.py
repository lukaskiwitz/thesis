import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import seaborn as sns
import random

from pandas import DataFrame

from my_debug import message
from scipy import spatial
from scipy.optimize import curve_fit
from scipy.constants import N_A


try:
    from thesis.scripts.patrick.niches.new_AC import AC_1_1, AC_1_2, AC_2, AC_2_1, AC_2_2
except ModuleNotFoundError:
    exec(open('/home/brunner/thesis/thesis/scripts/patrick/niches/new_AC.py').read())


#load the df

base_path = "/extra/brunner/thesis/kinetic/standard/8_big_scan/"
global_df = pd.read_hdf(base_path + "global_df.h5", mode="r")
cell_df = pd.read_hdf(base_path + "cell_df.h5", mode="r")


verbose = False
save_avg_hist_df = True
no_of_runs = 1
dt = 3600
cell_cell_distance = 20
max_distance = 13
pSTAT5 = True


frac = 0.25

# base_path = "/extra/brunner/thesis/kinetic/q_fraction_medium_pos_g_scan_multi/"
if pSTAT5 == True:
    runs_hist_df = pd.read_hdf(base_path + str(no_of_runs) + '_runs_' + str(frac) + '_Tsec_fraction_pSTAT5_hist_df.h5', mode="r")
    avg_hist_df = pd.read_hdf(base_path + str(no_of_runs) + '_runs_' + str(frac) + '_Tsec_fraction_pSTAT5_avg_hist_df.h5', mode="r")
else:
    runs_hist_df = pd.read_hdf(base_path + str(no_of_runs) + '_runs_' + str(frac) + '_Tsec_fraction_hist_df.h5', mode="r")
    avg_hist_df = pd.read_hdf(base_path + str(no_of_runs) + '_runs_' + str(frac) + '_Tsec_fraction_avg_hist_df.h5', mode="r")
mean_sec_dist = float(np.loadtxt(base_path + str(no_of_runs) + "_runs_" + str(frac) + "_Tsec_fraction_avg_sec_dists.txt"))/cell_cell_distance

fractioned_cell_df = cell_df.loc[(np.abs(cell_df["IL-2_Tsec_fraction"] - frac) < 0.0001)]

if pSTAT5 == True:
    pSTAT5_sec_histograms = pd.read_hdf(base_path + "pSTAT5_hist_dfs/sec_pSTAT5_" + str(frac) + "_Tsec_fraction_histograms.h5", mode="r")
sec_histograms = pd.read_hdf(base_path + "hist_dfs/sec_" + str(frac) + "_Tsec_fraction_histograms.h5", mode="r")

message("loaded runs_hist_df from " + base_path)

gammas = [1e-2,1e2]
timepoints = list(sec_histograms[gammas[0]].columns)[1:-1]
AC_over_gamma: DataFrame = pd.DataFrame(columns=gammas, index=["1_2"])
AC_1_2_autocorr_g = [[] for x in gammas]
unavg_EC50_niche_over_sv = [[] for x in gammas]

#############################################################################################################################################

message("calculating EC50")
EC50_list = [[] for gamma in gammas]
for g,gamma in enumerate(gammas):
    # Thresholding with EC50
    mean_R_list = []
    for t, time in enumerate(timepoints):
        # for niche_timepoint in [runs_hist_df.index[-2]]:
        niche = 0
        # for i in range(len(runs_hist_df.loc[niche_timepoint, (run,gamma)][:cut_off])):
        mean_R = fractioned_cell_df.loc[
            (np.abs(fractioned_cell_df["time"] - float(time)) < 0.0001) & (fractioned_cell_df["IL-2_gamma"] == gamma) & (
                        fractioned_cell_df["type_name"] == "Treg"), "IL-2_R"].mean()
        mean_R_list.append(mean_R)
        mean_R = mean_R * N_A ** -1 * 1e9  # into M
        R_k_max = 4e-11  # max pSTAT5 EC50 in M
        R_k_min = 2e-12  # min pSTAT55 EC50 in M
        R_k = 1e4 * N_A ** -1 * 1e9  # Receptors at half pSTAT5 maximum
        R_N = 1.2  # hill coefficient
        EC50 = (R_k_max * R_k ** R_N + R_k_min * mean_R ** R_N) / (R_k ** R_N + mean_R ** R_N)
        EC50 *= 1e12  # into pM
        # EC50 = 21.0
        EC50_list[g].append(EC50)




# calculate AC niches:

def func(x, a, b, c):
    return a * np.exp(-x / b) + c
def func2(x, a, b):
    return a - x / b


if pSTAT5 == True:
    # AC_over_gamma[gamma]["1_1"] = AC_1_1(pSTAT5_sec_histograms, gammas, timepoints, cell_cell_distance, max_distance, func)[0][0]
    AC_1_2_niche, AC_1_2_autocorr = AC_1_2(pSTAT5_sec_histograms, gammas, timepoints, cell_cell_distance, max_distance, func2, plot=False, linear=True, p0=None)
    for g, gamma in enumerate(gammas):
        AC_over_gamma[gamma]["1_2"] = AC_1_2_niche[g]
        AC_1_2_autocorr_g[g].append(AC_1_2_autocorr[g])
        # AC_over_gamma[gamma]["2"] = AC_2(sec_histograms, gammas, timepoints, cell_cell_distance, max_distance, func, plot=False)[0][0]
        # AC_over_gamma[gamma]["2_1"] = AC_2_1(sec_histograms, gammas, timepoints, cell_cell_distance, max_distance, func)[0][0]
        # AC_over_gamma[gamma]["2_2"] = AC_2_2(sec_histograms, gammas, timepoints, cell_cell_distance, max_distance, func, plot=False)[0][0]

message("calculating niche_avg_g")
niche_avg_g = [[] for g in gammas]
niche_effect_g = [[] for g in gammas]
for g,gamma in enumerate(gammas):
    for t, time in enumerate(timepoints):
        sec_niches = []
        concentrations = []
        for i in range(len(sec_histograms[gamma][str(time)])):
            for j in range(len(sec_histograms[gamma][str(time)].iloc[i]) + 1):
                try:
                    # if runs_hist_df.loc[niche_timepoint, (run, gamma)][i] - runs_hist_df.loc[niche_timepoint, (run, gamma)][i + 1] < slope_threshold * runs_hist_df.loc[niche_timepoint, (run, gamma)][i + 1]:
                    if sec_histograms[gamma][str(time)].iloc[i][j] < EC50_list[g][t]:
                        sec_niches.append(j)
                        break
                except IndexError:
                    sec_niches.append(j)
            for k in range(sec_niches[-1]):
                concentrations.append(sec_histograms[gamma][str(time)].iloc[i][k])
        niche_avg_g[g].append(np.mean(sec_niches))
        # calculate niche effect: take the timepoints niche and get the average concentration within that niche

        # print(concentrations)
        if no_of_runs > 1:
            try:
                global_df = pd.read_hdf(base_path + "run" + str(run) + "/global_df.h5", mode="r")
                average_concentration = global_df.loc[(np.abs(global_df["time"] - float(time)) < 0.0001) & (
                        global_df["IL-2_gamma"] == gamma), "surf_c"].values[0] * 1e3

            except IndexError:
                try:
                    global_df = pd.read_hdf(base_path_pos + "run" + str(run) + "/global_df.h5", mode="r")
                    average_concentration = global_df.loc[(np.abs(global_df["time"] - float(time)) < 0.0001) & (
                            global_df["IL-2_gamma"] == gamma), "surf_c"].values[0] * 1e3

                except IndexError:
                    global_df = pd.read_hdf(base_path_neg + "run" + str(run) + "/global_df.h5", mode="r")
                    average_concentration = global_df.loc[(np.abs(global_df["time"] - float(time)) < 0.0001) & (
                            global_df["IL-2_gamma"] == gamma), "surf_c"].values[0] * 1e3
        else:
            try:
                global_df = pd.read_hdf(base_path + "/global_df.h5", mode="r")
                average_concentration = global_df.loc[(np.abs(global_df["time"]- float(time)) < 0.0001) & (
                        global_df["IL-2_gamma"] == gamma), "surf_c"].values[0] * 1e3

            except IndexError:
                try:
                    global_df = pd.read_hdf(base_path_pos + "/global_df.h5", mode="r")
                    average_concentration = global_df.loc[(np.abs(global_df["time"] - float(time)) < 0.0001) & (
                            global_df["IL-2_gamma"] == gamma), "surf_c"].values[0] * 1e3

                except IndexError:
                    global_df = pd.read_hdf(base_path_neg + "/global_df.h5", mode="r")
                    average_concentration = global_df.loc[(np.abs(global_df["time"] - float(time)) < 0.0001) & (
                            global_df["IL-2_gamma"] == gamma), "surf_c"].values[0] * 1e3

        niche_effect_g[g].append(np.mean(concentrations) / average_concentration)

#############################################################################################################################################

# load average secretory distance and norm to cell-cell-distance
mean_sec_dist = float(np.loadtxt(base_path + str(no_of_runs) + "_runs_0.25_Tsec_fraction_avg_sec_dists.txt"))/cell_cell_distance

niche_timepoints_index = 15
niche_timepoint = runs_hist_df.index[niche_timepoints_index]
cut_off = None



ylim = (None,None)
fig = plt.figure(figsize=(12,12))
# plt.figure(figsize=(12,10))
# sns.set_context("talk", font_scale=1.1, rc={"lines.linewidth": 2.5})
# sns.set(rc={'figure.figsize':(8,8)})
sns.set_style("ticks")
sns.set_context("talk", font_scale=1, rc={"lines.linewidth": 3, 'lines.markersize': 5})
plt.subplots_adjust(wspace=.4, hspace=0.6)


#############################################################################################################################################

# Plotting

palette = sns.color_palette("bwr", len(gammas))
palette = [i for i in reversed(palette)]
for g,gamma in enumerate(gammas):
    fig = plt.figure(figsize=(20,6))
    sns.set_style("ticks")
    sns.set_context("talk", font_scale=1.3, rc={"lines.linewidth": 4, 'lines.markersize': 6})
    plt.subplots_adjust(wspace=.25, hspace=0.6)
    a_x = 1
    a_y = 2

    xscale = "linear"
    yscale = "linear"

    if np.isnan(mean_sec_dist) == True:
        mean_sec_dist = 0

    variable = "t"
    if gamma < 1:
        gamma_text = "negative feedback"
    if gamma > 1:
        gamma_text = "positive feedback"
    if gamma == 1:
        gamma_text = "no feedback"
    # fig = plt.figure()
    # sns.set_style("ticks")
    # sns.set_context("talk", font_scale=1.5, rc={"lines.linewidth": 5})
    # fig, ax = plt.subplots(figsize=(6, 5))
    # fig.add_subplot(a_x, a_y, g+1)

    x = niche_effect_g
    y = niche_avg_g


    for i in range(len(x)):
        # if gammas[i] > 1:
        #     y[i] += np.random.normal()/30
        # else:
        #     myColor = "black"
        plt.scatter(np.mean(x[i]), np.mean(y[i]), color=palette[i])
        # plt.errorbar(x[i], y[i], xerr=[[np.abs(x)] for x in niche_effect_error[i]], yerr=[[np.abs(x)] for x in niche_size_error[i]], linestyle='', c=palette[i])
    # plt.scatter(216.35519872163854, 1.6372549019607845, color="black")
    # plt.errorbar(216.35519872163854, 1.6372549019607845,xerr=13.160380344055845, yerr=0.22201019379779285, color="black")
    sns.set_context(rc={"lines.linewidth": 2})


    plt.xlabel("niche effect ($\%$)")
    plt.ylabel("niche size")

    # plt.title(str(round(slope_threshold*100,2)) + "%")
    # plt.ylim(ylim)
plt.show()
