import matplotlib.pyplot as plt
import numpy as np

import sys
sys.path.insert(0,'/home/brunner/thesis/thesis/main')
from my_debug import message

import pandas as pd
import seaborn as sns

from pandas import DataFrame

from scipy import spatial
from scipy.optimize import curve_fit
from scipy.constants import N_A
from scipy.stats import chisquare


try:
    from new_AC import AC_1_2, EC50, histogram_plotting, niche_avg
except ModuleNotFoundError:
    exec(open('/home/brunner/thesis/thesis/scripts/patrick/niches/new_AC.py').read())


#load the cell_df
prefix = "" #"/run/user/6205/gvfs/sftp:host=compute3"
# base_path = prefix + "/extra2/brunner/thesis/kinetic/standard/8_big_scan/"
base_path = r"/home/brunner/thesis/thesis/scripts/patrick/yukawa_analytical/test/"
global_df = pd.read_hdf(base_path + "global_df.h5", mode="r")
cell_df = pd.read_hdf(base_path + "cell_df.h5", mode="r")


verbose = False
save_avg_hist_df = True
no_of_runs = 1
dt = 1
cell_cell_distance = 20
max_distance = 13
pSTAT5 = False

frac = 0.25

offset = 2
hist_cutoff = 11

saving_string = "_Tsec_fraction"
# base_path = "/extra/brunner/thesis/kinetic/q_fraction_medium_pos_g_scan_multi/"
if pSTAT5 == True:
    runs_hist_df = pd.read_hdf(base_path + str(no_of_runs) + '_runs_' + str(frac) + saving_string + '_pSTAT5_hist_df.h5', mode="r")
    avg_hist_df = pd.read_hdf(base_path + str(no_of_runs) + '_runs_' + str(frac) + saving_string + '_pSTAT5_avg_hist_df.h5', mode="r")
else:
    runs_hist_df = pd.read_hdf(base_path + str(no_of_runs) + '_runs_' + str(frac) + saving_string + '_hist_df.h5', mode="r")
    avg_hist_df = pd.read_hdf(base_path + str(no_of_runs) + '_runs_' + str(frac) + saving_string + '_avg_hist_df.h5', mode="r")
mean_sec_dist = float(np.loadtxt(base_path + str(no_of_runs) + "_runs_" + str(frac) + saving_string + "_avg_sec_dists.txt"))/cell_cell_distance

fractioned_cell_df = cell_df.loc[(np.abs(cell_df["IL-2_Tsec_fraction"] - frac) < 0.0001)]

if pSTAT5 == True:
    pSTAT5_sec_histograms = pd.read_hdf(base_path + "pSTAT5_hist_dfs/sec_pSTAT5_" + str(frac) + saving_string + "_histograms.h5", mode="r")
sec_histograms = pd.read_hdf(base_path + "hist_dfs/sec_" + str(frac) + saving_string + "_histograms.h5", mode="r")

message("loaded runs_hist_df from " + base_path)

# gammas = sorted(global_df["IL-2_gamma"].unique())[:10:4] + sorted(global_df["IL-2_gamma"].unique())[11::4]
gammas = [0.01, 0.02, 0.1, 10, 50, 100]
timepoints = list(sec_histograms[gammas[0]].columns)[:-7]
# timepoints = [timepoints[idx] for idx,x in enumerate(timepoints) if idx%2 == 0]
# gammas = [runs_hist_df.columns[x][1] for x in range(len(runs_hist_df.columns))]

AC_over_gamma: DataFrame = pd.DataFrame(columns=gammas, index=["1_2"])
EC50_niches_over_gamma = [[] for f in gammas]
AC_1_2_autocorr_g = [[] for x in gammas]
unavg_EC50_niche_over_sv = [[] for x in gammas]

#############################################################################################################################################


message("calculating EC50_niches_over_gamma")
tmp_df = fractioned_cell_df
EC50_list = EC50(gammas,timepoints,fractioned_cell_df)


# calculate AC niches:
def func(x, a, b, c):
    return a * np.exp(-x / b) + c
def func2(x, a, b):
    return a - x / b

if pSTAT5 == True:
    AC_1_2_niche, AC_1_2_autocorr = AC_1_2(pSTAT5_sec_histograms, gammas, timepoints, cell_cell_distance, max_distance, func2, plot=False, linear=True, p0=None)
    for g, gamma in enumerate(gammas):
        AC_over_gamma[gamma]["1_2"] = AC_1_2_niche[g]
        AC_1_2_autocorr_g[g].append(AC_1_2_autocorr[g])

message("calculating niche_avg_g")
niche_avg_g = niche_avg(gammas,timepoints,fractioned_cell_df, offset,sec_histograms,EC50_list)

test_secs = sec_histograms.copy()
for g, gamma in enumerate(gammas):
    for t, time in enumerate(timepoints):
        sec_cells_ids = test_secs[gamma][str(time)].index
        for i in sec_cells_ids:
            for j in range(len(test_secs[gamma][str(time)].loc[i])):
                test_secs[gamma][str(time)].loc[i][j] *= 1/3

test = niche_avg(gammas[4:5],timepoints,fractioned_cell_df, offset,test_secs,EC50_list)

#############################################################################################################################################


# Plotting
message("plotting")
AC_over_gamma_T = AC_over_gamma.T


# timepoints = list(sec_histograms[gammas[0]].columns)[:-1]
# timepoints = [timepoints[idx] for idx,x in enumerate(timepoints) if idx%2 == 0]
timepoints = np.array([float(x) for x in timepoints])
# timepoints = timepoints[1:]


def histogram_plotting(avg_hist_df, gamma, time, color, label, p0=[0,0,0], exp_fit=False, hist_cutoff = -1):
    plt.plot(avg_hist_df[gamma].loc[np.abs(avg_hist_df.index - float(time)) < 0.001].iloc[0], label=label, color=color, marker="o")
    if exp_fit == True:
        data = avg_hist_df[gamma].loc[np.abs(avg_hist_df.index - float(time)) < 0.001].iloc[0][:hist_cutoff]
        data = data[~np.isnan(data)]
        def func(x,a,b,c):
            return a/x * np.exp(-x / b) + c
        myRange = np.arange(1, len(data)+1, 1)
        popt, pcov = curve_fit(func, myRange, data, maxfev=100000, p0=p0)
        plotting_range = np.arange(1, len(data)+1, 1)
        plt.plot(plotting_range-1, func(plotting_range, *popt), color="black", lw=1)
        return popt, chisquare(data, f_exp=func(myRange, *popt))
    else:
        return None


for g,gamma in enumerate(gammas[4:5]):
    fig = plt.figure(figsize=(20, 6))
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

    fig.add_subplot(a_x, a_y, 2)

    if pSTAT5 == True:
        plt.plot(timepoints/3600, AC_over_gamma_T["1_2"][gamma][:], "-o", label="AC niche")
    # plt.plot(timepoints/3600, EC50_niches_over_gamma[g], "-o", label="EC50 Threshold", color="black")
    plt.plot(timepoints/3600, niche_avg_g[g], "-o", label="Niche size by EC50", color="green")
    global_df["time"] = global_df["time"]/3600
    # global_df["niche_size"] /= cell_cell_distance
    try:
        niche_size_by_function = global_df.loc[global_df["IL-2_gamma"] == gamma].sort_values(by="time")[:-1]["niche_size"].values/cell_cell_distance
        plt.plot(timepoints/3600, niche_size_by_function[:len(timepoints)], "-o", label="Niche size by function", color="blue")
    except KeyError:
        pass
    plt.axhline(mean_sec_dist/2, color="black")

    plt.ylabel("niche size (c-c-d)")
    plt.xlabel("time (h)")
    # plt.ylim((-0.1,8.5))
    plt.title(gamma_text + ", f =" + str(frac))
    plt.legend()


    fig.add_subplot(a_x, a_y, 1)
    # timepoints = [timepoints[idx] for idx,x in enumerate(timepoints) if idx%4 == 0]
    # timepoints = [timepoints[idx] for idx in [0,2,5,7,11]]

    palette = sns.color_palette("Greys", len(timepoints))
    # palette = palette[len(timepoints):]
    # palette = ["blue","orange","green"]
    # palette = ["red"]

    perr_list = []
    avg_hist_df = pd.read_hdf(base_path + str(no_of_runs) + '_runs_' + str(frac) + saving_string + '_avg_hist_df.h5', mode="r")
    p0 = [4.48959398e+03, 6.24930170e-01, 1.21459107e+01]
    # p0 = [15,-1e5, 130]
    for t,time in enumerate(timepoints[1:]):
        #popt, error =
        try:
            popt, error = histogram_plotting(avg_hist_df,gamma,float(time),palette[t], label=str(round(float(time)/3600,1)) + "h", p0=p0, exp_fit = False, hist_cutoff = hist_cutoff+1)
        except:
            p0 = [p0[0]*10, p0[1], p0[2]]
        # print(popt)
        # p_value = error[1]
        # perr_list.append(p_value)
        # p0 = popt
    plt.axvline(mean_sec_dist/2, color="black")


    plt.title(gamma_text + ", Histograms" + " f =" + str(frac))
    plt.xticks([0,2,4,6,8,10,12])
    plt.xlabel("c-c-d")
    plt.ylabel("Concentration (pM)")
    # plt.legend()
    plt.savefig("/home/brunner/Documents/Current work/04122020/kin_analytical_histograms_f_" + str(frac) + "_g_" + str(gamma) + ".png", bbox_inches='tight')
    plt.show()

    # plt.plot(timepoints/3600,perr_list)
    # plt.title("Fitting p-value")
    # plt.xlabel("time (h)")
    # plt.ylabel("p-value")
    # plt.ylim((0,1.05))
    # # plt.savefig("/home/brunner/Documents/Current work/13112020/kin_saturation_fit_goodness_f_" + str(frac) + "_g_" + str(gamma) + ".png", bbox_inches='tight')
    # plt.savefig("/home/brunner/Documents/Current work/04122020/test" + ".png", bbox_inches='tight')
    #
    # plt.show()



# Plot a few hists


# sns.set_style("ticks")
# sns.set_context("talk", font_scale=1.5, rc={"lines.linewidth": 2.5, 'lines.markersize': 5})
# fig,ax = plt.subplots(figsize=(6,5))
#
# palette = sns.color_palette("bwr", len(np.array(gammas).flatten()))
# # neg = palette[:len(palette)//2]
# # pos = palette[len(palette)//2:]
# # palette = [neg,pos]
# if len(palette) == 2:
#     palette = ["Blue", "Red"]
#
# for g,gamma in enumerate(gammas):
#     for time in timepoints:
#         plt.plot(avg_hist_df[gamma].loc[np.abs(avg_hist_df.index - float(time)) < 0.001].iloc[0], label="test", color=palette[g], marker="o")
#
# test_mean = [avg_hist_df[0.1].loc[np.abs(avg_hist_df.index - float(1144.8792)) < 0.001].iloc[0],avg_hist_df[10.0].loc[np.abs(avg_hist_df.index - float(1144.8792)) < 0.001].iloc[0]]
# np.mean(test_mean,axis=0)
# plt.plot(np.mean(test_mean,axis=0), label="test", color="Black", marker="o")
# plt.axvline(mean_sec_dist/2, color="black")
# plt.ylabel("Concentration (pM)")
# plt.xlabel("cell-cell-distance")
#
# plt.savefig("/home/brunner/Documents/Current work/04122020/kin_stan_histograms_f_" + str(frac) + "_g_" + str(gamma) + ".png", bbox_inches='tight')
# # plt.show()

def func(x, a, b, c):
    return a / x * np.exp(-x / b) + c

r_range = np.linspace(1e-10,1e-8,100)
plt.plot(func(r_range,1e-10,1,0))
plt.xlabel("r")
plt.ylabel("$\phi$")
plt.savefig("/home/brunner/Documents/Current work/04122020/yukawa" + ".png", bbox_inches='tight')

plt.show()