import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import seaborn as sns
import random
try:
    import fenics as fcs
except RuntimeError:
    import os
    os.environ['PATH'] = '/home/brunner/anaconda3/envs/Lukas2/bin:/home/brunner/.local/bin:/home/brunner/anaconda3/condabin:/usr/local/bin:/usr/bin:/bin:/usr/local/games:/usr/games:/opt/puppetlabs/bin'
    import fenics as fcs

from pandas import DataFrame

import sys
sys.path.insert(0,'/home/brunner/thesis/thesis/main')
from my_debug import message
from scipy import spatial
from scipy.optimize import curve_fit
from scipy.constants import N_A
import matplotlib.lines as mlines


try:
    from thesis.scripts.patrick.niches.new_AC import AC_1_2, EC50, histogram_plotting, niche_avg, niche_effect_func
except ModuleNotFoundError:
    exec(open('/home/brunner/thesis/thesis/scripts/patrick/niches/new_AC.py').read())


#load the df

# base_path = r"/home/brunner/thesis/thesis/scripts/patrick/yukawa_analytical/test/"
prefix = "" #"/run/user/6205/gvfs/sftp:host=compute3"
base_path = prefix + "/extra2/brunner/thesis/kinetic/standard/8_big_scan/"
ext_cache = prefix + "/extra2/brunner/thesis/kinetic/standard/3D_ext_cache/"
global_df = pd.read_hdf(base_path + "global_df.h5", mode="r")
cell_df = pd.read_hdf(base_path + "cell_df.h5", mode="r")


verbose = False
save_avg_hist_df = True
no_of_runs = 1
dt = 3600
cell_cell_distance = 20
max_distance = 13
pSTAT5 = False


frac = 0.25
gammas = [x for x in sorted(list(global_df["IL-2_gamma"].unique()))]
offset = 2
niche_timepoint = -2
hist_cutoff = 11

# base_path = "/extra/brunner/thesis/kinetic/q_fraction_medium_pos_g_scan_multi/"
if pSTAT5 == True:
    runs_hist_df = pd.read_hdf(base_path + str(no_of_runs) + '_runs_' + str(frac) + '_Tsec_fraction_pSTAT5_hist_df.h5', mode="r")
    # avg_hist_df = pd.read_hdf(base_path + str(no_of_runs) + '_runs_' + str(frac) + '_Tsec_fraction_pSTAT5_avg_hist_df.h5', mode="r")
else:
    runs_hist_df = pd.read_hdf(base_path + str(no_of_runs) + '_runs_' + str(frac) + '_Tsec_fraction_hist_df.h5', mode="r")
    # avg_hist_df = pd.read_hdf(base_path + str(no_of_runs) + '_runs_' + str(frac) + '_Tsec_fraction_avg_hist_df.h5', mode="r")
mean_sec_dist = float(np.loadtxt(base_path + str(no_of_runs) + "_runs_" + str(frac) + "_Tsec_fraction_avg_sec_dists.txt"))/cell_cell_distance

fractioned_cell_df = cell_df.loc[(np.abs(cell_df["IL-2_Tsec_fraction"] - frac) < 0.0001)]

if pSTAT5 == True:
    pSTAT5_sec_histograms = pd.read_hdf(base_path + "pSTAT5_hist_dfs/sec_pSTAT5_" + str(frac) + "_Tsec_fraction_histograms.h5", mode="r")
sec_histograms = pd.read_hdf(base_path + "hist_dfs/sec_" + str(frac) + "_Tsec_fraction_histograms.h5", mode="r")

message("loaded runs_hist_df from " + base_path)


timepoints = [list(sec_histograms[gammas[0]].columns)[1:-1][niche_timepoint]]
time = timepoints[0]
AC_over_gamma: DataFrame = pd.DataFrame(columns=gammas, index=["1_2"])
AC_1_2_autocorr_g = [[] for x in gammas]
unavg_EC50_niche_over_sv = [[] for x in gammas]

#############################################################################################################################################

message("calculating EC50")
EC50_list = EC50(gammas,timepoints,fractioned_cell_df)

# calculate AC niches:
def func(x, a, b, c):
    return a * np.exp(-x / b) + c
def func2(x, a, b):
    return a - x / b

if pSTAT5 == True:
    # AC_over_gamma[gamma]["1_1"] = AC_1_1(pSTAT5_sec_histograms, gammas, timepoints, cell_cell_distance, max_distance, func)[0][0]
    AC_1_2_niche, AC_1_2_autocorr = AC_1_2(pSTAT5_sec_histograms, gammas, [time], cell_cell_distance, max_distance, func2, plot=False, linear=True, p0=None)
    for g, gamma in enumerate(gammas):
        AC_over_gamma[gamma]["1_2"] = AC_1_2_niche[g]
        AC_1_2_autocorr_g[g].append(AC_1_2_autocorr[g])
        # AC_over_gamma[gamma]["2"] = AC_2(sec_histograms, gammas, timepoints, cell_cell_distance, max_distance, func, plot=False)[0][0]
        # AC_over_gamma[gamma]["2_1"] = AC_2_1(sec_histograms, gammas, timepoints, cell_cell_distance, max_distance, func)[0][0]
        # AC_over_gamma[gamma]["2_2"] = AC_2_2(sec_histograms, gammas, timepoints, cell_cell_distance, max_distance, func, plot=False)[0][0]

message("calculating niche_avg_g")
niche_avg_g, niche_effect_g = niche_effect_func(gammas,time,frac,fractioned_cell_df,offset,sec_histograms,EC50_list, base_path, ext_cache)

#############################################################################################################################################

# load average secretory distance and norm to cell-cell-distance

# niche_timepoints_index = 15
# niche_timepoint = runs_hist_df.index[niche_timepoints_index]
# cut_off = None

#############################################################################################################################################
# run '/home/brunner/thesis/thesis/scripts/patrick/niches/plot_ne_gammas.py'

# Plotting
sns.set_style("ticks")
sns.set_context("talk", font_scale=1.5, rc={"lines.linewidth": 2.5, 'lines.markersize': 5})
fig,ax = plt.subplots(figsize=(6,5))

ylim = (None,None)
# fig = plt.figure(figsize=(12,12))
# sns.set_style("ticks")
# sns.set_context("talk", font_scale=1.3, rc={"lines.linewidth": 4, 'lines.markersize': 6})

# plt.subplots_adjust(wspace=.4, hspace=0.6)


pos_palette = sns.color_palette("Reds", len(gammas))
neg_palette = sns.color_palette("Blues", len(gammas))
new_handles = []
# palette = [i for i in reversed(palette)]
palette = sns.color_palette("bwr", len(gammas))

for g,gamma in enumerate(gammas):
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

    x = niche_effect_g[g]
    y = niche_avg_g[g]
    # if gamma > 1:
    #     halfed_palette = [i for i in reversed(pos_palette)] #[int((0.5*g-1)*len(palette))-1*g:int((1-g)*0.5*len(palette))-1*g]
    # else:
    #     halfed_palette = neg_palette
    # if g == 0:
    #     halfed_palette = [i for i in reversed(halfed_palette)]

    for i in range(len(x)):
        # if gammas[i] > 1:
        #     y[i] += np.random.normal()/30
        # else:
        #     myColor = "black"
        plt.scatter(x[i]*100, y[i], color=palette[g])
        # plt.errorbar(x[i], y[i], xerr=[[np.abs(x)] for x in niche_effect_error[i]], yerr=[[np.abs(x)] for x in niche_size_error[i]], linestyle='', c=palette[i])
    # plt.scatter(216.35519872163854, 1.6372549019607845, color="black")
    # plt.errorbar(216.35519872163854, 1.6372549019607845,xerr=13.160380344055845, yerr=0.22201019379779285, color="black")
    # sns.set_context(rc={"lines.linewidth": 2})
    # norm = plt.Normalize(np.min([float(a) for a in timepoints])/3600, np.max([float(a) for a in timepoints])/3600)
    # sm = plt.cm.ScalarMappable(cmap="bwr", norm=norm)
    # sm.set_array([])

    plt.xlabel("niche effect ($\%$)")
    plt.ylabel("niche size")

    # plt.title(str(round(slope_threshold*100,2)) + "%")
    # plt.ylim(ylim)

new_handles.append(mlines.Line2D([], [], color=palette[0], label="positive feedback"))
new_handles.append(mlines.Line2D([], [], color=palette[-1], label="negative feedback"))
new_handles.append(mlines.Line2D([], [], color='black', label='Black line'))
new_labels = ["neg. feedback", "pos. feedback", "no feedback"]
# plt.legend(handles=new_handles, labels=new_labels)
# plt.title("pSTAT5 = " + str(pSTAT5))
# if len(gammas) == 1:
#     plt.axes().get_legend().remove()
#     plt.axes().figure.colorbar(sm)
# fig.savefig(r"/home/brunner/Documents/Current work/04122020/" + "ne_gamma_f_" + str(frac) + ".png", bbox_inches='tight')
plt.show()
