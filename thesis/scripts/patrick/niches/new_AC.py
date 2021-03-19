import matplotlib.pyplot as plt
import numpy as np
import pandas as pd

from thesis.main.my_debug import message
from scipy import spatial
from scipy.signal import correlate
from scipy.optimize import curve_fit
from scipy.constants import N_A
from scipy.stats import chisquare
try:
    import fenics as fcs
except RuntimeError:
    import os
    os.environ['PATH'] = '/home/brunner/anaconda3/envs/Lukas2/bin:/home/brunner/.local/bin:/home/brunner/anaconda3/condabin:/usr/local/bin:/usr/bin:/bin:/usr/local/games:/usr/games:/opt/puppetlabs/bin'
    import fenics as fcs

import warnings
from scipy.optimize import OptimizeWarning

warnings.simplefilter("error", OptimizeWarning)

# cell_cell_distance = 20
# max_distance = 13
# pSTAT5 = True
#
# base_path = "/extra/brunner/thesis/kinetic/kin_large_saturation/"
#
# cell_df = pd.read_hdf(base_path + "cell_df.h5", mode="r")
# if pSTAT5 == True:
#     pSTAT5_sec_histograms = pd.read_hdf(base_path + "pSTAT5_hist_dfs/sec_histograms.h5", mode="r")
# sec_histograms = pd.read_hdf(base_path + "hist_dfs/sec_histograms.h5", mode="r")
#
# gammas = list(sorted(set([x[0] for x in sec_histograms.columns]))) #only unique values from columns
# gammas = [1/10.]
# timepoints = [str(y) for y in list(sorted(set([float(x[1]) for x in sec_histograms.columns])))]# cast into float to sort, then back to str
#
# timepoints = [timepoints[x] for x in [3,4,5,7,11,15,19]]
# float_timepoints = [float(x) for x in timepoints]
# sec_indices = sec_histograms.index.values


####### Version 1.1 ##########
def AC_1_1(pSTAT5_sec_histograms, gammas, timepoints, cell_cell_distance, max_distance, func):
    AC_niche_over_gamma_1_1 = [[] for gamma in gammas]
    AC_over_gamma_1_1 = [[] for gamma in gammas]
    sec_indices = pSTAT5_sec_histograms.index.values
    for g, gamma in enumerate(gammas):
        for t,time in enumerate(timepoints):
            tmp_niche = []
            tmp_auto_corr = []
            for sec in sec_indices:
                hist = [x for x in pSTAT5_sec_histograms.loc[sec, (gamma, time)] if str(x) != "nan"]
                auto_corr_1_1 = []
                for r_dash in range(len(hist) - 1):

                    A = correlate(hist, hist)
                    B = A[len(A) // 2:]
                    # tmp_auto_corr.append(hist[0]*hist[r_dash])
                    auto_corr_1_1.append(B[r_dash])
                try:
                    popt, pcov = curve_fit(func, np.arange(len(auto_corr_1_1)), auto_corr_1_1, maxfev=10000)
                    if popt[1] > 0:
                        tmp_niche.append(popt[1])
                except OptimizeWarning:
                    tmp_niche.append(np.nan)

                tmp_auto_corr.append(auto_corr_1_1)
            AC_niche_over_gamma_1_1[g].append(np.nanmean(tmp_niche))
            length = max(map(len, tmp_auto_corr))
            # make lists equal length with nan, then put into array
            tmp_auto_corr2 = np.array([xi + [np.nan] * (length - len(xi)) for xi in tmp_auto_corr])

            AC_over_gamma_1_1[g].append([np.nanmean(x) for x in tmp_auto_corr2.T])
    return AC_niche_over_gamma_1_1


####### Version 1.2 ##########
def AC_1_2(pSTAT5_sec_histograms, gammas, timepoints, cell_cell_distance, max_distance, func, plot=True, linear=False, p0 = None):
    AC_niche_over_gamma_1_2 = [[] for gamma in gammas]
    AC_over_gamma_1_2 = [[] for gamma in gammas]
    sec_indices = pSTAT5_sec_histograms.index.values
    for g, gamma in enumerate(gammas):
        for t,time in enumerate(timepoints):
            tmp_niche = []
            auto_corr_1_2 = []
            r_dash = 0
            for r_dash in range(max_distance):
                tmp_auto_corr = []
                for sec in sec_indices:
                    try:
                        hist = [x for x in pSTAT5_sec_histograms.loc[sec, (gamma, time)] if str(x) != "nan"][:-2]
                        A = correlate(hist,hist)
                        B = A[len(A) // 2:]
                        # tmp_auto_corr.append(hist[0]*hist[r_dash])
                        tmp_auto_corr.append(B[r_dash ])
                    except IndexError:
                        pass
                if tmp_auto_corr: #if not empty
                    auto_corr_1_2.append(np.sum(tmp_auto_corr))
            if linear == True:
                auto_corr_1_2 = np.log(auto_corr_1_2)
                for k_idx,k in enumerate(auto_corr_1_2):
                    if str(k) == "-inf":
                        auto_corr_1_2[k_idx] = 0
            try:
                popt = [None]
                popt, pcov = curve_fit(func, np.arange(len(auto_corr_1_2)), auto_corr_1_2, maxfev=10000, p0=p0)
                if popt[1] > 0:
                    tmp_niche.append(popt[1])
            except OptimizeWarning:
                message("Optimize Warning! Nan placed")
                tmp_niche.append(np.nan)

            AC_over_gamma_1_2[g].append(auto_corr_1_2)
            AC_niche_over_gamma_1_2[g].append(np.nanmean(tmp_niche))
            if plot == True and popt[0] != None:
                plt.plot(AC_over_gamma_1_2[g][t], label=str(gamma) + " " + str(round(float(time) / 3600, 1)) + "h")
                plt.plot(np.arange(len(auto_corr_1_2)), func(np.arange(len(auto_corr_1_2)), *popt), color="red", lw=1)
                plt.title("1_2")
                plt.xlabel("lag (c-c-d)")
                plt.ylabel("AC")
    if plot == True:
        plt.show()
    return AC_niche_over_gamma_1_2, AC_over_gamma_1_2


####### Version 2 ##########
def AC_2(sec_histograms, gammas, timepoints, cell_cell_distance, max_distance, func, plot=True):
    AC_niche_over_gamma_2 = [[] for gamma in gammas]
    AC_over_gamma_2 = [[] for gamma in gammas]
    sec_indices = sec_histograms.index.values
    for g, gamma in enumerate(gammas):
        for t,time in enumerate(timepoints):
            tmp_niche = []
            auto_corr_2 = []
            r_dash = 0
            for r_dash in range(max_distance):
                tmp_auto_corr = []
                for sec in sec_indices:
                    try:
                        hist = [x for x in sec_histograms.loc[sec, (gamma, time)] if str(x) != "nan"]
                        # A = correlate(hist,hist)
                        # B = A[len(A) // 2:]
                        tmp_auto_corr.append(hist[0]*hist[r_dash])
                        # tmp_auto_corr.append(B[r_dash ])
                    except IndexError:
                        pass
                if tmp_auto_corr: #if not empty
                    auto_corr_2.append(np.sum(tmp_auto_corr))
            try:
                popt, pcov = curve_fit(func, np.arange(len(auto_corr_2)), auto_corr_2, maxfev=10000)
                if popt[1] > 0:
                    tmp_niche.append(popt[1])
            except OptimizeWarning:
                tmp_niche.append(np.nan)

            AC_over_gamma_2[g].append(auto_corr_2)
            AC_niche_over_gamma_2[g].append(np.nanmean(tmp_niche))
            if plot == True:
                plt.plot(AC_over_gamma_2[g][t], label=str(gamma) + " " + str(round(float(time) / 3600, 1)) + "h")
                if "popt" in locals():
                    plt.plot(np.arange(len(auto_corr_2)), func(np.arange(len(auto_corr_2)), *popt), color="red", lw=1)
                plt.title("1")
                plt.xlabel("lag (c-c-d)")
                plt.ylabel("AC")
    if plot == True:
        plt.show()
    return AC_niche_over_gamma_2


####### Version 2.1 ##########
def AC_2_1(sec_histograms, gammas, timepoints, cell_cell_distance, max_distance, func):
    AC_niche_over_gamma_2_1 = [[] for gamma in gammas]
    AC_over_gamma_2_1 = [[] for gamma in gammas]
    sec_indices = sec_histograms.index.values
    for g, gamma in enumerate(gammas):
        for t,time in enumerate(timepoints):
            tmp_niche = []
            tmp_auto_corr = []
            for sec in sec_indices:
                hist = [x for x in sec_histograms.loc[sec, (gamma, time)] if str(x) != "nan"]
                auto_corr_2_1 = []
                for r_dash in range(len(hist) - 1):
                    # R = 0
                    # for r in range(len(hist)):
                    #     try:
                    #         R += hist[r + r_dash] * hist[r]
                    #     except IndexError:
                    #         pass
                    A = correlate(hist, hist)
                    B = A[len(A) // 2:]
                    # tmp_auto_corr.append(hist[0]*hist[r_dash])
                    auto_corr_2_1.append(B[r_dash])
                try:
                    popt, pcov = curve_fit(func, np.arange(len(auto_corr_2_1)), auto_corr_2_1, maxfev=10000)
                    if popt[1] > 0:
                        tmp_niche.append(popt[1])
                except OptimizeWarning:
                    tmp_niche.append(np.nan)

                tmp_auto_corr.append(auto_corr_2_1)
            AC_niche_over_gamma_2_1[g].append(np.nanmean(tmp_niche))
            length = max(map(len, tmp_auto_corr))
            # make lists equal length with nan, then put into array
            tmp_auto_corr2 = np.array([xi + [np.nan] * (length - len(xi)) for xi in tmp_auto_corr])

            AC_over_gamma_2_1[g].append([np.nanmean(x) for x in tmp_auto_corr2.T])
    return AC_niche_over_gamma_2_1



####### Version 2.2 ##########
def AC_2_2(sec_histograms, gammas, timepoints, cell_cell_distance, max_distance, func, plot=True):

    AC_niche_over_gamma_2_2 = [[] for gamma in gammas]
    AC_over_gamma_2_2 = [[] for gamma in gammas]
    sec_indices = sec_histograms.index.values
    for g, gamma in enumerate(gammas):
        for t,time in enumerate(timepoints):
            tmp_niche = []
            tmp_auto_corr = []
            for sec in sec_indices:
                hist = [x for x in sec_histograms.loc[sec, (gamma, time)] if str(x) != "nan"]
                auto_corr_2_2 = []
                for r_dash in range(len(hist) - 1):
                    # R = 0
                    # for r in range(len(hist)):
                    #     try:
                    #         R += hist[r + r_dash] * hist[r]
                    #     except IndexError:
                    #         pass
                    A = correlate(hist, hist)
                    B = A[len(A) // 2:]
                    # tmp_auto_corr.append(hist[0]*hist[r_dash])
                    auto_corr_2_2.append(B[r_dash])

                tmp_auto_corr.append(auto_corr_2_2)

            length = max(map(len, tmp_auto_corr))
            # make lists equal length with nan, then put into array
            tmp_auto_corr2 = np.array([xi + [np.nan] * (length - len(xi)) for xi in tmp_auto_corr])
            tmp_auto_corr2 = [np.nanmean(x) for x in tmp_auto_corr2.T]
            try:
                popt, pcov = curve_fit(func, np.arange(len(tmp_auto_corr2)), tmp_auto_corr2, maxfev=10000)
                if popt[1] > 0:
                    AC_niche_over_gamma_2_2[g].append(popt[1])
            except OptimizeWarning:
                AC_niche_over_gamma_2_2[g].append(np.nan)
            if not AC_niche_over_gamma_2_2[g]:
                AC_niche_over_gamma_2_2[g] = [np.nan]

            AC_over_gamma_2_2[g].append(tmp_auto_corr2)
            if plot == True:
                plt.plot(AC_over_gamma_2_2[g][t], label=str(gamma) + " " + str(round(float(time) / 3600, 1)) + "h")
                plt.plot(np.arange(len(auto_corr_2_2)), func(np.arange(len(auto_corr_2_2)), *popt), color="red", lw=1)
                plt.title("2_2")
                plt.xlabel("lag (c-c-d)")
                plt.ylabel("AC")
    if plot == True:
        plt.show()

    return AC_niche_over_gamma_2_2

def EC50(gammas,timepoints,fractioned_cell_df):
    EC50_list = [[] for gamma in gammas]
    for g, gamma in enumerate(gammas):
        # Thresholding with EC50
        mean_R_list = []
        for t, time in enumerate(timepoints):
            # for niche_timepoint in [runs_hist_df.index[-2]]:
            niche = 0
            # for i in range(len(runs_hist_df.loc[niche_timepoint, (run,gamma)][:cut_off])):
            mean_R = fractioned_cell_df.loc[
                (np.abs(fractioned_cell_df["time"] - float(time)) < 0.0001) & (
                            fractioned_cell_df["IL-2_gamma"] == gamma) & (
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
    return(EC50_list)


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

def niche_avg(gammas,timepoints,fractioned_cell_df,offset,sec_histograms,EC50_list):
    niche_avg_g = [[] for g in gammas]

    tmp_cell_df = fractioned_cell_df.loc[(np.isclose(fractioned_cell_df["time"], float(timepoints[0]))) & (
                fractioned_cell_df["IL-2_gamma"] == gammas[0])]

    x, y, z = (tmp_cell_df["x"].unique(), tmp_cell_df["y"].unique(), tmp_cell_df["z"].unique())

    offset_cells_ids = tmp_cell_df.loc[((tmp_cell_df["x"] < x[offset]) | \
                                        (tmp_cell_df["x"] > x[-offset - 1])) | \
                                       ((tmp_cell_df["y"] < y[offset]) | \
                                        (tmp_cell_df["y"] > y[-offset - 1])) | \
                                       (tmp_cell_df["z"] < z[offset]) | \
                                       (tmp_cell_df["z"] > z[-offset - 1]), "id"].unique()

    for g, gamma in enumerate(gammas):
        for t, time in enumerate(timepoints):
            sec_niches = []
            sec_cells_ids = sec_histograms[gamma][str(time)].index
            sec_cells_ids = [x for x in sec_cells_ids if x not in offset_cells_ids]
            for i in sec_cells_ids:
                for j in range(len(sec_histograms[gamma][str(time)].loc[i]) + 1):
                    try:
                        # if runs_hist_df.loc[niche_timepoint, (run, gamma)][i] - runs_hist_df.loc[niche_timepoint, (run, gamma)][i + 1] < slope_threshold * runs_hist_df.loc[niche_timepoint, (run, gamma)][i + 1]:
                        if sec_histograms[gamma][str(time)].loc[i][j] < EC50_list[g][t]:
                            sec_niches.append(j)
                            break
                    except IndexError:
                        sec_niches.append(j)
            niche_avg_g[g].append(np.mean(sec_niches))
    return niche_avg_g

def niche_effect_func(gammas,time,frac,fractioned_cell_df,offset,sec_histograms,EC50_list, base_path, ext_cache):
    niche_avg_g = [[] for g in gammas]
    niche_effect_g = [[] for g in gammas]

    tmp_cell_df = fractioned_cell_df.loc[(np.isclose(fractioned_cell_df["time"], float(time))) & (
                fractioned_cell_df["IL-2_gamma"] == gammas[0])]

    x, y, z = (tmp_cell_df["x"].unique(), tmp_cell_df["y"].unique(), tmp_cell_df["z"].unique())

    offset_cells_ids = tmp_cell_df.loc[((tmp_cell_df["x"] < x[offset]) | \
                                        (tmp_cell_df["x"] > x[-offset - 1])) | \
                                       ((tmp_cell_df["y"] < y[offset]) | \
                                        (tmp_cell_df["y"] > y[-offset - 1])) | \
                                       (tmp_cell_df["z"] < z[offset]) | \
                                       (tmp_cell_df["z"] > z[-offset - 1]), "id"].unique()

    for g, gamma in enumerate(gammas):
        t = 0
        sec_niches = []
        concentrations = []
        sec_cells_ids = sec_histograms[gamma][str(time)].index
        sec_cells_ids = [x for x in sec_cells_ids if x not in offset_cells_ids]
        for i in sec_cells_ids:
            for j in range(len(sec_histograms[gamma][str(time)].loc[i]) + 1):
                try:
                    # if runs_hist_df.loc[niche_timepoint, (run, gamma)][i] - runs_hist_df.loc[niche_timepoint, (run, gamma)][i + 1] < slope_threshold * runs_hist_df.loc[niche_timepoint, (run, gamma)][i + 1]:
                    if sec_histograms[gamma][str(time)].loc[i][j] < EC50_list[g][t]:
                        sec_niches.append(j)
                        break
                except IndexError:
                    sec_niches.append(j)
            if sec_niches[-1] != 0:
                for k in range(sec_niches[-1]):
                    concentrations.append(sec_histograms[gamma][str(time)].loc[i][k])
            else:
                concentrations.append(sec_histograms[gamma][str(time)].loc[i][0])
        niche_avg_g[g].append(np.mean(sec_niches))
        # calculate niche effect: take the timepoints niche and get the average concentration within that niche

        # average over runs

        try:
            global_df = pd.read_hdf(base_path + "/global_df.h5", mode="r")
            average_concentration = global_df.loc[(np.abs(global_df["time"] - float(time)) < 0.0001) & (
                    global_df["IL-2_gamma"] == gamma), "surf_c"].values[0] * 1e3
        except:
            print("Error!")

        try:
            mesh_c_sum = global_df.loc[(np.abs(global_df["time"] - float(time)) < 0.0001) &
                                       (global_df["IL-2_gamma"] == gamma), "Concentration"].values * \
                         global_df["no_of_vertices"][0]
            niche_effect_g[g].append(np.sum(concentrations) / mesh_c_sum)
        except KeyError:
            mesh = fcs.Mesh()
            with fcs.XDMFFile(fcs.MPI.comm_world, ext_cache + "mesh.xdmf") as f:
                f.read(mesh)
            mesh_c_sum = global_df.loc[
                             (np.abs(global_df["time"] - float(time)) < 0.0001) & (global_df["IL-2_gamma"] == gamma) & (
                                         global_df["IL-2_Tsec_fraction"] == frac), "Concentration"].values * 1e3 * len(
                mesh.coordinates())
            niche_effect_g[g].append(np.sum(concentrations) / mesh_c_sum)
    return niche_avg_g, niche_effect_g