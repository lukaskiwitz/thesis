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

verbose = True

no_of_runs = 10
timepoints = [0,2,5,10,20]
# gammas = [1/9, 1/8, 1/7, 1/6, 1/5, 1/4, 1/3]
gammas = [3,4,5,6,7,8,9]
cell_cell_distance = 20
get_dataframes = []

#create multiindex df
array1 = np.array([[x]*len(gammas) for x in range(no_of_runs)]).flatten()
array2 = np.array([[g for g in gammas]*no_of_runs]).flatten()
df_columns = [array1,array2]

df_indices = timepoints.copy()
df_indices.append("sec_dists")

runs_hist_df = pd.DataFrame(columns=df_columns, index=df_indices)
base_path = "/extra/brunner/thesis/kinetic/q_fraction_medium_pos_g_scan_multi/"
# base_path = "/extra/brunner/thesis/kinetic/q_fraction_medium_neg_g_scan_multi/"
if load_runs_hist_df == True:
    try:
        runs_hist_df = pd.read_hdf(base_path + str(no_of_runs) + "_runs_hist_df.h5", mode="r")
        message("loaded hist_df from " + base_path)
        df_loaded = True
    except FileNotFoundError:
        message("hist_df not found, commencing calculation")
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
            # cell_df = cell_df.loc[cell_df["IL-2_fraction"] == 0.25]
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
        message("saving hist_df: " +  base_path + str(no_of_runs) + '_runs_hist_df.h5')
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
    temp_sec_dists.append(transposed_runs_hist_df["sec_dists"][i][gammas[0]])
temp_sec_dists = [item for sublist in temp_sec_dists for item in sublist]
sec_dists = np.mean(temp_sec_dists)
# Write averaged dfs to file
if save_avg_hist_df == True:
    message("saving avg_hist_df: " + base_path + str(no_of_runs) + '_runs_avg_hist_df.h5')
    np.savetxt(base_path + str(no_of_runs) + "_runs_avg_sec_dist.txt", [sec_dists])
    avg_hist_df.to_hdf(base_path + str(no_of_runs) + '_runs_avg_hist_df.h5', key="data", mode="w")


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