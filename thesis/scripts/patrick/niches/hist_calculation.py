import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import os
import seaborn as sns
import random
from thesis.main.my_debug import message
from scipy import spatial

save_runs_hist_df = True
save_avg_hist_df = True

verbose = False

no_of_runs = 1

dt = 3600 #1h
length = 50
max_T = dt * 400

myRange = np.arange(0,length)
def exp_func(x,a,b,c):
    return a*np.exp(b*x) + c
a = 2*dt
c = -a
b = np.log((max_T-c)/a)/(length-1)

timepoints = [0,1,2] #list(exp_func(myRange,a,b,c))

# a = [1/x for x in reversed(np.arange(10,110,10))]
# b = [x for x in np.arange(10,110,10)]
# gammas = np.concatenate([a,b])
gammas = [1]



#create columns and indices for df
array1 = np.array([[x]*len(gammas) for x in range(no_of_runs)]).flatten()
array2 = np.array([[g for g in gammas]*no_of_runs]).flatten()
df_columns = [array1,array2]
df_indices = timepoints.copy()
df_indices.append("sec_dists")

array1 = np.array([[g]*len(timepoints) for g in gammas]).flatten()
array2 = np.array([[t for t in [str(float(x)) for x in timepoints]]*len(gammas)]).flatten()
sec_df_columns = [array1,array2]

cell_cell_distance = 30
base_path = "/extra/brunner/thesis/static/static_c_c_d_" + str(cell_cell_distance) + "/"
if no_of_runs > 1:
    path = base_path + "run" + str(1) + "/"
else:
    path = base_path
cell_df = pd.read_hdf(path + 'cell_df.h5', mode="r")
full_cell_df = cell_df.copy()
#check if correct times are used
for entry in df_indices[:-1]:
    assert any(np.isclose(entry,cell_df.sort_values(by="time")["time"].unique())) == True

# fractions = np.arange(0.01,0.25,0.02)
scan_variables = np.arange(0.01,0.25,0.02)
df_string = "IL-2_Tsec_fraction"
saving_string = "_Tsec_fraction"
for pSTAT5 in [True, False]:
    for scan_v in scan_variables:
        if len(str(scan_v)) > 5:
            scan_v = round(scan_v,2)
        runs_hist_df = pd.DataFrame(columns=df_columns, index=df_indices)
        for run in range(no_of_runs):
            message("calculating run " + str(run) + " with pSTAT5 = " + str(pSTAT5))
            Tsec_ids = full_cell_df.loc[(full_cell_df["time"] == timepoints[0]) &
                                   (full_cell_df["type_name"] == "Tsec") &
                                   (np.abs(full_cell_df[df_string] - scan_v) < 0.01) &
                                   (full_cell_df["IL-2_gamma"] == gammas[0]), "id"].unique()

            sec_histograms = pd.DataFrame(columns=sec_df_columns, index=Tsec_ids)
            for gamma in gammas:
                message("calculating gamma = " + str(gamma))
                if no_of_runs > 1:
                    path = base_path + "run" + str(run) + "/"
                else:
                    path = base_path

                global_df = pd.read_hdf(path + 'global_df.h5', mode="r")
                cell_df = pd.read_hdf(path + 'cell_df.h5', mode="r")

                variable = "time"
                cell_df = cell_df.loc[np.abs(cell_df[df_string] - scan_v) < 0.0001]
                cell_df = cell_df.loc[cell_df["IL-2_gamma"] == gamma]

                #for which timepoints to calculate for
                # points = [cell_df[variable].unique()[i] for i in timepoints]
                # points = timepoints
                if variable == "time":
                    points = cell_df.sort_values(by="time")["time"].unique()
                else:
                    print("please choose points to iterate over")
                    exit()

                #create multiindex dfs
                hist_df = pd.DataFrame(columns=[str(x) for x in points], index=cell_df.loc[(cell_df[variable] == points[0]) , "id"].unique())

                # prepare df's to build kd_tree
                cell_df_split = cell_df.loc[cell_df[variable] == points[0]].copy()
                cell_points = np.array([cell_df_split["x"].values, cell_df_split["y"].values, cell_df_split["z"].values]).T
                kd_tree = spatial.KDTree(cell_points)

                max_distance = 0

                #get the distances from the KDTree at the given time points
                for entry in points: #cell_df["time"].unique():
                    if verbose == True:
                        message("calculating KDTree at " + variable + " = " + str(entry/dt) + "h")
                    # create containers
                    cell_df_split = cell_df.loc[cell_df[variable] == entry].copy()
                    # iterate over cells
                    for il2cells in cell_df.loc[(cell_df[variable] == points[0]), "id"]:

                        kd_tree_data = kd_tree.query(cell_points[int(il2cells)], len(cell_points))
                        hist_df[str(entry)][str(il2cells)] = [np.zeros(len(cell_points)), np.zeros(len(cell_points))]

                        #write concentrations/pSTAT5
                        if pSTAT5 == True:
                            hist_df[str(entry)][str(il2cells)][0][kd_tree_data[1]] = cell_df_split.iloc[kd_tree_data[1]]["pSTAT5_pSTAT5"].values
                        else:
                            hist_df[str(entry)][str(il2cells)][0][kd_tree_data[1]] = cell_df_split.iloc[kd_tree_data[1]]["IL-2_surf_c"].mul(1e3).values
                        #write distances
                        hist_df[str(entry)][str(il2cells)][1][kd_tree_data[1]] = kd_tree_data[0]
                        #get the max distance between two cells while we are at it
                        if max_distance < kd_tree_data[0].max():
                            max_distance = kd_tree_data[0].max()
                try:
                    if pSTAT5 == True:
                        os.mkdir(base_path + "pSTAT5_hist_dfs")
                    else:
                        os.mkdir(base_path + "hist_dfs")
                except FileExistsError:
                    pass
                if pSTAT5 == True:
                    save_name = base_path + "pSTAT5_hist_dfs/hist_df_" + str(round(gamma, 4)) + "_" + str(scan_v) + saving_string + ".h5"
                else:
                    save_name = base_path + "hist_dfs/hist_df_" + str(round(gamma,4)) + "_" + str(scan_v) + saving_string + ".h5"
                hist_df.to_hdf(save_name, key="data", mode="w")
                # plt.plot(hist_df["2.0"][0][0])
                # plt.show()

                ############################################################################################################################################

                # Histogram building
                hist_timepoints = [str(x) for x in points]
                #sec histograms

                message("calculating hists with " + df_string + " " + str(scan_v) + " for gamma " + str(
                    gamma) + " with pSTAT5 = " + str(pSTAT5))

                # iterate over timepoints:
                for t_idx, hist_timepoint in enumerate(hist_timepoints):

                    # calculate necessary bins
                    bins = np.arange(0, (max_distance // cell_cell_distance + 1) * cell_cell_distance,
                                     cell_cell_distance)
                    center = bins[:-1] / cell_cell_distance

                    # iterate over all secreting cells in hist_df to average over them
                    # grab the secreting cells average distance from each other while we are at it
                    sec_dists = []
                    cell_tracker = 0
                    no_of_cells = len(hist_df[hist_timepoint][Tsec_ids])

                    for k in hist_df[hist_timepoint][Tsec_ids].index:  # .loc[hist_df[plot_time_point].index == str(2)].dropna():

                        cell = hist_df[hist_timepoint][Tsec_ids][k]
                        # sort the arrays and corresponding concentrations by distance via zipping und unzipping
                        cell_tracker += 1

                        zipped = zip(cell[0], cell[1])
                        sorted_zipped = sorted(list(zipped), key=lambda x: x[1])

                        c, distances = zip(*sorted_zipped)
                        c = list(c)
                        distances = list(distances)

                        c_bin_array = [[] for i in range((len(bins) - 1))]

                        # since "c" and "distances" are sorted, we can use np.where's first and last argument as a range for the entries of the bin
                        # this way we only iterate over the bins, the iteration over the distances is done by numpy.where
                        for i, bin in enumerate(bins[1:]):
                            if i == 0:
                                current_range = np.where(distances < bins[1:][i])[0]
                            else:
                                current_range = np.where((distances < bins[1:][i]) & (distances > bins[1:][i - 1]))[
                                    0]

                            if len(current_range) == 0:
                                pass
                            elif len(current_range) == 1:
                                c_bin_array[i].extend([c[current_range[0]]])
                            else:
                                c_bin_array[i].extend(c[current_range[0]:current_range[-1]])

                        sec_histograms.loc[k, (gamma, hist_timepoint)] = [np.mean(k) for k in c_bin_array if k]

                    # Full hist calculation
                    assert float(hist_timepoint) in points
                    if verbose == True:
                        message("creating hist at time = " + str(float(hist_timepoint)/dt) + "h")

                    # calculate necessary bins
                    bins = np.arange(0,(max_distance//cell_cell_distance + 1)*cell_cell_distance, cell_cell_distance)
                    center = bins[:-1]/cell_cell_distance

                    # create bins container
                    c_bin_array = [[] for i in range((len(bins)-1))]

                    # iterate over all secreting cells in hist_df to average over them
                    # grab the secreting cells average distance from each other while we are at it
                    sec_dists = []

                    # no_of_cells = len(hist_df[hist_timepoint][cell_df.loc[(cell_df["type_name"] == "Tsec") & (cell_df["time"] == float(hist_timepoint)), "id"]])
                    for k,cell in enumerate(hist_df[hist_timepoint][cell_df.loc[(cell_df["type_name"] == "Tsec") & (cell_df["time"] == float(hist_timepoint)), "id"]]):  #.loc[hist_df[plot_time_point].index == str(2)].dropna():
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
                        ids = [int(x) for x in cell_df.loc[(cell_df["type_name"] == "Tsec") & (cell_df["time_index"] == 0), "id"].values]
                        for i in cell[1][ids]:
                            if i != 0.0:
                                sec_dists.append(i)

                    # calc mean of bins
                    hist = [np.mean(k) for k in c_bin_array]
                    # write histogram into df
                    runs_hist_df.loc[:, (run, gamma)][t_idx] = hist
                    runs_hist_df.loc[:, (run, gamma)]["sec_dists"] = sec_dists

            try:
                if pSTAT5 == True:
                    os.mkdir(base_path + "pSTAT5_hist_dfs")
                else:
                    os.mkdir(base_path + "hist_dfs")
            except FileExistsError:
                pass

            if pSTAT5 == True:
                sec_histograms.to_hdf(
                    base_path + "pSTAT5_hist_dfs/sec_pSTAT5" + "_" + str(scan_v) + saving_string + "_histograms.h5",
                    key="data", mode="w")
            else:
                sec_histograms.to_hdf(
                    base_path + "hist_dfs/sec" + "_" + str(scan_v) + saving_string + "_histograms.h5", key="data",
                    mode="w")

        # write runs_hist_df to file
        if save_runs_hist_df == True:
            if pSTAT5 == True:
                save_name = base_path + str(no_of_runs) + '_runs_' + str(scan_v) + saving_string + '_pSTAT5_hist_df.h5'
            else:
                save_name = base_path + str(no_of_runs) + '_runs_' + str(scan_v) + saving_string + '_hist_df.h5'
            message("saving hist_df: " +  save_name)
            runs_hist_df.to_hdf(save_name, key="data", mode="w")
        #############################################################################################################################################

        # Avergaging runs. columns = gamma values, indices = timepoints
        # Transpose df for easier handling
        if pSTAT5 == True:
            load_name = base_path + str(no_of_runs) + '_runs_' + str(scan_v) + saving_string + '_pSTAT5_hist_df.h5'
        else:
            load_name = base_path + str(no_of_runs) + '_runs_' + str(scan_v) + saving_string + '_hist_df.h5'
        runs_hist_df = pd.read_hdf(load_name, mode="r")

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

        # plt.plot(avg_hist_df[1][2])
        # plt.show()
        # Averaging secretory cells distances
        temp_sec_dists = []
        for i in range(transposed_runs_hist_df["sec_dists"].index[-1][0]+1):
            temp_sec_dists.append(transposed_runs_hist_df["sec_dists"][i][gammas[0]])
        temp_sec_dists = [item for sublist in temp_sec_dists for item in sublist]
        sec_dists = np.mean(temp_sec_dists)
        # Write averaged dfs to file
        if save_avg_hist_df == True:
            np.savetxt(base_path + str(no_of_runs) + '_runs_' + str(scan_v) + saving_string + '_avg_sec_dists.txt', [sec_dists])
            if pSTAT5 == True:
                save_name = base_path + str(no_of_runs) + '_runs_' + str(scan_v) + saving_string + '_pSTAT5_avg_hist_df.h5'
            else:
                save_name = base_path + str(no_of_runs) + '_runs_' + str(scan_v) + saving_string + '_avg_hist_df.h5'
            message("saving avg_hist_df: " +  save_name)
            avg_hist_df.to_hdf(save_name, key="data", mode="w")
    ##########################################################################################################################