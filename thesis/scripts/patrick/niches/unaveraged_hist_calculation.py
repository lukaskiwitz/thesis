import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import os
import seaborn as sns
import random
from thesis.main.my_debug import message
from scipy import spatial


cell_cell_distance = 30
base_path = "/extra/brunner/thesis/static/static_c_c_d_" + str(cell_cell_distance) + "/"

pSTAT5 = False

cell_df = pd.read_hdf(base_path + "cell_df.h5", mode="r")


# a = [1/x for x in reversed(np.arange(10,110,10))]
# b = [x for x in np.arange(10,110,10)]
# gammas = np.concatenate([a,b])

gammas = [1]

# timepoints = cell_df["time"].unique()
timepoints = [0.,1.,2.]
hist_timepoints = [str(x) for x in timepoints]


array1 = np.array([[g]*len(timepoints) for g in gammas]).flatten()
array2 = np.array([[t for t in hist_timepoints]*len(gammas)]).flatten()
df_columns = [array1,array2]

# fractions = [0.1,0.15,0.2]
# fractions = np.arange(0.01,0.25,0.02)
# fractions = [0.25]

scan_variables = np.arange(0.01,0.25,0.02)
df_string = "IL-2_Tsec_fraction"
saving_string = "_Tsec_fraction"
for pSTAT5 in [True,False]:
    for scan_v in scan_variables:
        if len(str(scan_v)) > 5:
            scan_v = round(scan_v, 2)
        Tsec_ids = cell_df.loc[(cell_df["time"] == timepoints[0]) &
                                (cell_df["type_name"] == "Tsec") &
                                (np.abs(cell_df[df_string] - scan_v) < 0.0001) &
                                (cell_df["IL-2_gamma"] == gammas[0]), "id"].unique()
        sec_histograms = pd.DataFrame(columns=df_columns, index=Tsec_ids)
        for g, gamma in enumerate(gammas):
            message("calculating hists with " + df_string + " " + str(scan_v) + " for gamma " + str(gamma) +" with pSTAT5 = " + str(pSTAT5))
            if pSTAT5 == True:
                hist_df = pd.read_hdf(base_path + "pSTAT5_hist_dfs/hist_df_" + str(round(gamma, 4)) + "_" + str(scan_v) + saving_string + ".h5", mode="r")
            else:
                hist_df = pd.read_hdf(base_path + "hist_dfs/hist_df_" + str(round(gamma,4)) + "_" + str(scan_v) + saving_string +  ".h5", mode="r")
            max_distance = np.max(hist_df["0.0"][0])
            # Histogram building

            # iterate over timepoints:
            for t_idx, hist_timepoint in enumerate(hist_timepoints):

                # calculate necessary bins
                bins = np.arange(0, (max_distance // cell_cell_distance + 1) * cell_cell_distance, cell_cell_distance)
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

                    c_bin_array = [[] for i in range((len(bins)-1))]

                    # since "c" and "distances" are sorted, we can use np.where's first and last argument as a range for the entries of the bin
                    # this way we only iterate over the bins, the iteration over the distances is done by numpy.where
                    for i, bin in enumerate(bins[1:]):
                        if i == 0:
                            current_range = np.where(distances < bins[1:][i])[0]
                        else:
                            current_range = np.where((distances < bins[1:][i]) & (distances > bins[1:][i - 1]))[0]

                        if len(current_range) == 0:
                            pass
                        elif len(current_range) == 1:
                            c_bin_array[i].extend([c[current_range[0]]])
                        else:
                            c_bin_array[i].extend(c[current_range[0]:current_range[-1]])

                    sec_histograms.loc[k, (gamma, hist_timepoint)] = [np.mean(k) for k in c_bin_array if k]
        try:
            if pSTAT5 == True:
                os.mkdir(base_path + "pSTAT5_hist_dfs")
            else:
                os.mkdir(base_path + "hist_dfs")
        except FileExistsError:
            pass

        if pSTAT5 == True:
            sec_histograms.to_hdf(base_path + "pSTAT5_hist_dfs/sec_pSTAT5" + "_" + str(scan_v) + saving_string +  "_histograms.h5", key="data", mode="w")
        else:
            sec_histograms.to_hdf(base_path + "hist_dfs/sec" + "_" + str(scan_v) + saving_string + "_histograms.h5", key="data", mode="w")