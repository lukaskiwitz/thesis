import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import seaborn as sns
import random
from my_debug import message

# def get_id_loc(cell_df, id):
#     return np.array([cell_df.loc[cell_df["id"]==str(id), "x"].unique()[0], cell_df.loc[cell_df["id"]==str(id), "y"].unique()[0], cell_df.loc[cell_df["id"]==str(id), "z"].unique()[0]])
# get the distance between two cells. Written in a way so it can be applied to the cell_df_timed dataframe. Speeds up computation significantly

#fast applyable way to get the distances, packed in a function
def get_distances(origin_id, target_id):
    dist = np.linalg.norm(
        np.array(
            [cell_df_timed.loc[cell_df_timed["id"] == str(target_id), "x"].values,
             cell_df_timed.loc[cell_df_timed["id"] == str(target_id), "y"].values,
             cell_df_timed.loc[cell_df_timed["id"] == str(target_id), "z"].values
             ])
        -
        np.array(
            [cell_df_timed.loc[cell_df_timed["id"] == str(origin_id), "x"].values,
             cell_df_timed.loc[cell_df_timed["id"] == str(origin_id), "y"].values,
             cell_df_timed.loc[cell_df_timed["id"] == str(origin_id), "z"].values
             ])
    )
    cell_df_timed.loc[(cell_df_timed["id"] == str(target_id)), "distance"] = dist

#load the cell_df
# from parameters_q_fraction import path_kinetic
# path = path_kinetic
# path = "/extra/brunner/thesis/kinetic/q_fraction/"
# path = "/extra/brunner/thesis/kinetic/q_fraction_wave/"
path = "/extra/brunner/thesis/kinetic/q_fraction_wave_gamma_scan/"

cell_df = pd.read_hdf(path + 'cell_df.h5', mode="r")
global_df = pd.read_hdf(path + 'global_df.h5', mode="r")

cell_df = cell_df.loc[cell_df["IL-2_gamma"] == 0.1]

#for which timepoints to calculate for
time_points = [0., 5., 50.] #cell_df["time"].unique()
threshold = 1e10 #15000 #if less than x Receptors

try:
    hist_df = pd.read_hdf(path + "hist_df_c_" + str(threshold) + ".h5", mode="r")
    normalising_df = pd.read_hdf(path + "normalising_df_c_" + ".h5", mode="r")
    print("Loaded in previously calculated files")
except FileNotFoundError:
    #create multiindex dfs
    hist_df = pd.DataFrame(columns=[str(x) for x in time_points], index=cell_df.loc[(cell_df["type_name"] == "changed") & (cell_df["time"] == 0) , "id"].values)
    normalising_df = pd.DataFrame(columns=[str(x) for x in time_points], index=cell_df.loc[(cell_df["type_name"] == "changed") & (cell_df["time"] == 0), "id"].values)

    # small_cell_sample = random.sample(cell_df.loc[(cell_df["type_name"] == "changed") & (cell_df["time"] == time_points[0]), "id"].values.tolist(), 1)
    # for desired time_points
    for time in time_points: #cell_df["time"].unique():
        message("time = " + str(time))
        # for desired cells/all the cells
        for il2cells in cell_df.loc[(cell_df["type_name"] == "changed") & (cell_df["time"] == 0), "id"]: #small_cell_sample: #
            # create containers
            cell_df["distance"] = None
            cell_df_timed = cell_df.loc[cell_df["time"] == time].copy()
            cell_df_timed["distance"] = None
            #calculate distances by applying
            cell_df_timed.apply(lambda row: get_distances(il2cells, row["id"]), axis=1)
            data = []
            distance = []
            # threshold check
            for index,row in cell_df_timed.iterrows():
                # if row["IL-2_surf_c"] > threshold:
                if row["IL-2_R"] < threshold:
                    data.append(row["IL-2_surf_c"]*1e3)
                    distance.append(row["distance"])
            # put into df. [0] = data = surf_c, [1] = distancewave_c01.py
            hist_df[str(time)][str(il2cells)] = [np.array(data).astype(np.float).copy(), np.array(distance).astype(np.float).copy()] #cell_df.loc[cell_df["time"] == time, "distance"]
            normalising_df[str(time)][str(il2cells)] = [cell_df_timed["IL-2_surf_c"].mul(1e3).to_numpy().astype(np.float).copy(), cell_df_timed["distance"].copy()]
    # save to file
    hist_df.to_hdf(path + "hist_df_c_01_large_" + str(threshold) + ".h5",key="df", mode="w")
    normalising_df.to_hdf(path + "normalising_df_c_" + ".h5", key="df", mode="w")


# Morans I implementation. w_i_j = 1 with w_i_i = 0, so N/W = 1

#get intervals
cell_cell_distance = 20
max_distance = np.array([x for x in hist_df["0.0"][0][1]]).flatten().max()
bins = np.arange(cell_cell_distance,(max_distance//cell_cell_distance)*20, cell_cell_distance)

hist_df_0 = hist_df["0.0"]
#compute average
for r in [40]:
    concentrations = []
    for cell in hist_df_0:
        for i,value in enumerate(cell[1]):
            if value < r:
                concentrations.append(cell[0][i])
    avg_c = sum(concentrations)/len(concentrations)

from scipy import signal

# test = signal.convolve([[1,2],[3,4]], [[1,2],[3,4]])