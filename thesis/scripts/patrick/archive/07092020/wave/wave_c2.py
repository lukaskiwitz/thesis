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
from parameters_q_fraction import path_kinetic
path = path_kinetic
# path = "/extra/brunner/thesis/kinetic/q_fraction/"
# path = "/extra/brunner/thesis/kinetic/q_fraction_wave/"
# path = "/extra/brunner/thesis/kinetic/q_fraction_wave_gamma_scan/"
path = "/extra/brunner/thesis/kinetic/q_fraction_small_lower_q/"

cell_df = pd.read_hdf(path + 'cell_df.h5', mode="r")
global_df = pd.read_hdf(path + 'global_df.h5', mode="r")

cell_df = cell_df.loc[cell_df["IL-2_gamma"] == 0.1]

#for which timepoints to calculate for
time_points = [0., 5., 19.] #cell_df["time"].unique()
threshold = 1e10 #15000 #if less than x Receptors

try:
    # print("loading name changed to test")
    hist_df = pd.read_hdf(path + "hist_df_c_011111_large_" + str(threshold) + ".h5", mode="r")
    normalising_df = pd.read_hdf(path + "normalising_df_c_" + ".h5", mode="r")
    print("Loaded in previously calculated files")
except FileNotFoundError:
    #create multiindex dfs
    hist_df = pd.DataFrame(columns=[str(x) for x in time_points], index=cell_df.loc[(cell_df["type_name"] == "changed") & (cell_df["time"] == 0) , "id"].values)
    normalising_df = pd.DataFrame(columns=[str(x) for x in time_points], index=cell_df.loc[(cell_df["type_name"] == "changed") & (cell_df["time"] == 0), "id"].values)
    # exit()
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
            # write calculated distances into cell_df
            # for index,row in cell_df_timed.iterrows():
            #     cell_df.loc[(cell_df["time"] == time) & (cell_df["id"] == row["id"]), "distance"] = row["distance"]
                # cell_df.loc[cell_df["time"] == time, "distance"] = cell_df_timed["distance"]
            # cell_df.to_hdf(path + "cell_df_distance_" + str(time_point) + ".h5", key="df", mode="w")

            data = []
            distance = []
            # threshold check
            for index,row in cell_df_timed.iterrows():
                # if row["IL-2_surf_c"] > threshold:
                if row["IL-2_R"] < threshold:
                    data.append(row["IL-2_surf_c"]*1e3)
                    distance.append(row["distance"])
            # put into df
            hist_df[str(time)][str(il2cells)] = [np.array(data).astype(np.float).copy(), np.array(distance).astype(np.float).copy()] #cell_df.loc[cell_df["time"] == time, "distance"]
            normalising_df[str(time)][str(il2cells)] = [cell_df_timed["IL-2_surf_c"].mul(1e3).to_numpy().astype(np.float).copy(), cell_df_timed["distance"].copy()]
    # save to file
    hist_df.to_hdf(path + "hist_df_c_2_large_" + str(threshold) + ".h5",key="df", mode="w")
    normalising_df.to_hdf(path + "normalising_df_c_" + ".h5", key="df", mode="w")


sns.set_style("ticks")
sns.set_context("talk", font_scale=1, rc={"lines.linewidth": 3})
# plot_time_point = "12.0"
# time_list = ["0.0", "5.0", "10.0", "99.0"]
# time_list = ["0.0", "5.0", "50.0"]
time_list = ["0.0", "5.0", "19.0"]
# time_list.reverse()
for plot_time_point in time_list:
    # calculate necessary bins
    cell_cell_distance = 20
    max_distance = np.array([x for x in normalising_df[plot_time_point]]).flatten().max()
    # max_distance = np.array([x for x in normalising_df[plot_time_point][10]]).max()
    bins = np.arange(0,(max_distance//cell_cell_distance)*20, cell_cell_distance)
    # create container for histogram computation
    c_bin_array = np.zeros((len(bins)-1, len(hist_df[plot_time_point][0][1])*10))
    c_bin_array[:] = np.NaN
    # iterate over cells distances, filling the bins
    for cell in hist_df[plot_time_point]:
        for i in range(len(cell[1])):
            bin_index = 0
            while bin_index <= len(bins)-2:
                if cell[1][i] <= bins[bin_index+1]:
                    for k, entry in enumerate(c_bin_array[bin_index]):
                        # print(entry)
                        if np.isnan(entry) == True:
                            c_bin_array[bin_index][k] = cell[0][i]
                            break
                    break
                else:
                    bin_index += 1
    # calc mean and std of bins
    hist = [np.nanmean(k) for k in c_bin_array]
    error = [np.nanstd(k, dtype=np.float64) for k in c_bin_array]
    # lolims = [0 for k in error]

#
    width = 0.9 * (bins[1] - bins[0])
    center = (bins[:-1] + bins[1:]) / 2

    # plt.errorbar(center[:13], hist[:13], yerr=error[:13], label=str(float(plot_time_point)*10)+"h", fmt='-o', capsize=8)
    plt.plot(center[:13], hist[:13], '-o', label=str(float(plot_time_point) * 10) + "h")
    plt.legend(loc=1)
    #     # plt.title("time = " + str(int((float(plot_time_point)+0)*10)) + "h")
    #     plt.title("Change over time")
    plt.xlabel("Cell-cell-distances to secreting cell")
    plt.ylabel("Mean c. in pM")
    # plt.ylim((0,380))
    plt.xticks(bins[:14] - 10, [int(x / 20) for x in bins][:14])
    # plt.xlim((0,240))
# fig.savefig("/home/brunner/Documents/Current work/19062020/" + "wave_c01" + ".png", bbox_inches='tight')
plt.show()