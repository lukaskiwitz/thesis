import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import seaborn as sns
import random

# def get_id_loc(cell_df, id):
#     return np.array([cell_df.loc[cell_df["id"]==str(id), "x"].unique()[0], cell_df.loc[cell_df["id"]==str(id), "y"].unique()[0], cell_df.loc[cell_df["id"]==str(id), "z"].unique()[0]])
# get the distance between two cells. Written in a way so it can be applied to the cell_df_timed dataframe. Speeds up computation significantly
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
from parameters_q_fraction import path_kinetic2
path = path_kinetic2
# path = "/extra/brunner/thesis/kinetic/q_fraction/"
# path = "/extra/brunner/thesis/kinetic/q_fraction_wave/"
cell_df = pd.read_hdf(path + 'cell_df.h5', mode="r")
global_df = pd.read_hdf(path + 'global_df.h5', mode="r")

#for which timepoints to calculate for
time_points = [0., 5., 50.] #cell_df["time"].unique()
threshold = 19000 #if less than x Receptors

try:
    # print("loading name changed to test")
    hist_df = pd.read_hdf(path + "hist_df11_" + str(threshold) + ".h5", mode="r")
    normalising_df = pd.read_hdf(path + "normalising_df" + ".h5", mode="r")
    print("Loaded in previously calculated files")
except FileNotFoundError:
    hist_df = pd.DataFrame(columns=[str(x) for x in time_points], index=cell_df.loc[(cell_df["type_name"] == "changed") & (cell_df["time"] == 0) , "id"].values)
    normalising_df = pd.DataFrame(columns=[str(x) for x in time_points], index=cell_df.loc[(cell_df["type_name"] == "changed") & (cell_df["time"] == 0), "id"].values)
    # exit()
    small_cell_sample = random.sample(cell_df.loc[(cell_df["type_name"] == "changed") & (cell_df["time"] == time_points[0]), "id"].values.tolist(), 10)
    for time in time_points: #cell_df["time"].unique():
        print("time = " + str(time))
        for il2cells in cell_df.loc[(cell_df["type_name"] == "changed") & (cell_df["time"] == 0), "id"]: #small_cell_sample: #
            cell_df["distance"] = None
            cell_df_timed = cell_df.loc[cell_df["time"] == time].copy()
            cell_df_timed["distance"] = None
            cell_df_timed.apply(lambda row: get_distances(il2cells, row["id"]), axis=1)

            for index,row in cell_df_timed.iterrows():
                cell_df.loc[(cell_df["time"] == time) & (cell_df["id"] == row["id"]), "distance"] = row["distance"]
                # cell_df.loc[cell_df["time"] == time, "distance"] = cell_df_timed["distance"]
            # cell_df.to_hdf(path + "cell_df_distance_" + str(time_point) + ".h5", key="df", mode="w")

            data = []
            for index,row in cell_df_timed.iterrows():
                # if row["IL-2_surf_c"] > threshold:
                if row["IL-2_R"] < threshold:
                    data.append(row["distance"])
            hist_df[str(time)][str(il2cells)] = np.array(data).astype(np.float) #cell_df.loc[cell_df["time"] == time, "distance"]
            normalising_df[str(time)][str(il2cells)] = cell_df_timed["distance"].to_numpy().astype(np.float)

    hist_df.to_hdf(path + "hist_df_" + str(threshold) + ".h5",key="df", mode="w")
    normalising_df.to_hdf(path + "normalising_df" + ".h5", key="df", mode="w")

# sns.set()
# sns.set_context("talk", font_scale=1, rc={"lines.linewidth": 0.2})
# sns.set(rc={"lines.linewidth":0.2})

# add R=20000 to t=0
# zero_df = cell_df.loc[cell_df["time"] == 0.0].copy()
# zero_df["IL-2_R"] = 20000.0
# cell_df["time"] = cell_df["time"].add(1)
# cell_df = pd.concat([zero_df, cell_df])
# fig = plt.figure()
# ax = sns.lineplot(x="time", y="IL-2_R", data=cell_df, estimator=None, units="id_id", legend=False)
# ax.set(xlabel="time in h", ylabel="each cells R") #, ylim=(0.275,0.45))
# plt.show()

sns.set_style("ticks")
sns.set_context("talk", font_scale=1, rc={"lines.linewidth": 2})
# plot_time_point = "12.0"
time_list = ["0.0", "5.0", "50.0"]
# time_list.reverse()
for plot_time_point in time_list:
    cell_cell_distance = 20
    max_distance = np.array([x for x in normalising_df[plot_time_point]]).flatten().max()
    # max_distance = np.array([x for x in normalising_df[plot_time_point][10]]).max()
    bins = np.arange(0,(max_distance//cell_cell_distance + 2)*20, cell_cell_distance)
    norm_hist, norm_bins = np.histogram(np.array([x for x in normalising_df[plot_time_point]]).flatten(), bins = bins)
    hist, bins = np.histogram(np.array([x for x in hist_df[plot_time_point]]).flatten(), bins = norm_bins)
    # norm_hist, norm_bins = np.histogram(np.array([x for x in normalising_df[plot_time_point][10]]), bins = bins)
    # hist, bins = np.histogram(np.array([x for x in hist_df[plot_time_point][10]]), bins = norm_bins)

    width = 0.9 * (bins[1] - bins[0])
    center = (bins[:-1] + bins[1:]) / 2
    #normalising for every cell
    hist_norm = hist/len(hist_df[str(time_points[0])])
    #normalising for every cell in relation to how many cells there are in each bin
    hist_norm = hist_norm/(norm_hist/len(normalising_df[str(time_points[0])]))
    #standard error calculation:
    error = np.zeros(len(hist_norm))
    # create list with single cells histograms, normed
    testlist = []
    for j, cellsDistances in enumerate(hist_df[plot_time_point]):
        testlist.append([])
        testlist[j] = np.histogram(hist_df[plot_time_point][j], bins=norm_bins)[0]/np.histogram(normalising_df[plot_time_point][j], bins = norm_bins)[0]
    # remove nans from 0/0 and find largest difference to mean -> std
    for i in range(len(error)):
        tempList = np.transpose(testlist)[i]
        where_are_NaNs = np.isnan(tempList)
        tempList[where_are_NaNs] = 0
        error[i] = np.abs(tempList - hist_norm[i]).std()
    # print(testlist)
    plt.rcParams.update({'lines.markeredgewidth': 2})
    plt.errorbar(center[:13], hist_norm[:13], yerr=error[:13], label=str(float(plot_time_point)*10)+"h", fmt='-o', capsize=10)
    plt.legend()
    # plt.title("time = " + str(int((float(plot_time_point)+0)*10)) + "h")
    plt.title("slice 2")
    plt.xlabel("Cell-cell-distances to secreting cell")
    plt.ylabel("Activated cells in percent of total")
    plt.ylim((0,1.05))
    plt.xticks(bins[:14]-10, [int(x/20) for x in bins][:14])
    # plt.xlim((0,240))
plt.show()
