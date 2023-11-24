import getpass
import numpy as np
import pandas as pd
import os
from thesis.main.my_debug import message

def EC50_calculation(E_max, E_min, k, N, R):
    return (E_max * k ** N + E_min * R ** N) / (k ** N + R ** N)

yscale = "linear"
xscale = "linear"

startingPoint = None
stoppingPoint = None

xlim = (None,None)
ylim = (None, None)

sampling_timepoint = -1

offset = 0

pSTAT5 = True

myRange = np.arange(0,2,1)
get_dataframes = []#[[]]*3

message("loading dataframe")
from parameters import path

base_path = path

for run in myRange:
    path = base_path + "run" + str(run) + "/"
    save_path = path
    message("loading dataframe from" + save_path)

    try:
        cell_df = pd.read_hdf(path + 'cell_df.h5', mode="r")
        message("read h5")
    except:
        cell_df = pd.read_pickle(path + 'cell_df.pkl')
        message("read pkl")

    cell_df["time"] = cell_df["time"].div(3600) #s
    cell_df["IL-2_surf_c"] = cell_df["IL-2_surf_c"].mul(1e3) #pM

    gammas_list = np.sort(cell_df["IL-2_gamma"].unique())
    gamma = gammas_list[0]

    scan_indices = np.sort(cell_df["scan_index"].unique())
    timepoints = cell_df["time"].unique()

    cell_df["pSTAT5"] = cell_df["misc_pSTAT5"]

    sliced_cell_df = cell_df.loc[(cell_df["type_name"] == "Th")]
    activation_df = pd.DataFrame(columns=scan_indices, index=timepoints)
    for f, scan_index in enumerate(scan_indices[:]):
        message("calculating for scan_index = " + str(scan_index))
        amount_of_cells = len(sliced_cell_df.loc[(sliced_cell_df["scan_index"] == scan_index), "id_id"].unique())
        activation_list = []
        for t, time in enumerate(timepoints):
            try:
                activation_list.append(len(sliced_cell_df.loc[(np.abs(sliced_cell_df["time"] - float(time)) < 0.0001) &
                                                   (sliced_cell_df["scan_index"] == scan_index) &
                                                   (sliced_cell_df["pSTAT5"] >= cell_df.loc[(cell_df["type_name"] == "Th"), "misc_Km_neg"].unique()[
                                                         0])]) / amount_of_cells)
            except ZeroDivisionError:
                activation_list.append(0)
        activation_df[scan_index] = activation_list
    activation_df.to_hdf(save_path + "activation_df.h5", "w")
    print(activation_df)



