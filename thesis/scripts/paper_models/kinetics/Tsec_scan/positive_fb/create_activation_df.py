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

myRange = [0] #np.arange(0,10,1)
get_dataframes = []#[[]]*3

message("loading dataframe")
from parameters import path

for run in np.arange(1):
    path = path + "/run" + str(run) + "/"
    save_path = path

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

    fraction_list = np.sort(cell_df["IL-2_Tsec_fraction"].unique())
    timepoints = cell_df["time"].unique()

    cell_df["pSTAT5"] = cell_df["misc_pSTAT5"]

    sliced_cell_df = cell_df.loc[(cell_df["type_name"] == "Th")]
    activation_df = pd.DataFrame(columns=fraction_list, index=timepoints)
    for f, frac in enumerate(fraction_list[:-1]):
        message("calculating for Tsec_fraction = " + str(frac))
        amount_of_cells = len(sliced_cell_df.loc[(sliced_cell_df["IL-2_Tsec_fraction"] == frac), "id_id"].unique())
        activation_list = []
        for t, time in enumerate(timepoints):
            activation_list.append(len(sliced_cell_df.loc[(np.abs(sliced_cell_df["time"] - float(time)) < 0.0001) &
                                                   (sliced_cell_df["IL-2_Tsec_fraction"] == frac) &
                                                   (sliced_cell_df["pSTAT5"] >= cell_df.loc[(cell_df["type_name"] == "Th"), "misc_Km_pos"].unique()[
                                                         0])]) / amount_of_cells)
        activation_df[frac] = activation_list
    activation_df.to_hdf(save_path + "activation_df.h5", "w")


