import getpass

import matplotlib.pyplot as plt
import seaborn as sns
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

sampling_timepoint = -2

offset = 0

pSTAT5 = True

from parameters import path

save_path = path
message("loading dataframe from" + save_path)

try:
    cell_df = pd.read_hdf(path + 'cell_df.h5', mode="r")
    message("read h5")
except:
    cell_df = pd.read_pickle(path + 'cell_df.pkl')
    message("read pkl")

# cell_df["time"] = cell_df["time"].div(3600) #s
cell_df["IL-2_surf_c"] = cell_df["IL-2_surf_c"].mul(1e3) #pM

# gammas_list = np.sort(cell_df["IL-2_gamma"].unique())
# gamma = gammas_list[0]

scan_indices = np.sort(cell_df["scan_index"].unique())

replicats = np.sort(cell_df["replicat_index"].unique())
timepoints = cell_df["time"].unique()

cell_df["pSTAT5"] = cell_df["IL-2_surf_c"] ** 3 / (
        (EC50_calculation(E_max=125e-12, E_min=0, k=860, N=1, R=cell_df["IL-2_R"]) * 1e12) ** 3 + cell_df[
    "IL-2_surf_c"] ** 3).values

activation_df = pd.DataFrame()
for replicat in replicats:
    sliced_cell_df = cell_df.loc[(cell_df["type_name"] == "Th") & (cell_df["replicat_index"] == replicat)]
    df_dict_list = []
    for f, scan_index in enumerate(scan_indices[:]):
        message("calculating for scan_index = " + str(scan_index))
        amount_of_cells = len(sliced_cell_df.loc[(sliced_cell_df["scan_index"] == scan_index), "id_id"].unique())

        for t, time in enumerate(timepoints):
            try:
                act_value = len(sliced_cell_df.loc[(np.abs(sliced_cell_df["time"] - float(time)) < 0.0001) &
                                                       (sliced_cell_df["scan_index"] == scan_index) &
                                                       (sliced_cell_df["pSTAT5"] >= cell_df.loc[(cell_df["type_name"] == "Th"), "misc_Km_pos"].unique()[
                                                             0])]) #/ amount_of_cells
            except ZeroDivisionError:
                act_value = 0
            df_dict_list.append({
                "scan_index": scan_index,
                "time": time,
                "activation": act_value
            })

    activation_df = activation_df.append(df_dict_list, ignore_index=True)
activation_df.to_hdf(save_path + "activation_df_" + str(replicat) + ".h5", "w")
for si in activation_df["scan_index"].unique():
    print("scan_index =", si)
    print(activation_df.loc[activation_df["scan_index"] == si, "activation"])

# sns.lineplot(np.sort(cell_df["IL-2_R"].unique())[1:], activation_df.values[0])
# plt.tight_layout()
# plt.show()


