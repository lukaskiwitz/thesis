import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import seaborn as sns
import os
import getpass
from parameters import model_name, scan_name
hdd = "extra2" if os.path.exists("/extra2") else "extra"

user = getpass.getuser()

path = f"/{hdd}/{user}/paper_models/kinetics/{model_name}/dataframes_{scan_name}_timeseries/"
os.listdir(path)

paths = []
raw_paths = os.listdir(path)
for entry in raw_paths:
    paths.append(f"{path}{entry}/")
paths = np.array(paths)[[os.path.isdir(path) for path in paths]] #check for only directories

paths = np.array(paths).flatten()
print(paths)

df_list = []

for idx, sub_path in enumerate(paths):
    print("loading", idx + 1)
    try:
        cell_df = pd.read_hdf(sub_path + 'cell_df.h5', mode="r")
    except:
        cell_df = pd.read_pickle(sub_path + "cell_df.pkl")
    global_df = pd.read_hdf(sub_path + 'global_df.h5', mode="r")
    cell_df = cell_df.drop(
        columns=list(set(cell_df.columns) - {"IL-2_surf_c", "IL-2_R", "IL-2_gamma", "IL-2_pSTAT5", "IL-2_q",
                                             "IL-2_Tsec_fraction", "id", "scan_index", "scan_name_scan_name",
                                             "scan_value", "time", "time_index", "replicat_index", "x", "y", "z",
                                             "type_name", "misc_pos_half_time", "misc_neg_half_time", "misc_sec_start",
                                             "misc_gamma", "clustering_strength", "fractions_Treg", "fractions_Tsec"}))
    # cell_df["replicat_index"] += idx//4
    df_list.append((cell_df, global_df))
    print(len(cell_df))
    # if len (df_list) == 0:
    #     df_list.append((cell_df, global_df))
    # else:
    #     cell_df.replicat_index += (df_list[-1][0].replicat_index.max() + 1)
    #     global_df.replicat_index += (df_list[-1][1].replicat_index.max() + 1)
    #
    #     df_list.append((cell_df, global_df))

# df_list = np.array(df_list, dtype=object)

# cell_df = pd.concat(df_list[:,0])
# global_df = pd.concat(df_list[:,1])
cell_df = pd.concat([x[0] for x in df_list])
global_df = pd.concat([x[1] for x in df_list])

name = ""
print(f"saving combined dfs to .h5 at {path}")
try:
    global_df.to_hdf(path + "global_df_combined.h5", "w")
    cell_df.to_hdf(path + "cell_df_combined.h5", "w")
except:
    print("Saving the cell_df to .h5 failed, falling back to pickling...")
    global_df.to_pickle(path + "global_df_combined.pkl")
    cell_df.to_pickle(path + "cell_df_combined.pkl")
