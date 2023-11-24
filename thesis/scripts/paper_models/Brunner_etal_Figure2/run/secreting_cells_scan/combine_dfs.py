import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import seaborn as sns
import os
import getpass

hdd = "extra2" if os.path.exists("/extra2") else "extra"

user = getpass.getuser()


paths = []
repeats = [0]
for repeat in repeats:
    raw_paths = [
        f"pos_treg_{repeat}_sec_cs_strength_scanftreg_0_fth_0.95_fsec_0.05",
        f"pos_treg_{repeat}_sec_cs_strength_scanftreg_0_fth_0.975_fsec_0.025",
        f"pos_treg_{repeat}_sec_cs_strength_scanftreg_0_fth_0.99_fsec_0.01",
        f"pos_treg_{repeat}_sec_cs_strength_scanftreg_0_fth_0.9_fsec_0.1",
        f"pos_treg_{repeat}_treg_cs_strength_scanfth_0.75_ftreg_0.2_fsec_0.05",
        f"pos_treg_{repeat}_treg_cs_strength_scanfth_0.925_ftreg_0.025_fsec_0.05",
        f"pos_treg_{repeat}_treg_cs_strength_scanfth_0.95_ftreg_0_fsec_0.05",
        f"pos_treg_{repeat}_treg_cs_strength_scanfth_0.9_ftreg_0.05_fsec_0.05"]
    for entry in raw_paths:
        paths.append(f"/{hdd}/{user}/Tsec_clustering/300_3D/dataframes_pos_treg_{repeat}_timeseries/{entry}/")


paths = np.array(paths).flatten()
print(paths)

df_list = []

for idx, path in enumerate(paths):
    print("loading", idx + 1)
    cell_df = pd.read_hdf(path + 'cell_df.h5', mode="r")
    global_df = pd.read_hdf(path + 'global_df.h5', mode="r")
    cell_df = cell_df.drop(
        columns=list(set(cell_df.columns) - {"IL-2_surf_c", "IL-2_R", "IL-2_gamma", "IL-2_pSTAT5", "IL-2_q",
                                             "IL-2_Tsec_fraction", "id", "scan_index", "scan_name_scan_name",
                                             "scan_value", "time", "time_index", "replicat_index", "x", "y", "z",
                                             "type_name", "misc_pos_half_time", "misc_neg_half_time", "misc_sec_start",
                                             "clustering_strength", "fractions_Treg", "fractions_Tsec"}))
    cell_df["replicat_index"] += idx//4
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
path_0 = f"/extra2/brunner/Tsec_clustering/300_3D/dataframes_pos_treg_{repeats[0]}_timeseries/{name}/"
print(f"saving combined dfs to .h5 at {path_0}")
try:
    global_df.to_hdf(path_0 + "global_df_combined.h5", "w")
    cell_df.to_hdf(path_0 + "cell_df_combined.h5", "w")
except:
    print("Saving the cell_df to .h5 failed, falling back to pickling...")
    global_df.to_pickle(path_0 + "global_df_combined.pkl")
    cell_df.to_pickle(path_0 + "cell_df_combined.pkl")
