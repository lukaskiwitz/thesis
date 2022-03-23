import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import seaborn as sns
import os
import getpass

hdd = "extra2" if os.path.exists("/extra2") else "extra"

user = getpass.getuser()


model_name = "saturated"
name = "Tsec_scan_0.05_standard"
path_1 = "/{extra}/{u}/paper_models/statics/{mn}/{n}/".format(u=user, mn=model_name, n=name, extra = hdd)
name = "Tsec_scan_0.05_standard_2"
path_2 = "/{extra}/{u}/paper_models/statics/{mn}/{n}/".format(u=user, mn=model_name, n=name, extra = hdd)
name = "Tsec_scan_0.05_standard_3"
path_3 = "/{extra}/{u}/paper_models/statics/{mn}/{n}/".format(u=user, mn=model_name, n=name, extra = hdd)

df_list = []

for idx, path in enumerate([path_1, path_2, path_3]):
    print("loading", idx + 1)
    cell_df = pd.read_hdf(path + 'cell_df.h5', mode="r")
    global_df = pd.read_hdf(path + 'global_df.h5', mode="r")
    if len (df_list) == 0:
        df_list.append((cell_df, global_df))
    else:
        cell_df.replicat_index += (df_list[-1][0].replicat_index.max() + 1)
        global_df.replicat_index += (df_list[-1][1].replicat_index.max() + 1)

        df_list.append((cell_df, global_df))

df_list = np.array(df_list)

cell_df = pd.concat(df_list[:,0])
global_df = pd.concat(df_list[:,1])

global_df.to_hdf(path_3 + "global_df_combined.h5", "w")
cell_df.to_hdf(path_3 + "cell_df_combined.h5", "w")
