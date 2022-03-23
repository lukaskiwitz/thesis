import numpy as np
import pandas as pd
import os
import getpass

hdd = "extra2" if os.path.exists("/extra2") else "extra"
user = getpass.getuser()

'''
used to combine dataframes from multiple runs with different feedback strengths into one. Allows for easier plotting, optional.
'''

model_name = "Tsec_scan_25"
name = ""
path_1 = "/{extra}/{u}/paper_models/kinetics/{mn}/{n}/".format(u=user, mn=model_name, n=name, extra = hdd)
model_name = "Tsec_scan_50"
name = ""
path_2 = "/{extra}/{u}/paper_models/kinetics/{mn}/{n}/".format(u=user, mn=model_name, n=name, extra = hdd)

df_list = []

for idx, path in enumerate([path_1, path_2]):
    print("loading", idx + 1)
    cell_df = pd.read_hdf(path + 'cell_df.h5', mode="r")
    global_df = pd.read_hdf(path + 'global_df.h5', mode="r")
    if len (df_list) == 0:
        df_list.append((cell_df, global_df))
    else:
        cell_df.replicat_index += (df_list[-1][0].replicat_index.max() + 1)
        global_df.replicat_index += (df_list[-1][1].replicat_index.max() + 1)

        df_list.append((cell_df, global_df))

df_list = np.array(df_list, dtype=object)

cell_df = pd.concat(df_list[:,0])
global_df = pd.concat(df_list[:,1])
print("saving combined dfs to .h5")
global_df.to_hdf(path_2 + "global_df_combined.h5", "w")
cell_df.to_hdf(path_2 + "cell_df_combined.h5", "w")
