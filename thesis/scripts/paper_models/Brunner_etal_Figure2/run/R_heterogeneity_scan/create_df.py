import glob
import os
import re
import signal

import matplotlib.pyplot as plt
import numpy as np

from thesis.main.MyPlotter import Plotter
from thesis.main.my_debug import message
from parameters import path
from parameters import scan_name
collect_path = "/".join(path.split("/")[:-2])

# IMGPATH = os.path.join(collect_path+"_"+scan_name,"images/")
# DF_COLLECTION_PATH = os.path.join(collect_path+"_"+scan_name,"collected_dfs/")
#
# os.makedirs(IMGPATH,exist_ok=True)
# os.makedirs(DF_COLLECTION_PATH,exist_ok=True)
#


path_list = glob.glob(os.path.join(
    "/".join(path.split("/")[:-2])
    , scan_name + "*"))

all_paths = []


SS_FOLDER_PATH = os.path.join(collect_path, "dataframes_"+scan_name + "_steady_state")
FULL_FOLDER_PATH = os.path.join(collect_path, "dataframes_"+scan_name + "_timeseries")

os.makedirs(SS_FOLDER_PATH, exist_ok=True)
os.makedirs(FULL_FOLDER_PATH, exist_ok=True)

for path in path_list:

    refine_paths = glob.glob(os.path.join(path, "refine_*"))
    refine_paths = [re.search(r"refine_(?P<n>\d)", rp) for rp in refine_paths]
    refine_paths = [rp.string for rp in refine_paths if (int(rp.groupdict()["n"]) in [0, 1, 2, 3, 4, 5])]
    all_paths = all_paths + refine_paths

    try:
        plotter = Plotter(refine_paths)
        scan_folder_name = path.split("/")[-1]

        full_folder_path = os.path.join(FULL_FOLDER_PATH, scan_folder_name)
        ss_folder_path = os.path.join(SS_FOLDER_PATH, scan_folder_name)

        # print(full_folder_path)
        # print(ss_folder_path)
        os.makedirs(full_folder_path, exist_ok=True)
        os.makedirs(ss_folder_path, exist_ok=True)

        plotter.cell_df.to_hdf(os.path.join(full_folder_path, "cell_df.h5"), key="df", mode="w")
        try:
            plotter.global_df.to_hdf(os.path.join(full_folder_path, "global_df.h5"), key="df", mode="w")
        except Exception as e:
            print(e)

        df = plotter.cell_df
        df = df.loc[(df["time"] == plotter.t_max) | (df[plotter.time_index_key] == 0)]
        df.to_hdf(os.path.join(ss_folder_path, "cell_df.h5"), key="df", mode="w")

    except Exception as e:print(e)




all_paths = []

for path in path_list:
    refine_paths = glob.glob(os.path.join(path, "refine_*"))
    refine_paths = [re.search(r"refine_(?P<n>\d)", rp) for rp in refine_paths]
    refine_paths = [rp.string for rp in refine_paths if (int(rp.groupdict()["n"]) in [0, 1, 2, 3, 4, 5])]
    all_paths = all_paths + refine_paths
pass
print(FULL_FOLDER_PATH)
print(SS_FOLDER_PATH)

plotter = Plotter(all_paths, load_dataframes = {
            "global":True,
            "cell":True,
        })

import pandas as pd

def handler(signum, frame):
    raise Exception("This exception shouldnt appear, something went wrong...")
signal.signal(signal.SIGALRM, handler)

def get_steady_state(df):
    def f(df):
        return df.loc[(df.time_index == df.time_index.max()) | (df.time_index == 0)]

    ind = pd.MultiIndex.from_frame(
        df.loc[:,
        [plotter.model_index_key, plotter.scan_index_key, plotter.scan_name_key, plotter.replicat_index_key]]
    )
    df = df.set_index(ind, drop=True)
    df = df.drop(columns=df.index.names)

    signal.alarm(60*1)
    try:
        t_max = df.groupby(axis=0, level=df.index.names, group_keys=False).apply(f)
    except Exception:
        print("group_by exceeded alloted computation time, simplifying...")
        t_max = df.loc[(df.time_index == df.time_index.max()) | (df.time_index == 0)]
    finally:
        signal.alarm(0)
    t_max = t_max.reset_index()
    return t_max

# scan_folder_name = path.split("/")[-1]
# plotter.cell_df.to_hdf(os.path.join(FULL_FOLDER_PATH, "cell_df.h5"), key="df", mode="w")
# plotter.global_df.to_hdf(os.path.join(FULL_FOLDER_PATH, "global_df.h5"), key="df", mode="w")


df = plotter.cell_df
df = get_steady_state(df)
df.to_hdf(os.path.join(SS_FOLDER_PATH, "cell_df.h5"), key="df", mode="w")

df = plotter.global_df
df = get_steady_state(df)
df.to_hdf(os.path.join(SS_FOLDER_PATH, "global_df.h5"), key="df", mode="w")

