import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import seaborn as sns
import getpass
import KDEpy
import os
import multiprocessing as mp
from skimage.measure import label
import glob
from ext_kde import compute_kde

user = getpass.getuser()
path = "/home/{u}/hauser_ext_kde/".format(u=user)
IMGPATH = path
os.makedirs(IMGPATH,exist_ok = True)

files = []
for file_path in glob.glob(path+"*.csv"):
    with open(file_path) as f:

        df = pd.read_csv(f,names=["x","y"])

        type_name = file_path.split("/")[-1].split(".")[0].split("_")[0]

        series = pd.Series([type_name]*len(df))
        series.index = df.index

        df["type_name"] = series

        series = pd.Series([0] * len(df))
        series.index = df.index

        df["time_index"] = series

        series = pd.Series([0] * len(df))
        series.index = df.index

        df["z"] = series

        files.append(df)


cells = pd.concat(files)

# r,l,kernel_list = compute_kde(cells, bw = 50, lag=0)

