import glob
import os

import matplotlib.pyplot as plt
import numpy as np

from parameters import path

autocorr_path = os.path.join(path, "autocorr")
d = []
for file in glob.glob(autocorr_path + "/" + "*.npy"):
    c = np.load(file)
    si = file.split("/")[-1].split(".")[0].split("_")[1]
    ti = file.split("/")[-1].split(".")[0].split("_")[2]
    id = file.split("/")[-1].split(".")[0].split("_")[3]
    d.append({
        "scan_index": si,
        "time_index": ti,
        "id": id,
        "auto": c
    })
    plt.contourf(c, levels=100)
    plt.show()
