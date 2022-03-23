import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import seaborn as sns
import os
from thesis.scripts.paper_models.utilities.plot_helper import activation_to_global_df
from plotting_rc import rc_ticks

sns.set_theme(context = "talk", style = "ticks", rc = rc_ticks)

fig, ax = plt.subplots()

def EC50_calculation(E_max, E_min, k, N, R):
    return (E_max * k ** N + E_min * R ** N) / (k ** N + R ** N)

hdd = "/extra2" if os.path.exists("/extra2") else "/extra"
path = hdd + "/brunner/paper_models/prelim_boxed_static/distance_scan/"

fig_path = "/home/brunner/Documents/Current work/2022_03_04/"

global_df = pd.read_hdf(path + 'global_df.h5', mode="r")
cell_df = pd.read_hdf(path + 'cell_df.h5', mode="r")
global_df = global_df.loc[(global_df["model_name"] == "pde_model")]
cell_df = cell_df.loc[(cell_df["model_name"] == "pde_model")]

global_df = activation_to_global_df(cell_df, global_df)
#%%
linear_uptake = False
# uptake = "patrick_saturation"

scan_names = ["abs_R", "abs_R", "ratio", "sec_q", "D", "distance", "kd"]
x_tick_labels = ["R", "k$_{endo}$", "ratio", "q", "D", "c-c-d", "$\eta$"]
# scan_names = ["sec_q"]
scan_values = [0.1, 10]
ylim = (-1.05, 2.1)

scan_measure = "CV"
results = []
error = []

for name in scan_names:
    for value in scan_values:
        standard = global_df.loc[(global_df["numeric_linear"] == linear_uptake) &
                                     (global_df["scan_name_scan_name"] == name) &
                                     (global_df["scan_value"] == 1), scan_measure].mean()

        a_mean = global_df.loc[(global_df["numeric_linear"] == linear_uptake) &
                                     (global_df["scan_name_scan_name"] == name) &
                                     (global_df["scan_value"] == value), scan_measure].mean()
        results.append((a_mean - standard))

        a = global_df.loc[(global_df["numeric_linear"] == linear_uptake) &
                                     (global_df["scan_name_scan_name"] == name) &
                                     (global_df["scan_value"] == value), scan_measure]
        error.append((a - standard).std())

x_ticks = [0.5 + 3* x for x in range(len(scan_names))]
x = np.array([[3*x, 3*x + 1] for x in range(len(scan_names))]).flatten()

bars = np.array(results)
ax.bar(x=x[::2] , height= bars[::2], edgecolor = "black", color="black", yerr=error[::2])
ax.bar(x=x[1::2] , height= bars[1::2], edgecolor = "grey", color="grey", yerr=error[1::2])
ax.set_xticks(x_ticks)
ax.set_xticklabels(x_tick_labels)
ax.set_ylabel("change in cv")
ax.set_ylim(ylim)

import matplotlib.patches as mpatches

white_patch = mpatches.Patch(facecolor='white', edgecolor="white", alpha=1)
black_patch = mpatches.Patch(facecolor='black', edgecolor="black", alpha=1)
grey_patch = mpatches.Patch(facecolor='grey', edgecolor="grey", alpha=1)

plt.legend(handles = [white_patch, black_patch, grey_patch], labels = ["fold change", "0.1", "10"])

fig.savefig(fig_path + "test.pdf", bbox_inches='tight')
plt.tight_layout()
plt.show()