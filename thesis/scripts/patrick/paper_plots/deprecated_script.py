import numpy as np

import os.path

import sys

##
import pandas as pd
import subprocess
import matplotlib.pyplot as plt
import seaborn as sns

def averageRuns(basePath, myRange, groupVariable):
    get_dataframes = []  # [[]]*3
    for j in range(myRange):
        path = basePath + "run" + str(j) + "/"
        get_dataframes.append([])
        print("reading in run", j)

        global_df = pd.read_hdf(path + 'global_df.h5', mode="r")
        cell_df = pd.read_hdf(path + 'cell_df.h5', mode="r")

        try:
            global_df["surf_c_std"]
        except KeyError:
            global_df["surf_c_std"] = np.zeros(len(global_df["time"]))
            for scan_index in global_df["scan_index"].unique():
                for sigma in global_df[groupVariable].unique():
                    global_df.loc[(global_df["scan_index"] == scan_index) & (global_df[groupVariable] == sigma), "surf_c_std"] = \
                        cell_df.loc[(cell_df["scan_index"] == scan_index) & (cell_df[groupVariable] == sigma), "IL-2_surf_c"].std()
            global_df.to_hdf(path + 'global_df.h5', key="data", mode="w")

        get_dataframes[j] = [global_df, cell_df]

    global_df = pd.concat((get_dataframes[x][0] for x in range(len(get_dataframes)))).groupby([groupVariable],as_index=False).mean()
    global_df.reset_index(inplace=True)
    # global_df["IL-2_D"] = get_dataframes[0][0]["IL-2_D"]
    # cell_df = pd.concat((get_dataframes[x][1] for x in range(len(get_dataframes)))).groupby([groupVariable], as_index=False).mean()
    # global_df["IL-2_D"] = get_dataframes[0][1]["IL-2_D"]
    return global_df, cell_df

def averageRunsFraction(basePath, myRange, groupVariable):
    get_dataframes = []  # [[]]*3
    for j in range(myRange):
        # path = basePath + "run" + str(j) + "/"
        path= basePath
        get_dataframes.append([])
        print("reading in run", j)

        # global_df = pd.read_hdf(path + 'global_df.h5', mode="r")
        # cell_df = pd.read_hdf(path + 'cell_df.h5', mode="r")

        global_df = pd.read_hdf(path + 'global_df.h5', mode="r")
        cell_df = pd.read_hdf(path + 'cell_df.h5', mode="r")



        if groupVariable == "IL-2_fraction":
            global_df = global_df.groupby([groupVariable], as_index=False).mean()
            try:
                global_df["surf_c_std"]
            except KeyError:
                global_df["surf_c_std"] = np.zeros(len(global_df["time"]))
                for scan_index in global_df["scan_index"].unique():
                    for var in global_df[groupVariable].unique():
                        global_df.loc[(global_df[groupVariable] == round(var,7)), "surf_c_std"] = cell_df.loc[cell_df[groupVariable] == round(var,7), "IL-2_surf_c"].std()
                # global_df.to_hdf(path + 'global_df.h5', key="data", mode="w")
        else:
            try:
                global_df["surf_c_std"]
            except KeyError:
                global_df["surf_c_std"] = np.zeros(len(global_df["time"]))
                for scan_index in global_df["scan_index"].unique():
                    for var in global_df[groupVariable].unique():
                        global_df.loc[
                            (global_df["scan_index"] == scan_index) & (global_df[groupVariable] == var), "surf_c_std"] = \
                            cell_df.loc[(cell_df["scan_index"] == scan_index) & (
                                        cell_df[groupVariable] == var), "IL-2_surf_c"].std()
        get_dataframes[j] = [global_df, cell_df]

    global_df = pd.concat((get_dataframes[x][0] for x in range(len(get_dataframes)))).groupby([groupVariable],as_index=False).mean()
    global_df.reset_index(inplace=True)
    # global_df["IL-2_D"] = get_dataframes[0][0]["IL-2_D"]
    # cell_df = pd.concat((get_dataframes[x][1] for x in range(len(get_dataframes)))).groupby([groupVariable], as_index=False).mean()
    # global_df["IL-2_D"] = get_dataframes[0][1]["IL-2_D"]
    return global_df, cell_df


# path = "/extra/brunner/thesis/static/R_lognorm/"
path = "/extra/brunner/thesis/static/q_fraction/"

# global_df = pd.read_hdf(path + 'global_df.h5', mode="r")
# cell_df = pd.read_hdf(path + 'cell_df.h5', mode="r")

sns.set(rc={'figure.figsize':(5,5)})
sns.set_style("ticks")
sns.set_context("talk", font_scale=1, rc={"lines.linewidth": 5})

xscale = "linear"
yscale = "linear"

data_start = None
data_end = None

group_variables = ["IL-2_fraction", None]
# group_variables = ["IL-2_sigma", None]

# global_df = pd.read_hdf(path + 'global_df.h5', mode="r")
# cell_df = pd.read_hdf(path + 'cell_df.h5', mode="r")

global_df, cell_df = averageRunsFraction(path, 1, group_variables[0])

# if group_variables[0] == "IL-2_fraction":
#     global_df[group_variables[1]] = ["$%s$" % round(x, 3) for x in global_df[group_variables[1]]]

# if "IL-2_sigma" in group_variables:
#     global_df["IL-2_sigma"] = global_df["IL-2_sigma"].div(20000)

global_df["IL-2_fraction"] = (1 - global_df["IL-2_fraction"][7:])

fig = plt.figure()
global_df["surf_c_std_norm"] = global_df["surf_c_std"]/global_df["surf_c"]
ax1 = sns.lineplot(x=group_variables[0],
                   y="surf_c_std_norm",
                   data=global_df[data_start:data_end],
                   label="fraction",
                   legend=False,
                   color="firebrick")
ax1.set_ylabel("c$_v$", color="black")
ax1.set(xticks = [0, 0.5, 1], yscale=yscale, xscale=xscale)

# ax4 = sns.lineplot(x=group_variables[0], y="surf_c",
#                    data=global_df[data_start:data_end],
#                    label="surface c.",
#                    legend=False)
# ax4.set(yscale=yscale, xscale=xscale)
# if group_variables[0] == "IL-2_sigma":
#     ax1.set(yticks=[1e-3, 1e-2, 1e-1], yticklabels = ["$10^{-3}$", "$10^{-2}$", "$10^{-1}$"], ylabel = "nM", xlabel = "Heterogeneity", ylim=(4.5e-4, 1e-1))
# else:
#     ax1.set(xticks = [0, 0.5, 1], ylim=(None,None), xlabel = "Heterogeneity")

path = "/extra/brunner/para_handling/static/R_lognorm/"
group_variables = ["IL-2_sigma", None]
global_df, cell_df = averageRuns(path, 14, group_variables[0])
global_df["IL-2_sigma"] = global_df["IL-2_sigma"].div(20000)

data_end = -1

global_df["surf_c_std_norm"] = global_df["surf_c_std"]/global_df["surf_c"]
ax2 = sns.lineplot(x=group_variables[0],
                   y="surf_c_std_norm",
                   data=global_df[data_start:data_end],
                   label="sigma",
                   legend=False,
                   color="green")
ax2.set_ylabel("c$_v$", color="black")
ax2.set(yscale=yscale, xscale=xscale, ylim=(-0.1,None), xlabel = "Heterogeneity")
ax2.set_xticks([0, 0.5, 1])
plt.legend()
fig.savefig("/home/brunner/Documents/Current work/05062020/" + "static_c_v" + ".png", bbox_inches='tight')
plt.show()