import getpass
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import seaborn as sns
import os

saving_string =r"/home/brunner/Documents/Current work/2021_10_15/" + "R_log_surf_c_SD_test" + ".png"

# define which static sim you want to plot. "IL-2_Tsec_fraction" is the secreting cells scan, "IL-2_sigma" the R_lognorm scan
# scan_variable = "IL-2_Tsec_fraction"
scan_variable = "IL-2_sigma"
# run range
myRange = np.arange(0,1,1)

# plotting parameters
sns.set_style("ticks")
sns.set_context("talk", font_scale=2, rc={"lines.linewidth": 2.5})
fig,ax = plt.subplots(figsize=(7,6))

plot_legend = True

yscale = "linear"
xscale = "linear"
# define data start and end point
startingPoint = None
stoppingPoint = None
# maximum fraction of secreting cells
max_frac = 1

xlim = (None, None)
ylim = (None, None)

# defines a outer layer of N cells to be ignored in plotting. Used to further limit unwanted boundary effects.
offset = 0
# which cell type to plot in which colour. The second colour is the SD colour
cell_types = ["Th"]
y_colours = ["black", "orange"]

y_variable = "surf_c"
title = "Surface concentration"

########################################################################################################################
########################################################################################################################
##                                              Plotting                                                              ##
########################################################################################################################
########################################################################################################################

# load runs, apply offset, merge dataframes
get_dataframes = []#[[]]*3
for idx,value in enumerate(myRange):
    get_dataframes.append([])
    print("loading run", value)
    # path = "/extra/brunner/thesis/kinetic/standard/3_old/run" + str(value) + "/"
    hdd = "/extra2" if os.path.exists("/extra2") else "/extra"
    if scan_variable == "IL-2_Tsec_fraction":
        scan_path = "/extra2/brunner/paper_models/statics/Tsec_scan/"
    elif scan_variable == "IL-2_sigma":
        scan_path = "/extra2/brunner/paper_models/statics/R_lognorm/"

    path = scan_path + "run" + str(value) + "/"

    global_df = pd.read_hdf(path + 'global_df.h5', mode="r")
    cell_df = pd.read_hdf(path + 'cell_df.h5', mode="r")

    x, y, z = (cell_df["x"].unique(), cell_df["y"].unique(), cell_df["z"].unique())
    if offset != 0:
        try:
            offset_cells_ids = cell_df.loc[((cell_df["x"] < x[offset]) | \
                                            (cell_df["x"] > x[-offset - 1])) | \
                                           ((cell_df["y"] < y[offset]) | \
                                            (cell_df["y"] > y[-offset - 1])) | \
                                           (cell_df["z"] < z[offset]) | \
                                           (cell_df["z"] > z[-offset - 1]), "id"].unique()
        except IndexError:
            offset_cells_ids = cell_df.loc[((cell_df["x"] < x[offset]) | \
                                            (cell_df["x"] > x[-offset - 1])) | \
                                           ((cell_df["y"] < y[offset]) | \
                                            (cell_df["y"] > y[-offset - 1])), "id"].unique()
    else:
        offset_cells_ids = []
    cells_ids = cell_df["id"]
    cells_ids = [x for x in cells_ids if x not in offset_cells_ids]

    cell_df = cell_df[cell_df["id"].isin(cells_ids)]

    cell_df["run"] = idx
    global_df["run"] = idx

    get_dataframes[idx] = [global_df, cell_df.sort_values(by="time")[startingPoint:stoppingPoint]]

global_df = pd.concat(
    (get_dataframes[x][0] for x in range(len(get_dataframes))))  # .groupby([x_axis], as_index=False).mean()
global_df.reset_index(inplace=True)
cell_df = pd.concat(
    (get_dataframes[x][1] for x in range(len(get_dataframes))))  # .groupby(["sigma"], as_index=False).mean()
cell_df.reset_index(inplace=True)

# adjust units
global_df["surf_c_std_norm"] = global_df["surf_c_std"]/global_df["surf_c"]
global_df["surf_c"] *= 1e3
global_df["surf_c_std"] *= 1e3

global_df["time"] = global_df["time"].div(36000)

global_df["time"] = global_df["time"].div(3600/10)
global_df = global_df.sort_values(by="run")

cell_df["time"] = cell_df["time"].div(3600)
cell_df["IL-2_surf_c"] = cell_df["IL-2_surf_c"].mul(1e3)


########################################################################################################################
print("plotting")
########################################################################################################################


if scan_variable == "IL-2_Tsec_fraction":
    cell_df = cell_df.loc[cell_df[scan_variable] < max_frac]
    global_df = global_df.loc[global_df[scan_variable] < max_frac]

for c, cell_type in enumerate(cell_types):

    tmp_df = cell_df.loc[cell_df["type_name"] == cell_type]
    surf_c = []
    for s, sv in enumerate(np.sort(tmp_df[scan_variable].unique())):
        surf_c.append(tmp_df.loc[tmp_df[scan_variable] == sv, "IL-2_surf_c"].mean())
    sns.lineplot(np.sort(tmp_df[scan_variable].unique()), surf_c,
                 color=y_colours[0], ci=0, label="mean", legend=plot_legend)


    tmp_df = cell_df.loc[cell_df["type_name"] == cell_type]
    surf_c_std_norm = []
    for s, sv in enumerate(np.sort(tmp_df[scan_variable].unique())):
        surf_c_std_norm.append(tmp_df.loc[tmp_df[scan_variable] == sv, "IL-2_surf_c"].std())
    sns.lineplot(np.sort(tmp_df[scan_variable].unique()), surf_c_std_norm,
                 color=y_colours[1], ci=0, label="SD", legend=plot_legend)



for c, cell_type in enumerate(cell_types):

    bars_df = cell_df.loc[cell_df["type_name"] == cell_type]
    ste = []
    mean = []
    for s, scan in enumerate(np.sort(bars_df[scan_variable].unique())):
        ste.append(bars_df.loc[np.abs(bars_df[scan_variable] - scan) < 0.0001, "IL-2_" + y_variable].std() / np.sqrt(len(global_df["run"].unique())))
        mean.append(bars_df.loc[np.abs(bars_df[scan_variable] - scan) < 0.0001,  "IL-2_" + y_variable].mean())
    plt.fill_between(np.sort(bars_df[scan_variable].unique()),
                     np.clip(np.array(mean) - np.array(ste), 0, None),
                     np.clip(np.array(mean) + np.array(ste), 0, None), alpha=0.15, color=y_colours[0])

    ste = []
    mean = []
    for s, scan in enumerate(np.sort(global_df[scan_variable].unique())):
        sliced_df = cell_df.loc[(cell_df["type_name"] == cell_type) & (cell_df[scan_variable] == scan)]
        cvs = []
        for r, run in enumerate(np.sort(global_df["run"].unique())):
            cvs.append(sliced_df.loc[(sliced_df["run"] == run), "IL-2_surf_c"].std())

        ste.append(np.std(cvs) / np.sqrt(len(global_df["run"].unique())))
        mean.append(np.mean(cvs))
    plt.fill_between(np.sort(global_df[scan_variable].unique()),
                     np.clip(np.array(mean) - np.array(ste), 0, None),
                     np.clip(np.array(mean) + np.array(ste), 0, None), alpha=0.15, color=y_colours[1])



if scan_variable == "IL-2_Tsec_fraction":
    xlabel = "fraction of sec. cells"
elif scan_variable == "IL-2_sigma":
    xlabel = "Receptor heterogeneity"


ax.set(ylabel="pM", xlabel=xlabel, yscale=yscale, xscale=xscale, ylim=ylim, xlim=xlim, title=title)

fig.savefig(saving_string, bbox_inches='tight')
plt.tight_layout()
plt.show()


