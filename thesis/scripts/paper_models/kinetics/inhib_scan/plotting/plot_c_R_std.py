import getpass
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import seaborn as sns
import os

save_figure = True
hdd = "/extra2" if os.path.exists("/extra2") else "/extra"
user = getpass.getuser()
model_name = "inhib_scan"
saving_string = "/{extra}/{u}/paper_models/kinetics/{mn}/inhib_scan_c_R_std.png".format(u=user, mn=model_name, extra = hdd)

# run range
myRange = np.arange(0,1,1)
# plot positive or negative fb
name = "positive"
# name = "negative"

# which scan_index to plot
scan_index = 0

# plotting parameters
sns.set(rc={'figure.figsize': (23, 6)})
sns.set_style("ticks")
sns.set_context("talk", font_scale=1.5, rc={"lines.linewidth": 2.5})
# fig, ax = plt.subplots(figsize=(14, 6))
fig = plt.figure()
plt.subplots_adjust(wspace=.3)

show_legend = False

yscale = "log"
c_yscale = "log"
xscale = "linear"

# define data start and end point
startingPoint = None
stoppingPoint = None
# maximum time in h
time_limit = 200

xlim = (None,None) # lims for concentration c
ylim = (None, None) # lims for concentration c
ylimR = (1e1, 1e5) # ylim for receptors
ylimSD = (None, None) # ylim for SD

# whether to average or plot the SD
plot_every_cell = True
plot_std_area = False
# which cells to show
show_Tsecs = False
show_Tregs = False
show_Ths = True

plot_ODE = False
ODE_path = "/home/brunner/Documents/Current work/2021_07_16/ODE_kin_saturation/"

# defines a outer layer of N cells to be ignored in plotting. Used to further limit unwanted boundary effects.
offset = 0

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

    path = "/{extra}/{u}/paper_models/kinetics/{mn}/{n}/".format(u=user, n=name, mn=model_name, extra=hdd)
    path += "run" + str(idx) + "/"
    try:
        cell_df = pd.read_hdf(path + 'cell_df.h5', mode="r")
        print("read h5")
    except:
        cell_df = pd.read_pickle(path + 'cell_df.pkl')
        print("read pkl")
    # global_df = pd.read_hdf(path + 'global_df.h5', mode="r")

    if offset != 0:
        x, y, z = (cell_df["x"].unique(), cell_df["y"].unique(), cell_df["z"].unique())
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
        cells_ids = cell_df["id"]
        cells_ids = [x for x in cells_ids if x not in offset_cells_ids]

        cell_df = cell_df[cell_df["id"].isin(cells_ids)]
    else:
        offset_cells_ids = []

    try:
        cell_df["id_id"]
    except KeyError:
        cell_df["id_id"] = cell_df["cell_index"]
        cell_df["type_name"] = cell_df["cell_type"]
        cell_df["IL-2_surf_c"] = cell_df["IL2"] * 1e9
        cell_df["IL-2_R"] = cell_df["R"]
        cell_df["IL-2_gamma"] = cell_df["gamma"]


    if show_Tsecs == False:
        cell_df = cell_df.loc[cell_df["type_name"] != "Tsec"]
    if show_Tregs == False:
        cell_df = cell_df.loc[cell_df["type_name"] != "Treg"]
    if show_Ths == False:
        cell_df = cell_df.loc[cell_df["type_name"] != "Th"]


    cell_df  = cell_df.loc[cell_df["time"] < time_limit*3600]
    # global_df  = global_df.loc[global_df["time"] < time_limit*3600]
    get_dataframes[idx] = [cell_df.sort_values(by="time")[startingPoint:stoppingPoint]]



x_axis = "time"
my_hue = "scan_index"

if plot_every_cell == True:
    if len(get_dataframes) != 1:
        cell_df = pd.concat((get_dataframes[x][0] for x in range(len(get_dataframes)))).groupby(["time_index", my_hue, "id_id"], as_index=False).mean()
    else:
        cell_df = get_dataframes[0][0]
else:
    cell_df = pd.concat((get_dataframes[x][0] for x in range(
        len(get_dataframes))))  # .groupby(["time_index", "IL-2_gamma", "id_id"], as_index=False).mean()
cell_df.reset_index(inplace=True)

if plot_ODE == True:
    ODE_cell_df = pd.read_hdf(ODE_path + "cell_df.h5", mode="r")
    ODE_cell_df["time"] /= 3600
    ODE_cell_df = ODE_cell_df.loc[ODE_cell_df["time"] <= time_limit]
    ODE_cell_df["IL-2_surf_c"] = ODE_cell_df["IL-2_surf_c"] * 1e3

cell_df["time"] = cell_df["time"].div(3600)
cell_df["IL-2_surf_c"] = cell_df["IL-2_surf_c"].mul(1e3)

#%%
# scan_list = [10] #np.sort(cell_df[my_hue].unique())

for Tsec_scan_index in [scan_index]: # np.sort(cell_df[my_hue].unique())[::2]:
    scan_list = [Tsec_scan_index]
    print("plotting")

    sns.set(rc={'figure.figsize': (23, 6)})
    sns.set_style("ticks")
    sns.set_context("talk", font_scale=1.5, rc={"lines.linewidth": 2.5})
    # fig, ax = plt.subplots(figsize=(14, 6))
    fig = plt.figure()
    plt.subplots_adjust(wspace=.3)
    flat_scan_list = scan_list# [item for sublist in gammas_list for item in sublist]
    if all(np.array(flat_scan_list) < 1):
        cmap = "Blues"
    elif all(np.array(flat_scan_list) > 1):
        cmap = "Reds"
    else:
        cmap = "bwr"
    palette = sns.color_palette(cmap, len(flat_scan_list))
    ODE_palette = sns.color_palette("Greys", len(flat_scan_list))
    if cmap == "Blues":
        palette = [x for x in reversed(palette)]
    palette = ["blue"]

    timepoints = cell_df["time"].unique()[1:]

    fig.add_subplot(1,3, 1)
    for g, scan in enumerate(flat_scan_list):
        if scan < 1:
            label = "negative fb."
        elif scan > 1:
            label = "positive fb."
        else:
            label = ""

        if plot_every_cell == True:
            sns.set("talk", font_scale=1.5, rc={"lines.linewidth": 0.1})
            if show_Ths == True:
                ax_1 = sns.lineplot(x="time", y="IL-2_surf_c", data=cell_df.loc[(cell_df[my_hue] == scan) & (cell_df["type_name"] == "Th")].sort_values(by="time")[startingPoint:stoppingPoint], estimator=None,
                              units="id_id", color=palette[g])
            if show_Tregs == True:
                ax_1 = sns.lineplot(x="time", y="IL-2_surf_c", data=cell_df.loc[(cell_df[my_hue] == scan) & (cell_df["type_name"] == "Treg")].sort_values(by="time")[startingPoint:stoppingPoint], estimator=None,
                              units="id_id", color="black")

            if plot_ODE == True:
                sns.set("talk", font_scale=1.5, rc={"lines.linewidth": 1.5})
                ax_1 = sns.lineplot(x=x_axis, y="IL-2_surf_c",
                                    data=ODE_cell_df.loc[np.abs(ODE_cell_df[my_hue] - scan) < 0.001].sort_values(by="time")[
                                         startingPoint:stoppingPoint]
                                    , color=ODE_palette[g], alpha=0.3, label=label, ci=None, legend=show_legend)
        elif plot_ODE == True:
            ax_1 = sns.lineplot(x=x_axis, y="IL-2_surf_c",
                                data=cell_df.loc[cell_df[my_hue] == scan].sort_values(by="time")[
                                     startingPoint:stoppingPoint]
                                , color=palette[g], label=label, ci=None, legend=show_legend)
            ax_1 = sns.lineplot(x=x_axis, y="IL-2_surf_c",
                                data=ODE_cell_df.loc[np.abs(ODE_cell_df[my_hue] - scan) < 0.001].sort_values(by="time")[
                                     startingPoint:stoppingPoint]
                                , color=ODE_palette[g], alpha=0.3, label=label, ci=None, legend=show_legend)
        else:
            ax_1 = sns.lineplot(x=x_axis, y="IL-2_surf_c", data=cell_df.loc[cell_df[my_hue] == scan].sort_values(by="time")[startingPoint:stoppingPoint]
                         , color=palette[g], label=label, ci=0, legend=show_legend)
            # plt.fill_between(timepoints, 5,20, alpha=0.3)
            if plot_std_area == True:
                bars_df = cell_df.loc[cell_df[my_hue] == scan].sort_values(by="time")[startingPoint:stoppingPoint]
                std = []
                mean = []
                for t, time in enumerate(timepoints/3600):
                    std.append(bars_df.loc[np.abs(bars_df["time"] - time) < 0.0001,"IL-2_surf_c"].std())
                    mean.append(bars_df.loc[np.abs(bars_df["time"] - time) < 0.0001, "IL-2_surf_c"].mean())
                plt.fill_between(timepoints/3600, np.clip(np.array(mean) - np.array(std), 0, None), np.clip(np.array(mean) + np.array(std), 0, None), alpha=0.15, color=palette[g])

    sns.set_style("ticks")
    fig.add_subplot(1,3, 2)
    for g,scan in enumerate(flat_scan_list):
        if scan < 1:
            label = "negative fb."
        elif scan > 1:
            label = "positive fb."
        else:
            label = ""


        if plot_ODE == True:
            ax_2 = sns.lineplot(x=x_axis, y="IL-2_R",
                                data=cell_df.loc[cell_df[my_hue] == scan].sort_values(by="time")[
                                     startingPoint:stoppingPoint]
                                , color=palette[g], label=label, ci="sd", legend=show_legend)
            ax_2 = sns.lineplot(x=x_axis, y="IL-2_R",
                                data=ODE_cell_df.loc[np.abs(ODE_cell_df[my_hue] - scan) < 0.001].sort_values(by="time")[
                                     startingPoint:stoppingPoint]
                                , color=ODE_palette[g], alpha=0.3, label=label, ci=None, legend=show_legend)
        else:
            temp_cell_df = cell_df.loc[cell_df[my_hue] == scan].sort_values(by="time")[startingPoint:stoppingPoint]
            SD = []
            for t, time in enumerate(np.sort(temp_cell_df["time"].unique())):
                SD.append(temp_cell_df.loc[temp_cell_df["time_index"] == t + 1, "IL-2_surf_c"].values.std())
            SD_df = pd.DataFrame(columns=["time", "SD"])
            SD_df["SD"] = SD
            SD_df["time"] = np.sort(temp_cell_df["time"].unique())

            if SD_df["SD"].sum() > 1e-3:
                ax_2  = sns.lineplot(x="time", y="SD", data=SD_df
                             , color=palette[g], label=label, ci=0, legend=show_legend)

            # ax_2 = sns.lineplot(x=x_axis, y="surf_c_std", data=global_df.loc[global_df[my_hue] == scan].sort_values(by="time")[startingPoint:stoppingPoint]
            #              , color=palette[g], label=label, ci=0, legend=show_legend)
            # plt.fill_between(timepoints, 5,20, alpha=0.3)


    sns.set_style("ticks")
    fig.add_subplot(1,3, 3)
    for g,scan in enumerate(flat_scan_list):
        if scan < 1:
            label = "negative fb."
        elif scan > 1:
            label = "positive fb."
        else:
            label = ""

        if plot_every_cell == True:
            sns.set("talk", font_scale=1.5, rc={"lines.linewidth": 0.1})
            if show_Ths == True:
                ax_3 = sns.lineplot(x="time", y="IL-2_R", data=cell_df.loc[(cell_df[my_hue] == scan) & (cell_df["type_name"] == "Th")].sort_values(by="time")[startingPoint:stoppingPoint], estimator=None,
                              units="id_id", color=palette[g])
            if show_Tregs == True:
                ax_3 = sns.lineplot(x="time", y="IL-2_R", data=cell_df.loc[(cell_df[my_hue] == scan) & (cell_df["type_name"] == "Treg")].sort_values(by="time")[startingPoint:stoppingPoint], estimator=None,
                              units="id_id", color="black")
            if plot_ODE == True:
                sns.set("talk", font_scale=1.5, rc={"lines.linewidth": 1.5})
                ax_3 = sns.lineplot(x=x_axis, y="IL-2_R",
                                data=ODE_cell_df.loc[np.abs(ODE_cell_df[my_hue] - scan) < 0.001].sort_values(by="time")[
                                     startingPoint:stoppingPoint]
                                , color=ODE_palette[g], alpha=0.3, label=label, ci=None, legend=show_legend)

        elif plot_ODE == True:
            ax_3 = sns.lineplot(x=x_axis, y="IL-2_R",
                                data=cell_df.loc[cell_df[my_hue] == scan].sort_values(by="time")[
                                     startingPoint:stoppingPoint]
                                , color=palette[g], label=label, ci=None, legend=show_legend)
            ax_3 = sns.lineplot(x=x_axis, y="IL-2_R",
                                data=ODE_cell_df.loc[np.abs(ODE_cell_df[my_hue] - scan) < 0.001].sort_values(by="time")[
                                     startingPoint:stoppingPoint]
                                , color=ODE_palette[g], alpha=0.3, label=label, ci=None, legend=show_legend)
        else:
            ax_3 = sns.lineplot(x=x_axis, y="IL-2_R", data=cell_df.loc[cell_df[my_hue] == scan].sort_values(by="time")[startingPoint:stoppingPoint]
                         , color=palette[g], label=label, ci=0, legend=show_legend)
            # plt.fill_between(timepoints, 5,20, alpha=0.3)
            if plot_std_area == True:
                bars_df = cell_df.loc[cell_df[my_hue] == scan].sort_values(by="time")[startingPoint:stoppingPoint]
                std = []
                mean = []
                for t, time in enumerate(timepoints/3600):
                    std.append(bars_df.loc[np.abs(bars_df["time"] - time) < 0.0001,"IL-2_R"].std())
                    mean.append(bars_df.loc[np.abs(bars_df["time"] - time) < 0.0001, "IL-2_R"].mean())
                plt.fill_between(timepoints/3600, np.clip(np.array(mean) - np.array(std), 0, None), np.clip(np.array(mean) + np.array(std), 0, None), alpha=0.15, color=palette[g])


    # import matplotlib.lines as mlines
    #
    # white_line = mlines.Line2D([], [], color='white', alpha=0)
    # blue_line = mlines.Line2D([], [], color='blue')
    # lightblue_line = mlines.Line2D([], [], color='blue', alpha=0.3)
    # red_line = mlines.Line2D([], [], color='red')
    # lightred_line = mlines.Line2D([], [], color='red', alpha=0.3)
    # plt.legend([white_line, blue_line, lightblue_line, white_line, red_line, lightred_line],
    #                   ["Negative fb.", "Spatial", "ODE", "Positive fb.", "Spatial", "ODE"], loc='upper center',
    #                   bbox_to_anchor=(0.54, 0.4), ncol=2, prop={'size': 17})  # (loc=(0.385,0.2))
    # plt.legend(handles=new_handles, labels=new_labels)

    ax_1.set(ylabel="pM", xlabel="time (h)", yscale=c_yscale, xscale=xscale, ylim=ylim, xlim=xlim, title="Surf. c.")
    try:
        ax_2.set(ylabel="pM", xlabel="time (h)", yscale=yscale, xscale=xscale, ylim=ylimSD, xlim=xlim, title="SD")
    except:
        pass
    ax_3.set(ylabel="Molecules", xlabel="time (h)", yscale=yscale, xscale=xscale, ylim=ylimR, xlim=xlim, title="Receptors")
    # ax.set(ylabel="(pM)", xlabel="time (h)", yscale=yscale, xscale=xscale, ylim=ylim, xlim=xlim, title="surface c.")

    if save_figure == True:
        fig.savefig(saving_string[:-4] + "_" + str(Tsec_scan_index) + saving_string[-4:], bbox_inches='tight')
    plt.tight_layout()
    plt.show()
