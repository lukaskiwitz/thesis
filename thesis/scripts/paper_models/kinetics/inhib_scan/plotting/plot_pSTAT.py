import getpass
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import seaborn as sns
import os
from thesis.main.my_debug import message

def load_runs(path, no_of_runs):
    myRange = np.arange(0,no_of_runs,1)
    get_dataframes = []#[[]]*3
    for idx,value in enumerate(myRange):
        get_dataframes.append([])
        print("loading run" + str(value))
        load_path = path + "/run" + str(value) + "/"

        cell_df = pd.read_hdf(load_path + 'cell_df.h5', mode="r")
        cell_df["run"] = idx
        get_dataframes[idx] = [cell_df[startingPoint:stoppingPoint]]

    cell_df = pd.concat((get_dataframes[x][0] for x in range(len(get_dataframes))))#.groupby(["sigma"], as_index=False).mean()
    cell_df.reset_index(inplace=True)
    return cell_df

def EC50_calculation(E_max, E_min, k, N, R):
    return (E_max * k ** N + E_min * R ** N) / (k ** N + R ** N)

save_plot = True

hdd = "extra2" if os.path.exists("/extra2") else "extra"
user = getpass.getuser()
model_name = "inhib_scan"
saving_string = "/{extra}/{u}/paper_models/kinetics/{mn}/pSTAT5.png".format(u=user, mn=model_name, extra = hdd)

sns.set(rc={'figure.figsize': (7, 6)})
sns.set_style("ticks")
sns.set_context("talk", font_scale=1.5, rc={"lines.linewidth": 2.5})
fig = plt.figure()

plot_legend = True

yscale = "log"
xscale = "linear"
# define data start and end point
startingPoint = None
stoppingPoint = None
# plotting limits
xlim = (None,None)
ylim = (None, None)

plot_every_cell = True
show_Tsecs = False
show_Tregs = False
show_Ths = True

plot_std_area = True
# define which scan to plot
x_axis = "time"
my_hue = "scan_index"
scan = 0
# the first time index does not have pSTAT calculated yet.
skip_first_time_index = True


########################################################################################################################
########################################################################################################################
##                                              Plotting                                                              ##
########################################################################################################################
########################################################################################################################

pos_path = "/{extra}/{u}/paper_models/kinetics/{mn}/{n}/".format(u=user, n="positive", mn=model_name, extra = hdd)
neg_path = "/{extra}/{u}/paper_models/kinetics/{mn}/{n}/".format(u=user, n="positive", mn=model_name, extra = hdd)

pos_cell_df = load_runs(pos_path, 1)
neg_cell_df = load_runs(neg_path, 1)

pos_cell_df["time"] = pos_cell_df["time"].div(3600)
neg_cell_df["time"] = neg_cell_df["time"].div(3600)
pos_cell_df["IL-2_surf_c"] = pos_cell_df["IL-2_surf_c"].mul(1e3)
neg_cell_df["IL-2_surf_c"] = neg_cell_df["IL-2_surf_c"].mul(1e3)

if skip_first_time_index == True:
    pos_cell_df = pos_cell_df.loc[pos_cell_df["time_index"] != 0]
    neg_cell_df = neg_cell_df.loc[neg_cell_df["time_index"] != 0]

print("plotting")
for cell_df in [pos_cell_df, neg_cell_df]:
    try:
        cell_df["pSTAT5"] = cell_df["misc_pSTAT5"]
    except:
        cell_df["pSTAT5"] = cell_df["IL-2_surf_c"] ** 3 / (
                (EC50_calculation(E_max=40e-12, E_min=0, k=860, N=0.55, R=cell_df["IL-2_R"]) * 1e12) ** 3 + cell_df[
                        "IL-2_surf_c"] ** 3).values


    sns.set(rc={'figure.figsize': (7, 6)})
    sns.set_style("ticks")
    sns.set_context("talk", font_scale=1.5, rc={"lines.linewidth": 2.5})
    # fig, ax = plt.subplots(figsize=(14, 6))
    fig = plt.figure()
    plt.subplots_adjust(wspace=.3)

    palette = ["Grey"]

    timepoints = cell_df["time"].unique()[1:]

    sns.set_style("ticks")

    if cell_df["IL-2_gamma"].unique()[0] < 1:
        label = "negative fb."
    elif cell_df["IL-2_gamma"].unique()[0] > 1:
        label = "positive fb."
    else:
        label = ""

    if plot_every_cell == True:
        sns.set("talk", font_scale=1.5, rc={"lines.linewidth": 0.1})
        if show_Ths == True:
            ax_3 = sns.lineplot(x="time", y="pSTAT5", data=cell_df.loc[(cell_df["scan_index"] == scan)].sort_values(by="time")[startingPoint:stoppingPoint], estimator=None,
                          units="id_id", color=palette[0])
        if show_Tregs == True:
            ax_3 = sns.lineplot(x="time", y="pSTAT5", data=cell_df.loc[(cell_df["scan_index"] == scan) & (cell_df["type_name"] == "Treg")].sort_values(by="time")[startingPoint:stoppingPoint], estimator=None,
                          units="id_id", color="black")

    else:
        ax_3 = sns.lineplot(x=x_axis, y="pSTAT5", data=cell_df.loc[cell_df[my_hue] == scan].sort_values(by="time")[startingPoint:stoppingPoint]
                     , color=palette[0], label=label, ci=0, legend=plot_legend)
        # plt.fill_between(timepoints, 5,20, alpha=0.3)
        if plot_std_area == True:
            bars_df = cell_df.loc[cell_df[my_hue] == scan].sort_values(by="time")[startingPoint:stoppingPoint]
            std = []
            mean = []
            for t, time in enumerate(timepoints/3600):
                std.append(bars_df.loc[np.abs(bars_df["time"] - time) < 0.0001,"pSTAT5"].std())
                mean.append(bars_df.loc[np.abs(bars_df["time"] - time) < 0.0001, "pSTAT5"].mean())
            plt.fill_between(timepoints/3600, np.clip(np.array(mean) - np.array(std), 0, None), np.clip(np.array(mean) + np.array(std), 0, None), alpha=0.15, color=palette[0])

    ax_3.set(ylabel="pSTAT", xlabel="time (h)", yscale="linear", xscale=xscale, ylim=(-0.1, 1.1), xlim=xlim, title="pSTAT")
    sns.set("talk", font_scale=1.5, rc={"lines.linewidth": 2})
    plt.axhline(0.5, 0, 200, color="black")


    if save_plot == True:
        fig.savefig(saving_string, bbox_inches='tight')
    plt.tight_layout()
    plt.show()
