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

        cell_df = pd.read_hdf(load_path + 'activation_df.h5', mode="r")
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
model_name = "Tsec_scan"
saving_string = "/{extra}/{u}/paper_models/kinetics/{mn}/".format(u=user, mn=model_name, extra = hdd)

plot_ODE = False
ODE_path = "/home/brunner/Documents/Current work/2021_09_10/ODE_Tsec_scan/"

yscale = "log"
xscale = "linear"

# define data start and end point
startingPoint = None
stoppingPoint = None
# plotting limits
xlim = (None,0.86)
ylim = (None, None)

plot_every_cell = True
show_Tsecs = False
show_Tregs = False
show_Ths = True
save_figure = False

hdd = "extra2" if os.path.exists("/extra2") else "extra"
user = getpass.getuser()
model_name = "Tsec_scan_pSTAT"
saving_string = "/{extra}/{u}/paper_models/kinetics/{mn}/".format(u=user, mn=model_name, extra = hdd)


sns.set(rc={'figure.figsize': (7, 6)})
sns.set_style("ticks")
sns.set_context("talk", font_scale=1.5, rc={"lines.linewidth": 2.5})
fig = plt.figure()

pos_path = "/{extra}/{u}/paper_models/kinetics/{mn}/{n}/".format(u=user, n="positive", mn=model_name, extra = hdd)
neg_path = "/{extra}/{u}/paper_models/kinetics/{mn}/{n}/".format(u=user, n="negative", mn=model_name, extra = hdd)

pos_cell_df = load_runs(pos_path, 1)
neg_cell_df = load_runs(neg_path, 1)

if plot_ODE == True:
    ODE_pos_path = ODE_path + "pos_lin_10/"
    ODE_neg_path = ODE_path + "neg_lin_10/"
    ODE_pos_activation_df = load_runs(ODE_pos_path, 10)
    ODE_neg_activation_df = load_runs(ODE_neg_path, 10)


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
        if plot_ODE == True:
            sns.set("talk", font_scale=1.5, rc={"lines.linewidth": 1.5})
            ax_3 = sns.lineplot(x=x_axis, y="pSTAT5",
                            data=ODE_cell_df.loc[np.abs(ODE_cell_df["scan_index"] - scan) < 0.001].sort_values(by="time")[
                                 startingPoint:stoppingPoint]
                            , color=ODE_palette[0], alpha=0.3, label=label, ci=None, legend=show_legend)

    elif plot_ODE == True:
        ax_3 = sns.lineplot(x=x_axis, y="pSTAT5",
                            data=cell_df.loc[cell_df[my_hue] == scan].sort_values(by="time")[
                                 startingPoint:stoppingPoint]
                            , color=palette[g], label=label, ci=None, legend=show_legend)
        ax_3 = sns.lineplot(x=x_axis, y="pSTAT5",
                            data=ODE_cell_df.loc[np.abs(ODE_cell_df[my_hue] - scan) < 0.001].sort_values(by="time")[
                                 startingPoint:stoppingPoint]
                            , color=ODE_palette[g], alpha=0.3, label=label, ci=None, legend=show_legend)
    else:
        ax_3 = sns.lineplot(x=x_axis, y="pSTAT5", data=cell_df.loc[cell_df[my_hue] == scan].sort_values(by="time")[startingPoint:stoppingPoint]
                     , color=palette[g], label=label, ci=0, legend=show_legend)
        # plt.fill_between(timepoints, 5,20, alpha=0.3)
        if plot_std_area == True:
            bars_df = cell_df.loc[cell_df[my_hue] == scan].sort_values(by="time")[startingPoint:stoppingPoint]
            std = []
            mean = []
            for t, time in enumerate(timepoints/3600):
                std.append(bars_df.loc[np.abs(bars_df["time"] - time) < 0.0001,"pSTAT5"].std())
                mean.append(bars_df.loc[np.abs(bars_df["time"] - time) < 0.0001, "pSTAT5"].mean())
            plt.fill_between(timepoints/3600, np.clip(np.array(mean) - np.array(std), 0, None), np.clip(np.array(mean) + np.array(std), 0, None), alpha=0.15, color=palette[g])

    ax_3.set(ylabel="pSTAT", xlabel="time (h)", yscale="linear", xscale=xscale, ylim=(-0.1, 1.1), xlim=xlim, title="pSTAT")
    sns.set("talk", font_scale=1.5, rc={"lines.linewidth": 2})
    plt.axhline(0.5, 0, 200, color="black")


    if save_figure == True:
        fig.savefig(saving_string[:-4] + "_" + str(scan_list[0]) + saving_string[-4:], bbox_inches='tight')
    plt.tight_layout()
    plt.show()
