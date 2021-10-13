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

        activation_df = pd.read_hdf(load_path + 'activation_df.h5', mode="r")
        activation_df["run"] = idx
        get_dataframes[idx] = [activation_df[startingPoint:stoppingPoint]]

    pos_activation_df = pd.concat((get_dataframes[x][0] for x in range(len(get_dataframes))))#.groupby(["sigma"], as_index=False).mean()
    pos_activation_df.reset_index(inplace=True)
    return pos_activation_df

def EC50_calculation(E_max, E_min, k, N, R):
    return (E_max * k ** N + E_min * R ** N) / (k ** N + R ** N)

save_plot = True

hdd = "extra2" if os.path.exists("/extra2") else "extra"
user = getpass.getuser()
model_name = "Tsec_scan"
saving_string = "/{extra}/{u}/paper_models/kinetics/{mn}/".format(u=user, mn=model_name, extra = hdd)

plot_ODE = False
ODE_path = "/home/brunner/Documents/Current work/2021_09_10/ODE_Tsec_scan/"

yscale = "linear"
xscale = "linear"

# define data start and end point
startingPoint = None
stoppingPoint = None
# plotting limits
xlim = (None,0.86)
ylim = (None, None)
# time and fraction max limit
max_time = 100
max_frac = -1
# colours, first = positive, second = negative
colours = ["red", "blue"]
# opacity for the ODE
ODE_alpha = 0.25

sns.set(rc={'figure.figsize': (7, 5.5)})
sns.set_style("ticks")
sns.set_context("talk", font_scale=1.5, rc={"lines.linewidth": 2.5})
fig = plt.figure()

pos_path = "/{extra}/{u}/paper_models/kinetics/{mn}/{n}/".format(u=user, n="positive", mn=model_name, extra = hdd)
neg_path = "/{extra}/{u}/paper_models/kinetics/{mn}/{n}/".format(u=user, n="negative", mn=model_name, extra = hdd)

pos_activation_df = load_runs(pos_path, 1)
neg_activation_df = load_runs(neg_path, 1)

pos_activation_df = pos_activation_df[pos_activation_df.index < max_time]
neg_activation_df = neg_activation_df[neg_activation_df.index < max_time]

if plot_ODE == True:
    ODE_pos_path = ODE_path + "pos_lin_10/"
    ODE_neg_path = ODE_path + "neg_lin_10/"
    ODE_pos_activation_df = load_runs(ODE_pos_path, 10)
    ODE_neg_activation_df = load_runs(ODE_neg_path, 10)

    ODE_pos_activation_df = ODE_pos_activation_df[ODE_pos_activation_df.index < max_time]
    ODE_neg_activation_df = ODE_neg_activation_df[ODE_neg_activation_df.index < max_time]

for d,df in enumerate([pos_activation_df, neg_activation_df]):
    ste = []
    mean = []
    sliced_df = df.loc[df["index"] ==  df["index"].max()]
    for s, scan in enumerate(df.columns[1:max_frac]):
        ste.append(df[scan].std() / np.sqrt(len(df["run"].unique())))
        mean.append(df[scan].mean())
    plt.fill_between(np.sort(df.columns[1:max_frac].astype("float")), np.clip(np.array(mean) - np.array(ste), 0, 1),
                     np.clip(np.array(mean) + np.array(ste), 0, 1), alpha=0.15, color=colours[d])

    sns.scatterplot(x=np.sort(df.columns[1:max_frac].astype("float")), y=mean,
                    color=colours[d], s=30)
    sns.lineplot(x=np.sort(df.columns[1:max_frac].astype("float")), y=mean,
                    color=colours[d], linewidth = "2.5")

if plot_ODE == True:
    for d,df in enumerate([ODE_pos_activation_df, ODE_neg_activation_df]):
        ste = []
        mean = []
        sliced_df = df.loc[df["index"] == df["index"].max()]
        for s, scan in enumerate(df.columns[1:-1]):
            ste.append(df[scan].std() / np.sqrt(len(df["run"].unique())))
            mean.append(df[scan].mean())
        plt.fill_between(np.sort(df.columns[1:-1].astype("float")), np.clip(np.array(mean) - np.array(ste), 0, 1),
                         np.clip(np.array(mean) + np.array(ste), 0, 1), alpha=0.15, color=colours[d])

        sns.scatterplot(x=np.sort(df.columns[1:-1].astype("float")), y=mean,
                        color=colours[d], alpha = ODE_alpha, s=30)
        sns.lineplot(x=np.sort(df.columns[1:-1].astype("float")), y=mean,
                        color=colours[d], alpha = ODE_alpha, linewidth = "2.5")

import matplotlib.lines as mlines

white_line = mlines.Line2D([], [], color='white', alpha=1)
red_line = mlines.Line2D([], [], color='red')
blue_line = mlines.Line2D([], [], color='blue')
# red_line = mlines.Line2D([], [], color='red')
lightred_line = mlines.Line2D([], [], color='red', alpha=ODE_alpha)
lightblue_line = mlines.Line2D([], [], color='blue', alpha=ODE_alpha)
if plot_ODE == True:
    plt.legend([white_line, red_line, blue_line, white_line, lightred_line, lightblue_line],
               ["Spatial", "positive", "negative", "ODE", "positive", "negative"], loc='best',
               prop={'size': 17})
else:
    plt.legend([white_line, red_line, blue_line],
               ["Spatial", "positive", "negative"], loc='best',
               prop={'size': 17})

plt.ylabel("fraction of activated cells")
plt.xlabel("fraction of secreting cells")
plt.yscale(yscale)
plt.xscale(xscale)
plt.ylim(ylim)
plt.xlim(xlim)

plt.tight_layout()
if save_plot == True:
    fig.savefig(saving_string + "Tsec_scan.png", bbox_inches='tight')
plt.show()