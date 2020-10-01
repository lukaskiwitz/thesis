import getpass
import os

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import seaborn as sns

get_dataframes = []#[[]]*3
saving_dataframe = pd.DataFrame()

get_dataframes.append([])
# path = "/extra/brunner/10x10x10/R_lognorm/run" + str(j) + "/"
# path = "/extra/brunner/10x10x10/q_fraction_exact/run" + str(j) + "/"
user = getpass.getuser()
path = "/extra/brunner/thesis/kinetic/q_fraction_medium_g_0.1/"
# path = "/extra/brunner/thesis/kinetic/q_fraction_gamma_scan/"
# ext_cache="/extra/brunner/para_handling/kinetic/R_lognorm_ext_cache/"

T = np.arange(0, 200, 1)
dt = 3600

global_df = pd.read_hdf(path + 'global_df.h5', mode="r")
cell_df = pd.read_hdf(path + 'cell_df.h5', mode="r")

try:
    global_df["surf_c_std"]
except KeyError:
    print("calculating surf_c_std")
    global_df["surf_c_std"] = np.zeros(len(global_df["time"]))
    for scan_index in global_df["scan_index"].unique():
        for t in global_df["time"].unique():
            print("time:", t, "in scan", scan_index)
            global_df.loc[(global_df["scan_index"] == scan_index) & (global_df["time"] == t), "surf_c_std"] = cell_df.loc[(cell_df["scan_index"] == scan_index) & (cell_df["time"] == t), "IL-2_surf_c"].std()
    global_df.to_hdf(path + 'global_df.h5', key="data", mode="w")

# cell_stats_df = pd.read_hdf(path + 'cell_stats_dataframe_' + str(j) + '.h5', mode="r")


# sns.lineplot(x="time", y="Concentration", data=global_df.loc[global_df["IL-2_gamma"] == 2])
# plt.show()
# exit()
#8.93

get_dataframes[0] = [global_df, cell_df]#, cell_stats_df]


my_hue = "IL-2_gamma" #"D"
x_axis = "time"
averaging_values = np.sort(global_df[x_axis].unique())

global_df = pd.concat((get_dataframes[x][0] for x in range(len(get_dataframes))))#.groupby([x_axis], as_index=False).mean()
global_df.reset_index(inplace=True)
cell_df = pd.concat((get_dataframes[x][1] for x in range(len(get_dataframes))))#.groupby(["sigma"], as_index=False).mean()
# cell_df.reset_index(inplace=True)
# cell_stats_df = pd.concat((get_dataframes[x][2] for x in range(len(get_dataframes)))).groupby(["sigma"], as_index=False).mean()
# cell_stats_df.reset_index(inplace=True)


sns.set(rc={'figure.figsize':(11,8.27)})
fig = plt.figure()
plt.subplots_adjust(wspace=.3)
# sns.set_style("white")

# global_df = global_df.drop(0)
# global_df[:49]["time"] = global_df[:49]["time"] - 1

# subplot_style
a_x = 4
a_y = 3

#global_df["D"] /= ((10 ** 0.5 * 0.01) ** 2)

# global_df["IL-2_gamma"] = ["$%s$" % x for x in global_df["IL-2_gamma"]]
global_df["time"] = global_df["time"].mul(3600)
# global_df = global_df.loc[global_df["IL-2_gamma"] == 0.1]


print("averaging receptor levels")
R_mean_df_list = []
for j, value in enumerate(averaging_values):
    for my_hue_value in cell_df[my_hue].unique():
        a_cell_df = cell_df.loc[(cell_df[x_axis] == j) & (cell_df[my_hue] == my_hue_value)]
        R_mean_df = pd.DataFrame(columns=["scan_index", my_hue, x_axis, "R_mean"])
        # if j == 0:
        #     R_mean_df = R_mean_df.append(
        #         {"time": 0.0, "R_mean": 7955,
        #          "std": 11000.0,
        #          my_hue: my_hue_value}, ignore_index=True)
        # for i in range(len(averaging_values)):
        R_mean_df = R_mean_df.append(
            {x_axis: value+1, "R_mean": a_cell_df["IL-2_R"].mean(),"std": a_cell_df["IL-2_R"].std(),
             "D": a_cell_df["IL-2_D"].unique()[0],
             "scan_index": a_cell_df["scan_index"].unique()[0], my_hue: my_hue_value}, ignore_index=True)
        R_mean_df_list.append(R_mean_df)
R_mean_concat = pd.concat(R_mean_df_list[x] for x in range(len(R_mean_df_list)))
R_mean_plotting = R_mean_concat


# global_df[my_hue] /= round(((10 ** 0.5 * 0.01) ** 2 * 0.1),6)
global_df[my_hue] = ["$%s$" % x for x in global_df[my_hue]]
global_df["time"] = global_df["time"].div(3600/10)
cell_df["time"] = cell_df["time"].mul(10)
R_mean_plotting["time"] = R_mean_plotting["time"].mul(10)

# cell_df[my_hue] = ["$%s$" % x for x in cell_df[my_hue]]

R_mean_plotting[my_hue] = ["$%s$" % x for x in R_mean_plotting[my_hue]]
# R_mean_plotting["t"] = R_mean_plotting["t"].div(3600)
# global_df["sigma"] = global_df["sigma"].div(20000)
# R_mean_plotting["sigma"] = R_mean_plotting["sigma"].div(20000)

path = "/extra/brunner/thesis/kinetic/q_fraction_medium_g_0.1/"


plot_D = 0
stoppingPoint = -30
startingPoint = None

yscale = "linear"
xscale = "linear"

yscaleR = "linear"

# ylim = (1.0, 1.55) #(0.01, 0.045)
ylim = (None, 33)
xlim = (-1,50)

ylimR = (None, None) #(1, 23000)

hue_order = [10.0,0.5,0.1]
hue_order = ["$%s$" % x for x in hue_order]
# hue_order = None
# fig = plt.figure()
sns.set_style("ticks")
# fig = plt.figure(figsize=(8,8))
# sns.set_context("talk", font_scale=1.1, rc={"lines.linewidth": 2.5})
# sns.set(rc={'figure.figsize':(8,12)})
# sns.set_style("ticks")
# sns.set(rc={'figure.figsize':(6,5)})
sns.set_style("ticks")
sns.set_context("talk", font_scale=1.5, rc={"lines.linewidth": 5})
fig,ax = plt.subplots(figsize=(6,5))



# fig.add_subplot(a_x,a_y, 8)
global_df["surf_c_std_norm"] = global_df["surf_c_std"]/global_df["surf_c"]
global_df.loc[(global_df["time_index"] == 0) & (global_df["IL-2_gamma"] == "$2.0$"), "surf_c_std_norm"] = global_df.loc[(global_df["time_index"] == 0) & (global_df["IL-2_gamma"] == "$0.5$"), "surf_c_std_norm"].values[0]
global_df.loc[(global_df["time_index"] == 0) & (global_df["IL-2_gamma"] == "$0.1$"), "surf_c_std_norm"] = global_df.loc[(global_df["time_index"] == 0) & (global_df["IL-2_gamma"] == "$0.5$"), "surf_c_std_norm"].values[0]
global_df.loc[(global_df["time_index"] == 0) & (global_df["IL-2_gamma"] == "$10.0$"), "surf_c_std_norm"] = global_df.loc[(global_df["time_index"] == 0) & (global_df["IL-2_gamma"] == "$0.5$"), "surf_c_std_norm"].values[0]
# insert data from same sim just with finer time steps to even out the overshooting. todo: finer timesteps for whole sim!
global_df.loc[(global_df["IL-2_gamma"] == '$10.0$') & (global_df["time"] == 10), "surf_c"] = 8.93*1e-3
global_df.loc[(global_df["IL-2_gamma"] == '$10.0$') & (global_df["time"] == 10), "surf_c_std_norm"] = 0.00751/(8.93*1e-3)
# print(global_df["surf_c_std_norm"])
global_df["surf_c"] *= 1e3

palette = sns.color_palette("bwr", 4)

# ax7 = sns.lineplot(x=x_axis, y="surf_c_std_norm", data=global_df.sort_values(by="time")[startingPoint:stoppingPoint], hue=my_hue, hue_order = hue_order)
# sns.lineplot(ax=ax, x=x_axis, y="surf_c_std_norm", data=global_df.loc[global_df["IL-2_gamma"] == '$10.0$'].sort_values(by="time")[startingPoint:stoppingPoint], color=palette[0])
# sns.lineplot(ax=ax, x=x_axis, y="surf_c_std_norm", data=global_df.loc[global_df["IL-2_gamma"] == '$2.0$'].sort_values(by="time")[startingPoint:stoppingPoint], color=palette[1])
# sns.lineplot(ax=ax, x=x_axis, y="surf_c_std_norm", data=global_df.loc[global_df["IL-2_gamma"] == '$0.5$'].sort_values(by="time")[startingPoint:stoppingPoint], color=palette[2])
# sns.lineplot(ax=ax, x=x_axis, y="surf_c_std_norm", data=global_df.loc[global_df["IL-2_gamma"] == '$0.1$'].sort_values(by="time")[startingPoint:stoppingPoint], color=palette[3])

# plt.hlines(y=global_df.loc[(global_df["time_index"] == 0), "surf_c_std_norm"][0],xmin=0,xmax=600)
# ax.set(xlabel="time in h", ylabel="std./mean", yscale=yscale, xscale=xscale, ylim=ylim, xlim=xlim, xticks = [0, 50, 100,150,200])
#
# ax2 =ax.twinx()
sns.lineplot(x=x_axis, y="surf_c", data=global_df.loc[global_df["IL-2_gamma"] == '$10.0$'].sort_values(by="time")[startingPoint:stoppingPoint], color=palette[0])
sns.lineplot(x=x_axis, y="surf_c", data=global_df.loc[global_df["IL-2_gamma"] == '$2.0$'].sort_values(by="time")[startingPoint:stoppingPoint], color=palette[1])
sns.lineplot(x=x_axis, y="surf_c", data=global_df.loc[global_df["IL-2_gamma"] == '$0.5$'].sort_values(by="time")[startingPoint:stoppingPoint], color=palette[2])
sns.lineplot(x=x_axis, y="surf_c", data=global_df.loc[global_df["IL-2_gamma"] == '$0.1$'].sort_values(by="time")[startingPoint:stoppingPoint], color=palette[3])
plt.hlines(y=global_df.loc[(global_df["time_index"] == 0), "surf_c"][0],xmin=0,xmax=600)
ax.set(ylabel="surface c. (pM)", xlabel="time (h)", yscale=yscale, xscale=xscale, ylim=ylim, xlim=xlim, xticks = [0, 50, 100,150,200])
# [ax2.lines[i].set_linestyle("--") for i,line in enumerate(ax2.lines)]
# plt.ylim((None,10))
# ax2.spines['right'].set_color('white')




# plt.xticks([0, 50, 100,150])
# ax7.set(xlabel="time in h", ylabel="surface concentration in pM", yscale=yscale, xscale=xscale, ylim=ylim, xticks = [0, 50, 100,150,200])
# plt.legend(bbox_to_anchor=(1.04,1), loc="upper left")
# handles, labels = ax7.get_legend_handles_labels()
# labels[0] = "$\gamma$"

import matplotlib.lines as mlines

zero_line = mlines.Line2D([], [], color=palette[0], label='lightcoral line')
one_line = mlines.Line2D([], [], color=palette[1], label='red line')
black_line = mlines.Line2D([], [], color='black', label='Black line')
two_line = mlines.Line2D([], [], color=palette[2], label='lightsteelblue line')
three_line = mlines.Line2D([], [], color=palette[3], label='blue line')

grey_line = mlines.Line2D([], [], color="grey", label='grey line')
grey_dashed_line = mlines.Line2D([], [], color="grey", label='grey dashed line', linestyle="dashed")
white_line = mlines.Line2D([], [], color="white", label='white line')


# new_handles = [grey_line,  grey_dashed_line, white_line,three_line, black_line, zero_line]
# new_labels = ["c$_v$", "concentration", "", "neg. feedback", "no feedback", "pos. feedback"]

new_handles = [three_line, black_line, zero_line]
new_labels = ["neg. feedback", "no feedback", "pos. feedback"]

# new_handles = [three_line, two_line, black_line, one_line, zero_line]
# new_labels = ["strong negative feedback", "negative feedback", "no feedback", "positive feedback", "strong positive feedback"]
plt.legend(handles=new_handles, labels=new_labels, loc='upper right', fancybox=True)
fig.savefig("plots/" + "kin_q_c" + ".svg", bbox_inches='tight')
plt.show()
