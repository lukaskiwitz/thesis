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
path = "/extra/brunner/thesis/kinetic/q_fraction_gamma_scan/"
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
    # global_df.to_hdf(path + 'global_df.h5', key="data", mode="w")
path = "/extra/brunner/thesis/kinetic/isotropie1/"
global_df1 = pd.read_hdf(path + 'global_df.h5', mode="r")
cell_df1 = pd.read_hdf(path + 'cell_df.h5', mode="r")

try:
    global_df1["surf_c_std"]
except KeyError:
    print("calculating surf_c_std")
    global_df1["surf_c_std"] = np.zeros(len(global_df1["time"]))
    for scan_index in global_df1["scan_index"].unique():
        for t in global_df1["time"].unique():
            print("time:", t, "in scan", scan_index)
            global_df1.loc[(global_df1["scan_index"] == scan_index) & (global_df1["time"] == t), "surf_c_std"] = cell_df1.loc[(cell_df1["scan_index"] == scan_index) & (cell_df1["time"] == t), "IL-2_surf_c"].std()
    # global_df.to_hdf(path + 'global_df.h5', key="data", mode="w")
path = "/extra/brunner/thesis/kinetic/isotropie2/"
global_df2 = pd.read_hdf(path + 'global_df.h5', mode="r")
cell_df2 = pd.read_hdf(path + 'cell_df.h5', mode="r")

try:
    global_df2["surf_c_std"]
except KeyError:
    print("calculating surf_c_std")
    global_df2["surf_c_std"] = np.zeros(len(global_df2["time"]))
    for scan_index in global_df2["scan_index"].unique():
        for t in global_df2["time"].unique():
            print("time:", t, "in scan", scan_index)
            global_df2.loc[(global_df2["scan_index"] == scan_index) & (global_df2["time"] == t), "surf_c_std"] = cell_df2.loc[(cell_df2["scan_index"] == scan_index) & (cell_df2["time"] == t), "IL-2_surf_c"].std()
    # global_df.to_hdf(path + 'global_df.h5', key="data", mode="w")
path = "/extra/brunner/thesis/kinetic/isotropie3/"
global_df3 = pd.read_hdf(path + 'global_df.h5', mode="r")
cell_df3 = pd.read_hdf(path + 'cell_df.h5', mode="r")

try:
    global_df3["surf_c_std"]
except KeyError:
    print("calculating surf_c_std")
    global_df3["surf_c_std"] = np.zeros(len(global_df3["time"]))
    for scan_index in global_df3["scan_index"].unique():
        for t in global_df3["time"].unique():
            print("time:", t, "in scan", scan_index)
            global_df3.loc[(global_df3["scan_index"] == scan_index) & (global_df3["time"] == t), "surf_c_std"] = cell_df3.loc[(cell_df3["scan_index"] == scan_index) & (cell_df3["time"] == t), "IL-2_surf_c"].std()
    # global_df.to_hdf(path + 'global_df.h5', key="data", mode="w")

# cell_stats_df = pd.read_hdf(path + 'cell_stats_dataframe_' + str(j) + '.h5', mode="r")


# sns.lineplot(x="time", y="Concentration", data=global_df.loc[global_df["IL-2_gamma"] == 2])
# plt.show()
# exit()
get_dataframes1 = [[]]
get_dataframes2 = [[]]
get_dataframes3 = [[]]

get_dataframes[0] = [global_df, cell_df]#, cell_stats_df]
get_dataframes1[0] = [global_df1, cell_df1]
get_dataframes2[0] = [global_df2, cell_df2]
get_dataframes3[0] = [global_df3, cell_df3]

my_hue = "IL-2_gamma" #"D"
x_axis = "time"
averaging_values = np.sort(global_df[x_axis].unique())

global_df = pd.concat((get_dataframes[x][0] for x in range(len(get_dataframes))))#.groupby([x_axis], as_index=False).mean()
global_df.reset_index(inplace=True)
cell_df = pd.concat((get_dataframes[x][1] for x in range(len(get_dataframes))))#.groupby(["sigma"], as_index=False).mean()
# cell_df.reset_index(inplace=True)
# cell_stats_df = pd.concat((get_dataframes[x][2] for x in range(len(get_dataframes)))).groupby(["sigma"], as_index=False).mean()
# cell_stats_df.reset_index(inplace=True)
global_df1 = pd.concat((get_dataframes1[x][0] for x in range(len(get_dataframes1))))#.groupby([x_axis], as_index=False).mean()
global_df1.reset_index(inplace=True)
cell_df1 = pd.concat((get_dataframes1[x][1] for x in range(len(get_dataframes1))))
global_df2 = pd.concat((get_dataframes2[x][0] for x in range(len(get_dataframes2))))#.groupby([x_axis], as_index=False).mean()
global_df2.reset_index(inplace=True)
cell_df2 = pd.concat((get_dataframes2[x][1] for x in range(len(get_dataframes2))))
global_df3 = pd.concat((get_dataframes3[x][0] for x in range(len(get_dataframes3))))#.groupby([x_axis], as_index=False).mean()
global_df3.reset_index(inplace=True)
cell_df3 = pd.concat((get_dataframes3[x][1] for x in range(len(get_dataframes3))))

sns.set(rc={'figure.figsize':(11,8.27)})
fig = plt.figure()
plt.subplots_adjust(wspace=.3)
sns.set_style("white")

# subplot_style
a_x = 4
a_y = 3

#global_df["D"] /= ((10 ** 0.5 * 0.01) ** 2)

global_df["IL-2_gamma"] = ["$%s$" % x for x in global_df["IL-2_gamma"]]
global_df["time"] = global_df["time"].mul(10)
global_df = global_df.loc[global_df["IL-2_k_factor"] == 2]

global_df1["IL-2_gamma"] = ["$%s$" % x for x in global_df1["IL-2_gamma"]]
global_df1["time"] = global_df1["time"].mul(10)
global_df2["IL-2_gamma"] = ["$%s$" % x for x in global_df2["IL-2_gamma"]]
global_df2["time"] = global_df2["time"].mul(10)
global_df3["IL-2_gamma"] = ["$%s$" % x for x in global_df3["IL-2_gamma"]]
global_df3["time"] = global_df3["time"].mul(10)

print("averaging receptor levels")
R_mean_df_list = []
for j, value in enumerate(averaging_values):
    for my_hue_value in cell_df[my_hue].unique():
        a_cell_df = cell_df.loc[(cell_df[x_axis] == j) & (cell_df[my_hue] == my_hue_value)]
        R_mean_df = pd.DataFrame(columns=["scan_index", my_hue, x_axis, "R_mean"])
        # R_mean_df = R_mean_df.append(
        #     {"t": 0.0, "R_mean": p["R_l"]/(N_A**-1 * 1e9),
        #      "std": 10000.0,
        #      my_hue: value}, ignore_index=True)
        # for i in range(len(averaging_values)):
        R_mean_df = R_mean_df.append(
            {x_axis: value, "R_mean": a_cell_df["IL-2_R"].mean(),"std": a_cell_df["IL-2_R"].std(),
             "D": a_cell_df["IL-2_D"].unique()[0],
             "scan_index": a_cell_df["scan_index"].unique()[0], my_hue: my_hue_value}, ignore_index=True)
        R_mean_df_list.append(R_mean_df)
R_mean_concat = pd.concat(R_mean_df_list[x] for x in range(len(R_mean_df_list)))
R_mean_plotting = R_mean_concat
R_mean_plotting[my_hue] = ["$%s$" % x for x in R_mean_plotting[my_hue]]

R_mean_df_list1 = []
for j, value in enumerate(averaging_values):
    for my_hue_value in cell_df1[my_hue].unique():
        a_cell_df = cell_df1.loc[(cell_df1[x_axis] == j) & (cell_df1[my_hue] == my_hue_value)]
        R_mean_df1 = pd.DataFrame(columns=["scan_index", my_hue, x_axis, "R_mean"])
        # R_mean_df = R_mean_df.append(
        #     {"t": 0.0, "R_mean": p["R_l"]/(N_A**-1 * 1e9),
        #      "std": 10000.0,
        #      my_hue: value}, ignore_index=True)
        # for i in range(len(averaging_values)):
        R_mean_df1 = R_mean_df1.append(
            {x_axis: value, "R_mean": a_cell_df["IL-2_R"].mean(),"std": a_cell_df["IL-2_R"].std(),
             "D": a_cell_df["IL-2_D"].unique()[0],
             "scan_index": a_cell_df["scan_index"].unique()[0], my_hue: my_hue_value}, ignore_index=True)
        R_mean_df_list1.append(R_mean_df1)
R_mean_concat1 = pd.concat(R_mean_df_list1[x] for x in range(len(R_mean_df_list1)))
R_mean_plotting1 = R_mean_concat1
R_mean_plotting1[my_hue] = ["$%s$" % x for x in R_mean_plotting1[my_hue]]

R_mean_df_list2 = []
for j, value in enumerate(averaging_values):
    for my_hue_value in cell_df2[my_hue].unique():
        a_cell_df = cell_df2.loc[(cell_df2[x_axis] == j) & (cell_df2[my_hue] == my_hue_value)]
        R_mean_df2 = pd.DataFrame(columns=["scan_index", my_hue, x_axis, "R_mean"])
        # R_mean_df = R_mean_df.append(
        #     {"t": 0.0, "R_mean": p["R_l"]/(N_A**-1 * 1e9),
        #      "std": 10000.0,
        #      my_hue: value}, ignore_index=True)
        # for i in range(len(averaging_values)):
        R_mean_df2 = R_mean_df2.append(
            {x_axis: value, "R_mean": a_cell_df["IL-2_R"].mean(),"std": a_cell_df["IL-2_R"].std(),
             "D": a_cell_df["IL-2_D"].unique()[0],
             "scan_index": a_cell_df["scan_index"].unique()[0], my_hue: my_hue_value}, ignore_index=True)
        R_mean_df_list2.append(R_mean_df2)
R_mean_concat2 = pd.concat(R_mean_df_list2[x] for x in range(len(R_mean_df_list2)))
R_mean_plotting2 = R_mean_concat2
R_mean_plotting2[my_hue] = ["$%s$" % x for x in R_mean_plotting2[my_hue]]

R_mean_df_list3 = []
for j, value in enumerate(averaging_values):
    for my_hue_value in cell_df3[my_hue].unique():
        a_cell_df = cell_df3.loc[(cell_df3[x_axis] == j) & (cell_df3[my_hue] == my_hue_value)]
        R_mean_df3 = pd.DataFrame(columns=["scan_index", my_hue, x_axis, "R_mean"])
        # R_mean_df = R_mean_df.append(
        #     {"t": 0.0, "R_mean": p["R_l"]/(N_A**-1 * 1e9),
        #      "std": 10000.0,
        #      my_hue: value}, ignore_index=True)
        # for i in range(len(averaging_values)):
        R_mean_df3 = R_mean_df3.append(
            {x_axis: value, "R_mean": a_cell_df["IL-2_R"].mean(),"std": a_cell_df["IL-2_R"].std(),
             "D": a_cell_df["IL-2_D"].unique()[0],
             "scan_index": a_cell_df["scan_index"].unique()[0], my_hue: my_hue_value}, ignore_index=True)
        R_mean_df_list3.append(R_mean_df3)
R_mean_concat3 = pd.concat(R_mean_df_list3[x] for x in range(len(R_mean_df_list3)))
R_mean_plotting3 = R_mean_concat3
R_mean_plotting3[my_hue] = ["$%s$" % x for x in R_mean_plotting3[my_hue]]

# global_df[my_hue] /= round(((10 ** 0.5 * 0.01) ** 2 * 0.1),6)
# global_df[my_hue] = ["$%s$" % round(x,3) for x in global_df[my_hue]]
# global_df["time"] = global_df["time"].div(3600/10)
# cell_df["time"] = cell_df["time"].div(3600/10)

# cell_df[my_hue] = ["$%s$" % x for x in cell_df[my_hue]]


# R_mean_plotting["t"] = R_mean_plotting["t"].div(3600)
# global_df["sigma"] = global_df["sigma"].div(20000)
# R_mean_plotting["sigma"] = R_mean_plotting["sigma"].div(20000)

plot_D = 0
stoppingPoint = None
startingPoint = None

yscale = "linear"
xscale = "linear"

yscaleR = "linear"

ylim = (None, None) #(0.01, 0.045)
ylimR = (None, None) #(1, 23000)

# hue_order = [0.1,0.125,0.25,0.5,2.0,4.0,8.0,10.0]
hue_order = [0.1]
hue_order = ["$%s$" % x for x in hue_order]
# hue_order = None
# sns.set_style("ticks")
# plt.figure(figsize=(12,7.5))
# sns.set_context("talk", font_scale=1, rc={"lines.linewidth": 2.5})
# sns.set(rc={'figure.figsize':(5,5)})
# sns.set_style("ticks")
# sns.set_context("talk", font_scale=1, rc={"lines.linewidth": 5})
# fig = plt.figure()

# fig.add_subplot(a_x,a_y, 8)
# global_df["surf_c_std_norm"] = global_df["surf_c_std"]/global_df["surf_c"]
# print(global_df["surf_c_std_norm"])
# ax7 = sns.lineplot(x=x_axis, y="surf_c_std_norm", data=global_df.sort_values(by="time")[startingPoint:stoppingPoint], hue=my_hue, hue_order = hue_order)
# ax7.set(xlabel="time", ylabel="$c_v$", title="normed surf. std", yscale=yscale, xscale=xscale, ylim=ylim, xticks = [0, 50, 100,150,200])
# # plt.legend(bbox_to_anchor=(1.04,1), loc="upper left")
# plt.legend()
# plt.show()
# fig.savefig("/home/brunner/Documents/Current work/05062020/" + "kin_q_gamma" + ".png", bbox_inches='tight')
#
# exit()
global_df1 = (global_df1 + global_df2 + global_df3)
global_df1["Concentration"] = global_df1["Concentration"].div(3)
global_df1["SD"] = global_df1["SD"].div(3)
global_df1["Gradient"] = global_df1["Gradient"].div(3)
global_df1["surf_c"] = global_df1["surf_c"].div(3)
global_df1["surf_c_std"] = global_df1["surf_c_std"].div(3)
global_df1["time"] = global_df1["time"].div(3)

label1 = "averaged"

fig.add_subplot(a_x,a_y, 1)
#sns.lineplot(x="fraction", y="mean_surf_c_il2", data=saving_dataframe,hue=group_variables[1])
ax1  = sns.lineplot(x=x_axis, y="Concentration", data=global_df.loc[global_df["IL-2_gamma"] == "$0.1$"], label="3D", color="red")#,legend=False)
ax1  = sns.lineplot(x=x_axis, y="Concentration", data=global_df1.sort_values(by="time")[startingPoint:stoppingPoint],label=label1)#,legend=False)
# ax1  = sns.lineplot(x=x_axis, y="Concentration", data=global_df2.sort_values(by="time")[startingPoint:stoppingPoint],label=2)#,legend=False)
# ax1  = sns.lineplot(x=x_axis, y="Concentration", data=global_df3.sort_values(by="time")[startingPoint:stoppingPoint],label=3)#,legend=False)

ax1.set(xlabel="", ylabel="mean c. in nM", title="concentration", yscale=yscale, xscale=xscale, ylim=ylim)#, xticklabels=[]) #, ylim=(0.275,0.45))
ax1.set_xticklabels([])
#ax1.errorbar(np.linspace(0.00005,1.0,20), plotting_df[5:]["surf_c"], yerr=plotting_df[5:]["sd"], fmt='-o')

fig.add_subplot(a_x,a_y, 2)
#sns.lineplot(x="fraction", y="mean_surf_c_il2", data=saving_dataframe,hue=group_variables[1])
ax3  = sns.lineplot(x=x_axis, y="surf_c", data=global_df.loc[global_df["IL-2_gamma"] == "$0.1$"], label="3D", color="red",legend=False)
ax3  = sns.lineplot(x=x_axis, y="surf_c", data=global_df1.sort_values(by="time")[startingPoint:stoppingPoint],label=label1,legend=False)
# ax3  = sns.lineplot(x=x_axis, y="surf_c", data=global_df2.sort_values(by="time")[startingPoint:stoppingPoint],label=2,legend=False)
# ax3  = sns.lineplot(x=x_axis, y="surf_c", data=global_df3.sort_values(by="time")[startingPoint:stoppingPoint],label=3,legend=False)
# ax3  = sns.lineplot(x="t", y="surf_c_il2", data=cell_df.loc[cell_df["id"] == 1][:stoppingPoint],hue=my_hue, hue_order=hue_order, color="black")
#ax1.errorbar(np.linspace(0.00005,1.0,20), plotting_df[5:]["surf_c"], yerr=plotting_df[5:]["sd"], fmt='-o')
ax3.set(xlabel="", ylabel="mean surf.c. in nM", title="surface c.", yscale=yscale, xscale=xscale, ylim=ylim)#, xticklabels=[]) #, ylim=(0.275,0.45))
ax3.set_xticklabels([])

fig.add_subplot(a_x,a_y, 10)
global_df["Gradient"] = global_df["Gradient"]*1e5
global_df1["Gradient"] = global_df1["Gradient"]*1e5
global_df2["Gradient"] = global_df2["Gradient"]*1e5
global_df3["Gradient"] = global_df3["Gradient"]*1e5
ax2  = sns.lineplot(x=x_axis, y="Gradient", data=global_df.loc[global_df["IL-2_gamma"] == "$0.1$"], label="3D", color="red",legend=False)
ax2  = sns.lineplot(x=x_axis, y="Gradient", data=global_df1.sort_values(by="time")[startingPoint:stoppingPoint],label=label1,legend=False)
# ax2  = sns.lineplot(x=x_axis, y="Gradient", data=global_df2.sort_values(by="time")[startingPoint:stoppingPoint],label=2,legend=False)
# ax2  = sns.lineplot(x=x_axis, y="Gradient", data=global_df3.sort_values(by="time")[startingPoint:stoppingPoint],label=3,legend=False)
ax2.set(xlabel="time", ylabel="mean gradient in nM/um", title="gradient", yscale=yscale, xscale=xscale)

fig.add_subplot(a_x,a_y, 7)
global_df["std_norm"] = global_df["SD"]/global_df["Concentration"]
global_df1["std_norm"] = global_df1["SD"]/global_df1["Concentration"]
global_df2["std_norm"] = global_df2["SD"]/global_df2["Concentration"]
global_df3["std_norm"] = global_df3["SD"]/global_df3["Concentration"]
ax5  = sns.lineplot(x=x_axis, y="std_norm", data=global_df.loc[global_df["IL-2_gamma"] == "$0.1$"], label="3D", color="red",legend=False)
ax5  = sns.lineplot(x=x_axis, y="std_norm", data=global_df1.sort_values(by="time")[startingPoint:stoppingPoint],label=label1,legend=False)
# ax5  = sns.lineplot(x=x_axis, y="std_norm", data=global_df2.sort_values(by="time")[startingPoint:stoppingPoint],label=2,legend=False)
# ax5  = sns.lineplot(x=x_axis, y="std_norm", data=global_df3.sort_values(by="time")[startingPoint:stoppingPoint],label=3,legend=False)
ax5.set(xlabel="", ylabel="normed std.", title="normed std", yscale=yscale, xscale=xscale, xticklabels=[])
ax5.set_xticklabels([])
# ax5.set_yticks(np.array([1e-2, 1e-1, 1e0,1e1]),minor=True)

fig.add_subplot(a_x,a_y, 4)
# global_df["SD"] = global_df["SD"]
ax4  = sns.lineplot(x=x_axis, y="SD", data=global_df.loc[global_df["IL-2_gamma"] == "$0.1$"], label="3D", color="red",legend=False)
ax4  = sns.lineplot(x=x_axis, y="SD", data=global_df1.sort_values(by="time")[startingPoint:stoppingPoint],label=label1,legend=False)
# ax4  = sns.lineplot(x=x_axis, y="SD", data=global_df2.sort_values(by="time")[startingPoint:stoppingPoint],label=2,legend=False)
# ax4  = sns.lineplot(x=x_axis, y="SD", data=global_df3.sort_values(by="time")[startingPoint:stoppingPoint],label=3,legend=False)
ax4.set(xlabel="", ylabel="std. in nM", title="std", xticklabels=[], yscale=yscale, xscale=xscale, ylim=ylim)
ax4.set_xticklabels([])

fig.add_subplot(a_x,a_y, 5)
#std["surf_c_std_norm"] = std.set_index("t")["surf_c_il2"]/global_df.set_index("timeIndex")["surf_c"]
ax6  = sns.lineplot(x=x_axis, y="surf_c_std", data=global_df.loc[global_df["IL-2_gamma"] == "$0.1$"], label="3D", color="red",legend=False)
ax6  = sns.lineplot(x=x_axis, y="surf_c_std", data=global_df1.sort_values(by="time")[startingPoint:stoppingPoint],label=label1,legend=False)
# ax6  = sns.lineplot(x=x_axis, y="surf_c_std", data=global_df2.sort_values(by="time")[startingPoint:stoppingPoint],label=2,legend=False)
# ax6  = sns.lineplot(x=x_axis, y="surf_c_std", data=global_df3.sort_values(by="time")[startingPoint:stoppingPoint],label=3,legend=False)
ax6.set(xlabel="", ylabel="surf.c. std. in nM", title="surf. std", xticklabels=[], yscale=yscale, xscale=xscale, ylim=ylim)
ax6.set_xticklabels([])
#
#
fig.add_subplot(a_x,a_y, 8)
global_df["surf_c_std_norm"] = global_df["surf_c_std"]/global_df["surf_c"]
global_df1["surf_c_std_norm"] = global_df1["surf_c_std"]/global_df1["surf_c"]
global_df2["surf_c_std_norm"] = global_df2["surf_c_std"]/global_df2["surf_c"]
global_df3["surf_c_std_norm"] = global_df3["surf_c_std"]/global_df3["surf_c"]
ax7  = sns.lineplot(x=x_axis, y="surf_c_std_norm", data=global_df.loc[global_df["IL-2_gamma"] == "$0.1$"], label="3D", color="red",legend=False)
ax7  = sns.lineplot(x=x_axis, y="surf_c_std_norm", data=global_df1.sort_values(by="time")[startingPoint:stoppingPoint],label=label1,legend=False)
# ax7  = sns.lineplot(x=x_axis, y="surf_c_std_norm", data=global_df2.sort_values(by="time")[startingPoint:stoppingPoint],label=2,legend=False)
# ax7  = sns.lineplot(x=x_axis, y="surf_c_std_norm", data=global_df3.sort_values(by="time")[startingPoint:stoppingPoint],label=3,legend=False)
ax7.set(xlabel="time", ylabel="surf.c. std./mean surf.c.", title="normed surf. std", yscale=yscale, xscale=xscale)

# R_mean_plotting = R_mean_plotting.sort_values(by="t")[:stoppingPoint]
# x_axis = "t"

fig.add_subplot(a_x, a_y,3)
ax8 = sns.lineplot(x=x_axis, y="R_mean", data=R_mean_plotting.loc[R_mean_plotting["IL-2_gamma"] == "$0.1$"], label="3D", color="red", legend=False)
ax8 = sns.lineplot(x=x_axis, y="R_mean", data=R_mean_plotting1, label=label1, legend=False)
# ax8 = sns.lineplot(x=x_axis, y="R_mean", data=R_mean_plotting2, label=2, legend=False)
# ax8 = sns.lineplot(x=x_axis, y="R_mean", data=R_mean_plotting3, label=3, legend=False)
ax8.set(xlabel="", ylabel="Receptors", title="R mean", xticklabels=[], xscale="linear", yscale=yscaleR, ylim=ylimR)
ax8.set_xticklabels([])

fig.add_subplot(a_x,a_y, 6)
# R_mean_plotting["R_std_norm"] = R_mean_plotting["std"]/R_mean_plotting["R_mean"]
ax9 = sns.lineplot(x=x_axis, y="std", data=R_mean_plotting.loc[R_mean_plotting["IL-2_gamma"] == "$0.1$"], label="3D", color="red", legend=False)
ax9 = sns.lineplot(x=x_axis, y="std", data=R_mean_plotting1, label=label1, legend=False)
# ax9 = sns.lineplot(x=x_axis, y="std", data=R_mean_plotting2, label=2, legend=False)
# ax9 = sns.lineplot(x=x_axis, y="std", data=R_mean_plotting3, label=3, legend=False)
ax9.set(xlabel="", ylabel="std. in Receptors", title="R std", xticklabels=[], xscale="linear", yscale=yscaleR, ylim=ylimR)
ax9.set_xticklabels([])

fig.add_subplot(a_x,a_y, 9)
R_mean_plotting["R_std_norm"] = R_mean_plotting["std"]/R_mean_plotting["R_mean"]
R_mean_plotting1["R_std_norm"] = R_mean_plotting1["std"]/R_mean_plotting1["R_mean"]
R_mean_plotting2["R_std_norm"] = R_mean_plotting2["std"]/R_mean_plotting2["R_mean"]
R_mean_plotting3["R_std_norm"] = R_mean_plotting3["std"]/R_mean_plotting3["R_mean"]
ax10 = sns.lineplot(x=x_axis, y="R_std_norm", data=R_mean_plotting.loc[R_mean_plotting["IL-2_gamma"] == "$0.1$"], label="3D", color="red", legend=False)
ax10 = sns.lineplot(x=x_axis, y="R_std_norm", data=R_mean_plotting1, label=label1, legend=False)
# ax10 = sns.lineplot(x=x_axis, y="R_std_norm", data=R_mean_plotting2, label=2, legend=False)
# ax10 = sns.lineplot(x=x_axis, y="R_std_norm", data=R_mean_plotting3, label=3, legend=False)
ax10.set(xlabel="time", ylabel="std./mean", title="normed R std", xscale="linear", yscale=yscaleR, ylim = ylimR)
plt.show()

# sns.set(rc={"lines.linewidth":0.1})
#
# fig = plt.figure()
# ax = sns.lineplot(x="time_index", y="IL-2_R", data=cell_df.loc[cell_df["scan_index"]==0], estimator=None, units="id", legend=False)
# ax.set(xlabel="t", ylabel="each cells R, gamma = 0.5") #, ylim=(0.275,0.45))
# plt.show()
#
# fig = plt.figure()
# ax = sns.lineplot(x="time_index", y="IL-2_R", data=cell_df.loc[cell_df["scan_index"]==1], estimator=None, units="id", legend=False)
# ax.set(xlabel="t", ylabel="each cells R, gamma = 2") #, ylim=(0.275,0.45))
# plt.show()
