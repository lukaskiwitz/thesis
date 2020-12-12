import getpass
import os

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import seaborn as sns

get_dataframes = [[]]#[[]]*3
saving_dataframe = pd.DataFrame()

# get_dataframes.append([])
# path = "/extra/brunner/10x10x10/R_lognorm/run" + str(j) + "/"
# path = "/extra/brunner/10x10x10/q_fraction_exact/run" + str(j) + "/"
user = getpass.getuser()
path = "/extra/brunner/thesis/kinetic/saturation/R_lognorm" \
       "/"
ext_cache="/extra/brunner/para_handling/kinetic/R_lognorm_ext_cache/"

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


get_dataframes[0] = [global_df, cell_df]#, cell_stats_df]

# global_df_2 = pd.read_hdf("/extra/brunner/para_handling/kinetic/R_lognorm_test_2/" + 'global_df.h5', mode="r")
# cell_df_2 = pd.read_hdf("/extra/brunner/para_handling/kinetic/R_lognorm_test_2/" + 'cell_df.h5', mode="r")
#
# get_dataframes[1] = [global_df_2, cell_df_2]

my_hue = "IL-2_gamma" #"D"
x_axis = "time"
averaging_values = np.sort(global_df[x_axis].unique())

global_df = pd.concat((get_dataframes[x][0] for x in range(len(get_dataframes))))#.groupby([x_axis], as_index=False).mean()
global_df.reset_index(inplace=True)
cell_df = pd.concat((get_dataframes[x][1] for x in range(len(get_dataframes))))#.groupby(["sigma"], as_index=False).mean()
# cell_df.reset_index(inplace=True)
# cell_stats_df = pd.concat((get_dataframes[x][2] for x in range(len(get_dataframes)))).groupby(["sigma"], as_index=False).mean()
# cell_stats_df.reset_index(inplace=True)

# sns.set(rc={'figure.figsize':(14,8.27)})
fig = plt.figure()
plt.subplots_adjust(wspace=.3)
# sns.set_style("white")

# subplot_style
a_x = 4
a_y = 3

#global_df["D"] /= ((10 ** 0.5 * 0.01) ** 2)

global_df["IL-2_gamma"] = ["$%s$" % x for x in global_df["IL-2_gamma"]]
global_df["time"] = global_df["time"].mul(3600)
global_df = global_df.loc[global_df["IL-2_k_factor"] == 2]


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


# global_df[my_hue] /= round(((10 ** 0.5 * 0.01) ** 2 * 0.1),6)
# global_df[my_hue] = ["$%s$" % round(x,3) for x in global_df[my_hue]]
global_df["time"] = global_df["time"].div(3600/10)
cell_df["time"] = cell_df["time"].mul(10)

# cell_df[my_hue] = ["$%s$" % x for x in cell_df[my_hue]]

R_mean_plotting[my_hue] = ["$%s$" % x for x in R_mean_plotting[my_hue]]
# R_mean_plotting["t"] = R_mean_plotting["t"].div(3600)
# global_df["sigma"] = global_df["sigma"].div(20000)
# R_mean_plotting["sigma"] = R_mean_plotting["sigma"].div(20000)

plot_D = 0
stoppingPoint = None
startingPoint = None

yscale = "linear"
xscale = "linear"

yscaleR = "log"

ylim = (None, None) #(0.003, 0.045)
ylimR = (None, None)

hue_order = None# ["$0.5$", "$2.0$"]
# hue_order = [0.1,0.125,0.25,0.5,2.0,4.0,8.0,10.0]
# hue_order = ["$%s$" % x for x in hue_order]
#
# sns.set(rc={'figure.figsize':(12,7.5)})
# sns.set_style("ticks")
# plt.figure(figsize=(12,7.5))
# sns.set_context("talk", font_scale=1, rc={"lines.linewidth": 2.5})
# fig = plt.figure()
# # fig.add_subplot(a_x,a_y, 8)
# global_df["surf_c_std_norm"] = global_df["surf_c_std"]/global_df["surf_c"]
# print(global_df["surf_c_std_norm"])
# ax7 = sns.lineplot(x=x_axis, y="surf_c_std_norm", data=global_df.sort_values(by="time")[startingPoint:stoppingPoint], hue=my_hue, hue_order = hue_order)
# ax7.set(xlabel="time", ylabel="surf.c. std./mean surf.c.", title="normed surf. std", yscale=yscale, xscale=xscale)
# plt.legend(bbox_to_anchor=(1.04,1), loc="upper left")
# plt.show()

fig.add_subplot(a_x,a_y, 1)
#sns.lineplot(x="fraction", y="mean_surf_c_il2", data=saving_dataframe,hue=group_variables[1])
ax1  = sns.lineplot(x=x_axis, y="Concentration", data=global_df.sort_values(by="time")[startingPoint:stoppingPoint],hue=my_hue, hue_order=hue_order)#,legend=False)
ax1.set(xlabel="", ylabel="mean c. in nM", title="concentration", yscale=yscale, xscale=xscale, ylim=ylim)#, xticklabels=[]) #, ylim=(0.275,0.45))
ax1.set_xticklabels([])
#ax1.errorbar(np.linspace(0.00005,1.0,20), plotting_df[5:]["surf_c"], yerr=plotting_df[5:]["sd"], fmt='-o')

fig.add_subplot(a_x,a_y, 2)
#sns.lineplot(x="fraction", y="mean_surf_c_il2", data=saving_dataframe,hue=group_variables[1])
ax3  = sns.lineplot(x=x_axis, y="surf_c", data=global_df.sort_values(by="time")[startingPoint:stoppingPoint],hue=my_hue, hue_order=hue_order,legend=False)
# ax3  = sns.lineplot(x="t", y="surf_c_il2", data=cell_df.loc[cell_df["id"] == 1][:stoppingPoint],hue=my_hue, hue_order=hue_order, color="black")
#ax1.errorbar(np.linspace(0.00005,1.0,20), plotting_df[5:]["surf_c"], yerr=plotting_df[5:]["sd"], fmt='-o')
ax3.set(xlabel="", ylabel="mean surf.c. in nM", title="surface c.", yscale=yscale, xscale=xscale, ylim=ylim)#, xticklabels=[]) #, ylim=(0.275,0.45))
ax3.set_xticklabels([])

fig.add_subplot(a_x,a_y, 10)
global_df["Gradient"] = global_df["Gradient"]*1e5
ax2 = sns.lineplot(x=x_axis, y="Gradient", hue=my_hue, hue_order=hue_order, data=global_df.sort_values(by="time")[startingPoint:stoppingPoint], legend=False)
ax2.set(xlabel="time", ylabel="mean gradient in nM/um", title="gradient", yscale=yscale, xscale=xscale)

fig.add_subplot(a_x,a_y, 7)
global_df["std_norm"] = global_df["SD"]/global_df["Concentration"]
ax5 = sns.lineplot(x=x_axis, y="std_norm", data=global_df.sort_values(by="time")[startingPoint:stoppingPoint], hue=my_hue, hue_order=hue_order, label="std_norm", legend=False)
ax5.set(xlabel="", ylabel="normed std.", title="normed std", yscale=yscale, xscale=xscale, ylim=(None,None), xticklabels=[])
ax5.set_xticklabels([])
# ax5.set_yticks(np.array([1e-2, 1e-1, 1e0,1e1]),minor=True)

fig.add_subplot(a_x,a_y, 4)
global_df["SD"] = global_df["SD"]
ax4 = sns.lineplot(x=x_axis, y="SD", data=global_df.sort_values(by="time")[startingPoint:stoppingPoint], hue=my_hue, hue_order=hue_order, label="std", legend=False)
ax4.set(xlabel="", ylabel="std. in nM", title="std", xticklabels=[], yscale=yscale, xscale=xscale, ylim=ylim)
ax4.set_xticklabels([])

fig.add_subplot(a_x,a_y, 5)
#std["surf_c_std_norm"] = std.set_index("t")["surf_c_il2"]/global_df.set_index("timeIndex")["surf_c"]
ax6 = sns.lineplot(x=x_axis, y="surf_c_std", data=global_df.sort_values(by="time")[startingPoint:stoppingPoint], hue=my_hue, hue_order=hue_order, label="surf_c_std", legend=False)
ax6.set(xlabel="", ylabel="surf.c. std. in nM", title="surf. std", xticklabels=[], yscale=yscale, xscale=xscale, ylim=ylim)
ax6.set_xticklabels([])
#
#
fig.add_subplot(a_x,a_y, 8)
global_df["surf_c_std_norm"] = global_df["surf_c_std"]/global_df["surf_c"]
ax7 = sns.lineplot(x=x_axis, y="surf_c_std_norm", data=global_df.sort_values(by="time")[startingPoint:stoppingPoint], hue=my_hue, hue_order=hue_order, label="surf_c_std", legend=False)
ax7.set(xlabel="time", ylabel="surf.c. std./mean surf.c.", title="normed surf. std", yscale=yscale, xscale=xscale, ylim=(None, None))

# R_mean_plotting = R_mean_plotting.sort_values(by="t")[:stoppingPoint]
# x_axis = "t"

fig.add_subplot(a_x, a_y,3)
ax8 = sns.lineplot(x=x_axis, y="R_mean", data=R_mean_plotting.sort_values(by=x_axis)[startingPoint:stoppingPoint],hue=my_hue, hue_order=hue_order, label="R_mean", legend=False)
ax8.set(xlabel="", ylabel="Receptors", title="R mean", xticklabels=[], xscale="linear", yscale=yscaleR, ylim=ylimR)
ax8.set_xticklabels([])

fig.add_subplot(a_x,a_y, 6)
# R_mean_plotting["R_std_norm"] = R_mean_plotting["std"]/R_mean_plotting["R_mean"]
ax9 = sns.lineplot(x=x_axis, y="std", data=R_mean_plotting.sort_values(by=x_axis)[startingPoint:stoppingPoint], hue=my_hue, hue_order=hue_order, label="R_std", legend=False)
ax9.set(xlabel="", ylabel="std. in Receptors", title="R std", xticklabels=[], xscale="linear", yscale=yscaleR, ylim=ylimR)
ax9.set_xticklabels([])

fig.add_subplot(a_x,a_y, 9)
R_mean_plotting["R_std_norm"] = R_mean_plotting["std"]/R_mean_plotting["R_mean"]
ax10 = sns.lineplot(x=x_axis, y="R_std_norm", data=R_mean_plotting.sort_values(by=x_axis)[startingPoint:stoppingPoint], hue=my_hue, hue_order=hue_order, label="R_std", legend=False)
ax10.set(xlabel="time", ylabel="std./mean", title="normed R std", xscale="linear", yscale=yscaleR, ylim = ylimR)
plt.show()
#
# [0.1,0.125,0.25,0.5,2,4,8,10]
# sns.set(rc={"lines.linewidth":0.1})
#
# fig = plt.figure()
# ax = sns.lineplot(x="time", y="IL-2_R", data=cell_df.loc[cell_df["scan_index"]==2][:133100], estimator=None, units="id", legend=False)
# ax.set(xlabel="t", ylabel="each cells R, gamma = 0.5") #, ylim=(0.275,0.45))
# plt.show()
#
# fig = plt.figure()
# ax = sns.lineplot(x="time", y="IL-2_R", data=cell_df.loc[cell_df["scan_index"]==5], estimator=None, units="id", legend=False)
# ax.set(xlabel="t", ylabel="each cells R, gamma = 2") #, ylim=(0.275,0.45))
# plt.show()