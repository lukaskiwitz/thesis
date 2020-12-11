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

# path = "/extra/brunner/thesis/kinetic/test_new_hill/"
path = "/extra/brunner/thesis/kinetic/halfed_hill/3/"
# ext_cache="/extra/brunner/para_handling/kinetic/R_lognorm_ext_cache/"

T = np.arange(0, 200, 1)
dt = 3600*5

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
# global_df["time"] = global_df["time"].mul(dt)
cell_df_1 = cell_df.copy()

global_df = global_df.loc[global_df["IL-2_Tsec_fraction"] == 0.25]
cell_df = cell_df.loc[cell_df["IL-2_Tsec_fraction"] == 0.25]


print("averaging receptor levels")
R_mean_df_list = []
for j, value in enumerate(averaging_values):
    for my_hue_value in cell_df[my_hue].unique():
        a_cell_df = cell_df.loc[(cell_df[x_axis] == value) & (cell_df[my_hue] == my_hue_value)]
        R_mean_df = pd.DataFrame(columns=["scan_index", my_hue, x_axis, "R_mean"])
        if j == 0:
            R_mean_df = R_mean_df.append(
                {"time": 0.0, "R_mean": 10000,
                 "std": 0.0,
                 my_hue: my_hue_value}, ignore_index=True)
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
global_df["time"] = global_df["time"].mul(1/3600)
# cell_df["time"] = cell_df["time"].mul(dt/3600)
R_mean_plotting["time"] = R_mean_plotting["time"].mul(1/3600)

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

# ylim = (1.0, 1.55) #(0.01, 0.045)
ylim = (None, None)
ylimR = (None, None) #(1, 23000)

hue_order = [1/10,10.]
hue_order = ["$%s$" % x for x in hue_order]
# hue_order = None

fig.add_subplot(a_x,a_y, 1)
#sns.lineplot(x="fraction", y="mean_surf_c_il2", data=saving_dataframe,hue=group_variables[1])
global_df["Concentration"] = global_df["Concentration"].mul(1e3)
ax1  = sns.lineplot(x=x_axis, y="Concentration", data=global_df.sort_values(by="time")[startingPoint:stoppingPoint],hue=my_hue, hue_order=hue_order)#,legend=False)
ax1.set(xlabel="", ylabel="mean c. in pM", title="concentration", yscale=yscale, xscale=xscale, ylim=ylim)#, xticklabels=[]) #, ylim=(0.275,0.45))
ax1.set_xticklabels([])
#ax1.errorbar(np.linspace(0.00005,1.0,20), plotting_df[5:]["surf_c"], yerr=plotting_df[5:]["sd"], fmt='-o')

fig.add_subplot(a_x,a_y, 2)
global_df["surf_c"] = global_df["surf_c"].mul(1e3)
#sns.lineplot(x="fraction", y="mean_surf_c_il2", data=saving_dataframe,hue=group_variables[1])
ax3  = sns.lineplot(x=x_axis, y="surf_c", data=global_df.sort_values(by="time")[startingPoint:stoppingPoint],hue=my_hue, hue_order=hue_order,legend=False)
# ax3  = sns.lineplot(x="t", y="surf_c_il2", data=cell_df.loc[cell_df["id"] == 1][:stoppingPoint],hue=my_hue, hue_order=hue_order, color="black")
#ax1.errorbar(np.linspace(0.00005,1.0,20), plotting_df[5:]["surf_c"], yerr=plotting_df[5:]["sd"], fmt='-o')
ax3.set(xlabel="", ylabel="mean surf.c. in pM", title="surface c.", yscale=yscale, xscale=xscale, ylim=ylim)#, xticklabels=[]) #, ylim=(0.275,0.45))
ax3.set_xticklabels([])

fig.add_subplot(a_x,a_y, 10)
global_df["Gradient"] = global_df["Gradient"]*1e5
ax2 = sns.lineplot(x=x_axis, y="Gradient", hue=my_hue, hue_order=hue_order, data=global_df.sort_values(by="time")[startingPoint:stoppingPoint], legend=False)
ax2.set(xlabel="time in h", ylabel="mean gradient in nM/um", title="gradient", yscale=yscale, xscale=xscale)

fig.add_subplot(a_x,a_y, 7)
global_df["std_norm"] = global_df["SD"].mul(1e3)/global_df["Concentration"]
ax5 = sns.lineplot(x=x_axis, y="std_norm", data=global_df.sort_values(by="time")[startingPoint:stoppingPoint], hue=my_hue, hue_order=hue_order, label="std_norm", legend=False)
ax5.set(xlabel="", ylabel="normed std.", title="normed std", yscale=yscale, xscale=xscale, ylim=(None,None), xticklabels=[])
ax5.set_xticklabels([])
# ax5.set_yticks(np.array([1e-2, 1e-1, 1e0,1e1]),minor=True)

fig.add_subplot(a_x,a_y, 4)
global_df["SD"] = global_df["SD"]*1e3
ax4 = sns.lineplot(x=x_axis, y="SD", data=global_df.sort_values(by="time")[startingPoint:stoppingPoint], hue=my_hue, hue_order=hue_order, label="std", legend=False)
ax4.set(xlabel="", ylabel="std. in pM", title="std", xticklabels=[], yscale=yscale, xscale=xscale, ylim=ylim)
ax4.set_xticklabels([])

fig.add_subplot(a_x,a_y, 5)
#std["surf_c_std_norm"] = std.set_index("t")["surf_c_il2"]/global_df.set_index("timeIndex")["surf_c"]
global_df["surf_c_std"] = global_df["surf_c_std"].mul(1e3)
ax6 = sns.lineplot(x=x_axis, y="surf_c_std", data=global_df.sort_values(by="time")[startingPoint:stoppingPoint], hue=my_hue, hue_order=hue_order, label="surf_c_std", legend=False)
ax6.set(xlabel="", ylabel="surf.c. std. in pM", title="surf. std", xticklabels=[], yscale=yscale, xscale=xscale, ylim=ylim)
ax6.set_xticklabels([])
#
#
fig.add_subplot(a_x,a_y, 8)
global_df["surf_c_std_norm"] = global_df["surf_c_std"]/global_df["surf_c"]
ax7 = sns.lineplot(x=x_axis, y="surf_c_std_norm", data=global_df.sort_values(by="time")[startingPoint:stoppingPoint], hue=my_hue, hue_order=hue_order, label="surf_c_std", legend=False)
ax7.set(xlabel="", ylabel="surf.c. std./mean surf.c.", title="normed surf. std", yscale=yscale, xscale=xscale, ylim=(None, None))
ax7.set_xticklabels([])

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
ax10.set(xlabel="time in h", ylabel="std./mean", title="normed R std", xscale="linear", yscale=yscaleR, ylim = ylimR)

fig.add_subplot(a_x,a_y, 11)
# R_mean_plotting["R_std_norm"] = R_mean_plotting["std"]/R_mean_plotting["R_mean"]
global_df["h1_norm"] = global_df["h1_norm"]*1e5
ax11 = sns.lineplot(x=x_axis, y="h1_norm", data=global_df.sort_values(by="time")[startingPoint:stoppingPoint], hue=my_hue, hue_order=hue_order, label="R_std", legend=False)
ax11.set(xlabel="time", ylabel="H1 norm", title="H1 norm", xscale="linear", yscale=yscale, ylim = ylim)

plt.show()



sns.set(font_scale=2, rc={"lines.linewidth":0.1})


cell_df_1["time"] = cell_df_1["time"].div(3600)

Tsec_fraction = 0.25
gamma = 1/10
tmp_df = cell_df_1.loc[(cell_df_1["IL-2_gamma"] == gamma) & (cell_df_1["IL-2_Tsec_fraction"] == Tsec_fraction)]

fig = plt.figure()
ax = sns.lineplot(x="time", y="IL-2_R", data=tmp_df, estimator=None, units="id", legend=False)
# ax.set(xlabel="t", ylabel="each cells R", title="gamma = " + str(cell_df_1.loc[cell_df["scan_index"]==cell_df["scan_index"].unique()[0], "IL-2_gamma"].values[0])) #, ylim=(0.275,0.45))
ax.set(xlabel="time (h)", ylabel="each cells R", title = "g = 1/" + str(int(1/gamma))) #, ylim=(0.275,0.45))
# fig.savefig(r"/home/brunner/Documents/Current work/02102020/" + "kin_q_R_neg_f_" + str(global_df["IL-2_Tsec_fraction"].unique()[0]) + ".svg", bbox_inches='tight')
plt.show()

gamma =  10
tmp_df = cell_df_1.loc[(cell_df_1["IL-2_gamma"] == gamma) & (cell_df_1["IL-2_Tsec_fraction"] == Tsec_fraction)]
fig = plt.figure()
ax = sns.lineplot(x="time", y="IL-2_R", data=tmp_df, estimator=None, units="id", legend=False)
# ax.set(xlabel="t", ylabel="each cells R", title="gamma = " + str(cell_df_1.loc[cell_df["scan_index"]==cell_df["scan_index"].unique()[1], "IL-2_gamma"].values[0])) #, ylim=(0.275,0.45))
ax.set(xlabel="time (h)", ylabel="each cells R", title="g = " + str(gamma)) #, ylim=(0.275,0.45))
# fig.savefig(r"/home/brunner/Documents/Current work/02102020/" + "kin_q_R_pos_f_" + str(global_df["IL-2_Tsec_fraction"].unique()[0]) + ".svg", bbox_inches='tight')

plt.show()
# del fig, ax
# fig = plt.figure()
# ax = sns.lineplot(x="time_index", y="IL-2_R", data=cell_df.loc[cell_df["scan_index"]==1], estimator=None, units="id", legend=False)
# ax.set(xlabel="t", ylabel="each cells R, gamma = 2") #, ylim=(0.275,0.45))
# plt.show()
