import getpass
import os

import matplotlib.pyplot as plt
import matplotlib as mpl
import numpy as np
import pandas as pd
import seaborn as sns

get_dataframes = []#[[]]*3
saving_dataframe = pd.DataFrame()

get_dataframes.append([])
user = getpass.getuser()
# path = "/extra/brunner/thesis/kinetic/standard/8_big_scan/"
path = r"/home/brunner/thesis/thesis/scripts/patrick/yukawa_analytical/test/"

# ext_cache="/extra/brunner/para_handling/kinetic/R_lognorm_ext_cache/"

T = np.arange(0, 200, 1)
dt = 3600


global_df = pd.read_hdf(path + 'global_df.h5', mode="r")

offset = 2
frac = 0.25
gammas_list = [sorted(global_df["IL-2_gamma"].unique())[:10],sorted(global_df["IL-2_gamma"].unique())[10:]]
sns.set_style("ticks")
# fig = plt.figure(figsize=(8,8))
# sns.set_context("talk", font_scale=1.1, rc={"lines.linewidth": 2.5})
# sns.set(rc={'figure.figsize':(8,12)})
# sns.set_style("ticks")
# sns.set(rc={'figure.figsize':(6,5)})

for gammas in gammas_list:
    for frac in [0.25]:
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
        # global_df = global_df.loc[global_df["IL-2_Tsec_fraction"] == frac]
        cell_df = cell_df.loc[cell_df["IL-2_Tsec_fraction"] == frac]

        # global_df = global_df.drop(0)
        # global_df[:49]["time"] = global_df[:49]["time"] - 1

        # subplot_style
        a_x = 4
        a_y = 3

        #global_df["D"] /= ((10 ** 0.5 * 0.01) ** 2)

        # global_df["IL-2_gamma"] = ["$%s$" % x for x in global_df["IL-2_gamma"]]
        global_df["time"] = global_df["time"].mul(3600)
        # global_df = global_df.loc[global_df["IL-2_gamma"] == 1/100]


        print("calculating c")
        x, y, z = (cell_df["x"].unique(), cell_df["y"].unique(), cell_df["z"].unique())

        offset_cells_ids = cell_df.loc[((cell_df["x"] < x[offset]) | \
                                     (cell_df["x"] > x[-offset - 1])) | \
                                    ((cell_df["y"] < y[offset]) | \
                                     (cell_df["y"] > y[-offset - 1])) | \
                                    (cell_df["z"] < z[offset]) | \
                                    (cell_df["z"] > z[-offset - 1]), "id"].unique()
        cells_ids = cell_df.loc[(cell_df[my_hue] == cell_df[my_hue].unique()[-1]) & (cell_df[x_axis] == averaging_values[-1]), "id"]
        cells_ids = [x for x in cells_ids if x not in offset_cells_ids]

        c_df = pd.DataFrame(index=range(len(averaging_values)), columns=["time"] + list(cell_df[my_hue].unique()))
        c_df["time"] = averaging_values
        for my_hue_value in cell_df[my_hue].unique(): #g
            c_over_t = []
            std = []
            time_index = []
            for v, value in enumerate(averaging_values): #t
                time_index.append(v)
                tmp_df = cell_df.loc[(cell_df[my_hue] == my_hue_value) & (cell_df[x_axis] == value)]
                tmp_values = tmp_df[tmp_df["id"].isin(cells_ids)]["IL-2_surf_c"].values
                c_over_t.append(np.mean(tmp_values))
            c_df[my_hue_value] = c_over_t
            c_df["time_index"] = time_index

global_df["surf_c_std_norm"] = global_df["surf_c_std"]/global_df["surf_c"]
global_df["surf_c"] *= 1e3

global_df["time"] = global_df["time"].div(36000)
c_df["time"] = c_df["time"].div(3600)


sns.set_style("ticks")
sns.set_context("talk", font_scale=1.5, rc={"lines.linewidth": 2.5})
fig,ax = plt.subplots(figsize=(7,5))


print("plotting")
global_df["time"] = global_df["time"].div(3600/10)
cell_df["time"] = cell_df["time"].mul(10)


plot_D = 0
stoppingPoint = None
startingPoint = None

yscale = "linear"
xscale = "linear"

yscaleR = "linear"

# ylim = (1.0, 1.55) #(0.01, 0.045)
ylim = (None, None)
xlim = (None,None)

ylimR = (None, None) #(1, 23000)

# hue_order = [10.0,0.5,0.1]
# hue_order = ["$%s$" % x for x in hue_order]
hue_order = None

palette = sns.color_palette("bwr", len(np.array(gammas_list).flatten()))
neg = palette[:len(palette)//2]
pos = palette[len(palette)//2:]
palette = [neg,pos]
norm = plt.Normalize(np.min(np.array(gammas_list).flatten()), np.max(np.array(gammas_list).flatten()))
sm = plt.cm.ScalarMappable(cmap="bwr", norm=norm)
sm.set_array([])
# plt.axes().get_legend().remove()
# ax.figure.colorbar(sm, label='feedback strength')
cbar = fig.colorbar(sm, label='feedback', ticks=[np.min(np.array(gammas_list).flatten()), np.max(np.array(gammas_list).flatten())])
cbar.ax.set_yticklabels(["neg.", 'pos.'])

for ga,gammas in enumerate(gammas_list):
    # if gammas[0] < 1:
    #     palette = sns.color_palette("Blues", len(gammas))
    #     norm = plt.Normalize(np.min(gammas), np.max(gammas))
    #     cmap = mpl.cm.Blues
    #     sm = plt.cm.ScalarMappable(cmap=cmap.reversed(), norm=norm)
    #     sm.set_array([])
    #     # plt.axes().get_legend().remove()
    #     ax.figure.colorbar(sm, label='feedback strength')
    #     palette = [x for x in reversed(palette)]
    # elif gammas[0] > 1:
    #     palette = sns.color_palette("Reds", len(gammas))
    #     norm = plt.Normalize(np.min(gammas), np.max(gammas))
    #     sm = plt.cm.ScalarMappable(cmap="Reds", norm=norm)
    #     sm.set_array([])
    #     # plt.axes().get_legend().remove()
    #     ax.figure.colorbar(sm, label='feedback strength')
    # else:
    #     raise ValueError
    if gammas[0] < 1:
        # halfed_palette = [x for x in reversed(palette[ga])]
        halfed_palette = palette[ga]
    else:
        halfed_palette = palette[ga]

    # test = []
    # ax7 = sns.lineplot(x=x_axis, y="surf_c_std_norm", data=global_df.sort_values(by="time")[startingPoint:stoppingPoint], hue=my_hue, hue_order = hue_order)
    for g,gamma in enumerate(gammas):
    #     test.append(c_df[gamma].values[-10])
    # plt.plot(gammas,test)
    # plt.show()
        temp_df = c_df.copy()
        temp_df[gamma] = temp_df[gamma].div(temp_df[gamma].values[0])
        # temp_df["surf_c"] *= 1e3
        if gamma < 1:
            sns.lineplot(x=x_axis, y=gamma, data=temp_df.sort_values(by="time")[startingPoint:stoppingPoint], color=halfed_palette[g])
        else:
            sns.lineplot(x=x_axis, y=gamma, data=temp_df.sort_values(by="time")[startingPoint:stoppingPoint],
                     color=halfed_palette[g])

    plt.hlines(y=temp_df.loc[(temp_df["time_index"] == 0), gamma].iloc[0],xmin=0,xmax=temp_df.sort_values(by="time")["time"][startingPoint:stoppingPoint].unique()[-1], color = "black", linestyles="dashed")

    import matplotlib.lines as mlines

    new_handles = []
    for g,gamma in enumerate(gammas):
        new_handles.append(mlines.Line2D([], [], color=halfed_palette[g], label=str(gamma)))

    black_line = mlines.Line2D([], [], color='black', label='Black line')
    new_handles.append(black_line)


    new_labels = ["neg. feedback", "pos. feedback", "no feedback"]
    # plt.legend(handles=new_handles, labels=new_labels)

ax.set(ylabel="fold change", xlabel="time (h)", yscale=yscale, xscale=xscale, ylim=ylim, xlim=xlim, title="surface c.")

# if gamma < 1:
#     saving_string =r"/home/brunner/Documents/Current work/04122020/" + "kin_satu_neg_c_f_" + str(cell_df["IL-2_Tsec_fraction"].unique()[0]) + ".png"
# elif gamma > 1:
#     saving_string =r"/home/brunner/Documents/Current work/04122020/" + "kin_satu_pos_c_f_" + str(cell_df["IL-2_Tsec_fraction"].unique()[0]) + ".png"

saving_string =r"/home/brunner/Documents/Current work/04122020/" + "kin_analytical_c_f_" + str(cell_df["IL-2_Tsec_fraction"].unique()[0]) + ".png"

fig.savefig(saving_string, bbox_inches='tight')
plt.show()
