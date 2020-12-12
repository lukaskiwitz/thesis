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
# path = "/extra/brunner/10x10x10/R_lognorm/run" + str(j) + "/"
# path = "/extra/brunner/10x10x10/q_fraction_exact/run" + str(j) + "/"
user = getpass.getuser()
path = "/extra/brunner/thesis/kinetic/standard/8_big_scan/"
# path = "/extra/brunner/thesis/kinetic/q_fraction_gamma_scan/"
# ext_cache="/extra/brunner/para_handling/kinetic/R_lognorm_ext_cache/"

T = np.arange(0, 200, 1)
dt = 3600
global_df = pd.read_hdf(path + 'global_df.h5', mode="r")


frac = 0.25
gammas_list = [sorted(global_df["IL-2_gamma"].unique())[:10],sorted(global_df["IL-2_gamma"].unique())[10:]]
for gammas in gammas_list:
    for frac in [0.25,0.5]:
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

        global_df = global_df.loc[global_df["IL-2_Tsec_fraction"] == frac]
        cell_df = cell_df.loc[cell_df["IL-2_Tsec_fraction"] == frac]

        sns.set(rc={'figure.figsize':(11,8.27)})
        fig = plt.figure()
        plt.subplots_adjust(wspace=.3)

        # subplot_style
        a_x = 4
        a_y = 3


        global_df["time"] = global_df["time"].div(3600)
        cell_df["time"] = cell_df["time"].mul(10)



        plot_D = 0
        stoppingPoint = None
        startingPoint = None

        yscale = "linear"
        xscale = "linear"

        yscaleR = "linear"

        # ylim = (1.0, 1.55) #(0.01, 0.045)
        # ylim = (0.65, 0.99)
        ylim = (0,4)
        xlim = (None,None)

        ylimR = (None, None) #(1, 23000)

        # hue_order = [10.0,0.5,0.1]
        # hue_order = ["$%s$" % x for x in hue_order]
        hue_order = None
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
        # global_df.loc[(global_df["time_index"] == 0) & (global_df["IL-2_gamma"] == "$2.0$"), "surf_c_std_norm"] = global_df.loc[(global_df["time_index"] == 0) & (global_df["IL-2_gamma"] == "$0.5$"), "surf_c_std_norm"].values[0]
        # global_df.loc[(global_df["time_index"] == 0) & (global_df["IL-2_gamma"] == "$0.1$"), "surf_c_std_norm"] = global_df.loc[(global_df["time_index"] == 0) & (global_df["IL-2_gamma"] == "$0.5$"), "surf_c_std_norm"].values[0]
        # global_df.loc[(global_df["time_index"] == 0) & (global_df["IL-2_gamma"] == "$10.0$"), "surf_c_std_norm"] = global_df.loc[(global_df["time_index"] == 0) & (global_df["IL-2_gamma"] == "$0.5$"), "surf_c_std_norm"].values[0]
        # insert data from same sim just with finer time steps to even out the overshooting. todo: finer timesteps for whole sim!
        # global_df.loc[(global_df["IL-2_gamma"] == '$10.0$') & (global_df["time"] == 10), "surf_c"] = 8.93*1e-3
        # global_df.loc[(global_df["IL-2_gamma"] == '$10.0$') & (global_df["time"] == 10), "surf_c_std_norm"] = 0.00751/(8.93*1e-3)
        # print(global_df["surf_c_std_norm"])
        global_df["surf_c"] *= 1e3

        if gammas[0] < 1:
            palette = sns.color_palette("Blues", len(gammas))
            norm = plt.Normalize(np.min(gammas), np.max(gammas))
            cmap = mpl.cm.Blues
            sm = plt.cm.ScalarMappable(cmap=cmap.reversed(), norm=norm)
            sm.set_array([])
            # plt.axes().get_legend().remove()
            ax.figure.colorbar(sm)
        elif gammas[0] > 1:
            palette = sns.color_palette("Reds", len(gammas))
            norm = plt.Normalize(np.min(gammas), np.max(gammas))
            sm = plt.cm.ScalarMappable(cmap="Reds", norm=norm)
            sm.set_array([])
            # plt.axes().get_legend().remove()
            ax.figure.colorbar(sm)
        else:
            raise ValueError


        for g,gamma in enumerate(gammas):
            temp_df = global_df.loc[global_df["IL-2_gamma"] == gamma].copy()
            temp_df["h1_norm"] = temp_df["h1_norm"].div(temp_df.loc[temp_df["time_index"] == 1, "h1_norm"].values[0])

            sns.lineplot(x=x_axis, y="h1_norm", data=temp_df.sort_values(by="time")[startingPoint:stoppingPoint], color=palette[g])

            # sns.lineplot(x=x_axis, y="h1_norm", data=global_df.loc[global_df["IL-2_gamma"] == gamma].sort_values(by="time")[startingPoint:stoppingPoint], color=palette[g])

        plt.hlines(y=temp_df.loc[(temp_df["time_index"] == 1), "h1_norm"].iloc[0],xmin=0,xmax=temp_df.loc[temp_df["IL-2_gamma"] == gamma].sort_values(by="time")["time"][startingPoint:stoppingPoint].unique()[-1])
        # plt.hlines(y=global_df.loc[(global_df["time_index"] == 1), "h1_norm"].iloc[0],xmin=0,xmax=global_df.loc[global_df["IL-2_gamma"] == gamma].sort_values(by="time")["time"].unique()[-1])
        ax.set(xlabel="time (h)", ylabel="fold change", yscale=yscale, xscale=xscale, ylim=ylim, xlim=xlim, title="H$_1$ norm")

        import matplotlib.lines as mlines

        new_handles = []
        for g,gamma in enumerate(gammas):
            new_handles.append(mlines.Line2D([], [], color=palette[g], label=str(gamma)))
        # zero_line = mlines.Line2D([], [], color=palette[0], label='lightcoral line')
        # one_line = mlines.Line2D([], [], color=palette[1], label='red line')
        black_line = mlines.Line2D([], [], color='black', label='Black line')
        new_handles.append(black_line)
        # two_line = mlines.Line2D([], [], color=palette[2], label='lightsteelblue line')
        # three_line = mlines.Line2D([], [], color=palette[3], label='blue line')
        #
        # grey_line = mlines.Line2D([], [], color="grey", label='grey line')
        # grey_dashed_line = mlines.Line2D([], [], color="grey", label='grey dashed line', linestyle="dashed")
        # white_line = mlines.Line2D([], [], color="white", label='white line')
        #
        #
        # new_handles = [three_line, black_line, zero_line]
        new_labels = ["neg. feedback", "pos. feedback", "no feedback"]
        # plt.legend(handles=new_handles, labels=new_labels,loc=(0.2, 0.44))

        if gamma < 1:
            saving_string =r"/home/brunner/Documents/Current work/04122020/" + "kin_stan_neg_h1_f_" + str(global_df["IL-2_Tsec_fraction"].unique()[0]) + ".png"
        elif gamma > 1:
            saving_string =r"/home/brunner/Documents/Current work/04122020/" + "kin_stan_pos_h1_f_" + str(global_df["IL-2_Tsec_fraction"].unique()[0]) + ".png"

        else:
            raise ValueError
        fig.savefig(saving_string, bbox_inches='tight')
        plt.show()
