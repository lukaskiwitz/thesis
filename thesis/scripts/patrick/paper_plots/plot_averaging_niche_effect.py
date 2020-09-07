import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import seaborn as sns
import random
from my_debug import message
from scipy import spatial

#load the df

no_of_runs = 10
base_path = "/extra/brunner/thesis/kinetic/q_fraction_medium_g_scan_multi/"
cell_cell_distance = 20
gammas = [1/9, 1/8, 1/7, 1/6, 1/5, 1/4, 1/3, 1/2, 2, 3, 4, 5, 6, 7, 8, 9, 10]
# gammas = [2, 3, 4, 5, 6, 7, 8, 9, 10]
# gammas = [0.1, 0.5, 1/3, 1/4, 1/5, 1/6, 1/7, 1/8, 1/9, 2,10]
# gammas = [0.1,0.5,2,10]

# avg_hist_df = pd.read_hdf(base_path + str(no_of_runs) + '_runs_avg_hist_df.h5', mode="r")
runs_hist_df = pd.read_hdf(base_path + str(no_of_runs) + "_runs_hist_df.h5", mode="r")
# print("SHORTENING base runs hist df")
# runs_hist_df = runs_hist_df.transpose()[:6*4].transpose()

base_path_pos = "/extra/brunner/thesis/kinetic/q_fraction_medium_pos_g_scan_multi/"
base_path_neg = "/extra/brunner/thesis/kinetic/q_fraction_medium_neg_g_scan_multi/"


runs_hist_df_pos = pd.read_hdf("/extra/brunner/thesis/kinetic/q_fraction_medium_pos_g_scan_multi/" + "10_runs_hist_df.h5", mode="r")
runs_hist_df_neg = pd.read_hdf("/extra/brunner/thesis/kinetic/q_fraction_medium_neg_g_scan_multi/" + "10_runs_hist_df.h5", mode="r")
message("loaded runs_hist_df's from " + base_path)

runs_hist_df = pd.concat([runs_hist_df, runs_hist_df_pos, runs_hist_df_neg], axis=1, sort=False)
print("Concatenated the different runs")
#############################################################################################################################################
# Average niche over runs

# load average secretory distance and norm to cell-cell-distance
mean_sec_dist = float(np.loadtxt(base_path + str(no_of_runs) + "_runs_avg_sec_dist.txt"))/cell_cell_distance
mean_sec_dist_pos = float(np.loadtxt("/extra/brunner/thesis/kinetic/q_fraction_medium_pos_g_scan_multi/" + "6_runs_avg_sec_dist.txt"))/cell_cell_distance
mean_sec_dist_neg = float(np.loadtxt("/extra/brunner/thesis/kinetic/q_fraction_medium_neg_g_scan_multi/" + "6_runs_avg_sec_dist.txt"))/cell_cell_distance

mean_sec_dist = np.mean([mean_sec_dist, mean_sec_dist_pos, mean_sec_dist_neg])

niche_timepoint = 20
cut_off = -1

slope_thresholds = [0.02] #np.arange(0.01,0.08,0.01)

ylim = (0.95,None)

for slope_threshold in slope_thresholds:
    slope_threshold = round(slope_threshold,2)

    niche_size_at_gamma = [[] for gamma in gammas]
    niche_effect_at_gamma = [[] for gamma in gammas]
    for g,gamma in enumerate(gammas):
        # niche_size_at_gamma.append([])
        for run in range(runs_hist_df.columns[-1][0]+1):
            niche = 0
            for i in range(len(runs_hist_df.loc[niche_timepoint, (run,gamma)][:cut_off])):
                try:
                    if runs_hist_df.loc[niche_timepoint, (run,gamma)][i] - runs_hist_df.loc[niche_timepoint, (run,gamma)][i+1] < slope_threshold * runs_hist_df.loc[niche_timepoint, (run,gamma)][i+1]:
                        niche_size_at_gamma[g].append(i)
                        niche = i
                        break
                    # if runs_hist_df.loc[niche_timepoint, (run,gamma)][i+1] - runs_hist_df.loc[niche_timepoint, (run,gamma)][i] < slope_threshold* runs_hist_df.loc[niche_timepoint, (run,gamma)][i]:
                    #     niche_size_at_gamma[g].append(i+1)
                    #     niche = i
                    #     break
                except IndexError:
                    niche_size_at_gamma[g].append(i)
                    niche = i
            # calculate niche effect: take the timepoints niche and get the average concentration within that niche
            concentrations = []
            for k in range(niche):
                concentrations.append(runs_hist_df.loc[niche_timepoint, (run,gamma)][k])
            # print(concentrations)
            try:
                global_df = pd.read_hdf(base_path + "run" + str(run) + "/global_df.h5", mode="r")
                average_concentration = global_df.loc[(global_df["time"] == niche_timepoint) & (
                            global_df["IL-2_gamma"] == gamma), "surf_c"].values[0] * 1e3

            except IndexError:
                try:
                    global_df = pd.read_hdf(base_path_pos + "run" + str(run) + "/global_df.h5", mode="r")
                    average_concentration = global_df.loc[(global_df["time"] == niche_timepoint) & (
                                global_df["IL-2_gamma"] == gamma), "surf_c"].values[0] * 1e3

                except IndexError:
                    global_df = pd.read_hdf(base_path_neg + "run" + str(run) + "/global_df.h5", mode="r")
                    average_concentration = global_df.loc[(global_df["time"] == niche_timepoint) & (
                                global_df["IL-2_gamma"] == gamma), "surf_c"].values[0] * 1e3

            niche_effect_at_gamma[g].append(np.mean(concentrations)/average_concentration)

    #get errors:
    niche_size_error = [[] for x in range(len(niche_size_at_gamma))]
    avg_niche_size = [np.mean(x) for x in niche_size_at_gamma]

    import scipy.stats as st

    # print(st.t.interval(0.95, len(niche_size_at_gamma) - 1, loc=np.mean(niche_size_at_gamma), scale=st.sem(niche_size_at_gamma[0])))


    for i in range(len(niche_size_at_gamma)):
        # niche_size_error[i].append(np.min(np.array(niche_size_at_gamma[i]) - avg_niche_size[i]))
        # niche_size_error[i].append(np.max(np.array(niche_size_at_gamma[i]) - avg_niche_size[i]))
        # or confidence interval with p=0.95
        # conf_int = st.t.interval(0.95, len(niche_size_at_gamma[i]) - 1, loc=np.mean(niche_size_at_gamma[i]), scale=st.sem(niche_size_at_gamma[i]))
        # niche_size_error[i] = [conf_int[0], conf_int[1]]
        # or std
        niche_size_error[i].append(np.std(niche_size_at_gamma[i])/np.sqrt(no_of_runs))
        niche_size_error[i].append(np.std(niche_size_at_gamma[i])/np.sqrt(no_of_runs))

    niche_size_at_gamma = [np.mean(x) for x in niche_size_at_gamma]

    niche_effect_at_gamma = np.array(niche_effect_at_gamma)*100
    niche_effect_error = [[] for x in range(len(niche_effect_at_gamma))]
    avg_niche_effect = [np.mean(x) for x in niche_effect_at_gamma]
    for i in range(len(niche_effect_at_gamma)):
        # niche_effect_error[i].append(np.min(np.array(niche_effect_at_gamma[i]) - avg_niche_effect[i]))
        # niche_effect_error[i].append(np.max(np.array(niche_effect_at_gamma[i]) - avg_niche_effect[i]))
        # or confidence interval with p=0.95
        # conf_int = st.t.interval(0.95, len(niche_effect_at_gamma[i]) - 1, loc=np.mean(niche_effect_at_gamma[i]), scale=st.sem(niche_effect_at_gamma[i]))
        # niche_effect_error[i] = [conf_int[0], conf_int[1]]
        # or std
        niche_effect_error[i].append(np.std(niche_effect_at_gamma[i])/np.sqrt(no_of_runs))
        niche_effect_error[i].append(np.std(niche_effect_at_gamma[i])/np.sqrt(no_of_runs))

    niche_effect_at_gamma = [np.mean(x) for x in niche_effect_at_gamma]
    #############################################################################################################################################

    # Plotting

    # fig = plt.figure()
    sns.set_style("ticks")
    sns.set_context("talk", font_scale=1.5, rc={"lines.linewidth": 5})
    fig, ax = plt.subplots(figsize=(6, 5))

    palette = sns.color_palette("bwr", len(gammas)+1)
    palette = [i for i in reversed(palette)]
    #
    x = niche_effect_at_gamma
    y = niche_size_at_gamma


    for i in range(len(x)):
        if gammas[i] < 1:
            myColor = "red"
        else:
            myColor = "black"
        plt.scatter(x[i], y[i], color=palette[i])
        plt.errorbar(x[i], y[i], xerr=[[np.abs(x)] for x in niche_effect_error[i]], yerr=[[np.abs(x)] for x in niche_size_error[i]], linestyle='', c=palette[i])
    plt.scatter(216.35519872163854, 1.6372549019607845, color="black")
    plt.errorbar(216.35519872163854, 1.6372549019607845,xerr=13.160380344055845, yerr=0.22201019379779285, color="black")
    sns.set_context(rc={"lines.linewidth": 2})
    # plt.axhline(mean_sec_dist/2, color="black")
    # plt.show()
        # exit()
    plt.xlabel("niche effect ($\%$)")
    plt.ylabel("niche size")
    # plt.title(str(int(slope_threshold*100)) + "%")
    plt.ylim(ylim)
    plt.yticks([1,2,3])

    # import matplotlib.lines as mlines
    #
    # zero_line = mlines.Line2D([], [], color=palette[0], label='lightcoral line')
    # one_line = mlines.Line2D([], [], color=palette[1], label='red line')
    # two_line = mlines.Line2D([], [], color=palette[2], label='lightsteelblue line')
    # three_line = mlines.Line2D([], [], color=palette[-1], label='blue line')
    #
    # grey_line = mlines.Line2D([], [], color="grey", label='grey line')
    # grey_dashed_line = mlines.Line2D([], [], color="grey", label='grey dashed line', linestyle="dashed")
    # white_line = mlines.Line2D([], [], color="white", label='white line')
    #
    #
    # new_handles = [three_line, zero_line]
    # new_labels = ["neg. feedback", "pos. feedback"]
    # # new_handles = [three_line, two_line, black_line, one_line, zero_line]
    # # new_labels = ["strong negative feedback", "negative feedback", "no feedback", "positive feedback", "strong positive feedback"]
    # plt.legend(handles=new_handles, labels=new_labels, loc='upper right', fancybox=True)
    fig.savefig("plots/" + "averaged_niche_effect_" + str(slope_threshold) + ".svg", bbox_inches='tight')
    plt.show()
