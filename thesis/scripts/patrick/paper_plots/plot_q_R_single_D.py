# -*- coding: utf-8 -*-

# import fenics as fcs

import numpy as np

import os.path

import sys

##
import pandas as pd
import subprocess
import matplotlib.pyplot as plt
import seaborn as sns
import xml.etree.ElementTree as ET

group_variables = ["IL-2_sigma", None]  # hue over the second one
# group_variables = ["IL-2_fraction", "IL-2_D"]
# subplot_style
D_to_iter_over = ['$10.0$']  # ['$1.0$', '$10.0$', '$100.0$']

common_y_axis = False
save_plots = True

sns.set(rc={'figure.figsize':(5,5)})
sns.set_style("ticks")
sns.set_context("talk", font_scale=1, rc={"lines.linewidth": 5})
xscale = "linear"
yscale = "linear"
if "IL-2_fraction" in group_variables:
    data_start = 1
    data_end = -1
else:
    data_start = None
    data_end = None

# plot_variables = ["std", "concentration"]  # ["surf_c_std"]#, "std", "gradient", "concentration", "surf_c"]
plot_variables = ["surf_c", "surf_c_std"]  # ["surf_c_std"]#, "std", "gradient", "concentration", "surf_c"]
# plot_variables = ["gradient"]

get_dataframes = []  # [[]]*3
saving_dataframe = pd.DataFrame()
if group_variables[0] == "IL-2_fraction":
    myRangeInt = 1
else:
    myRangeInt = 30
for j in range(myRangeInt):
    get_dataframes.append([])
    print("reading in run", j)
    if group_variables[0] == "IL-2_fraction":
        # path = "/extra/brunner/10x10x10/q_fraction_exact/run" + str(j) + "/"
        path = "/extra/brunner/thesis/static/q_fraction/"
        ext_cache = "/extra/brunner/10x10x10/q_fraction/run_ext_cache/"
    elif group_variables[0] == "IL-2_sigma":
        # path = "/extra/brunner/10x10x10/R_lognorm/run" + str(j) + "/"
        path = "/extra/brunner/para_handling/static/R_lognorm/run" + str(j) + "/"
        ext_cache = "/extra/brunner/10x10x10/R_lognorm/run_ext_cache/"
    else:
        print("Unknown grouping variable")
        exit()

    T = range(1)
    dt = 1

    global_df = pd.read_hdf(path + 'global_df.h5', mode="r")
    cell_df = pd.read_hdf(path + 'cell_df.h5', mode="r")
    # cell_stats_df = pd.read_hdf(path + 'cell_stats_dataframe_' + str(j) + '.h5', mode="r")

    # mySeries = pd.Series(np.ones(cell_df.shape[0]) * j).apply(lambda x: int(x))
    # cell_df["rep_index"] = mySeries

    '''
    column_name = "scan_index"
    saving_dataframe[column_name] = cell_df["scan_index"]

    column_name = "D"
    saving_dataframe[column_name] = cell_df["D"]

    column_name = "fraction"
    saving_dataframe[column_name] = cell_df["fraction"]

    column_name = "run_" + str(j) + "_surf_c_il2"
    saving_dataframe[column_name] = cell_df["surf_c_il2"]
    '''
    # add the surf_c std to global_df
    # temp_df = cell_df.groupby(["scan_index", "D"], as_index=False).std()["surf_c_il2"]
    # global_df_1["surf_c_std"] = cell_df.groupby(["scan_index", "D"], as_index=False).std()["surf_c_il2"]

    get_dataframes[j] = [global_df, cell_df]

'''
std = cell_stats_df.loc[:, (slice(None), "std")]
std.columns = std.columns.droplevel(1)
std.reset_index(inplace=True)
'''

concat_df = pd.concat((get_dataframes[x][0] for x in range(len(get_dataframes))))
# concat_df = pd.concat([get_dataframes[0][1], get_dataframes[1][1]])


# means = concat_df.groupby(["scan_index"], as_index=False).describe()
# means.columns = means.columns.droplevel(1)
# means.reset_index(inplace=True)

if "IL-2_sigma" in group_variables:
    plotting_df = concat_df.groupby(["IL-2_sigma"], as_index=False).mean()
    plotting_df["IL-2_D"] = "$10.0$"
else:
    plotting_df = concat_df.groupby(group_variables, as_index=False).mean()
plotting_df.reset_index(inplace=True)

global_df = global_df.groupby([group_variables[0]], as_index=False).mean()
try:
    global_df["surf_c_std"]
except KeyError:
    global_df["surf_c_std"] = np.zeros(len(global_df["time"]))
    for scan_index in global_df["scan_index"].unique():
        for var in global_df[group_variables[0]].unique():
            global_df.loc[(global_df[group_variables[0]] == round(var, 7)), "surf_c_std"] = cell_df.loc[
                cell_df[group_variables[0]] == round(var, 7), "IL-2_surf_c"].std()
if "IL-2_fraction" in group_variables:
    plotting_df = global_df
'''
saving_dataframe["mean_surf_c_il2"] =  saving_dataframe.drop(["scan_index", "D", "fraction"], axis=1).mean(axis=1)

si_to_frac_dict = {}
for scan_index in global_df["scanIndex"]:
    si_to_frac_dict[scan_index] = float(global_df.loc[global_df["scanIndex"] == scan_index]["fraction"])


for i in range(len(saving_dataframe)):
    saving_dataframe["fraction"].loc[i] = si_to_frac_dict[saving_dataframe["scan_index"][i]]
'''

######################### PLOTS #################################
# plotting_df[group_variables[1]] /= round(((10 ** 0.5 * 0.01) ** 2 * 0.1), 6)
if group_variables[0] == "IL-2_fraction":
    plotting_df[group_variables[1]] = ["$%s$" % round(x, 3) for x in plotting_df[group_variables[1]]]

if "IL-2_sigma" in group_variables:
    plotting_df["IL-2_sigma"] = plotting_df["IL-2_sigma"].div(20000)
# sns.set(rc={'figure.figsize':(x_size,y_size)})
# fig = plt.figure()
# plt.subplots_adjust(wspace=.4)



for single_D in D_to_iter_over:
    ax1 = ax2 = ax3 = ax4 = ax5 = None
    fig = plt.figure()
    if "surf_c_std" in plot_variables:
        # plotting_df["surf_c_std"] = 1 * plotting_df["surf_c_std"]  # /plotting_df["surf_c"]
        if len(plot_variables) == 1:
            ax1 = sns.lineplot(x=group_variables[0], y="surf_c_std",
                               data=plotting_df[data_start:data_end].loc[plotting_df[data_start:data_end]["IL-2_D"] == single_D], label="surface c. std",
                               legend=False)
            ax1.set_ylabel("standard deviation in nM", color="black")
            ax1.set(yscale=yscale, xscale=xscale)
        else:
            plotting_df["surf_c_std"] *= 1e3
            ax1 = sns.lineplot(x=group_variables[0], y="surf_c_std",
                               data=plotting_df[data_start:data_end].loc[plotting_df[data_start:data_end]["IL-2_D"] == single_D], label="surface concentration std",
                               legend=False)
            ax1.set(yscale=yscale, xscale=xscale)
            # ax1.set(yticks = [])
    if "std" in plot_variables:
        # plotting_df["std_norm"] = plotting_df["SD"]  # /plotting_df["concentration"]
        if len(plot_variables) == 1:
            ax2 = sns.lineplot(x=group_variables[0], y="std_norm",
                               data=plotting_df[data_start:data_end].loc[plotting_df[data_start:data_end]["IL-2_D"] == single_D], label="concentration std", legend=False,
                               color="orange")
            ax2.set_ylabel("standard deviation in nM", color="black")
            ax2.set(yscale=yscale, xscale=xscale)
        else:
            plotting_df["SD"] *= 1e3
            ax2 = sns.lineplot(x=group_variables[0], y="SD",
                               data=plotting_df[data_start:data_end].loc[plotting_df[data_start:data_end]["IL-2_D"] == single_D], label="concentration std", legend=False)
            ax2.set(yscale=yscale, xscale=xscale)
            # ax2.set_ylabel("standard deviation in nM", color="black")
    if "concentration" in plot_variables:
        if len(plot_variables) == 1:
            ax3 = sns.lineplot(x=group_variables[0], y="Concentration",
                               data=plotting_df[data_start:data_end].loc[plotting_df[data_start:data_end]["IL-2_D"] == single_D], label="average concentration",
                               legend=False)  # ,legend=False)
            ax3.set(ylabel="mean concentration in nM", yscale=yscale, xscale=xscale)  # , ylim=(0.275,0.45))
        else:
            plotting_df["Concentration"] *= 1e3
            ax3 = sns.lineplot(x=group_variables[0], y="Concentration",
                               data=plotting_df[data_start:data_end].loc[plotting_df[data_start:data_end]["IL-2_D"] == single_D], label="average concentration",
                               legend=False)  # ,legend=False)
            ax3.set(yscale=yscale, xscale=xscale)
            if group_variables[0] == "IL-2_sigma":
                # pass
                ax3.set(ylabel = "pM", xlabel = "Heterogeneity", xlim=(-0.06,1.08))
            else:
                # ax3.set(xticks = [0.5, 1], ylim=(2e-3, 1e-1), xlabel = "fraction of sec. cells", ylabel = "nM")
                ax3.set(xticks=[0,0.5, 1], xlabel="fraction of sec. cells", ylabel="pM")

    if "surf_c" in plot_variables:
        if len(plot_variables) == 1:
            ax4 = sns.lineplot(x=group_variables[0], y="surf_c",
                               data=plotting_df[data_start:data_end].loc[plotting_df[data_start:data_end]["IL-2_D"] == single_D], label="surface c.",
                               legend=False)
            # ax1.errorbar(np.linspace(0.00005,1.0,20), plotting_df[5:]["surf_c"], yerr=plotting_df[5:]["sd"], fmt='-o')
            ax4.set(ylabel="mean surf. c. in nM")  # , ylim=(0.275,0.45))
        else:
            plotting_df["surf_c"] *= 1e3
            ax4 = sns.lineplot(x=group_variables[0], y="surf_c",
                               data=plotting_df[data_start:data_end].loc[plotting_df[data_start:data_end]["IL-2_D"] == single_D], label="surface concentration",
                               legend=False)
            ax4.set(yscale=yscale, xscale=xscale)
            if group_variables[0] == "IL-2_sigma":
                # pass
                ax4.set(ylabel = "pM", xlabel = "Heterogeneity", xlim=(-0.06,1.08))
            else:
                # ax4.set(xticks = [0.5, 1], ylim=(2e-3, 1e-1), xlabel = "fraction of sec. cells", ylabel = "nM")
                ax4.set(xticks=[0,0.5, 1], xlabel="fraction of sec. cells", ylabel="pM")
    # if len(plot_variables) > 1:
    #     if ax1:
    #         ax1.set_ylabel("concentration in nM")
    #     elif ax2:
    #         ax2.set_ylabel("concentration in nM")
    #     elif ax3:
    #         ax3.set_ylabel("concentration in nM")
    #     elif ax4:
    #         ax4.set_ylabel("concentration in nM")
    # if common_y_axis == True:
    #     if group_variables[0] == "fraction":
    #         if ax1:
    #             ax1.set(ylim=(0, 0.0088))
    #         elif ax2:
    #             ax2.set(ylim=(0, 0.0088))
    #     else:
    #         if ax1:
    #             ax1.set(ylim=(0.0023, 0.005))
    #         elif ax2:
    #             ax2.set(ylim=(0, 0.005))
    #
    # if "gradient" in plot_variables:
    #     # twin object for two different y-axis on the sample plot
    #     if ax1:
    #         ax5 = ax1.twinx()
    #     elif ax2:
    #         ax5 = ax2.twinx()
    #     if ax5:
    #         sns.lineplot(x=group_variables[0], y="gradient", data=plotting_df[data_start:].loc[plotting_df[data_start:]["D"] == single_D],
    #                      label="mean gradient", ax=ax5, color="green", legend=False)
    #     else:
    #         ax5 = sns.lineplot(x=group_variables[0], y="gradient",
    #                            data=plotting_df[5:].loc[plotting_df[data_start:]["D"] == single_D],
    #                            label="mean gradient", color="green", legend=False)
    #
    #     if common_y_axis == True:
    #         if group_variables[0] == "fraction":
    #             ax5.set(ylim=(0, 0.00252))
    #         else:
    #             ax5.set(ylim=(0, 0.0017))
    #     # ax1.get_yaxis().set_visible(False)
    #     if ax1 or ax2:
    #         ax5.set_ylabel("mean gradient in nM/um", color="green")
    #     else:
    #         ax5.set_ylabel("mean gradient in nM/um")

    if single_D == '$10.0$':
        if len(plot_variables) == 1:
            plt.legend()
        else:
            if ax1:
                handles, labels = ax1.get_legend_handles_labels()
                if group_variables[0] == "IL-2_sigma":
                    ax1.figure.legend(reversed(handles), reversed(labels), loc='upper center', bbox_to_anchor=(0.56, 1.115), fancybox=True)  # (loc=(0.385,0.2))
                else:
                    ax1.figure.legend(reversed(handles), reversed(labels), loc='upper center', bbox_to_anchor=(0.56, 1.115), fancybox=True)  # (loc=(0.385,0.2))
            elif ax2:
                handles, labels = ax2.get_legend_handles_labels()
                if group_variables[0] == "IL-2_sigma":
                    ax2.figure.legend(reversed(handles), reversed(labels), loc='upper center', bbox_to_anchor=(0.56, 1.115), fancybox=True)
                else:
                    ax2.figure.legend(reversed(handles), reversed(labels), loc='upper center', bbox_to_anchor=(0.56, 1.115), fancybox=True)
            elif ax3:
                ax3.figure.legend()
            elif ax4:
                ax4.figure.legend()
            elif ax5:
                ax5.figure.legend(loc=(0.55, 0.7))
    else:
        if len(D_to_iter_over) == 1:
            if ax1:
                handles, labels = ax1.get_legend_handles_labels()
                if group_variables[0] == "IL-2_sigma":
                    ax1.figure.legend(reversed(handles), reversed(labels), loc=(0.4, 0.64))  # (loc=(0.385,0.2))
                else:
                    ax1.figure.legend(reversed(handles), reversed(labels), loc=(0.2, 0.64))  # (loc=(0.385,0.2))
            elif ax2:
                handles, labels = ax2.get_legend_handles_labels()
                if group_variables[0] == "IL-2_sigma":
                    ax2.figure.legend(reversed(handles), reversed(labels),loc=(0.45, 0.64))
                else:
                    ax2.figure.legend(reversed(handles), reversed(labels), loc=(0.2, 0.64))
                # ax2.figure.legend(
            elif ax3:
                ax3.figure.legend()
            elif ax4:
                ax4.figure.legend()
            elif ax5:
                ax5.figure.legend(loc=(0.55, 0.7))

    if save_plots == True:
        if group_variables[0] == "IL-2_fraction":
                temp_string = ""
                for var in plot_variables:
                    temp_string += "_" + str(var)
                fig.savefig("plots/q_fraction_" + single_D[1:-1] + temp_string + ".png", bbox_inches='tight')

        elif group_variables[0] == "IL-2_sigma":
                temp_string = ""
                for var in plot_variables:
                    temp_string += "_" + str(var)
                fig.savefig("plots/R_lognorm_" + single_D[1:-1] + temp_string + ".png", bbox_inches='tight')
    # plt.show()
