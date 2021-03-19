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
# D_to_iter_over = ['$10.0$']  # ['$1.0$', '$10.0$', '$100.0$']

common_y_axis = False
save_plots = False

sns.set(rc={'figure.figsize':(6,5)})
sns.set_style("ticks")
sns.set_context("talk", font_scale=1.5, rc={"lines.linewidth": 5})
xscale = "linear"
yscale = "linear"

error_band = True
error_band_style = "band"

if "IL-2_fraction" in group_variables:
    data_start = None
    data_end = None
else:
    data_start = None
    data_end = None

# plot_variables = ["std", "concentration"]  # ["surf_c_std"]#, "std", "gradient", "concentration", "surf_c"]
plot_variables = ["surf_c", "surf_c_std"]  # ["surf_c_std"]#, "std", "gradient", "concentration", "surf_c"]
# plot_variables = ["gradient"]

get_dataframes = []  # [[]]*3
saving_dataframe = pd.DataFrame()
if group_variables[0] == "IL-2_fraction":
    myRangeInt = 5
else:
    myRangeInt = 1
for j in range(myRangeInt):
    get_dataframes.append([])
    print("reading in run", j)
    if group_variables[0] == "IL-2_fraction":
        # path = "/extra/brunner/10x10x10/q_fraction_exact/run" + str(j) + "/"
        path = "/extra/brunner/thesis/static/q_fraction_multi_medium/run" + str(j) + "/"
        # path = "/extra/brunner/thesis/static/q_fraction/"
        ext_cache = "/extra/brunner/10x10x10/q_fraction/run_ext_cache/"
    elif group_variables[0] == "IL-2_sigma":
        # path = "/extra/brunner/10x10x10/R_lognorm/run" + str(j) + "/"
        path = "/extra/brunner/thesis/static/saturation/R_lognorm/"
        ext_cache = "/extra/brunner/10x10x10/R_lognorm/run_ext_cache/"
    else:
        print("Unknown grouping variable")
        exit()

    T = range(1)
    dt = 1

    global_df = pd.read_hdf(path + 'global_df.h5', mode="r")
    cell_df = pd.read_hdf(path + 'cell_df.h5', mode="r")

    global_df["surf_c_std"] = np.zeros(len(global_df["time"]))
    for scan_index in global_df["scan_index"].unique():
        for var in global_df[group_variables[0]].unique():
            global_df.loc[(global_df[group_variables[0]] == var), "surf_c_std"] = cell_df.loc[
                cell_df[group_variables[0]] == var, "IL-2_surf_c"].std()

    get_dataframes[j] = [global_df[data_start:data_end], cell_df]


plotting_df = pd.concat((get_dataframes[x][0] for x in range(len(get_dataframes))))
cell_concat_df = pd.concat((get_dataframes[x][1] for x in range(len(get_dataframes))))
if "IL-2_sigma" in group_variables:
    plotting_df["IL-2_D"] = "$10.0$"

if error_band == False:
    if "IL-2_sigma" in group_variables:
        plotting_df = cell_concat_df.groupby(["IL-2_sigma"], as_index=False).mean()
    else:
        plotting_df = cell_concat_df.groupby(group_variables, as_index=False).mean()
    plotting_df.reset_index(inplace=True)

    global_df = pd.concat((get_dataframes[x][0] for x in range(len(get_dataframes)))).groupby([group_variables[0]], as_index=False).mean()
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

######################### PLOTS #################################
# format entries
if group_variables[0] == "IL-2_fraction":
    plotting_df[group_variables[1]] = ["$%s$" % round(x, 3) for x in plotting_df[group_variables[1]]]

if "IL-2_sigma" in group_variables:
    plotting_df["IL-2_sigma"] = plotting_df["IL-2_sigma"].div(20000)


# do the plotting depending on what is to be plotted

ax1 = ax2 = ax3 = ax4 = ax5 = None
fig = plt.figure()
if "surf_c_std" in plot_variables:
    # plotting_df["surf_c_std"] = 1 * plotting_df["surf_c_std"]  # /plotting_df["surf_c"]
    if len(plot_variables) == 1:
        ax1 = sns.lineplot(x=group_variables[0], y="surf_c_std",
                           data=plotting_df[data_start:data_end].loc[plotting_df[data_start:data_end]], label="surface c. std.",
                           legend=False, color="black", err_style=error_band_style)
        ax1.set_ylabel("standard deviation in nM", color="black")
        ax1.set(yscale=yscale, xscale=xscale)
    else:
        plotting_df["surf_c_std"] *= 1e3
        ax1 = sns.lineplot(x=group_variables[0], y="surf_c_std",
                           data=plotting_df[data_start:data_end].loc[plotting_df[data_start:data_end]], label="std. surf.",
                           legend=False, color="black", err_style=error_band_style)
        ax1.set(yscale=yscale, xscale=xscale)
        # ax1.set(yticks = [])
if "std" in plot_variables:
    # plotting_df["std_norm"] = plotting_df["SD"]  # /plotting_df["concentration"]
    if len(plot_variables) == 1:
        ax2 = sns.lineplot(x=group_variables[0], y="std_norm",
                           data=plotting_df[data_start:data_end].loc[plotting_df[data_start:data_end]], label="concentration std", legend=False,
                           color="orange", err_style=error_band_style)
        ax2.set_ylabel("standard deviation in nM", color="black")
        ax2.set(yscale=yscale, xscale=xscale)
    else:
        plotting_df["SD"] *= 1e3
        ax2 = sns.lineplot(x=group_variables[0], y="SD",
                           data=plotting_df[data_start:data_end], label="concentration std", legend=False, err_style=error_band_style)
        ax2.set(yscale=yscale, xscale=xscale)
        # ax2.set_ylabel("standard deviation in nM", color="black")
if "concentration" in plot_variables:
    if len(plot_variables) == 1:
        ax3 = sns.lineplot(x=group_variables[0], y="Concentration",
                           data=plotting_df[data_start:data_end].loc[plotting_df[data_start:data_end]], label="avg. conc.",
                           legend=False, err_style=error_band_style)  # ,legend=False)
        ax3.set(ylabel="mean concentration in nM", yscale=yscale, xscale=xscale)  # , ylim=(0.275,0.45))
    else:
        plotting_df["Concentration"] *= 1e3
        ax3 = sns.lineplot(x=group_variables[0], y="Concentration",
                           data=plotting_df[data_start:data_end].loc[plotting_df[data_start:data_end]], label="avg. conc.",
                           legend=False, err_style=error_band_style)  # ,legend=False)
        ax3.set(yscale=yscale, xscale=xscale)
        if group_variables[0] == "IL-2_sigma":
            # pass
            ax3.set(ylabel = "cytokine conc. [pM]", xlabel = "Receptor expression heterogeneity", xlim=(-0.06,1.08))
        else:
            # ax3.set(xticks = [0.5, 1], ylim=(2e-3, 1e-1), xlabel = "fraction of sec. cells", ylabel = "nM")
            ax3.set(xticks=[0,0.5, 1], xlabel="fraction of sec. cells", ylabel="pM")

if "surf_c" in plot_variables:
    if len(plot_variables) == 1:
        ax4 = sns.lineplot(x=group_variables[0], y="surf_c",
                           data=plotting_df[data_start:data_end].loc[plotting_df[data_start:data_end]], label="avg. surf.",
                           legend=False, err_style=error_band_style)
        # ax1.errorbar(np.linspace(0.00005,1.0,20), plotting_df[5:]["surf_c"], yerr=plotting_df[5:]["sd"], fmt='-o')
        ax4.set(ylabel="mean surf. c. in nM")  # , ylim=(0.275,0.45))
    else:
        plotting_df["surf_c"] *= 1e3
        ax4 = sns.lineplot(x=group_variables[0], y="surf_c",
                           data=plotting_df[data_start:data_end].loc[plotting_df[data_start:data_end]], label="avg. surf.",
                           legend=False, color="darkgray", err_style=error_band_style)
        ax4.set(yscale=yscale, xscale=xscale, ylim=(-3,73))
        if group_variables[0] == "IL-2_sigma":
            # pass
            ax4.set(ylabel = "cytokine conc. [pM]", xlabel = "Receptor expression heterogeneity", xlim=(-0.06,1.08))
        else:
            # ax4.set(xticks = [0.5, 1], ylim=(2e-3, 1e-1), xlabel = "fraction of sec. cells", ylabel = "nM")
            ax4.set(xticks=[0,0.5, 1], xlabel="fraction of sec. cells", ylabel="cytokine (pM)")

# create and place the legend

if len(plot_variables) == 1:
    plt.legend()
else:
    if ax1: #surf_c_std
        handles, labels = ax1.get_legend_handles_labels()
        if group_variables[0] == "IL-2_sigma":
            ax1.figure.legend(reversed(handles), reversed(labels), loc='upper left', bbox_to_anchor=(0.34, 0.93), fancybox=True)  # (loc=(0.385,0.2))
        else:
            ax1.figure.legend(handles, labels, loc='upper left', bbox_to_anchor=(0.45, 1), fancybox=True)  # (loc=(0.385,0.2))
    elif ax2: #std
        handles, labels = ax2.get_legend_handles_labels()
        if group_variables[0] == "IL-2_sigma":
            ax2.figure.legend(reversed(handles), reversed(labels), loc='upper left', bbox_to_anchor=(0.56, 0.9), fancybox=True)
        else:
            ax2.figure.legend(reversed(handles), reversed(labels), loc='upper center', bbox_to_anchor=(0.56, 1.115), fancybox=True)

# save the plots

if save_plots == True:
    if group_variables[0] == "IL-2_fraction":
            temp_string = ""
            for var in plot_variables:
                temp_string += "_" + str(var)
            fig.savefig("plots/q_fraction_" + temp_string + ".svg", bbox_inches='tight')

    elif group_variables[0] == "IL-2_sigma":
            temp_string = ""
            for var in plot_variables:
                temp_string += "_" + str(var)
            fig.savefig("plots/R_lognorm_" + temp_string + ".png", bbox_inches='tight')
plt.show()
