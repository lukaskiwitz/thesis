import getpass
import matplotlib.pyplot as plt
import numpy as np
import seaborn as sns
import os

from scipy.interpolate import UnivariateSpline
from plotting_rc import rc_ticks
from thesis.scripts.paper_models.utilities.plot_helper import my_load_df, my_interpolation, EC50_calculation

save_plot = True

hdd = "extra2" if os.path.exists("/extra2") else "extra"
user = getpass.getuser()
model_name = "Tsec_scan_higher_EC50_N"
saving_string = r"/home/brunner/Documents/Current work/2022_03_04/"

name = ""
# which cvs to calculate, which color the feedback line should have and what the labels should be
calculate_cv_of = [ "IL-2_R", "IL-2_surf_c"]
colors = ["blue", "grey"]
labels = ["receptors", "surface conc."]


yscale = "log"
xscale = "linear"

# define data start and end point
startingPoint = None
stoppingPoint = None
# plotting limits
x_lims = [(-0.01, 0.45), (-0.01, 0.45)]
y_lims = [(0.095, 10.5), (0.095, 10.5)]
# time and fraction max limit
time_limit = 49
max_frac = -1
# if interpolation is desired
interpolate = True
sf = 0.003


path = "/{extra}/{u}/paper_models/kinetics/{mn}/{n}/".format(u=user, n=name, mn=model_name, extra=hdd)
spatial_cell_df, global_df =  my_load_df(path, offset=0, run_range=[0], custom_ending = "")
spatial_cell_df["IL-2_surf_c"] *= 1e3
spatial_cell_df  = spatial_cell_df.loc[spatial_cell_df["time"] < time_limit*3600]

# define which scan to plot
x_axis = "IL-2_Tsec_fraction"

# the first time index does not have pSTAT calculated yet.
skip_first_time_index = True
show_Ths = True
show_Tregs = False

plot_legend = False

########################################################################################################################
########################################################################################################################
##                                              Plotting                                                              ##
########################################################################################################################
########################################################################################################################
#%%


if skip_first_time_index == True:
    spatial_cell_df = spatial_cell_df.loc[spatial_cell_df["time_index"] != 0]
    # neg_cell_df = neg_cell_df.loc[neg_cell_df["time_index"] != 0]

print("plotting")
max_cell_df = spatial_cell_df.loc[spatial_cell_df["time"] == spatial_cell_df["time"].max()]
min_cell_df = spatial_cell_df.loc[spatial_cell_df["time"] == spatial_cell_df["time"].min()]
# can also be done for various other dfs
for cv_index, cv in enumerate(calculate_cv_of):
    sns.set_theme(context="talk", style="ticks", rc=rc_ticks)
    fig, ax = plt.subplots()

    for c, cell_df in enumerate([max_cell_df, min_cell_df]):
        try:
            cell_df["pSTAT5"] = cell_df["IL-2_pSTAT5"]
        except:
            cell_df["pSTAT5"] = cell_df["IL-2_surf_c"] ** 3 / (
                    (EC50_calculation(E_max=125e-12, E_min=0, k=860, N=1.5, R=cell_df["IL-2_R"]) * 1e12) ** 3 + cell_df[
                            "IL-2_surf_c"] ** 3).values

        cell_df = cell_df.loc[cell_df[x_axis] > 0.0005]
        tmp_calc_cvs = [cv, x_axis]
        rep_cvs = []
        for rep in cell_df["replicat_index"].unique():
            cvs = []
            tmp_df = cell_df.loc[cell_df["replicat_index"] == rep, tmp_calc_cvs]
            for x, axs in enumerate(cell_df[x_axis].unique()):
                timed_df = tmp_df.loc[tmp_df[x_axis] == axs]
                result = np.std(timed_df.to_numpy().T, axis=1) / np.mean(timed_df.to_numpy().T, axis=1)
                cvs.append(result[:-1])
            cvs = np.array(cvs)
            rep_cvs.append(cvs)
        plot_cvs = np.mean(np.array(rep_cvs).T, axis=2)
        plot_stds = np.std(np.array(rep_cvs).T, axis=2)

        # split into positive and negative stds
        for e, entry in enumerate(plot_cvs):
            pos_stds = []
            neg_stds = []
            for i, x in enumerate(np.array(rep_cvs).T[e]):
                pos_stds.append(np.std(x[np.where(x > plot_cvs[e][i])]))
                neg_stds.append(np.std(x[np.where(x <= plot_cvs[e][i])]))
            # to make it standard error
            pos_stds = np.array(pos_stds)/np.sqrt(len(cell_df["replicat_index"].unique()))
            neg_stds = np.array(neg_stds)/np.sqrt(len(cell_df["replicat_index"].unique()))
            sorting = np.argsort(cell_df[x_axis].unique())
            if interpolate == False:
                plt.plot(cell_df[x_axis].unique()[sorting], plot_cvs[e][sorting], color=colors[e], label="feedback")
                plt.fill_between(cell_df[x_axis].unique(), plot_cvs[e] - neg_stds,
                                 plot_cvs[e] + pos_stds , alpha=0.3, color=colors[c])
            else:
                my_interpolation(cell_df[x_axis].unique()[sorting], plot_cvs[e][sorting], smoothing_factor = 0.01, plot=True,
                             label="feedback", fill=True, neg_std=neg_stds, pos_std=pos_stds, color=colors[c])


    plt.ylabel("CV of " + labels[cv_index])
    plt.xlabel("fraction of secreting cells")
    plt.yscale(yscale)
    plt.xscale(xscale)
    plt.ylim(y_lims[cv_index])
    plt.xlim(x_lims[cv_index])

    plt.legend()
    plt.yticks([0.1, 1, 10], [0.1, 1, 10])
    plt.xticks([0, 0.2, 0.4], ["*", "0.2", "0.4"])

    if save_plot == True:
        fig.savefig(saving_string + f"{cv}_over_{x_axis}.pdf", bbox_inches='tight')
    plt.tight_layout()
    plt.show()
