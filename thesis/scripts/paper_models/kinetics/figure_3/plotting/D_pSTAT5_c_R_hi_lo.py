import pandas as pd
import getpass
import matplotlib.pyplot as plt
import numpy as np
import seaborn as sns
import os
from plotting_rc import rc_ticks

from thesis.scripts.paper_models.utilities.plot_helper import my_load_df, EC50_calculation

save_plot = True
hdd = "/extra2" if os.path.exists("/extra2") else "/extra"
user = getpass.getuser()
model_name = "satu_test"
name = "q_ramp"
path = "/{extra}/{u}/paper_models/kinetics/{mn}/{n}/".format(u=user, n=name, mn=model_name, extra=hdd)
saving_string =r"/home/brunner/Documents/Current work/2022_03_04/"
if not os.path.exists(saving_string):
    os.mkdir(saving_string)

# which metrics to plot and their plotting parameters
metrics = ["IL-2_surf_c", "IL-2_R", "IL-2_pSTAT5"]
y_scales = ["linear", "log", "linear"]
x_scales = ["linear", "linear", "linear"]

y_lims = [(-1, 25), (1.1e1, 3.1e5), (-0.1, 1.1)]
x_lims = [(None,None), (None,None), (None,None)]

x_labels = ["time (days)", "time (days)", "time (days)"]
y_labels = ["surface conc. (pM)", "receptors", "pSTAT5"]

# choose a scan and replicat index. For these plots the sec_cells_scan was used, from it a standard setup was chosen and
# plotted here.
scan_index = 0
replicat_index = 0

# plotting parameters
colours = ["grey" , "red"]
plot_ODE = False
show_legend = False


# define data start and end point
startingPoint = None
stoppingPoint = None
# maximum time in h
time_limit = 49


# which cells to show
show_Tsecs = False
show_Tregs = False
show_Ths = True

skip_first_time_index = True

# defines a outer layer of N cells to be ignored in plotting. Used to further limit unwanted boundary effects.
offset = 0

########################################################################################################################
########################################################################################################################
##                                              Plotting                                                              ##
########################################################################################################################
########################################################################################################################

# load runs, apply offset
cell_df, global_df =  my_load_df(path, offset=offset, run_range=[0], custom_ending = "")

# if plotting ODE this is necessary
try:
    cell_df["id_id"]
except KeyError:
    cell_df["id_id"] = cell_df["cell_index"]
    cell_df["type_name"] = cell_df["cell_type"]
    cell_df["IL-2_surf_c"] = cell_df["IL2"] * 1e9
    cell_df["IL-2_R"] = cell_df["R"]
    cell_df["IL-2_gamma"] = cell_df["gamma"]



x_axis = "time"

if plot_ODE == True:
    ODE_cell_df = pd.read_hdf(ODE_path + "cell_df.h5", mode="r")
    ODE_cell_df["time"] /= 3600
    ODE_cell_df = ODE_cell_df.loc[ODE_cell_df["time"] <= time_limit]
    ODE_cell_df["IL-2_surf_c"] = ODE_cell_df["IL-2_surf_c"] * 1e3

cell_df["time"] = cell_df["time"].div(3600 * 24)
cell_df["IL-2_surf_c"] = cell_df["IL-2_surf_c"].mul(1e3)
cell_df = cell_df.loc[cell_df["replicat_index"] == replicat_index]

if skip_first_time_index == True:
    cell_df = cell_df.loc[cell_df["time_index"] != 0]

try:
    cell_df["pSTAT5"] = cell_df["IL-2_pSTAT5"]
except:
    cell_df["pSTAT5"] = cell_df["IL-2_surf_c"] ** 3 / (
                (EC50_calculation(E_max=125e-12, E_min=0, k=860, N=1, R=cell_df["IL-2_R"]) * 1e12) ** 3 + cell_df[
                        "IL-2_surf_c"] ** 3).values


#%%
print("plotting")

timepoints = cell_df["time"].unique()[1:]
label = ""

cell_df["activated"] = "pSTAT5$^-$"
act_ids = cell_df.loc[(cell_df["scan_index"] == scan_index) & (cell_df["time"] == cell_df["time"].max()) & (cell_df["IL-2_R"] > 1e4) & (cell_df["type_name"] == "Th"), "id"]
cell_df.loc[cell_df["id"].isin(act_ids), "activated"] = "pSTAT5$^+$"

for m,metric in enumerate(metrics):
    sns.set_theme(context="talk", style="ticks", rc=rc_ticks)
    fig, ax = plt.subplots()
    ax = sns.lineplot(x="time", y=metric,
                        data=cell_df.loc[(cell_df["scan_index"] == scan_index) & (cell_df["type_name"] == "Th")].sort_values(by="time"), palette=colours, ci="sd", hue="activated")


    ax.set(yscale=y_scales[m], xscale=x_scales[m], ylim=y_lims[m], xlim=x_lims[m], xlabel=x_labels[m], ylabel=y_labels[m])

    plt.legend()

    if save_plot == True:
        fig.savefig(saving_string + f"{metric}_over_{x_axis}.pdf", bbox_inches='tight')
    plt.tight_layout()
    plt.show()