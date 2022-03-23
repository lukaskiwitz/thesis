import getpass
import matplotlib.pyplot as plt
import numpy as np
import seaborn as sns
import os
from plotting_rc import rc_ticks

sns.set_theme(context = "talk", style = "ticks", rc = rc_ticks)
fig,ax = plt.subplots()

from thesis.scripts.paper_models.utilities.plot_helper import my_load_df

def EC50_calculation(E_max, E_min, k, N, R):
    return (E_max * k ** N + E_min * R ** N) / (k ** N + R ** N)

save_plot = True

hdd = "extra2" if os.path.exists("/extra2") else "extra"
user = getpass.getuser()

model_name = "satu_test"
name = "q_ramp"
saving_string = "/{extra}/{u}/paper_models/kinetics/{mn}/".format(u=user, mn=model_name, extra = hdd)
if not os.path.exists(saving_string):
    os.mkdir(saving_string)
filename = "act_over_time_every_cell.svg"

# choose a scan and replicat index. For these plots the sec_cells_scan was used, from it a standard setup was chosen and
# plotted here.
scan_index = 0
replicat_index = 0

yscale = "linear"
xscale = "linear"

# define data start and end point
startingPoint = None
stoppingPoint = None
# plotting limits
xlim = (None,None)
ylim = (-0.1, 1.1)
# time and fraction max limit
max_time = 100
max_frac = -1
# colour
colours = ["black"]

# if line smoothing is desired, sfs are the smoothing factors
c_interpolate = False
c_sf = 0


path = "/{extra}/{u}/paper_models/kinetics/{mn}/{n}/".format(u=user, n=name, mn=model_name, extra=hdd)
cell_df, global_df =  my_load_df(path, offset=0, run_range=[0], custom_ending = "")
cell_df["IL-2_surf_c"] *= 1e3

# plot only a random fraction of the cells due to visibility
picked_fraction = 1/4
cell_ids = cell_df["id_id"].unique()
picked_cell_ids = np.random.choice(cell_ids, int(len(cell_ids) * picked_fraction), replace=False)
cell_df = cell_df.loc[cell_df["id_id"].isin(picked_cell_ids)]

# define which x_axis to plot
x_axis = "time"


cell_df = cell_df.loc[(cell_df["scan_index"] == scan_index)]

cell_df = cell_df.loc[cell_df["replicat_index"] == replicat_index]
print(cell_df["misc_EC50_N"].unique())

# the first time index does not have pSTAT5 calculated yet.
skip_first_time_index = True
show_Ths = True

plot_legend = False

########################################################################################################################
########################################################################################################################
##                                              Plotting                                                              ##
########################################################################################################################
########################################################################################################################


if skip_first_time_index == True:
    cell_df = cell_df.loc[cell_df["time_index"] > 0]

print("plotting")
for cell_df in [cell_df]:
    try:
        cell_df["pSTAT5"] = cell_df["IL-2_pSTAT5"]
    except:
        cell_df["pSTAT5"] = cell_df["IL-2_surf_c"] ** 3 / (
                (EC50_calculation(E_max=125e-12, E_min=0, k=860, N=1.5, R=cell_df["IL-2_R"]) * 1e12) ** 3 + cell_df[
            "IL-2_surf_c"] ** 3).values

    cell_df["time"] /= (3600 * 24)

    timepoints = cell_df["time"].unique()[1:]

    sns.set(rc={"lines.linewidth": 0.2})
    if show_Ths == True:
        ax_3 = sns.lineplot(x="time", y="pSTAT5", data=cell_df.loc[(cell_df["type_name"] == "Th")].sort_values(by="time")[startingPoint:stoppingPoint], estimator=None,
                            units="id_id", color=colours[0])


    plt.ylabel("pSTAT5")
    plt.xlabel("time (days)")
    plt.yscale(yscale)
    plt.xscale(xscale)
    plt.ylim(ylim)
    plt.xlim(xlim)

    plt.axhline(0.5, 0, 200, color="red")

    if save_plot == True:
        fig.savefig(saving_string + filename, bbox_inches='tight')
    plt.tight_layout()
    plt.show()