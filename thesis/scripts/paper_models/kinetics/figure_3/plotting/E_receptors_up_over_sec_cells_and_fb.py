import getpass
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import seaborn as sns
import os

sns.set_theme(context="talk", style="ticks", rc=rc_ticks)
fig, ax = plt.subplots()

save_plot = True

hdd = "extra2" if os.path.exists("/extra2") else "extra"
user = getpass.getuser()
path = "/{extra}/{u}/paper_models/kinetics/feedback_scan/".format(u=user, extra=hdd)
saving_string =r"/home/brunner/Documents/Current work/2022_03_04/CD25_up_over_sec_cells_and_fbs"

name = ""

yscale = "linear"
xscale = "linear"

# define data start and end point
startingPoint = None
stoppingPoint = None
# plotting limits
xlim = (0.001, 0.55)
ylim = (-2, 101)
# time and fraction max limit
max_time = 100
max_frac = -1
# colours, first = positive, second = negative
colours = ["black", "blue"]
# opacity for the ODE
ODE_alpha = 0.25
# if interpolation is desired
interpolate = True
sf = 0.003



spatial_cell_df, global_df = my_load_df(path, offset=0, run_range=[0], custom_ending = "_combined")
spatial_cell_df["IL-2_surf_c"] *= 1e3


path = f"/{hdd}/brunner/paper_models/ODE/saturated/kinetics/"
ODE_cell_df =  pd.read_hdf(path + "cell_df_combined" + ".h5", mode="r")
ODE_cell_df["IL-2_surf_c"] *= 1e3

plot_every_cell = True
plot_std_area = True
# define which scan to plot
x_axis = "time"
my_hue = "scan_index"
scan = 0
# the first time index does not have pSTAT calculated yet.
show_Ths = True
show_Tregs = False

plot_legend = False

colors = sns.color_palette("rocket", 5)
no_fb_color = colors[0]
colors = colors[1:]

########################################################################################################################
########################################################################################################################
##                                              Plotting                                                              ##
########################################################################################################################
########################################################################################################################

#%%

print("plotting")
plt.axhline(0, 0.01, 0.8, color=no_fb_color, label=0)
for c,cell_df in enumerate([ODE_cell_df]):
    cell_df = cell_df.loc[cell_df["time"] == cell_df["time"].max()]
    try:
        cell_df["R_up_factor"] = cell_df["IL-2_R"] / [item[-1] for item in cell_df["R_start_R_start"].values]
    except IndexError:
        cell_df["R_up_factor"] = cell_df["IL-2_R"] / cell_df["R_start_R_start"]

    rep_act = []
    # desired amount of runs to average over. Chosen at 10 because RD-system was also 10
    for rep in range(10):
        act = []
        for f, frac in enumerate(cell_df["fractions_Tsec"].unique()):
            frac_df = cell_df.loc[(cell_df["IL-2_gamma"] == 50)]
            # this rep + 10 weirdness was necessary because I had to stich multiple dataframes together. Can be removed if everything is one run
            frac_df = frac_df.loc[(frac_df["replicat_index"] == rep + 10) & (frac_df["fractions_Tsec"] == frac) & (frac_df["type_name"] == "Th")]
            act.append(len(frac_df.loc[(frac_df["R_up_factor"] > (0.02 * frac_df["IL-2_gamma"].unique()[0] + 3))]) / len(frac_df))
        rep_act.append(act)

    plot_act = np.mean(np.array(rep_act).T, axis=1)

    pos_stds = []
    neg_stds = []
    for i, x in enumerate(np.array(rep_act).T):
        pos_stds.append(np.std(x[np.where(x > plot_act[i])]))
        neg_stds.append(np.std(x[np.where(x <= plot_act[i])]))
    pos_stds = np.array(pos_stds) * 100
    neg_stds = np.array(neg_stds) * 100
    plot_act *= 100
    # sorting = np.where(cell_df["fractions_Tsec"].unique() >= 0.01)
    if interpolate == False:
        plt.plot(cell_df["fractions_Tsec"].unique(), plot_act, "--", color="grey", label="well-mixed")
        plt.fill_between(cell_df["fractions_Tsec"].unique(), plot_act - neg_stds / 2,
                         plot_act + pos_stds / 2, alpha=0.3, color="grey")
    else:
        my_interpolation(cell_df["fractions_Tsec"].unique(), plot_act, smoothing_factor = sf, plot=True, label="well-mixed", fill=True, pos_std=pos_stds, neg_std=neg_stds, color="grey")


for g, gamma in enumerate([25, 50, 75, 100]):
    print(f"gamma = {gamma}")
    for c,cell_df in enumerate([spatial_cell_df]):
        c += 1
        cell_df = cell_df.loc[(cell_df["time"] == cell_df["time"].max()) ]
        try:
            cell_df["R_up_factor"] = cell_df["IL-2_R"] / [item[-1] for item in cell_df["R_start_R_start"].values]
        except IndexError:
            cell_df["R_up_factor"] = cell_df["IL-2_R"] / cell_df["R_start_R_start"]

        act = []
        ste = []
        for f, frac in enumerate(cell_df["fractions_Tsec"].unique()):
            frac_df = cell_df.loc[(cell_df["fractions_Tsec"] == frac) &
                                  (cell_df["type_name"] == "Th")]
            frac_df = frac_df.loc[
                (frac_df["IL-2_gamma"] == np.abs(frac_df["IL-2_gamma"].unique() - gamma).min() + gamma)]
            # scaling the receptor upregulation definition with feedback strength
            act.append(
                len(frac_df.loc[(frac_df["R_up_factor"] > (0.02 * frac_df["IL-2_gamma"].unique()[0] + 3))]) / len(
                    frac_df))

            ste.append(len(frac_df.loc[(frac_df["R_up_factor"] > (0.02 * frac_df["IL-2_gamma"].unique()[0] + 3)), "replicat_index"].unique()))

        xs = cell_df["fractions_Tsec"].unique()
        sorting = np.argsort(xs)
        act = np.array(act) * 100
        ste = np.array(ste)
        if interpolate == False:
            plt.plot(xs[sorting], act[sorting], label=gamma, color=colors[g])
            plt.fill_between(xs[sorting], act[sorting] - act[sorting] / ste,
                             act[sorting] + act[sorting] / ste, alpha=0.3, color=colors[g])
        else:
            my_interpolation(xs[sorting], act[sorting], smoothing_factor = sf, plot=True, label=gamma, fill=True, std= act[sorting] / ste, color=colors[g])


plt.ylabel("%CD25$^+$")
plt.xlabel("fraction of secreting cells")
plt.yscale(yscale)
plt.xscale(xscale)
plt.ylim(ylim)
plt.xlim(xlim)

plt.legend(title="feedback strength")

if save_plot == True:
    fig.savefig(saving_string + ".pdf", bbox_inches='tight')
plt.tight_layout()
plt.show()