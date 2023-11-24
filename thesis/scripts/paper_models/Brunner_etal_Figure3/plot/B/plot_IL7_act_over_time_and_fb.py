import getpass
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import seaborn as sns
import os
from thesis.scripts.paper_models.utilities.plotting_rc import rc_ticks
from scipy.interpolate import UnivariateSpline
from thesis.scripts.paper_models.utilities.plot_helper import my_load_df
from thesis.main.my_debug import message

def EC50_calculation(E_max, E_min, k, N, R):
    return (E_max * k ** N + E_min * R ** N) / (k ** N + R ** N)

save_plot = True
hdd = "extra2" if os.path.exists("/extra2") else "extra"
user = getpass.getuser()

prefix = "neg"
# prefix = "pos"

if prefix == "pos":
    model_name = "feedback_scan_5"
    name = "dataframes_positive_0.176_time_till_4_days_timeseries"
    fraction = 0
    ODE_path = f"/{hdd}/brunner/paper_models/ODE/saturated/kinetics/gamma_10/"
elif prefix == "neg":
    model_name = "feedback_scan_4"
    name = "dataframes_negative_timeseries"
    fraction = 4
    ODE_path = f"/{hdd}/brunner/paper_models/ODE/saturated/kinetics/gamma_0.1/"


saving_string =r"/home/brunner/Documents/Current work/2023_11_17/"
if not os.path.exists(saving_string):
    os.mkdir(saving_string)

filename = prefix + "_act_over_time_and_fb"

ODE_cell_df =  pd.read_hdf(ODE_path + "cell_df" + ".h5", mode="r")
ODE_cell_df["IL-2_surf_c"] *= 1e3

base_path = "/{extra}/{u}/paper_models/kinetics/{mn}/{n}/".format(u=user, mn=model_name, extra=hdd, n=name)

spatial_cell_df, global_df =  my_load_df(base_path, offset=0, custom_ending = "_combined")
try:
    spatial_cell_df["IL-2_gamma"]
except KeyError:
    spatial_cell_df["IL-2_gamma"] = spatial_cell_df["misc_gamma"]
    spatial_cell_df["IL-2_Tsec_fraction"] = spatial_cell_df["fractions_Tsec"]

spatial_cell_df["IL-2_surf_c"] *= 1e3
spatial_cell_df = spatial_cell_df.loc[spatial_cell_df["time_index"] != 0] #first time index does not contain pSTAT5 values
#%%
print("plotting")
f_cell_df = spatial_cell_df.loc[(spatial_cell_df["IL-2_Tsec_fraction"] == np.sort(spatial_cell_df["IL-2_Tsec_fraction"].unique())[fraction])]
f_ODE_cell_df = ODE_cell_df.loc[(ODE_cell_df["IL-2_Tsec_fraction"] == (0.1 if prefix == "pos" else 0.064))]

rc_ticks['figure.figsize'] = (1.67475 * 1.15, 1.386 * 0.5)
sns.set_theme(context="talk", style="ticks", rc=rc_ticks)
fig, ax = plt.subplots()
amount_of_reps = 20
for c,cell_df in enumerate([f_cell_df, f_ODE_cell_df]):
    gammas = cell_df["IL-2_gamma"].unique()
    gammas = np.sort(gammas[np.where(gammas != 1.)])
    cell_df["time"] /= 3600

    if prefix == "neg":
        cmap_name = "Blues"
        color = "blue"
        reverse = True
        gammas = gammas[::-1] # reversed
    else:
        cmap_name = "Reds"
        color = "red"
        reverse = False
        gamma_range = np.around(np.logspace(-2, np.log10(12), 30), 3)
        gamma_range = gamma_range[gamma_range > 1]
        gammas = gamma_range
    if c == 1:
        color = "black"
        reverse = False
        gammas = cell_df["IL-2_gamma"].unique()

    try:
        cell_df["pSTAT5"] = cell_df["IL-2_pSTAT5"]
    except:
        cell_df["pSTAT5"] = cell_df["IL-2_surf_c"] ** 3 / (
                (EC50_calculation(E_max=125e-12, E_min=0, k=860, N=1.5, R=cell_df["IL-2_R"]) * 1e12) ** 3 +
                cell_df[
                    "IL-2_surf_c"] ** 3).values

    alphas = np.logspace(-1, 0, len(gammas) + 1) if c != 1 else [1] #1 for ODE
    for g, gamma in enumerate(gammas):
        frac_df = cell_df.loc[(cell_df["IL-2_gamma"] == gamma)]

        rep_act = []
        for rep in np.random.choice(frac_df["replicat_index"].unique(), amount_of_reps):
            act = []
            rep_df = frac_df.loc[(frac_df["replicat_index"] == rep) & (frac_df["type_name"] == "Th")]
            if len(rep_df["time"].unique()) == len(frac_df["time"].unique()):
                for t, time in enumerate(rep_df["time"].unique()):
                    timed_df = rep_df.loc[(rep_df["time"] == time)]
                    act.append(len(timed_df.loc[(timed_df["pSTAT5"] > 0.5)]) / len(timed_df) * 100)
                rep_act.append(act)
        sns.lineplot(x = frac_df["time"].unique(), y = np.mean(rep_act, axis = 0), label=gamma, color=color, alpha=alphas[g], legend=False)
        plt.fill_between(frac_df["time"].unique(), np.mean(rep_act, axis = 0) - np.std(rep_act, axis = 0)/np.sqrt(amount_of_reps),
                         np.mean(rep_act, axis = 0) + np.std(rep_act, axis = 0)/np.sqrt(amount_of_reps), alpha=0.3 * alphas[g], color=color, linewidth=0)

        plt.ylabel(r"%pSTAT5$^+$")
        plt.xlabel("time (d)")
        plt.yscale("linear")
        plt.xscale("linear")
        plt.ylim(0, 100)
        plt.yticks([0, 50, 100])
        if prefix == "pos":
            plt.xticks([0, 48, 96], [0,2,4])
            plt.xlim((0, 4 * 24))
        else:
            plt.xticks([0, 4 * 24, 8 * 24], [0,4,8])
            plt.xlim((0, 8 * 24))

norm = plt.Normalize(np.min(f_cell_df["IL-2_gamma"].unique()), np.max(f_cell_df["IL-2_gamma"].unique()))
sm = plt.cm.ScalarMappable(cmap=cmap_name + "_r" if reverse == True else cmap_name, norm=norm)
sm.set_array([])
cbar = fig.colorbar(sm, label="feedback", ticks=[np.min(f_cell_df["IL-2_gamma"].unique()), np.max(f_cell_df["IL-2_gamma"].unique())])
if prefix == "neg":
    cbar.ax.set_yticklabels([round(np.max(f_cell_df["IL-2_gamma"].unique()),5), round(np.min(f_cell_df["IL-2_gamma"].unique()),5)])
else:
    cbar.ax.set_yticklabels([round(np.min(f_cell_df["IL-2_gamma"].unique()),5), round(np.max(f_cell_df["IL-2_gamma"].unique()),5)])

if save_plot == True:
    fig.savefig(saving_string + filename + ".pdf", bbox_inches='tight', transparent=True)
plt.tight_layout()
plt.show()