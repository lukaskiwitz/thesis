import getpass
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import seaborn as sns
import os
from thesis.scripts.paper_models.utilities.plotting_rc import rc_ticks
from thesis.scripts.paper_models.utilities.plot_helper import my_load_df, ashmans_D, wangs_BI

def EC50_calculation(E_max, E_min, k, N, R):
    return (E_max * k ** N + E_min * R ** N) / (k ** N + R ** N)

save_plot = True

hdd = "extra2" if os.path.exists("/extra2") else "extra"
user = getpass.getuser()

setups = [
    "pos",
    "neg",
]

for prefix in setups:
    if prefix == "pos":
        model_name = "feedback_scan_for_time_plot/"
        name = "dataframes_positive_0.1_timeseries"  # michaelis
    elif prefix == "neg":
        model_name = "feedback_scan_4"
        name = "dataframes_negative_for_Fig3C_act_plot_2_steady_state" #michaelis

    saving_string = f"/home/brunner/Documents/Current work/2023_11_17/"

    path = "/{extra}/{u}/paper_models/kinetics/{mn}/{n}/".format(u=user, n=name, mn=model_name, extra=hdd)
    spatial_cell_df, global_df =  my_load_df(path, offset=0, custom_ending = "_combined")
    spatial_cell_df["IL-2_surf_c"] *= 1e3
    try:
        spatial_cell_df["IL-2_gamma"]
    except KeyError:
        spatial_cell_df["IL-2_gamma"] = spatial_cell_df["misc_gamma"]
        spatial_cell_df["IL-2_Tsec_fraction"] = spatial_cell_df["fractions_Tsec"]
    print(spatial_cell_df["IL-2_gamma"].unique())

    if prefix == "pos":
        spatial_cell_df = spatial_cell_df.loc[(spatial_cell_df["IL-2_Tsec_fraction"] == 0.1)]
    else:
        spatial_cell_df = spatial_cell_df.loc[(spatial_cell_df["IL-2_Tsec_fraction"] == 0.05)]

    #%%
    labels = ["feedback", "ODE"]
    plot_x = []
    plot_ashmans = []
    plot_wangs = []
    acts = []
    for c,cell_df in enumerate([spatial_cell_df]):
        try:
            cell_df["pSTAT5"] = cell_df["IL-2_pSTAT5"]
        except:
            print(f"calculating own activation for c = {c}")
            cell_df["pSTAT5"] = cell_df["IL-2_surf_c"] ** 3 / (
                    (EC50_calculation(E_max=125e-12, E_min=0, k=860, N=1.5, R=cell_df["IL-2_R"]) * 1e12) ** 3 + cell_df[
                            "IL-2_surf_c"] ** 3).values
        cell_df = cell_df.loc[cell_df["time"] == cell_df["time"].max()]
        cell_df = cell_df.loc[cell_df["type_name"] == "Th"]

        ashmans = np.zeros((len(cell_df["IL-2_gamma"].unique()),
                        len(cell_df["replicat_index"].unique())))  # amount_of_components / amount_of_Tsecs
        wangs = np.zeros((len(cell_df["IL-2_gamma"].unique()),
                        len(cell_df["replicat_index"].unique())))  # amount_of_components / amount_of_Tsecs
        ashmans[:] = np.nan
        wangs[:] = np.nan

        for g, gamma in enumerate(np.sort(cell_df["IL-2_gamma"].unique())):
            gamma_df = cell_df.loc[(cell_df["IL-2_gamma"] == gamma)]
            print("gamma:", gamma)
            for r, rep in enumerate(np.sort(gamma_df["replicat_index"].unique())):
                rep_df = gamma_df.loc[(gamma_df["replicat_index"] == rep)]
                test = rep_df["pSTAT5"].values
                ashmans[g][r] = ashmans_D(rep_df["pSTAT5"].values)
                wangs[g][r] = wangs_BI(rep_df["pSTAT5"].values)

        plot_x.append(np.sort(cell_df["IL-2_gamma"].unique()))
        plot_ashmans.append(ashmans)
        plot_wangs.append(wangs)

    #%%
    print("plotting")
    rc_ticks["figure.figsize"] = (1.67475 * 1.27, 0.43)
    sns.set_theme(context="talk", style="ticks", rc=rc_ticks)
    fig, ax = plt.subplots()
    ylabels = ["Ashmans D", "Wangs index"]
    for v, value in enumerate([plot_ashmans]):
        for e, entry in enumerate(value): #cell_dfs
            alphas = np.logspace(-1, 0, len(entry))

            x = plot_x[e]
            y = np.mean(entry, axis=1)
            std = np.std(entry, axis=1)/np.sqrt(len(entry))
            if prefix == "neg":
                colour = "Blue"
                x = 1 / x
                x = x[::-1]
                y = y[::-1]
            elif prefix == "pos":
                colour = "Red"
            plt.fill_between(x,
                         np.clip(y - std, 0, None),
                         np.clip(y + std, 0, None), color=colour if e != 1 else "black",
                             alpha=0.3, linewidth=0.0)
            plt.plot(x, y, label=labels[e], color=colour if e != 1 else "black", alpha=1, marker=".")
        plt.axhline(2, color="black", linestyle="--") if v == 0 else None
        plt.ylabel(ylabels[v])
        plt.ylabel("bimodality", labelpad=0.3)
        plt.xlabel("feedback fold change")
        plt.yscale("log")
        plt.xscale("log")
        plt.xlim(1,100)
        plt.ylim((1, 30))

        if save_plot == True:
            fig.savefig(saving_string + f"{prefix}_{ylabels[v]}_bimodality_over_fb.pdf", bbox_inches='tight', transparent=True)
        plt.tight_layout()
        plt.show()