import getpass
import matplotlib.pyplot as plt
import numpy as np
import seaborn as sns
import os
from thesis.scripts.paper_models.utilities.plotting_rc import rc_ticks
from thesis.scripts.paper_models.utilities.plot_helper import my_load_df, myDBscan, EC50_calculation

save_plot = True
hdd = "/extra2" if os.path.exists("/extra2") else "/extra"
user = getpass.getuser()
saving_string =r"/home/brunner/Documents/Current work/2023_11_17/"
if not os.path.exists(saving_string):
    os.mkdir(saving_string)

setups = [
    "pos",
    "neg",
]

for prefix in setups:
    if prefix == "neg":
        model_name = "feedback_scan_4"
        name = "dataframes_negative_for_Fig3C_act_plot_2_steady_state"
        Tsec_fraction = 0
        color = "blue"
    elif prefix == "pos":
        model_name = "feedback_scan_4/"
        name = "dataframes_positive_for_Fig3C_act_plot_steady_state"
        Tsec_fraction = 0
        color = "red"

    path = "/{extra}/{u}/paper_models/kinetics/{mn}/{n}/".format(u=user, mn=model_name, n=name, extra=hdd)
    cell_df, _ = my_load_df(path, offset=0, custom_ending = "_combined_2" if prefix == "neg" else "_combined")
    cell_df["IL-2_surf_c"] *= 1e3
    try:
        cell_df["IL-2_gamma"]
    except KeyError:
        cell_df["IL-2_gamma"] = cell_df["misc_gamma"]
        cell_df["IL-2_Tsec_fraction"] = cell_df["fractions_Tsec"]
    cell_df = cell_df.loc[
        (cell_df["IL-2_Tsec_fraction"] == np.sort(cell_df["IL-2_Tsec_fraction"].unique())[Tsec_fraction])]


    cell_df = cell_df.loc[(cell_df["time_index"] == cell_df["time_index"].max())]
    try:
        cell_df["pSTAT5"] = cell_df["IL-2_pSTAT5"]
    except:
        print(f"calculating own activation for c = {c}")
        cell_df["pSTAT5"] = cell_df["IL-2_surf_c"] ** 3 / (
                (EC50_calculation(E_max=125e-12, E_min=0, k=860, N=1.5, R=cell_df["IL-2_R"]) * 1e12) ** 3 + cell_df[
                        "IL-2_surf_c"] ** 3).values

    signaling_range = np.empty((len(cell_df["IL-2_gamma"].unique()), len(cell_df["replicat_index"].unique())))
    signaling_range[:] = np.nan
    for g, gamma in enumerate(np.sort(cell_df["IL-2_gamma"].unique())):
        gamma_df = cell_df.loc[(cell_df["IL-2_gamma"] == gamma)]
        print("gamma:", gamma)
        for r, rep in enumerate(np.sort(gamma_df["replicat_index"].unique())):
            rep_df = gamma_df.loc[(gamma_df["replicat_index"] == rep)]

            DBres, sliced_df, DBcoords = myDBscan(rep_df, 20, with_Tsecs=True)
            coordinates = sliced_df[["x", "y", "z"]].values
            from scipy.spatial import distance_matrix

            dm = distance_matrix(coordinates, coordinates)
            niche_dists = []
            for niche_index in np.where(np.unique(DBres, return_counts=True)[1] > 1)[0]:
                niche_df = sliced_df.iloc[np.where(DBres == niche_index)]
                Tsec_ids = niche_df.loc[niche_df["type_name"] == "Tsec"].index.values
                if len(Tsec_ids) > 0:
                    for entry in niche_df.index:
                        if entry not in Tsec_ids:
                            niche_dists.append(np.min(dm[entry][Tsec_ids]))
            if len(niche_dists) < 0:
                assert False
            signaling_range[g][r] = np.quantile(niche_dists, 0.95) if len(niche_dists) > 0 else np.nan

    plot_x = np.sort(cell_df["IL-2_gamma"].unique())
    signaling_range /= 20 #standard c-c-d

    #%%
    print("plotting")
    rc_ticks['figure.figsize'] = [1.67475 * 1.25, 0.65]
    sns.set_theme(context="talk", style="ticks", rc=rc_ticks)
    fig, ax = plt.subplots()

    y = np.nanmean(signaling_range, axis=1)
    std = np.nanstd(signaling_range, axis=1)/np.sqrt(np.sum(~np.isnan(signaling_range[0])))
    if prefix == "neg":
        plot_x = 1 / plot_x
        plot_x = plot_x[::-1]
        y = y[::-1]
        std = std[::-1]
    plt.plot(plot_x, y, color=color, marker=".")
    plt.fill_between(np.unique(plot_x), np.clip(y - std, 0, None), np.clip(y + std, 0, None), color=color, alpha=0.3, linewidth=0.0)

    plt.ylabel(r"signaling" "\n" "range (c-c-d)")
    plt.xlabel("feedback fold change")
    plt.yscale("linear")
    plt.xscale("log")
    plt.ylim(0.95, 2.4)
    plt.yticks([1, 1.7, 2.4])
    plt.xlim(1,100)

    if save_plot == True:
        fig.savefig(saving_string + prefix +  f"_signaling_range_over_fb.pdf", bbox_inches='tight', transparent=True)
    plt.tight_layout()
    plt.show()
