import getpass
import matplotlib.pyplot as plt
import numpy as np
import seaborn as sns
import os
import pandas as pd
from scipy.spatial import distance_matrix
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
        model_name = "feedback_tsec_scan_for_ns_ne_plot" #michaelis
        name = "dataframes_negative_steady_state"
        my_hue = "scan_name_scan_name"
    elif prefix == "pos":
        model_name = "feedback_tsec_scan_for_ns_ne_plot" #michaelis
        name = "dataframes_positive_steady_state"
        my_hue = "scan_name_scan_name"

    path = "/{extra}/{u}/paper_models/kinetics/{mn}/{n}/".format(u=user, mn=model_name, n=name, extra=hdd)
    spatial_cell_df, _ = my_load_df(path, offset=0, custom_ending = "_combined")
    spatial_cell_df["IL-2_Tsec_fraction"] = spatial_cell_df["fractions_Tsec"]

    #%%
    epsilon = 20
    for c, cell_df in enumerate([spatial_cell_df]):
        try:
            cell_df["pSTAT5"] = cell_df["IL-2_pSTAT5"]
        except:
            print("pSTAT5 loading failed, recalculating")
            cell_df["IL-2_surf_c"] *= 1e3
            cell_df["pSTAT5"] = cell_df["IL-2_surf_c"] ** 3 / (
                    (EC50_calculation(E_max=125e-12, E_min=0, k=860, N=1.5, R=cell_df["IL-2_R"]) * 1e12) ** 3 + cell_df[
                "IL-2_surf_c"] ** 3).values

        if my_hue != "time":
            cell_df = cell_df.loc[(cell_df["time_index"] == cell_df["time_index"].max())]

        df_dict_list = []
        for Tsec_fraction in [0,1,2,3,4]:
            if Tsec_fraction != None:
                frac_df = cell_df.loc[(cell_df["IL-2_Tsec_fraction"] == np.sort(cell_df["IL-2_Tsec_fraction"].unique())[Tsec_fraction])]

            niche_concentrations = [[[] for x in frac_df[my_hue].unique()] for y in frac_df["replicat_index"].unique()]
            niche_pSTAT5 = [[[] for x in frac_df[my_hue].unique()] for y in frac_df["replicat_index"].unique()]

            niche_score = np.zeros((len(frac_df["replicat_index"].unique()),
                                    len(frac_df[my_hue].unique())))  # amount_of_components / amount_of_Tsecs
            niche_effect = np.zeros((len(frac_df["replicat_index"].unique()),
                                    len(frac_df[my_hue].unique())))  # amount_of_components / amount_of_Tsecs

            gammas = np.zeros((len(frac_df["replicat_index"].unique()),
                                    len(frac_df[my_hue].unique())))

            for r, rep in enumerate(frac_df["replicat_index"].unique()):
                rep_df = frac_df.loc[frac_df["replicat_index"] == rep]
                for idx, value in enumerate(np.sort(rep_df[my_hue].unique())):
                    gamma_df = rep_df.loc[(rep_df[my_hue] == value)]
                    coordinates = gamma_df[["x", "y", "z"]].values
                    from scipy.spatial import distance_matrix

                    dist_m = distance_matrix(coordinates, coordinates)
                    Tsec_ids = gamma_df.loc[gamma_df["type_name"] == "Tsec", "id"].values

                    DBres, sliced_df, DBcoords = myDBscan(gamma_df, epsilon, with_Tsecs=True)
                    sliced_df["cluster"] = DBres

                    mean_score = np.sum(~np.isnan(np.where(np.unique(DBres, return_counts=True)[1] > 1))) / len(
                                            gamma_df.loc[gamma_df["type_name"] == "Tsec"]) if len(
                                            gamma_df.loc[gamma_df["type_name"] == "Tsec"]) != 0 else 0
                    niche_score[r][idx] = mean_score
                    try:
                        gammas[r][idx] = gamma_df["misc_gamma"].unique()[gamma_df["misc_gamma"].unique() != 1][0]
                    except KeyError:
                        gammas[r][idx] = gamma_df["IL-2_gamma"].unique()[gamma_df["IL-2_gamma"].unique() != 1][0]

                    effects = []
                    for u in np.unique(DBres):
                        if len(np.where(DBres == u)[0]) >= 1:
                            ids = sliced_df.loc[sliced_df["cluster"] == u, "id"].values

                            in_cluster_values = sliced_df.loc[sliced_df["cluster"] == u, "pSTAT5"].values
                            outside_cluster_values = gamma_df.loc[
                                (gamma_df["pSTAT5"] < 0.5) & (gamma_df["type_name"] != "Tsec"), "pSTAT5"].values

                            if not all(np.isnan(x) for x in np.unique(in_cluster_values)) and len(
                                    outside_cluster_values) > 0:
                                effects.append(np.nanmean(in_cluster_values) / np.nanmean(outside_cluster_values))
                    mean_effect = np.nanmean(effects) if len(effects) > 0 else None
                    niche_effect[r][idx] = mean_effect

                    coordinates = sliced_df[["x", "y", "z"]].values
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
                    signaling_range = np.quantile(niche_dists, 0.95) if len(niche_dists) > 0 else np.nan

                    df_dict_list.append({"IL-2_Tsec_fraction" : np.sort(cell_df["IL-2_Tsec_fraction"].unique())[Tsec_fraction],
                                         "replicat_index": rep, "hue_value": value, "niche_effect": mean_effect,
                                         "niche_score": mean_score, "gamma": gammas[r][idx],
                                         "signaling_range": signaling_range / rep_df.geometry_distance.unique()[0]})
    df = pd.DataFrame(df_dict_list)

    #%%
    palette = ["mistyrose", "red"] if prefix == "pos" else ["lavender", "blue"]
    gamma_df = df.loc[(df["gamma"] == df["gamma"].min()) | (df["gamma"] == df["gamma"].max())]
    gamma_df["IL-2_Tsec_fraction"] *= 100

    rc_ticks['figure.figsize'] = (0.93, 1.25)
    sns.set_theme(context="talk", style="ticks", rc=rc_ticks)
    fig, ax = plt.subplots()

    gammas = np.sort(gamma_df.gamma.unique()) if prefix == "pos" else np.sort(gamma_df.gamma.unique())[::-1]

    g = sns.barplot(data=gamma_df, x = "IL-2_Tsec_fraction", y="niche_score", hue="gamma", capsize=.08, errorbar="sd",
                ax=ax, errwidth=0.5, linewidth=0.5, edgecolor="black", dodge=True, hue_order=gammas, palette=palette, width=0.6, alpha=1)
    g.legend_.remove()
    ax.set(xlabel="secreting cells (%)", ylabel="niche score")
    plt.xticks(rotation=-45, ha="center")
    if prefix == "pos":
        plt.ylim(0, 0.4)
        plt.yticks([0, 0.2, 0.4])
    else:
        plt.ylim(0, 0.8)
        plt.yticks([0, 0.4, 0.8])
    plt.locator_params(axis='y', nbins=3)
    fig.savefig(saving_string + f"/{prefix}_niche_score_over_secreting_cells" + ".pdf", bbox_inches='tight',
                transparent=True)
    plt.tight_layout()
    plt.show()

    #%%
    sns.set_theme(context="talk", style="ticks", rc=rc_ticks)
    fig, ax = plt.subplots()

    g = sns.barplot(data=gamma_df, x = "IL-2_Tsec_fraction", y="niche_effect", hue="gamma", capsize=.08, errorbar="sd",
                ax=ax, errwidth=0.5, linewidth=0.5, edgecolor="black", dodge=True, hue_order=gammas, palette=palette, width=0.6, alpha=1)
    g.legend_.remove()
    ax.set(xlabel="secreting cells (%)", ylabel="niche effect")
    plt.xticks(rotation=-45, ha="center")
    plt.locator_params(axis='y', nbins=3)
    if prefix == "pos":
        plt.ylim(0, 1600)
        plt.yticks([0, 800, 1600])
    else:
        plt.ylim(0, 40)
        plt.yticks([0, 20, 40])
    fig.savefig(saving_string + f"/{prefix}_niche_effect_over_secreting_cells" + ".pdf", bbox_inches='tight',
                transparent=True)
    plt.tight_layout()
    plt.show()

    #%%
    rc_ticks['figure.figsize'] = (0.93, 1.25)
    sns.set_theme(context="talk", style="ticks", rc=rc_ticks)
    fig, ax = plt.subplots()

    gammas = np.sort(gamma_df.gamma.unique()) if prefix == "pos" else np.sort(gamma_df.gamma.unique())[::-1]

    g = sns.barplot(data=gamma_df, x = "IL-2_Tsec_fraction", y="signaling_range", hue="gamma", capsize=.08, errorbar="sd",
                ax=ax, errwidth=0.5, linewidth=0.5, edgecolor="black", dodge=True, hue_order=gammas, palette=palette, width=0.6, alpha=1) #alpha=0.8 for grey
    g.legend_.remove()
    ax.set(xlabel="secreting cells (%)", ylabel=r"signaling range (c-c-d)")
    plt.xticks(rotation=-45, ha="center")
    plt.locator_params(axis='y', nbins=3)

    if prefix == "pos":
        plt.ylim(0, 2)
        plt.yticks([0, 1, 2])
    else:
        plt.ylim(0, 3)
        plt.yticks([0, 1.5, 3])

    fig.savefig(saving_string + f"/{prefix}_signaling_range_over_secreting_cells" + ".pdf", bbox_inches='tight',
                transparent=True)
    plt.tight_layout()
    plt.show()
