import getpass
import matplotlib.pyplot as plt
import numpy as np
import seaborn as sns
import os
from scipy.optimize import curve_fit
from scipy.stats import spearmanr
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
        model_name = "feedback_scan" #c4
        name = "dataframes_negative_for_FigS5C_2_steady_state"  # menten
        my_hue = "scan_name_scan_name"
        Tsec_fraction = 1
    elif prefix == "pos":
        model_name = "feedback_scan_4" #c4
        name = "dataframes_positive_10_steady_state"
        Tsec_fraction = 1
        my_hue = "scan_name_scan_name"

    path = "/{extra}/{u}/paper_models/kinetics/{mn}/{n}/".format(u=user, mn=model_name, n=name, extra=hdd)
    spatial_cell_df, _ = my_load_df(path, offset=0, custom_ending = "_combined_2" if prefix == "neg" else "_combined")
    spatial_cell_df["IL-2_Tsec_fraction"] = spatial_cell_df["fractions_Tsec"]

    epsilon = 20
    results = []
    for c, cell_df in enumerate([spatial_cell_df]):

        cell_df["IL-2_surf_c"] *= 1e3
        try:
            cell_df["pSTAT5"] = cell_df["IL-2_pSTAT5"]
        except:
            print("pSTAT5 loading failed, recalculating")
            cell_df["pSTAT5"] = cell_df["IL-2_surf_c"] ** 3 / (
                    (EC50_calculation(E_max=125e-12, E_min=0, k=860, N=1.5, R=cell_df["IL-2_R"]) * 1e12) ** 3 + cell_df[
                "IL-2_surf_c"] ** 3).values

        if my_hue != "time":
            cell_df = cell_df.loc[(cell_df["time_index"] == cell_df["time_index"].max())]
        if Tsec_fraction != None:
            cell_df = cell_df.loc[(cell_df["IL-2_Tsec_fraction"] == np.sort(cell_df["IL-2_Tsec_fraction"].unique())[Tsec_fraction])]

        niche_concentrations = [[[] for x in cell_df[my_hue].unique()] for y in cell_df["replicat_index"].unique()]
        niche_pSTAT5 = [[[] for x in cell_df[my_hue].unique()] for y in cell_df["replicat_index"].unique()]

        niche_score = np.zeros((len(cell_df["replicat_index"].unique()),
                                len(cell_df[my_hue].unique())))
        niche_effect = np.zeros((len(cell_df["replicat_index"].unique()),
                                len(cell_df[my_hue].unique())))

        gammas = np.zeros((len(cell_df["replicat_index"].unique()),
                                len(cell_df[my_hue].unique())))

        for r, rep in enumerate(cell_df["replicat_index"].unique()):
            print(rep)
            rep_df = cell_df.loc[cell_df["replicat_index"] == rep]
            for idx, value in enumerate(np.sort(rep_df[my_hue].unique())[1:]):
                frac_df = rep_df.loc[(rep_df[my_hue] == value)]
                Tsec_ids = frac_df.loc[frac_df["type_name"] == "Tsec", "id"].values

                DBres, sliced_df, DBcoords = myDBscan(frac_df, epsilon, with_Tsecs=True)
                sliced_df["cluster"] = DBres

                niche_score[r][idx] = np.sum(~np.isnan(np.where(np.unique(DBres, return_counts=True)[1] > 1))) / len(
                                        frac_df.loc[frac_df["type_name"] == "Tsec"]) if len(
                                        frac_df.loc[frac_df["type_name"] == "Tsec"]) != 0 else 0
                try:
                    gammas[r][idx] = frac_df["misc_gamma"].unique()[frac_df["misc_gamma"].unique() != 1][0]
                except KeyError:
                    gammas[r][idx] = frac_df["IL-2_gamma"].unique()[frac_df["IL-2_gamma"].unique() != 1][0]

                effects = []
                for u in np.unique(DBres):
                    if len(np.where(DBres == u)[0]) >= 1:
                        ids = sliced_df.loc[sliced_df["cluster"] == u, "id"].values

                        in_cluster_values = sliced_df.loc[sliced_df["cluster"] == u, "pSTAT5"].values
                        outside_cluster_values = frac_df.loc[
                            (frac_df["pSTAT5"] < 0.5) & (frac_df["type_name"] != "Tsec"), "pSTAT5"].values

                        if len(in_cluster_values[~np.isnan(in_cluster_values)]) > 0 and len(outside_cluster_values[~np.isnan(outside_cluster_values)]) > 0:
                            effects.append(np.nanmean(in_cluster_values) / np.nanmean(outside_cluster_values))
                niche_effect[r][idx] = np.nanmean(effects) if len(effects) > 0 else None
        results.append([niche_score, niche_effect])

    #%%
    print("plotting")
    factor = 1 if prefix == "pos" else 1.0
    rc_ticks['figure.figsize'] = [1.67475, 1.386]
    sns.set_theme(context="talk", style="ticks", rc=rc_ticks)
    fig, ax = plt.subplots()
    for e, entry in enumerate([spatial_cell_df]):
        niche_score = results[e][0]
        niche_effect = results[e][1]

        if prefix == "pos":
            cmap_name = "Red"
            reverse = False
        elif prefix == "neg":
            cmap_name = "Blue"
            reverse = True
        no_fb_colors = sns.color_palette("Greys", len(np.nanmean(niche_score, axis=0)))
        import matplotlib as mpl
        for i in range(len(niche_score)):
            scatterplot = plt.scatter(niche_score[i][np.argsort(gammas[0])] + np.random.normal(0, 0.002, len(niche_score[i])),
                        niche_effect[i][np.argsort(gammas[0])], c=cmap_name, s = 7, linewidth=0, norm=mpl.colors.LogNorm(),
                        alpha = np.logspace(0,-0.9, len(gammas[0])) if reverse == True else np.logspace(-0.9,0, len(gammas[0])))
        def func(x,a,b):
            return a + b*x
        no_nans = ~np.isnan(niche_effect.flatten())
        no_zeros = np.where(niche_effect.flatten()[no_nans] != 0)[0]
        popt, pcov = curve_fit(func, niche_score.flatten()[no_nans][no_zeros], niche_effect.flatten()[no_nans][no_zeros])
        if prefix == "neg":
            if Tsec_fraction == 0:
                base = [0.3, func(0.3, *popt)]
                dx, dy = [0.6, func(0.6, *popt)]
            elif Tsec_fraction == 1:
                base = [0.03, func(0.03, *popt)]
                dx, dy = [0.25, func(0.25, *popt)]
            elif Tsec_fraction == 2:
                popt, pcov = curve_fit(func, niche_score.flatten()[no_nans],
                                       niche_effect.flatten()[no_nans])
                base = [0.001, func(0.001, *popt)]
                dx, dy = [0.024, func(0.024, *popt)]

            xoffset = 0
            yoffset = 0
            plt.annotate("", xytext=(dx + xoffset, dy + yoffset), xy=[base[0] + xoffset, base[1] + yoffset],
                         arrowprops={"linewidth": 1.5, "width": 1, "headwidth": 5, "headlength": 5}) #white line around arrow, controlled via linewidth
            plt.annotate("", xytext=(dx + xoffset, dy + yoffset), xy=[base[0] + xoffset, base[1] + yoffset],
                         arrowprops={"linewidth": 0, "width": 1, "headwidth": 5, "headlength": 5, "color": "black"})
        else:
            if Tsec_fraction == 0:
                base = [0.035, func(0.035, *popt)]
                dx, dy = [0.235, func(0.235, *popt)]
            elif Tsec_fraction == 1:
                base = [0.1, func(0.1, *popt)]
                dx, dy = [0.35, func(0.35, *popt)]
            elif Tsec_fraction == 2:
                base = [0.25, func(0.25, *popt)]
                dx, dy = [0.3, func(0.3, *popt)]
            xoffset = 0
            yoffset = 0
            plt.annotate("", xy=(dx + xoffset, dy + yoffset), xytext=[base[0] + xoffset, base[1] + yoffset],
                         arrowprops={"linewidth": 1.5, "width": 1, "headwidth": 5, "headlength": 5}) #white line around arrow, controlled via linewidth
            plt.annotate("", xy=(dx + xoffset, dy + yoffset), xytext=[base[0] + xoffset, base[1] + yoffset],
                         arrowprops={"linewidth": 0, "width": 1, "headwidth": 5, "headlength": 5, "color": "black"})
    plt.xlabel("niche score")
    plt.ylabel("niche effect")
    plt.yscale("linear")
    plt.xscale("linear")
    if prefix == "neg":
        if Tsec_fraction == 0:
            plt.xlim((0.15, 0.75))
            plt.ylim((10, 40))
            plt.yticks([10, 25, 40])
        elif Tsec_fraction == 1:
            plt.xlim((0., 0.3))
            plt.ylim((6, 16))
            plt.yticks([6,11, 16])
        elif Tsec_fraction == 2:
            plt.xlim((0, 0.03))
            plt.ylim((0, 8))
            plt.yticks([0, 4, 8])
    else:
        if Tsec_fraction == 0:
            plt.xlim((0.02, 0.25))
            plt.ylim((250, 1750))
            plt.yticks([250, 1000, 1750])
            plt.yticks()
        elif Tsec_fraction == 1:
            plt.xlim((0.05, 0.4))
            plt.ylim((0, 500))
            plt.yticks([0,250, 500])
        elif Tsec_fraction == 2:
            plt.xlim((0.14, 0.4))
            plt.ylim((0, 340))
            plt.yticks([0, 170, 340])

    x = niche_score.flatten()
    y = niche_effect.flatten()
    res = spearmanr(x, y)
    rs = round(res[0], 3)
    plt.text(0.01, 0.9, r"r$_{\rm s}$ = " + str(rs), ha="left", fontsize=6, transform=plt.gca().transAxes,
             bbox={"pad": 0.5, "color": "white", "alpha": 0})

    if save_plot == True:
        fig.savefig(saving_string + f"_{prefix}" + "niche_effect_over_niche_score.pdf", bbox_inches='tight', transparent=True)
    plt.tight_layout()
    plt.show()