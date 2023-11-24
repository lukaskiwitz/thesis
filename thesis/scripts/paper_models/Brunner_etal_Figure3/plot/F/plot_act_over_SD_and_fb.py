import matplotlib.pyplot as plt
import numpy as np
import seaborn as sns
import os
import getpass
from thesis.scripts.paper_models.utilities.plotting_rc import rc_ticks
import pandas as pd

from thesis.scripts.paper_models.utilities.plot_helper import my_load_df, my_interpolation

hdd = "extra2" if os.path.exists("/extra2") else "extra"
user = getpass.getuser()
# prefix = "pos"
prefix = "neg"
# prefix = "ODE_pos"
# if prefix == "pos":
#     model_name = "Tsec_scan_6"
#     name = "dataframes_positive_2_steady_state" #michaelis
#     xlim = (1, 50)
# elif prefix == "neg":
#     model_name = "Tsec_scan_6"
#     # name = "gamma_0.5"  # compute4
#     # model_name = "feedback_scan_4"
#     name = "dataframes_negative_steady_state" # menten
#     xlim = (1, 20)
if prefix == "pos":
    model_name = "feedback_scan_4"
    name = "dataframes_positive_for_Fig3D_pos_steady_state" #michaelis
elif prefix == "neg":
    model_name = "feedback_scan_4"
    name = "dataframes_negative_for_Fig3D_neg_steady_state"  # compute4

saving_string = f"/home/brunner/Documents/Current work/2023_11_03/{prefix}_correlate_act_with_SD_over_fb.pdf"

path = "/{extra}/{u}/paper_models/kinetics/{mn}/{n}/".format(u=user, n=name, mn=model_name, extra=hdd)
spatial_cell_df, global_df =  my_load_df(path, offset=0, custom_ending = "_combined")
spatial_cell_df["IL-2_surf_c"] *= 1e3
try:
    spatial_cell_df["IL-2_gamma"]
except KeyError:
    spatial_cell_df["IL-2_gamma"] = spatial_cell_df["misc_gamma"]
    spatial_cell_df["IL-2_Tsec_fraction"] = spatial_cell_df["fractions_Tsec"]

ss_df = spatial_cell_df.loc[(spatial_cell_df["time_index"] == spatial_cell_df["time_index"].max())]
########################################################################################################################
print("plotting")
########################################################################################################################
#%%
def EC50_calculation(E_max, E_min, k, N, R):
    return (E_max * k ** N + E_min * R ** N) / (k ** N + R ** N)

dict_list = []
for c,cell_df in enumerate([ss_df]):
    try:
        cell_df["pSTAT5"] = cell_df["IL-2_pSTAT5"]
    except:
        print(f"calculating own activation for c = {c}")
        cell_df["pSTAT5"] = cell_df["IL-2_surf_c"] ** 3 / (
                (EC50_calculation(E_max=125e-12, E_min=0, k=860, N=1.5, R=cell_df["IL-2_R"]) * 1e12) ** 3 + cell_df[
                        "IL-2_surf_c"] ** 3).values
    # cell_df = cell_df.loc[cell_df["time"] == cell_df["time"].max()]
    cell_df = cell_df.loc[cell_df["type_name"] == "Th"]

    for g, gamma in enumerate(np.sort(cell_df["IL-2_gamma"].unique())[:-9 if prefix == "pos" else -1]):
        gamma_df = cell_df.loc[(cell_df["IL-2_gamma"] == gamma)]
        print("gamma:", gamma)
        for f, frac in enumerate(np.sort(gamma_df["IL-2_Tsec_fraction"].unique())):
            frac_df = gamma_df.loc[(gamma_df["IL-2_Tsec_fraction"] == frac)]
            for r, rep in enumerate(np.sort(frac_df["replicat_index"].unique())):
                rep_df = frac_df.loc[(frac_df["replicat_index"] == rep)]
                act = len(rep_df.loc[(rep_df["pSTAT5"] >= 0.5)]) / len(rep_df)
                SD = rep_df["IL-2_surf_c"].std()
                dict_list.append({"gamma": gamma, "frac": frac, "activation": act * 100, "SD": SD, "replicat_index": rep})
#%%
data_df = pd.DataFrame(dict_list)
for frac in data_df["frac"].unique():
    subset_df = data_df.loc[data_df["frac"] == frac]
    # rc_ticks['figure.figsize'] = (1.67475 * 1.05, 1.386 * 1.06 / 2)
    rc_ticks['figure.figsize'] = (1.67475 * 0.84, 1.386 * 0.5)
    # rc_ticks['figure.figsize'] = (5,5)
    sns.set_theme(context="talk", style="ticks", rc=rc_ticks)
    fig, ax = plt.subplots()


    color = "red" if prefix == "pos" else "blue"
    from matplotlib.colors import to_rgba
    from copy import deepcopy
    rgba_color = list(to_rgba(color))
    palette = []
    if prefix == "pos":
        alphas = np.logspace(-0.9,-0.4, len(subset_df["gamma"].unique()))
    else:
        alphas = np.logspace(-1,0, len(subset_df["gamma"].unique()))
    for e,entry in enumerate(np.sort(subset_df["gamma"].unique())):
        alphaed_color = deepcopy(rgba_color)
        alphaed_color[-1] = alphas[e]
        palette.append(alphaed_color)
    if prefix == "neg":
        palette = palette[::-1]
    sns.scatterplot(data=subset_df, x = "SD", y = "activation", hue = "gamma", palette=palette, legend=False, s=5, linewidth=0.075, zorder=-10)
    # plt.title(frac)
    # plt.ylabel("%pSTAT$^+$")
    plt.xlabel(r"spatial s.d. (pM)")
    plt.xscale("linear")
    plt.yscale("linear")
    plt.ylim(-4, 104)
    plt.xscale("log")
    if prefix == "pos":
        plt.xlim(3, 10)
        plt.xlabel("")
    else:
        plt.xlim(3, 150)
    plt.xlim(2, 150)
    plt.ylabel("%pSTAT$^+$")
    # yticks = plt.yticks()[0]
    # yticks = yticks[(yticks >= -4) & (yticks <= 104)]
    # plt.xticks([0,5,10,15])
    # plt.yticks([0,50,100])
    # plt.yticks(yticks,["" for x in range(len(yticks))])
    x = subset_df["SD"].values
    y = subset_df.activation.values
    res = spearmanr(x, y)
    rs = round(res[0], 3)
    plt.text(0.01, 0.76, r"r$_{\rm s}$ = " + str(rs), ha="left",fontsize=6,transform = plt.gca().transAxes, bbox={"pad": 0.5, "color": "white", "alpha": 0})
    ax.set_rasterization_zorder(0)
    fig.savefig(saving_string, bbox_inches='tight', transparent=True)
    plt.tight_layout()
    plt.show()

    from scipy.stats import pearsonr, spearmanr


    print(res)

