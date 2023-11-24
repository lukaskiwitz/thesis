import getpass
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import seaborn as sns
import os
from thesis.scripts.paper_models.utilities.plotting_rc import rc_ticks
from thesis.scripts.paper_models.utilities.plot_helper import my_load_df, EC50_calculation

save_plot = True

saving_string = r"/home/brunner/Documents/Current work/2023_11_03/pSTAT+_over_dist_to_tsec.pdf"

hdd = "extra2" if os.path.exists("/extra2") else "extra"
user = getpass.getuser()

base_path = "/extra2/brunner/paper_models/kinetics/Figure_1C/dataframes_act_over_Tsec_large_steady_state/"
big_c_df, global_df =  my_load_df(base_path, offset=0, custom_ending = "_combined")

try:
    big_c_df["IL-2_Tsec_fraction"].unique()
except KeyError:
    big_c_df["IL-2_Tsec_fraction"] = big_c_df["fractions_Tsec"]
try:
    big_c_df["IL-2_gamma"].unique()
except KeyError:
    big_c_df["IL-2_gamma"] = big_c_df["misc_gamma"]
if len(big_c_df["IL-2_gamma"].unique()) > 1:
    big_c_df = big_c_df.loc[big_c_df["IL-2_gamma"] == 12.]

big_c_df["IL-2_surf_c"] *= 1e3
big_c_df["pSTAT5"] = big_c_df["IL-2_surf_c"] ** 3 / (
        (EC50_calculation(E_max=125e-12, E_min=0, k=860, N=1.5, R=big_c_df["IL-2_R"]) * 1e12) ** 3 +
        big_c_df["IL-2_surf_c"] ** 3).values

#%%
df_dicts = []
for f, frac in enumerate([0.05, 0.1, 0.15, 0.2, 0.3, 0.4]):
    print(frac)
    frac_df = big_c_df.loc[np.abs(big_c_df["IL-2_Tsec_fraction"] - frac) < 5e-3]
    for r, replicat_index in enumerate(frac_df["replicat_index"].unique()):
        rep_df = frac_df.loc[(frac_df["replicat_index"] == replicat_index)]
        rep_df["IL-2_surf_c"] *= 1e3

        from scipy.spatial import distance_matrix
        coordinates = rep_df[["x", "y", "z"]].values
        dm = distance_matrix(coordinates, coordinates)
        Tsec_ids = rep_df.loc[rep_df["type_name"] == "Tsec", "id"].values
        surf_c_over_dist_tsec = np.zeros((len(rep_df["id"]), 2))
        for index, id in enumerate(rep_df["id"]):
            distance = np.min(dm[id][Tsec_ids][np.argwhere(dm[id][Tsec_ids] > 0)]) if id not in Tsec_ids else np.nan
            act = 1 if rep_df.loc[rep_df["id"] == id, "pSTAT5"].values[0] >= 0.5 else 0
            df_dicts.append({"IL-2_Tsec_fraction": frac, "replicat_index": replicat_index, "id": id, "distance": distance,
                             "activated": act})

#%%
rc_ticks['figure.figsize'] = (1.67475 * 0.56 * 1.12, 1.386 * 0.95 / 2 * 1.1)
color = "#aa0000ff"
sns.set_theme(context="talk", style="ticks", rc=rc_ticks)
fig, ax = plt.subplots()
plot_df = pd.DataFrame(df_dicts)
plot_df["distance"] /= 20 #c-c-d
alphas = np.logspace(-1, 0, len(plot_df["IL-2_Tsec_fraction"].unique()))
distances = np.sort(plot_df.distance.unique())

for f, frac in enumerate(plot_df["IL-2_Tsec_fraction"].unique()):
    frac_df = plot_df.loc[(plot_df["IL-2_Tsec_fraction"] == frac)]
    rep_y = []
    for r, rep in enumerate(frac_df.replicat_index.unique()):
        rep_df = frac_df.loc[frac_df.replicat_index == rep]
        x = []
        y = []
        for e, entry in enumerate(distances):
            x.append(entry)
            bin_df = rep_df.loc[(rep_df.distance == entry)]
            y.append(len(bin_df[bin_df["activated"] == 1]) / len(bin_df) * 100 if len(bin_df) > 0 else np.nan)
        x = np.array(x)
        rep_y.append(y)
    mean = np.nanmean(rep_y, axis=0)
    std = np.nanstd(rep_y, axis=0)/np.sqrt(len(plot_df.replicat_index.unique()))
    plt.plot(x, mean, color=color, alpha=alphas[f], marker="o", markersize=3, markeredgewidth=0, linewidth=0.5)
    plt.fill_between(x, np.clip(mean - std, 0, 100), mean + std, color = color, alpha=alphas[f]/2, linewidth=0)

plt.xlabel(r"distance to nearest" "\n" r"T$_{\rm sec}$ (c-c-d)")
plt.ylabel(r"%pSTAT$^+$")
plt.xlim((1, 3))
plt.ylim((-0.5, 100))

fig.savefig(saving_string, bbox_inches='tight', transparent=True)
plt.tight_layout()
plt.show()