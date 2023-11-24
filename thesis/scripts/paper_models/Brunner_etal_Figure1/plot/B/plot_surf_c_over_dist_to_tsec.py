import getpass
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import seaborn as sns
import os
from thesis.scripts.paper_models.utilities.plotting_rc import rc_ticks
from thesis.scripts.paper_models.utilities.plot_helper import my_load_df

save_plot = True

saving_string = r"/home/brunner/Documents/Current work/2023_11_03/surf_c_over_dist_to_tsec.pdf"

hdd = "extra2" if os.path.exists("/extra2") else "extra"
user = getpass.getuser()

base_path = "/extra2/brunner/paper_models/kinetics/Figure_1C/dataframes_act_over_Tsec_large_steady_state/" #michaelis
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
            distance = np.min(dm[id][Tsec_ids][np.argwhere(dm[id][Tsec_ids] > 0)]) if id not in Tsec_ids else 0
            conc = rep_df.loc[rep_df["id"] == id, "IL-2_surf_c"].values
            assert len(conc) == 1
            df_dicts.append({"IL-2_Tsec_fraction": frac, "replicat_index": replicat_index, "id": id, "distance": distance,
                             "IL-2_surf_c": conc[0]})

#%%
rc_ticks['figure.figsize'] = (1.67475 * 0.56 * 1.12, 1.386 * 0.95 / 2 * 1.1)
color = "#aa0000ff"
sns.set_theme(context="talk", style="ticks", rc=rc_ticks)
fig, ax = plt.subplots()
plot_df = pd.DataFrame(df_dicts)
plot_df["distance"] /= 20
alphas = np.logspace(-1, 0, len(plot_df["IL-2_Tsec_fraction"].unique()))

for f, frac in enumerate(plot_df["IL-2_Tsec_fraction"].unique()):
    if frac == 0.3:
        plot_df = plot_df.loc[(plot_df.distance < 2.4)] # filter the ones where there is only one occurence
    if frac == 0.4:
        plot_df = plot_df.loc[(plot_df.distance < 1.8)] # filter the ones where there is only one occurence
    sns.lineplot(data=plot_df.loc[(plot_df["IL-2_Tsec_fraction"] == frac)], x = "distance", y = "IL-2_surf_c",
                 color=color, alpha=alphas[f], err_kws={"linewidth": 0}, marker="o", markersize=3, markeredgewidth=0,
                 errorbar=("se"), linewidth=0.5)

plt.xlabel(r"")
plt.ylabel(r"surface conc." "\n" "avg. (pM)")
plt.xlim((1, 3))
plt.ylim((-0.5, 25))
plt.yticks([0, 10, 20])
plt.xticks()
xticks = plt.xticks()[0]
plt.xticks(xticks,["" for x in range(len(xticks))])
fig.savefig(saving_string, bbox_inches='tight', transparent=True)
plt.tight_layout()
plt.show()