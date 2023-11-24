import matplotlib.pyplot as plt
from mpl_toolkits.axes_grid1.inset_locator import inset_axes
import numpy as np
import pandas as pd
import seaborn as sns
from scipy.spatial import distance_matrix
from thesis.scripts.paper_models.utilities.plotting_rc import rc_ticks
from thesis.scripts.paper_models.utilities.plot_helper import myDBscan

saving_path = f"/home/brunner/Documents/Current work/2023_11_17/"

path = "/extra2/brunner/paper_models/kinetics/il7_boundary_bridson/dataframes_Tsec_clustering_scan_steady_state/"
try:
    cell_df = pd.read_hdf(path + "cell_df" + ".h5", mode="r")
except:
    cell_df = pd.read_pickle(path + 'cell_df' + '.pkl')

treg_scan_names = cell_df.loc[cell_df.type_name == "Treg", "scan_name_scan_name"].unique()
spliced_df = cell_df.loc[~cell_df.scan_name_scan_name.isin(treg_scan_names)]

spliced_df["pSTAT5"] = spliced_df["IL-2_pSTAT5"]
#%% calculate metrics
dict_list = []
hue = "fractions_Tsec"
for f, frac in enumerate(np.sort(spliced_df[hue].unique())):
    scan_df = spliced_df.loc[(spliced_df[hue] == frac)]
    scan_df = scan_df.loc[(scan_df["time_index"] == scan_df["time_index"].max())]
    act_cells = np.zeros((len(scan_df["replicat_index"].unique()),
                             len(scan_df["scan_value"].unique())))
    act_cells[:] = None
    for r, rep in enumerate(scan_df["replicat_index"].unique()):
        rep_df = scan_df.loc[(scan_df["replicat_index"] == rep)]
        for s, scan in enumerate(np.sort(rep_df["scan_value"].unique())):
            df = rep_df.loc[(rep_df["scan_value"] == scan) & (rep_df["type_name"] == "Th")]
            if len(df) > 0:
                act_cells = len(df.loc[df["IL-2_pSTAT5"] > 0.5])/len(df) * 100
            else:
                act_cells = None

            niche_df = rep_df.loc[(rep_df["scan_value"] == scan)]
            DBres, sliced_df, DBcoords = myDBscan(niche_df, 20, with_Tsecs=True)
            sliced_df["cluster"] = DBres
            niche_score = np.sum(~np.isnan(np.where(np.unique(DBres, return_counts=True)[1] > 1))) / len(
                                    niche_df.loc[niche_df["type_name"] == "Tsec"]) if len(
                                    niche_df.loc[niche_df["type_name"] == "Tsec"]) != 0 else np.nan
            effects = []
            for u in np.unique(DBres):
                if len(np.where(DBres == u)[0]) > 1:
                    ids = sliced_df.loc[sliced_df["cluster"] == u, "id"].values

                    in_cluster_values = sliced_df.loc[sliced_df["cluster"] == u, "pSTAT5"].values
                    outside_cluster_values = niche_df.loc[
                        (niche_df["IL-2_pSTAT5"] < 0.5) & (niche_df["type_name"] == "Th"), "pSTAT5"].values

                    if len(in_cluster_values) > 0 and len(outside_cluster_values) > 0:
                        if ~np.isnan(in_cluster_values).all() and ~np.isnan(outside_cluster_values).all():
                            effects.append(np.nanmean(in_cluster_values) / np.nanmean(outside_cluster_values))
            niche_effect = np.nanmean(effects) if len(effects) > 0 else np.nan

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

            dict_list.append({"scan_value": scan, "frac": frac, "replicat_index": rep, "act_cells": act_cells,
                              "niche_score": niche_score, "niche_effect": niche_effect, "surf_c": df["IL-2_surf_c"].mean()*1e3,
                              "signaling_range": signaling_range/df.geometry_distance.unique()[0]})
plot_df = pd.DataFrame(dict_list)
#%%
factor = 0.89
rc_ticks['figure.figsize'] = [1.67475 * factor, 1.386 * factor]
sns.set_theme(context="talk", style="ticks", rc=rc_ticks)
fig, ax = plt.subplots()

fracs  = [0.01, 0.02, 0.03, 0.04]
alphas = np.logspace(-1.05, 0, len(fracs))
linestyles = ["-" for x in range(len(fracs))]
for f, frac in enumerate(fracs):
    color = "#aa0000ff"
    sns.lineplot(data=plot_df.loc[plot_df.frac == frac], x = "scan_value", y = "surf_c", errorbar="se", color=color,
                 alpha=alphas[f], linestyle=linestyles[f], err_kws={"linewidth":0, "alpha":alphas[f]*0.3}, ax=ax)

ax.set(yscale="log", ylabel=r"surface conc. s.d. (pM)", xlabel=r"$\varphi$(FRC)", ylim=(1e-1,1e2),
       xticks=[0, 0.5, 1], xticklabels=[0, 0.5, 1], yticks=[1e-1, 1e0, 1e1, 1e2], yticklabels=[r"10$^{-1}$", r"10$^{0}$", r"10$^{1}$", r"10$^{2}$"], xlim=(0,1))

# plt.ylim(1e0, 4e4)
fig.savefig(saving_path + "FRC_surface_concentration_over_phi" + ".pdf", bbox_inches='tight', transparent=True)
plt.tight_layout()
plt.show()

#%%
rc_ticks['figure.figsize'] = [1.81, 1.553]
sns.set_theme(context="talk", style="ticks", rc=rc_ticks)

fig, ax = plt.subplots()
for f, frac in enumerate(fracs):
    color = "#aa0000ff"
    sns.lineplot(data=plot_df.loc[plot_df.frac == frac], x = "scan_value", y = "act_cells", errorbar="sd", color=color,
                 alpha=alphas[f], linestyle=linestyles[f], err_kws={"linewidth":0, "alpha":alphas[f]*0.3}, ax=ax)

ax.set(ylabel=r"% pSTAT5$^+$", xlabel=r"$\varphi$(FRC)", ylim=(0,20),
       xticks=[0, 0.5, 1], xticklabels=[0, 0.5, 1], yticks=[0, 10, 20], xlim=(0,1))

plt.tight_layout()
fig.savefig(saving_path + "FRC_activation_over_phi" + ".pdf", bbox_inches='tight', transparent=True)
plt.show()
#%%
factor = 1
rc_ticks['figure.figsize'] = [1.67475 * factor, 1.386 * factor]
sns.set_theme(context="talk", style="ticks", rc=rc_ticks)
fig, ax = plt.subplots()

color = "#aa0000ff"
for f, frac in enumerate(fracs):
    sns.lineplot(data=plot_df.loc[plot_df.frac == frac], x = "scan_value", y = "niche_score", errorbar="se", color=color,
                 alpha=alphas[f], linestyle=linestyles[f], err_kws={"linewidth":0, "alpha":alphas[f]*0.3})

plt.ylabel(r"niche score")
plt.xlabel(r"$\varphi$(FRC)")
plt.xticks([0, 0.5, 1])
plt.yscale("linear")
plt.xlim(0, 1)
fig.savefig(saving_path + "FRC_niche_score_over_phi" + ".pdf", bbox_inches='tight', transparent=True)
plt.tight_layout()
plt.show()

#%%
factor = 1
rc_ticks['figure.figsize'] = [1.67475 * factor, 1.386 * factor]
sns.set_theme(context="talk", style="ticks", rc=rc_ticks)
fig, ax = plt.subplots()

color = "#aa0000ff"
for f, frac in enumerate(fracs):
    color = "#aa0000ff"
    if np.sum(~np.isnan(plot_df.loc[plot_df.frac == frac, "niche_effect"].values)) > 0:
        sns.lineplot(data=plot_df.loc[plot_df.frac == frac], x="scan_value", y="niche_effect", errorbar="se", color=color,
                     alpha=alphas[f], linestyle=linestyles[f], err_kws={"linewidth":0, "alpha":alphas[f]*0.3})
plt.ylabel(r"niche effect")
plt.xlabel(r"$\varphi$(FRC)")
plt.xticks([0, 0.5, 1])
plt.yscale("log")
plt.ylim(9, 1.1e2)
plt.yticks([9, 10, 30, 100], ["9", "", "30", "100"])
plt.xlim(0, 1)
fig.savefig(saving_path + "FRC_niche_effect_over_phi" + ".pdf", bbox_inches='tight', transparent=True)
plt.tight_layout()
plt.show()

#%%
factor = 0.89
rc_ticks['figure.figsize'] = [1.67475 * factor, 1.386 * factor]
sns.set_theme(context="talk", style="ticks", rc=rc_ticks)
fig, ax = plt.subplots()

for f, frac in enumerate(fracs):

    color = "#aa0000ff"
    sns.lineplot(data=plot_df.loc[plot_df.frac == frac], x="scan_value", y="signaling_range", errorbar="se", color=color,
                 alpha=alphas[f], linestyle=linestyles[f], err_kws={"linewidth":0, "alpha":alphas[f]*0.3})
plt.ylabel(r"signaling range (c-c-d)")
plt.xlabel(r"$\varphi$(FRC)")
plt.xticks([0, 0.5, 1])
plt.yscale("linear")
plt.ylim(1,2.6)
plt.yticks([1, 1.8, 2.6])
plt.xlim(0, 1)
fig.savefig(saving_path + "FRC_signaling_range_over_phi" + ".pdf", bbox_inches='tight', transparent=True)
plt.tight_layout()
plt.show()