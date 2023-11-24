import matplotlib.pyplot as plt
from mpl_toolkits.axes_grid1.inset_locator import inset_axes
import numpy as np
import pandas as pd
import seaborn as sns
from scipy.spatial import distance_matrix
from thesis.scripts.paper_models.utilities.plotting_rc import rc_ticks
from thesis.scripts.paper_models.utilities.plot_helper import myDBscan

saving_path = f"/home/brunner/Documents/Current work/2023_11_17/"

path = "/extra2/brunner/clustering_both/300_3D_fine/dataframes_Fig5C_more_reps_2_steady_state/"

try:
    cell_df = pd.read_hdf(path + "cell_df_combined" + ".h5", mode="r")
except:
    cell_df = pd.read_pickle(path + 'cell_df_combined' + '.pkl')


treg_scan_names = cell_df.loc[cell_df.type_name == "Treg", "scan_name_scan_name"].unique()
spliced_df = cell_df.loc[cell_df.scan_name_scan_name.isin(treg_scan_names)]

spliced_df["pSTAT5"] = spliced_df["IL-2_pSTAT5"]
#%% calculate metrics
dict_list = []
hue = "fractions_Treg"
for f, frac in enumerate(np.sort(spliced_df[hue].unique())):
    scan_df = spliced_df.loc[(spliced_df[hue] == frac) & (spliced_df["time_index"] == spliced_df["time_index"].max())]
    if frac == 0:
        scan_df = scan_df.loc[scan_df.scan_value == 1]
    act_cells = np.zeros((len(scan_df["replicat_index"].unique()),
                             len(scan_df["scan_value"].unique())))
    Treg_0_act_cells = np.zeros((len(scan_df["replicat_index"].unique()),
                             len(scan_df["scan_value"].unique())))
    act_cells[:] = None
    Treg_0_act_cells[:] = None
    for r, rep in enumerate(scan_df["replicat_index"].unique()):
        rep_df = scan_df.loc[(scan_df["replicat_index"] == rep)]
        for s, scan in enumerate(np.sort(rep_df["scan_value"].unique())):
            df = rep_df.loc[(rep_df["scan_value"] == scan) & (rep_df["type_name"] == "Th")]
            if len(df) > 0:
                act_cells = len(df.loc[df["IL-2_pSTAT5"] > 0.5])/len(df) * 100
            else:
                act_cells = None

            niche_df = rep_df.loc[(rep_df["scan_value"] == scan) & (rep_df["type_name"] != "Treg")]
            DBres, sliced_df, DBcoords = myDBscan(niche_df, niche_df.geometry_distance.unique()[0] + 8, with_Tsecs=True) # = 20, doesnt really matter
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

            DBres, sliced_df, DBcoords = myDBscan(rep_df.loc[(rep_df["scan_value"] == scan) & (rep_df["type_name"] != "Treg")], niche_df.geometry_distance.unique()[0] + 8, with_Tsecs=True)
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
                 alpha=alphas[f], linestyle=linestyles[f], err_kws={"linewidth":0, "alpha":alphas[f]*0.3})
plt.ylabel(r"surface conc. s.d. (pM)")
plt.xlabel(r"$\varphi$(T$_{\rm reg}$)")
plt.yscale("linear")
plt.xticks([0, 0.5, 1], ["0", "0.5", "1"])
plt.xlim(0, 1)
plt.yticks([0,0.3, 0.6])
fig.savefig(saving_path + "Treg_surface_concentration_over_phi" + ".pdf", bbox_inches='tight', transparent=True)
plt.tight_layout()
plt.show()

#%%
rc_ticks['figure.figsize'] = [1.795, 1.553]
sns.set_theme(context="talk", style="ticks", rc=rc_ticks)

fig, ax = plt.subplots()
for f, frac in enumerate(fracs):
    color = "#aa0000ff"
    sns.lineplot(data=plot_df.loc[plot_df.frac == frac], x = "scan_value", y = "act_cells", errorbar="sd", color=color,
                 alpha=alphas[f], linestyle=linestyles[f], err_kws={"linewidth":0, "alpha":alphas[f]*0.3}, ax=ax)
ax.set(ylabel=r"% pSTAT5$^+$", xlabel=r"$\varphi$(T$_{\rm reg}$)", ylim=(0,3),
       xticks=[0, 0.5, 1], xticklabels=[0, 0.5, 1], yticks=[0, 1.5, 3], xlim=(0,1))

plt.tight_layout()
fig.savefig(saving_path + "Treg_activation_over_phi" + ".pdf", bbox_inches='tight', transparent=True)
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
plt.xlabel(r"$\varphi$(T$_{\rm reg}$)")
plt.xticks([0, 0.5, 1])
plt.yscale("linear")
plt.xlim(0, 1)
fig.savefig(saving_path + "Treg_niche_score_over_phi" + ".pdf", bbox_inches='tight', transparent=True)
plt.tight_layout()
plt.show()

#%%
factor = 1
rc_ticks['figure.figsize'] = [1.67475 * factor, 1.386 * factor]
sns.set_theme(context="talk", style="ticks", rc=rc_ticks)
fig, ax = plt.subplots()

color = "#aa0000ff"
for f, frac in enumerate(fracs):
    frac_df = plot_df.loc[plot_df.frac == frac]
    if np.sum(~np.isnan(frac_df["niche_effect"].values)) > 0:
        scan_values = frac_df["scan_value"].unique()
        niche_effects = []
        niche_effect_se = []
        for s, sv in enumerate(scan_values):
            niche_effects.append(frac_df.loc[frac_df.scan_value == sv, "niche_effect"].mean())
            niche_effect_se.append(frac_df.loc[frac_df.scan_value == sv, "niche_effect"].std()/np.sqrt(len(frac_df.replicat_index.unique())))
        niche_effects = np.array(niche_effects)
        niche_effect_se = np.array(niche_effect_se)
        plt.plot(scan_values, niche_effects, color=color, alpha=alphas[f], linestyle=linestyles[f])
        plt.fill_between(scan_values, np.clip(niche_effects - niche_effect_se, 0, None), np.clip(niche_effects + niche_effect_se, 0, None),
                         color=color, alpha=0.3 * alphas[f], linewidth=0.0)
plt.ylabel(r"niche effect")
plt.xlabel(r"$\varphi$(T$_{\rm reg}$)")
plt.xticks([0, 0.5, 1])
plt.yscale("log")
plt.yticks([2e3, 1e4, 4e4], [r"2x10$^3$", r"10$^4$",r"4x10$^4$"])
plt.xlim(0, 1)
fig.savefig(saving_path + "Treg_niche_effect_over_phi" + ".pdf", bbox_inches='tight', transparent=True)
plt.tight_layout()
plt.show()

#%%
factor = 0.89
rc_ticks['figure.figsize'] = [1.67475 * factor, 1.386 * factor]
sns.set_theme(context="talk", style="ticks", rc=rc_ticks)
fig, ax = plt.subplots()

for f, frac in enumerate(fracs):
    frac_df = plot_df.loc[plot_df.frac == frac]
    color = "#aa0000ff"
    scan_values = frac_df["scan_value"].unique()
    SRs = []
    SR_se = []
    for s, sv in enumerate(scan_values):
        SRs.append(frac_df.loc[frac_df.scan_value == sv, "signaling_range"].mean())
        SR_se.append(frac_df.loc[frac_df.scan_value == sv, "signaling_range"].std()/np.sqrt(len(frac_df.replicat_index.unique())))
    SRs = np.array(SRs)
    SR_se = np.array(SR_se)
    plt.plot(scan_values, SRs, color=color, alpha=alphas[f], linestyle=linestyles[f])
    plt.fill_between(scan_values, np.clip(SRs - SR_se, 0, None), np.clip(SRs + SR_se, 0, None),
                     color=color, alpha=0.3 * alphas[f], linewidth=0.0)
plt.ylabel(r"signaling range (c-c-d)")
plt.xlabel(r"$\varphi$(T$_{\rm reg}$)")
plt.xticks([0, 0.5, 1])
plt.yscale("linear")
plt.xlim(0, 1)
plt.yticks([1.1, 1.5 ,1.9])
fig.savefig(saving_path + "Treg_signaling_range_over_phi" + ".pdf", bbox_inches='tight', transparent=True)
plt.tight_layout()
plt.show()