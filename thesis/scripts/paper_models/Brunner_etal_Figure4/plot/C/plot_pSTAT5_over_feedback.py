import getpass
import matplotlib.pyplot as plt
import numpy as np
import seaborn as sns
import os
from thesis.scripts.paper_models.utilities.plotting_rc import rc_ticks
from thesis.scripts.paper_models.utilities.plot_helper import my_load_df, my_interpolation, myDBscan, EC50_calculation

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

plot_inside = []
plot_outside = []
plot_x = []
for prefix in setups:
    if prefix == "pos":
        model_name = "feedback_scan_4/"
        name = "dataframes_positive_for_Fig3C_act_plot_steady_state"
        Tsec_fraction = 0
        color = "red"
    elif prefix == "neg":
        model_name = "feedback_scan_4"
        name = "dataframes_negative_for_Fig3C_act_plot_2_steady_state"
        Tsec_fraction = 0
        color = "blue"

    path = "/{extra}/{u}/paper_models/kinetics/{mn}/{n}/".format(u=user, n=name, mn=model_name, extra=hdd)
    cell_df, global_df =  my_load_df(path, offset=0, custom_ending = "_combined")
    cell_df["IL-2_surf_c"] *= 1e3
    try:
        cell_df["IL-2_gamma"]
    except KeyError:
        cell_df["IL-2_gamma"] = cell_df["misc_gamma"]
        cell_df["IL-2_Tsec_fraction"] = cell_df["fractions_Tsec"]
    cell_df = cell_df.loc[(cell_df["IL-2_Tsec_fraction"] == np.sort(cell_df["IL-2_Tsec_fraction"].unique())[Tsec_fraction])]

    cell_df = cell_df.loc[(cell_df["time_index"] == cell_df["time_index"].max())]
    try:
        cell_df["pSTAT5"] = cell_df["IL-2_pSTAT5"]
    except:
        print(f"calculating own activation for c = {c}")
        cell_df["pSTAT5"] = cell_df["IL-2_surf_c"] ** 3 / (
                (EC50_calculation(E_max=125e-12, E_min=0, k=860, N=1.5, R=cell_df["IL-2_R"]) * 1e12) ** 3 + cell_df[
                        "IL-2_surf_c"] ** 3).values
    niche_effect = np.zeros((len(cell_df.replicat_index.unique()),
                             len(cell_df["IL-2_gamma"].unique())))
    outside_vals = np.zeros((len(cell_df.replicat_index.unique()),
                             len(cell_df["IL-2_gamma"].unique())))
    inside_vals = np.zeros((len(cell_df.replicat_index.unique()),
                            len(cell_df["IL-2_gamma"].unique())))

    niche_effect[:] = np.nan
    outside_vals[:] = np.nan
    inside_vals[:] = np.nan
    for g, gamma in enumerate(np.sort(cell_df["IL-2_gamma"].unique())):
        gamma_df = cell_df.loc[(cell_df["IL-2_gamma"] == gamma)]
        print("gamma:", gamma)
        for r, replicat in enumerate(np.sort(gamma_df["replicat_index"].unique())[:20]):
            rep_df = gamma_df.loc[(gamma_df["replicat_index"] == replicat)]
            DBres, sliced_df, DBcoords = myDBscan(rep_df, 20, with_Tsecs=True)
            sliced_df["cluster"] = DBres
            effects = []
            in_vals = []
            out_vals = []
            for u in np.unique(DBres):
                if len(np.where(DBres == u)[0]) > 1:
                    ids = sliced_df.loc[sliced_df["cluster"] == u, "id"].values
                    in_cluster_values = sliced_df.loc[(sliced_df["cluster"] == u) & (sliced_df.type_name == "Th"), "pSTAT5"].values
                    outside_cluster_values = rep_df.loc[
                        (rep_df["IL-2_pSTAT5"] < 0.5) & (rep_df["type_name"] == "Th"), "pSTAT5"].values

                    if len(in_cluster_values) > 0 and len(outside_cluster_values) > 0:
                        if ~np.isnan(in_cluster_values).all() and ~np.isnan(outside_cluster_values).all():
                            effects.append(np.nanmean(in_cluster_values) / np.nanmean(outside_cluster_values))
                            in_vals.append(in_cluster_values.astype(float))
                            out_vals.append(outside_cluster_values.astype(float))
            niche_effect[r][g] = np.nanmean(effects) if len(effects) > 0 else np.nan
            inside_vals[r][g] = np.nanmean(np.concatenate(in_vals) if len(in_vals) > 0 else np.nan)
            outside_vals[r][g] = np.nanmean(np.concatenate(out_vals) if len(out_vals) > 0 else np.nan)

    plot_x.append(np.sort(cell_df["IL-2_gamma"].unique()))
    plot_inside.append(inside_vals)
    plot_outside.append(outside_vals)

#%%
print("plotting")
rc_ticks['figure.figsize'] = [1.67475 * 1.25, 0.6]
sns.set_theme(context="talk", style="ticks", rc=rc_ticks)
fig, ax = plt.subplots()

ylabel = r"pSTAT5"
colors = ["red", "blue"]
prefixes = ["pos", "neg"]
ylims = [[[0.75, 1], [0, 0.25]], [[0.78, 0.88], [0, 0.1]]]
for v, value in enumerate(plot_x):
    alphas = np.logspace(-1, 0, len(plot_inside[v]))
    inside = np.nanmean(plot_inside[v], axis=0)
    outside = np.nanmean(plot_outside[v], axis=0)
    std = np.nanstd(plot_inside[v], axis=0)/np.sqrt(len(plot_outside[v]))
    x = value
    if v == 1:
        x = 1 / x
        x = x[::-1]
        inside = inside[::-1]
        outside = outside[::-1]
        std = std[::-1]
    ax.plot(x, inside, color=colors[v])
    ax.fill_between(np.unique(x),
                 np.clip(inside - std, 0, None),
                 np.clip(inside + std, 0, None), color=colors[v],
                     alpha=0.3, linewidth=0.0)
    ax.set_ylabel(ylabel)
    ax.set_xlabel("")
    ax.set_xscale("log")
    ax.set_yscale("linear")
    ax.set_ylim([0.7, 1])
    ax.set_xlim([1, 100])
    ax.set_yticks([0.7, 1])
    ax.set_xticks([1, 10, 100], ["", "", ""])
if save_plot == True:
    fig.savefig(saving_string + f"inside_pSTAT5_inside_over_fb.pdf", bbox_inches='tight', transparent=True)
plt.tight_layout()
plt.show()

#%%
rc_ticks['figure.figsize'] = [1.67475 * 1.25, 0.6]
sns.set_theme(context="talk", style="ticks", rc=rc_ticks)
fig, ax = plt.subplots()
ylabel = r"pSTAT5"
colors = ["red", "blue"]
for v, value in enumerate(plot_x):
    alphas = np.logspace(-1, 0, len(plot_inside[v]))
    outside = np.nanmean(plot_outside[v], axis=0)
    std = np.nanstd(plot_outside[v], axis=0)/np.sqrt(len(plot_outside[v]))
    x = value
    if v == 1:
        x = 1 / x
        x = x[::-1]
        outside = outside[::-1]
        std = std[::-1]
    plt.plot(x, outside, color=colors[v])
    plt.fill_between(np.unique(x),
                 np.clip(outside - std, 0, None),
                 np.clip(outside + std, 0, None), color=colors[v],
                     alpha=0.3, linewidth=0.0)
    plt.ylabel(ylabel)
    plt.xlabel("feedback fold change")
    plt.yscale("linear")
    plt.xscale("log")
    plt.ylim(0, 0.1)
    plt.yticks([0, 0.1])
    plt.xlim(1,100)

if save_plot == True:
    fig.savefig(saving_string + f"outside_pSTAT5_outside_over_fb.pdf", bbox_inches='tight', transparent=True)
plt.tight_layout()
plt.show()