import getpass
import matplotlib.pyplot as plt
import numpy as np
import seaborn as sns
import os
from thesis.scripts.paper_models.utilities.plotting_rc import rc_ticks
from thesis.scripts.paper_models.utilities.plot_helper import my_load_df, ashmans_D, wangs_BI, my_interpolation

def EC50_calculation(E_max, E_min, k, N, R):
    return (E_max * k ** N + E_min * R ** N) / (k ** N + R ** N)

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
        name = "dataframes_negative_for_Fig3C_act_plot_2_timeseries"

        path = "/{extra}/{u}/paper_models/kinetics/{mn}/{n}/".format(u=user, n=name, mn=model_name, extra=hdd)
        colours = ["blue", "grey"]
        Tsec_fraction = 0
        bins = [25]
        color = "blue"
    elif prefix == "pos":
        model_name = "feedback_scan_for_time_plot/"
        name = "dataframes_positive_0.1_timeseries"

        path = "/{extra}/{u}/paper_models/kinetics/{mn}/{n}/".format(u=user, n=name, mn=model_name, extra=hdd)
        colours = ["red", "grey"]
        Tsec_fraction = 0
        bins = [25]
        color = "red"

    cell_df, global_df = my_load_df(path, offset=0, custom_ending="_combined")
    cell_df["IL-2_surf_c"] *= 1e3
    try:
        cell_df["IL-2_gamma"]
    except KeyError:
        cell_df["IL-2_gamma"] = cell_df["misc_gamma"]
        cell_df["IL-2_Tsec_fraction"] = cell_df["fractions_Tsec"]
    cell_df = cell_df.loc[
        (cell_df["IL-2_Tsec_fraction"] == np.sort(cell_df["IL-2_Tsec_fraction"].unique())[Tsec_fraction])]

    #%%
    plot_x = []
    plot_times = []
    try:
        cell_df["pSTAT5"] = cell_df["IL-2_pSTAT5"]
    except:
        print("recalculating activation")
        cell_df["pSTAT5"] = cell_df["IL-2_surf_c"] ** 3 / (
                (EC50_calculation(E_max=125e-12, E_min=0, k=860, N=1.5, R=cell_df["IL-2_R"]) * 1e12) ** 3 + cell_df[
                        "IL-2_surf_c"] ** 3).values

    cell_df = cell_df.loc[cell_df["type_name"] == "Th"]
    print(cell_df["IL-2_Tsec_fraction"].unique())
    avg_times = np.zeros((len(cell_df["IL-2_gamma"].unique()),
                    len(cell_df["replicat_index"].unique())))
    avg_times[:] = np.nan

    for g, gamma in enumerate(np.sort(cell_df["IL-2_gamma"].unique())):
        gamma_df = cell_df.loc[(cell_df["IL-2_gamma"] == gamma)]
        for r, rep in enumerate(np.sort(gamma_df["replicat_index"].unique())):
            rep_df = gamma_df.loc[(gamma_df["replicat_index"] == rep)]
            ptable = rep_df.loc[rep_df["type_name"] != "Tsec"].pivot_table('pSTAT5', index='id', columns='time')
            times = rep_df["time"].unique() / 3600
            where_object = np.where(ptable > 0.5, 1, 0)
            cell_times = np.zeros((len(ptable), 2))

            for index, id in enumerate(ptable.index):
                res = times[np.where(where_object[index] == 1)[0]]
                if len(np.where(where_object[index] == 1)[0]) > 0:
                    if np.where(where_object[index] == 1)[0][-1] != len(times) - 1:
                        res = []
                cell_times[index][0] = id
                cell_times[index][1] = np.min(res) if len(res) > 0 else None
            avg_times[g][r] = np.nanmean(cell_times[:, 1])/24

    plot_x.append(np.sort(cell_df["IL-2_gamma"].unique()))
    plot_times.append(avg_times)

    #%%
    print("plotting")
    rc_ticks['figure.figsize'] = (1.67475 * 1.15, 1.386 * 0.5)
    ylabel = r"time to activation (d)"
    for e, entry in enumerate(plot_times):
        sns.set_theme(context="talk", style="ticks", rc=rc_ticks)
        fig, ax = plt.subplots()

        alphas = np.logspace(-1, 0, len(entry))
        y = np.nanmean(entry, axis=1)
        std = np.nanstd(entry, axis=1)
        x = plot_x[e]
        if prefix == "neg":
            x = 1 / x
            x = x[::-1]
            y = y[::-1]
            std = std[::-1]
        plt.plot(x, y, color=color, marker=".")
        plt.fill_between(np.unique(x),np.clip(y - std, 0, None),np.clip(y + std, 0, None), color=color,
                         alpha=0.3, linewidth=0.0)
    plt.ylabel(ylabel)
    plt.xlabel("feedback fold change")
    plt.yscale("linear")
    plt.xscale("log")
    plt.ylim((0, 3))
    plt.yticks([0,1.5,3])
    plt.xlim(1,100)

    if save_plot == True:
        fig.savefig(saving_string + prefix + f"_time_to_act_over_fb.pdf", bbox_inches='tight', transparent=True)
    plt.tight_layout()
    plt.show()