import getpass
import matplotlib.pyplot as plt
import seaborn as sns
import numpy as np
import os
from thesis.scripts.paper_models.utilities.plotting_rc import rc_ticks
from thesis.scripts.paper_models.utilities.plot_helper import my_load_df, my_interpolation, myDBscan, EC50_calculation

rc_ticks['figure.figsize'] = [1.67475 * 1.25, 1.386]
sns.set_theme(context="talk", style="ticks", rc=rc_ticks)
fig, ax = plt.subplots()

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
    if prefix == "pos":
        model_name = "Tsec_scan_for_ns_plot"
        name = "dataframes_gamma_100_steady_state/"
        color = "red"
    elif prefix == "neg":
        model_name = "Tsec_scan_for_ns_plot"
        name = "dataframes_gamma_0.5_steady_state/"
        color = "blue"

    path = "/{extra}/{u}/paper_models/kinetics/{mn}/{n}/".format(u=user, n=name, mn=model_name, extra=hdd)
    cell_df, global_df = my_load_df(path, offset=0, custom_ending="")
    cell_df["IL-2_Tsec_fraction"] = cell_df["fractions_Tsec"]
    cell_df["IL-2_surf_c"] *= 1e3
    try:
        cell_df["pSTAT5"] = cell_df["IL-2_pSTAT5"]
    except:
        print("pSTAT5 loading failed, recalculating")
        cell_df["pSTAT5"] = cell_df["IL-2_surf_c"] ** 3 / (
                (EC50_calculation(E_max=125e-12, E_min=0, k=860, N=1, R=cell_df["IL-2_R"]) * 1e12) ** 3 + cell_df[
                        "IL-2_surf_c"] ** 3).values

    cell_df = cell_df.loc[(cell_df["time_index"] == cell_df["time_index"].max())]
    replicats_to_be_plotted = len(cell_df["replicat_index"].unique())
    no_of_niches = np.zeros((replicats_to_be_plotted, len(cell_df["IL-2_Tsec_fraction"].unique()))) # amount_of_niches
    no_of_niches[:] = np.nan
    for r, replicat in enumerate(np.sort(cell_df["replicat_index"].unique())):
        for f, frac in enumerate(np.sort(cell_df["IL-2_Tsec_fraction"].unique())):
            frac_df = cell_df.loc[(cell_df["IL-2_Tsec_fraction"] == frac) & (cell_df["replicat_index"] == replicat)]
            DBres, sliced_df, DBcoords = myDBscan(frac_df, 20, with_Tsecs = True)

            no_of_niches[r][f] = np.sum(~np.isnan(np.where(np.unique(DBres, return_counts=True)[1] > 1)))

    x_axes = cell_df["IL-2_Tsec_fraction"].unique()
    plot_x = x_axes[x_axes.argsort()]
    unique_indices = np.unique(plot_x, return_index=True)[1]
    plot_x = plot_x[unique_indices] * 1000
    mean = np.mean(no_of_niches, axis=0)
    std = np.std(no_of_niches, axis=0)/np.sqrt(len(cell_df["replicat_index"].unique()))
    my_interpolation(plot_x, mean, smoothing_factor = 3e-3, color = color, fill = True, std = std)

plt.plot(np.linspace(0,0.4,1000) * 1e3, np.linspace(0,0.4,1000) * 1000, color="black")
plt.xlim((0, 200))
plt.ylim((0, 30))
plt.xticks([0, 100, 200])
plt.yticks([0, 15, 30])
plt.xlabel("no. of secreting cells")
plt.ylabel("no. of niches")
if save_plot == True:
    fig.savefig(saving_string + "no_of_niches_over_secreting_cells" + ".pdf", bbox_inches='tight', transparent=True)
plt.tight_layout()
plt.show()

