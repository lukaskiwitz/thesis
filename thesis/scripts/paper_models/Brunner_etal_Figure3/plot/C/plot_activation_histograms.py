import pandas as pd
import getpass
import matplotlib.pyplot as plt
import seaborn as sns
import os
from thesis.scripts.paper_models.utilities.plotting_rc import rc_ticks
from thesis.scripts.paper_models.utilities.plot_helper import my_load_df, EC50_calculation

save_plot = True
hdd = "/extra2" if os.path.exists("/extra2") else "/extra"
user = getpass.getuser()
saving_string =r"/home/brunner/Documents/Current work/2023_11_17/"
if not os.path.exists(saving_string):
    os.mkdir(saving_string)

setups = [
    "pos",
    "neg",
    "ODE_pos",
    "ODE_neg"
]

for prefix in setups:
    metrics = ["pSTAT5"]
    h_x_scales = ["linear"]

    if prefix == "neg":
        model_name = "feedback_scan_4"
        name = "dataframes_negative_for_Fig3C_act_plot_steady_state" #michaelis

        path = "/{extra}/{u}/paper_models/kinetics/{mn}/{n}/".format(u=user, n=name, mn=model_name, extra=hdd)
        colours = ["blue", "grey"]
        scan_index = 1
        bins = [25]
    elif prefix == "pos":
        model_name = "feedback_scan_for_time_plot/"
        name = "dataframes_positive_0.1_steady_state" #michaelis

        path = "/{extra}/{u}/paper_models/kinetics/{mn}/{n}/".format(u=user, n=name, mn=model_name, extra=hdd)
        colours = ["red" , "grey"]
        scan_index = 3
        bins = [25]
    elif prefix == "ODE_pos":
        model_name = "Tsec_scan_5"
        name = "gamma_100"
        path = "/{extra}/{u}/paper_models/ODE/saturated/kinetics/{mn}/{n}/".format(u=user, n=name, mn=model_name, extra=hdd)

        colours = ["red" , "grey"]
        scan_index = 4
        bins = [25]
    elif prefix == "ODE_neg":
        model_name = "Tsec_scan_5"
        name = "gamma_0.01"
        path = "/{extra}/{u}/paper_models/ODE/saturated/kinetics/{mn}/{n}/".format(u=user, n=name, mn=model_name, extra=hdd)

        colours = ["blue" , "grey"]
        scan_index = 4
        bins = [25]
    else:
        colours = ["blue", "grey"]
        scan_index = 0
        bins = [25]

    # defines a outer layer of N cells to be ignored in plotting. Used to further limit unwanted boundary effects.
    offset = 0
    try:
        cell_df, global_df =  my_load_df(path, offset=offset, custom_ending = "_combined")
    except FileNotFoundError:
        cell_df = pd.read_hdf(path + "cell_df" + ".h5", mode="r")
    try:
        cell_df["id_id"]
    except KeyError:
        cell_df["id_id"] = cell_df["id"]
    cell_df["IL-2_surf_c"] = cell_df["IL-2_surf_c"].mul(1e3)
    cell_df = cell_df.loc[cell_df["replicat_index"] == 0] #could also choose randomly
    cell_df = cell_df.loc[cell_df["time_index"] != 0]

    try:
        cell_df["pSTAT5"] = cell_df["IL-2_pSTAT5"]
    except:
        cell_df["pSTAT5"] = cell_df["IL-2_surf_c"] ** 3 / (
                    (EC50_calculation(E_max=125e-12, E_min=0, k=860, N=1.5, R=cell_df["IL-2_R"]) * 1e12) ** 3 + cell_df[
                            "IL-2_surf_c"] ** 3).values

    #%%
    print("plotting")
    if scan_index in cell_df["scan_index"].unique():
        hist_df = cell_df.loc[(cell_df["scan_index"] == scan_index) & (cell_df["type_name"] == "Th")]
    else:
        if prefix == "pos":
            hist_df = cell_df.loc[(cell_df["misc_gamma"] == 100.) & (cell_df["type_name"] == "Th") & (cell_df["fractions_Tsec"] == 0.1)]
            hist_df = hist_df.loc[hist_df["scan_name_scan_name"] == hist_df["scan_name_scan_name"].unique()[0]]
        elif prefix == "neg":
            hist_df = cell_df.loc[(cell_df["misc_gamma"] == 0.01) & (cell_df["type_name"] == "Th") & (cell_df["fractions_Tsec"] == 0.05)]
            hist_df = hist_df.loc[hist_df["scan_name_scan_name"] == hist_df["scan_name_scan_name"].unique()[0]]
        elif prefix == "ODE_pos":
            hist_df = cell_df.loc[(cell_df["misc_gamma"] == 100) & (cell_df["type_name"] == "Th") & (cell_df["fractions_Tsec"] == 0.101)]
            hist_df = hist_df.loc[hist_df["scan_name_scan_name"] == hist_df["scan_name_scan_name"].unique()[0]]
        elif prefix == "ODE_neg":
            hist_df = cell_df.loc[(cell_df["misc_gamma"] == 0.01) & (cell_df["type_name"] == "Th") & (cell_df["fractions_Tsec"] == 0.033)]
            hist_df = hist_df.loc[hist_df["scan_name_scan_name"] == hist_df["scan_name_scan_name"].unique()[0]]
    hist_df = hist_df.loc[(hist_df["time_index"] == hist_df["time_index"].max())]
    hist_df["activated"] = False
    act_ids = hist_df.loc[(hist_df["pSTAT5"] > 0.5), "id"]
    hist_df.loc[hist_df["id"].isin(act_ids), "activated"] = True

    rc_ticks["figure.figsize"] = (1.67475 * 1.27, 0.35)

    h_colours = [x for x in reversed(colours)]
    metrics = ["pSTAT5"]
    for m, metric in enumerate(metrics):
        sns.set_theme(context="talk", style="ticks", rc=rc_ticks)
        fig, ax = plt.subplots()
        ax_1 = sns.histplot(data=hist_df, x=metric, bins=bins[m], log_scale=False if h_x_scales[m] == "linear" else True, hue="activated",
                            palette=h_colours, legend=False, linewidth=0.3, alpha=1, binwidth=1/bins[0])
        ax_1.set(yscale="log")
        ax_1.set(xlabel="pSTAT5 normalized", ylabel="Count", ylim=(0.5, 1e3), xlim=(0, 1.0), xticks=[0,0.5,1],
                 xticklabels=[0, 0.5, 1], yticks=[1e0, 1e1, 1e2, 1e3], yticklabels = [r"10$^0$", "", "", r"10$^3$"])
        plt.axvline(0.5,0, 1, color="black", linewidth=0.5)
        plt.text(0.5, 0.63, "well-mixed" if "ODE" in prefix else "RD-system", ha="center",fontsize=6,transform = plt.gca().transAxes, bbox={"pad": 0.5, "color": "white"})
        plt.minorticks_off()
        if save_plot == True:
            fig.savefig(saving_string + prefix + f"_histogram_of_{metric}.pdf", bbox_inches='tight', transparent=True)
        plt.tight_layout()
        plt.show()
