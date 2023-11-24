import matplotlib.pyplot as plt
import numpy as np
import seaborn as sns
import os
import pandas as pd
from brokenaxes import brokenaxes
from thesis.scripts.paper_models.utilities.plotting_rc import rc_ticks
rc_ticks["xtick.labelsize"] = 7
rc_ticks["ytick.labelsize"] = 6
from thesis.scripts.paper_models.Brunner_etal_Figure2.plot.funcs_for_plotting import get_path, scale_dataframes
from thesis.scripts.paper_models.utilities.plot_helper import my_load_df
# bc = "standard"
bc = "saturated"

scan_variables = ["IL-2_sigma"]

# define the standard values (fold-change = 1)
standards = [1]

models = ["well_mixed", "loc. q and R"]

plot_legend = True

saving_string = r"/home/brunner/Documents/Current work/2023_09_22/" + "pSTAT5_{model}_{m}_static_{sv}_{bc}".format(model = "ODE" if models[0] == "well_mixed" else "spatial", m = "box", sv=scan_variables[0], bc = bc) + ".pdf"


replicat = "replicat_index"

# load runs, apply, merge dataframes
dataframes = [[] for x in models]
print("loading data")
hdd = "/extra2" if os.path.exists("/extra2") else "/extra"
for m,model in enumerate(models):
    for sv, scan_variable in enumerate(scan_variables):
        m_sv_path = get_path(bc, hdd, model, scan_variable)
        c_df, g_df = my_load_df(m_sv_path, custom_ending = "_combined")
        try:
            c_df["IL-2_sigma"]
        except KeyError:
            c_df["IL-2_sigma"] = c_df["misc_sigma"]
            try:
                g_df["IL-2_sigma"] = g_df["scan_value"]
            except:
                pass
        # scales the dataframes. Depends on your desired units. Currently everything is in pM
        c_df, g_df = scale_dataframes(c_df, g_df, model, scan_variable, sv, standards)
        dataframes[m].append([c_df, g_df])
#%%
########################################################################################################################
print("plotting")
########################################################################################################################
fc_to_plot = [1/10, 1, 10]

def EC50_calculation(E_max, E_min, k, N, R):
    return (E_max * k ** N + E_min * R ** N) / (k ** N + R ** N)

df_list = []
for m,model in enumerate(models):
    x = np.arange(0, len(scan_variables) * 2)
    for sv,scan_variable in enumerate(scan_variables):
        c_df = dataframes[m][sv][0]
        g_df = dataframes[m][sv][1]

        ste = []
        act = []
        surf_c = []
        cv = []
        for f, fc in enumerate(fc_to_plot):
            scan = c_df[scan_variable].unique()[(np.abs(c_df[scan_variable].unique() - (standards[sv] * fc))).argmin()]
            sliced_df = c_df.loc[(c_df["type_name"] == "Th") & (c_df[scan_variable] == scan)]
            sliced_df["pSTAT5"] = sliced_df["IL-2_surf_c"] ** 3 / (
                    (EC50_calculation(E_max=125e-12, E_min=0, k=860, N=0.8, R=sliced_df["IL-2_R"]) * 1e12) ** 3 +
                    sliced_df["IL-2_surf_c"] ** 3).values
            run_act = []
            try:
                for r, run in enumerate(np.sort(sliced_df[replicat].unique())):
                    run_df = sliced_df.loc[(sliced_df[replicat] == run)]
                    run_act.append(len(run_df.loc[(run_df["pSTAT5"] >= 0.5)].values)/len(run_df.values))
                    df_list.append({"model": model, "scan_variable": scan_variable, "fold_change": fc, "replicat": run, "activation": run_act[-1]})
            except KeyError: # ODE system does not have a gradient
                run_act = [0]

            ste.append(np.std(run_act) / np.sqrt(len(sliced_df[replicat].unique())))
            act.append(np.mean(run_act))
            surf_c.append(sliced_df["IL-2_surf_c"].mean())
            cv.append(sliced_df["IL-2_R"].std()/sliced_df["IL-2_R"].mean())
plot_df = pd.DataFrame(df_list)
#%%
from copy import deepcopy
copy_df = deepcopy(plot_df)
copy_df["activation"] *= 100

rc_ticks["figure.figsize"] = (0.55, 1.2)
rc_ticks['xtick.major.size'] = 2.5
rc_ticks['xtick.major.pad'] = 2.5
rc_ticks['axes.titlesize'] = 6
rc_ticks['xtick.labelsize'] = 5.5

sns.set_theme(context = "talk", style = "ticks", rc = rc_ticks)
fig,ax = plt.subplots()

from matplotlib.colors import to_rgba
base_colors = ["black", "black"]
base_alphas = np.logspace(np.log10(0.15),0,len(fc_to_plot))
model_colors = []
for bc in base_colors:
    model_colors.append([to_rgba(bc, alpha) for alpha in base_alphas])

for sv, scan in enumerate(copy_df.scan_variable.unique()):
    scan_df = copy_df.loc[(copy_df.scan_variable == scan)]
    xs = np.array(np.arange(len(scan_df.fold_change.unique())))
    for m, model in enumerate(scan_df.model.unique()):
        model_df = scan_df.loc[(scan_df.model == model)]
        xs += (np.max(xs) + 2)*m
        ys = model_df.groupby("fold_change").mean()["activation"].values
        yerr = model_df.groupby("fold_change").std()["activation"].values
        plt.bar(xs, ys, yerr=yerr, linewidth=0.2, width=0.8, color=model_colors[m], edgecolor="black", error_kw=dict(lw=0.3))

ax.set(ylabel=r"pSTAT$^+$ (%)", xlabel="", title="R hetero.", ylim=(-1/6,12))
first = np.arange(len(scan_df.fold_change.unique()))[int(np.floor(len(scan_df.fold_change.unique())/2))]
second = xs[int(np.floor(len(xs)/2))]
plt.xticks([first, second], [r"well-" "\n" r"mixed", "RD"])
plt.yscale("linear")
ax.set_yticks([0, 6, 12])
fig.savefig(saving_string, bbox_inches='tight', transparent=True)
plt.tight_layout()
plt.show()
