import matplotlib.pyplot as plt
import numpy as np
import seaborn as sns
import os
import pandas as pd
from copy import deepcopy
from thesis.scripts.paper_models.utilities.plotting_rc import rc_ticks
from thesis.scripts.paper_models.Brunner_etal_Figure2.plot.funcs_for_plotting import get_path, scale_dataframes
from thesis.scripts.paper_models.utilities.plot_helper import my_load_df

# bc = "standard"
bc = "saturated"

scan_variables = ["IL-2_Tsec_fraction", "IL-2_sigma", "IL-2_KD"]

saving_string = r"/home/brunner/Documents/Current work/2023_11_03/" + "SD_{m}_static_{sv}_{bc}".format(m = "box", sv=scan_variables[0], bc = bc) + ".pdf"

sv_styles = ['-','-.',':']
plot_names = ["sec. scan", "R scan", "KD scan"]
# define the standard values (fold-change = 0) for each of the scan_variables
standards = [0.05, 1, 7.4]

# which models to plot

# models = ["loc. q and R", "loc. q", "well_mixed"]
models = ["loc. q and R"]
# defines a outer layer of N cells to be ignored in plotting. Used to further limit unwanted boundary effects.


replicat = "replicat_index"

# load runs, apply , merge dataframes
dataframes = [[] for x in models]
print("loading data")
hdd = "/extra2" if os.path.exists("/extra2") else "/extra"
for m,model in enumerate(models):
    for sv, scan_variable in enumerate(scan_variables):
        m_sv_path = get_path(bc, hdd, model, scan_variable)
        print(m_sv_path)
        c_df, g_df = my_load_df(m_sv_path, custom_ending = "_combined")
        try:
            c_df["IL-2_Tsec_fraction"]
        except KeyError:
            c_df["IL-2_Tsec_fraction"] = c_df["fractions_Tsec"]
            try:
                g_df["IL-2_Tsec_fraction"] = g_df["fractions_Tsec"]
            except:
                pass
        try:
            c_df["IL-2_sigma"]
        except KeyError:
            c_df["IL-2_sigma"] = c_df["misc_sigma"]
            try:
                g_df["IL-2_sigma"] = g_df["scan_value"]
            except:
                pass
        c_df, g_df = scale_dataframes(c_df, g_df, model, scan_variable, sv, min_value, max_value, standards)
        dataframes[m].append([c_df, g_df])

#%%
########################################################################################################################
print("plotting")
########################################################################################################################
fc_to_plot = [1/5, 1, 5]
sns.set_theme(context = "talk", style = "ticks", rc = rc_ticks)
fig,ax = plt.subplots()
df_list = []
for m,model in enumerate(models):
    x = np.arange(0, len(scan_variables) * 2)
    for sv,scan_variable in enumerate(scan_variables):
        c_df = dataframes[m][sv][0]
        g_df = dataframes[m][sv][1]
        ste = []
        cv = []
        for f, fc in enumerate(fc_to_plot):
            scan = c_df[scan_variable].unique()[(np.abs(c_df[scan_variable].unique() - (standards[sv] * fc))).argmin()]
            print((np.abs(c_df[scan_variable].unique() - (standards[sv] * fc))).argmin())
            sliced_df = c_df.loc[(c_df["type_name"] == "Th") & (c_df[scan_variable] == scan)]
            run_cv = []
            mean = []
            std = []
            for r, run in enumerate(np.sort(sliced_df[replicat].unique())):
                std.append(sliced_df.loc[(sliced_df[replicat] == run), scan_measure].std())
                mean.append(sliced_df.loc[(sliced_df[replicat] == run), scan_measure].mean())
                run_cv.append(std[-1] / mean[-1])
                df_list.append({"model": model, "scan_variable": scan_variable, "fold_change": fc, "replicat": run,
                                "cv": run_cv[-1], "std": std[-1], "mean": mean[-1]})
            ste.append(np.std(run_cv))# / np.sqrt(len(g_df[replicat].unique())))
            cv.append(np.mean(run_cv))
plot_df = pd.DataFrame(df_list)
#%%
if len(fc_to_plot) == 3:
    model_colors = ["gainsboro", "grey", "black"]
    rc_ticks["figure.figsize"] = (1.7 * 0.57, 1.2)
else:
    base_color = "black"
    base_alphas = np.logspace(np.log10(0.15),0,len(fc_to_plot))
    from matplotlib.colors import to_rgba
    model_colors = [to_rgba(base_color, alpha) for alpha in base_alphas]
    rc_ticks["figure.figsize"] = (1.7, 1.2)

sns.set_theme(context = "talk", style = "ticks", rc = rc_ticks)
fig,ax = plt.subplots()

scaled_df = deepcopy(plot_df)
xs = np.array(np.arange(len(scaled_df.fold_change.unique())))
padding = np.max(xs)
xticks = []
for sv, scan in enumerate(scaled_df.scan_variable.unique()):
    scan_df = scaled_df.loc[(scaled_df.scan_variable == scan)]
    ys = scan_df.groupby("fold_change").mean()["std"].values
    yerr = scan_df.groupby("fold_change").std()["std"].values
    plt.bar(xs, ys, yerr=yerr, linewidth=0.2, width=0.8, color=model_colors, edgecolor="black", error_kw=dict(lw=0.3))
    xticks.append(np.mean(xs))
    xs += (padding + 2)

ax.set(ylabel=r"surface conc. s.d. (pM)", xlabel="", yticks=[0,7,14], ylim=(0,14))
plt.xticks(xticks, [r"f$_{\rm sec}$: secreting cells", r"$\sigma$: receptor heterog.", r"K$_{\rm D}$: saturation const."], rotation=-90, ha="center")
fig.savefig(saving_string, bbox_inches='tight', transparent=True)
plt.tight_layout()
plt.show()
