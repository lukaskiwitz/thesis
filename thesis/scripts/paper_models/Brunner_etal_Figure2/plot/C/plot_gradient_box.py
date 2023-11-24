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

saving_string = r"/home/brunner/Documents/Current work/2023_11_03/" + "gradient_{m}_static_{sv}_{bc}".format(m = "box", sv=scan_variables[0], bc = bc) + ".pdf"

sv_styles = ['-','-.',':']
plot_names = ["sec. scan", "R scan", "KD scan"]
# define the standard values (fold-change = 0) for each of the scan_variables
standards = [0.05, 1, 7.4]

# which models to plot

# models = ["loc. q and R", "loc. q", "well_mixed"]
models = ["loc. q and R"]
# defines a outer layer of N cells to be ignored in plotting. Used to further limit unwanted boundary effects.


replicat = "replicat_index"

# load runs, apply, merge dataframes
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
            g_df["IL-2_Tsec_fraction"] = g_df["fractions_Tsec"]
        if "ODE" not in m_sv_path:
            try:
                c_df["IL-2_sigma"]
            except KeyError:
                c_df["IL-2_sigma"] = c_df["misc_sigma"]
                try:
                    g_df["IL-2_sigma"] = g_df["scan_value"]
                except:
                    pass
            if model == "k_on_linear":
                c_df["IL-2_KD"] = c_df["IL-2_k_off"].unique()[0] / c_df["IL-2_k_on"] * 1e3
                g_df["IL-2_KD"] = c_df["IL-2_k_off"].unique()[0] / g_df["IL-2_k_on"] * 1e3
        # scales the dataframes. Depends on your desired units. Currently everything is in pM
        c_df, g_df = scale_dataframes(c_df, g_df, model, scan_variable, sv, standards)
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
    for sv, scan_variable in enumerate(scan_variables):
        c_df = dataframes[m][sv][0]
        g_df = dataframes[m][sv][1]
        ste = []
        grad = []
        for f, fc in enumerate(fc_to_plot):
            scan = g_df[scan_variable].unique()[(np.abs(g_df[scan_variable].unique() - (standards[sv] * fc))).argmin()]
            sliced_df = g_df.loc[(g_df[scan_variable] == scan)]
            run_grad = []
            mean = []
            try:
                for r, run in enumerate(np.sort(sliced_df[replicat].unique())):
                    run_grad.append(sliced_df.loc[(sliced_df[replicat] == run), "Gradient"].values) #/sliced_df.loc[(sliced_df[replicat] == run), "Concentration"].values)
                    mean.append(sliced_df.loc[(sliced_df[replicat] == run), "surf_c"].mean())
                    df_list.append({"model": model, "scan_variable": scan_variable, "fold_change": fc, "replicat": run,
                                    "gradient": run_grad[-1][0] * 1e3, "gradient_norm": run_grad[-1][0]/mean[-1] * 1e3})
            except KeyError as e: # ODE system does not have a gradient
                print("KeyError:", e)
                run_grad = [0]

            ste.append(np.std(run_grad) / np.sqrt(len(g_df[replicat].unique())))
            grad.append(np.mean(run_grad))
plot_df = pd.DataFrame(df_list)
assert len(plot_df) > 0
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
    ys = scan_df.groupby("fold_change").mean()["gradient"].values
    yerr = scan_df.groupby("fold_change").std()["gradient"].values
    plt.bar(xs, ys, yerr=yerr, linewidth=0.2, width=0.8, color=model_colors, edgecolor="black", error_kw=dict(lw=0.3))
    xticks.append(np.mean(xs))
    xs += (padding + 2)

ylim = (None, None)
if len(fc_to_plot) == 3:
    ylim = (0.1, None)
    ax.set(yticks=[0.1,0.3,0.5])
ax.set(ylabel=r"gradient (pM/µm)", xlabel="", ylim=ylim)
plt.xticks(xticks, [r"f$_{\rm sec}$: secreting cells", r"$\sigma$: receptor heterog.", r"K$_{\rm D}$: saturation const."], rotation=-90, ha="center")
fig.savefig(saving_string, bbox_inches='tight', transparent=True)
# plt.tight_layout()
plt.show()
