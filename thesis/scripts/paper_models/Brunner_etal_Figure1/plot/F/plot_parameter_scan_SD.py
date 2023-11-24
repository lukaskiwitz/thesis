import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import seaborn as sns
import os
from thesis.scripts.paper_models.utilities.plot_helper import activation_to_global_df
from thesis.scripts.paper_models.utilities.plotting_rc import rc_ticks

def EC50_calculation(E_max, E_min, k, N, R):
    return (E_max * k ** N + E_min * R ** N) / (k ** N + R ** N)

hdd = "/extra2" if os.path.exists("/extra2") else "/extra"
path = "/extra2/brunner/paper_models/boxed_static/parameter_scan_test/"

fig_path = "/home/brunner/Documents/Current work/2023_11_03/"  + "Fig1C_short_"

global_df = pd.read_hdf(path + 'global_df.h5', mode="r")
cell_df = pd.read_hdf(path + 'cell_df.h5', mode="r")

global_df = global_df.loc[(global_df["model_name"] == "pde_model")]
cell_df = cell_df.loc[(cell_df["model_name"] == "pde_model")]
cell_df["IL-2_surf_c"] *= 1e3

global_df = activation_to_global_df(cell_df, global_df)
#%%
print("plotting")
linear_uptake = False

scan_names = ["R", "sigma", "q", "T_sec", "KD"]
parameter_names = ["scan_index", "IL-2_sigma", "IL-2_q", "IL-2_Tsec_fraction", "IL-2_KD"]
labels = [r"R: receptors", r"$\sigma$: receptor heterog.", r" q: secretion rate", r"f$_{\rm sec}$: secreting cells", r"K$_{\rm D}$: saturation const."]
standards = [1.5e3, 1, 10., 0.1, 0.007437]

no_of_scan_points = 9
scan_values = np.logspace(-1,1,no_of_scan_points)
plotted_sv = scan_values[::4]

ylim = (5e-2, 6e1)

scan_measure = "SD"
cell_df["pSTAT5"] = cell_df["IL-2_surf_c"] ** 3 / (
        (EC50_calculation(E_max=125e-12, E_min=0, k=860, N=1.5, R=cell_df["IL-2_R"]) * 1e12) ** 3 + cell_df[
    "IL-2_surf_c"] ** 3).values

#%%
for s, sn in enumerate(scan_names):
    if sn == "kd":
        kd_values = scan_values
        for i, si in enumerate(cell_df.loc[(cell_df["scan_name_scan_name"] == sn), "scan_index"].unique()):
            cell_df.loc[(cell_df["scan_name_scan_name"] == sn) & (cell_df["scan_index"] == si), "scan_value"] = kd_values[i]
    if sn == "D":
        D_values = scan_values
        for i, si in enumerate(cell_df.loc[(cell_df["scan_name_scan_name"] == sn), "scan_index"].unique()):
            cell_df.loc[(cell_df["scan_name_scan_name"] == sn) & (cell_df["scan_index"] == si), "scan_value"] = D_values[i]
    elif sn in ["R", "sigma", "KD"]: #only on Ths
        for pn, para_names in enumerate(
                cell_df.loc[(cell_df["scan_name_scan_name"] == sn) & (cell_df.type_name == "Th"), parameter_names[s]].unique()):
            scan_index = cell_df.loc[(cell_df["scan_name_scan_name"] == sn) & (
                        np.abs(cell_df[parameter_names[s]] - para_names) < 1e-10), "scan_index"].unique()[0]
            cell_df.loc[(cell_df["scan_name_scan_name"] == sn) & (cell_df.scan_index == scan_index), "scan_value"] = scan_values[pn]
    elif sn in ["q"]: #only on Tsecs
        for pn, para_names in enumerate(
                cell_df.loc[(cell_df["scan_name_scan_name"] == sn) & (cell_df.type_name == "Tsec"), parameter_names[s]].unique()):
            scan_index = cell_df.loc[(cell_df["scan_name_scan_name"] == sn) & (
                        np.abs(cell_df[parameter_names[s]] - para_names) < 1e-10), "scan_index"].unique()[0]
            cell_df.loc[(cell_df["scan_name_scan_name"] == sn) & (cell_df.scan_index == scan_index), "scan_value"] = scan_values[pn]
    else:
        for pn, para_names in enumerate(cell_df.loc[(cell_df["scan_name_scan_name"] == sn), parameter_names[s]].unique()):
            cell_df.loc[(cell_df["scan_name_scan_name"] == sn) & (np.abs(cell_df[parameter_names[s]] - para_names) < 1e-10), "scan_value"] = scan_values[pn]
#%%
results = []
error = []
for n,name in enumerate(scan_names):
    name_df = cell_df.loc[(cell_df["scan_name_scan_name"] == name)]
    for value in plotted_sv:
        rep_standards = []
        rep_as = []
        for rep in name_df["replicat_index"].unique():
            if scan_measure == "pSTAT5":
                try:
                    standard = len(name_df.loc[(name_df["replicat_index"] == rep) & (name_df["scan_value"] == 1) & (
                                name_df["pSTAT5"] > 0.5)]) / \
                               len(name_df.loc[(name_df["replicat_index"] == rep) & (name_df["scan_value"] == 1)])
                    a_mean = len(name_df.loc[(name_df["replicat_index"] == rep) & (name_df["scan_value"] == value) & (
                            name_df["pSTAT5"] > 0.5)]) / \
                             len(name_df.loc[(name_df["replicat_index"] == rep) & (name_df["scan_value"] == value)])

                    rep_standards.append(standard)
                    rep_as.append(a_mean)
                except ZeroDivisionError:
                    pass
            elif scan_measure == "SD":
                standard = name_df.loc[(name_df["replicat_index"] == rep) & (name_df["scan_value"] == 1), "IL-2_surf_c"].std()
                rep_standards.append(standard)
                a_mean = name_df.loc[(name_df["replicat_index"] == rep) & (name_df["scan_value"] == value), "IL-2_surf_c"].std()
                rep_as.append(a_mean)
            elif scan_measure == "CV":
                standard = name_df.loc[(name_df["replicat_index"] == rep) & (np.abs(name_df["scan_value"] - 1) < 1e-12), "IL-2_surf_c"].std() / \
                           name_df.loc[(name_df["replicat_index"] == rep) & (np.abs(name_df["scan_value"] - 1) < 1e-12), "IL-2_surf_c"].mean()
                rep_standards.append(standard)
                a_mean = name_df.loc[(name_df["replicat_index"] == rep) & (np.abs(name_df["scan_value"] - value) < 1e-12), "IL-2_surf_c"].std() / \
                         name_df.loc[(name_df["replicat_index"] == rep) & (np.abs(name_df["scan_value"] - value) < 1e-12), "IL-2_surf_c"].mean()
                rep_as.append(a_mean)
            elif scan_measure == "Gradient":
                name_df = global_df.loc[(global_df["scan_name_scan_name"] == name)]
                name_df["Gradient"] *= 1e3
                standard = name_df.loc[(name_df["replicat_index"] == rep) & (name_df["scan_value"] == 1), scan_measure].mean()
                rep_standards.append(standard)
                a_mean = name_df.loc[(name_df["replicat_index"] == rep) & (name_df["scan_value"] == value), scan_measure].mean()
                rep_as.append(a_mean)
            else:
                standard = name_df.loc[
                    (name_df["replicat_index"] == rep) & (name_df["scan_value"] == 1), scan_measure].mean()
                rep_standards.append(standard)
                a_mean = name_df.loc[
                    (name_df["replicat_index"] == rep) & (name_df["scan_value"] == value), scan_measure].mean()
                rep_as.append(a_mean)
        results.append(np.nanmean(np.array(rep_as)))# / np.array(rep_standards)))
        error.append(np.nanstd(np.array(rep_as)))# / np.array(rep_standards)))
#%%
N = len(plotted_sv)

x = np.array([[N*(x + 0) + b for b in np.arange(N)] for x in range(len(scan_names))]).flatten()
np.array([[N*x + b for b in np.arange(N)] for x in range(len(scan_names))]).flatten()
x_list = []
x_ticks = []
for x in range(len(scan_names)):
    x_list.append(np.arange(x*1.25*N, x*1.25*N+N))
    x_ticks.append(np.mean(x_list[-1]))
bars = np.array(results).reshape(-1,N)
x_reshaped = np.array(x_list)
error_reshaped = np.array(error).reshape(-1,N)
#%%
rc_ticks['figure.figsize'] = (1.67475 * 1.32, 1.2 * 1.1) #(1.7, 1.2)
sns.set_theme(context = "talk", style = "ticks", rc = rc_ticks)
fig, ax = plt.subplots()
width = 0.75 if len(plotted_sv) > 3 else 1
alphas = np.logspace(-0.8,0, N)
if scan_measure == "surf. c.":
    for e,entry in enumerate(bars):
        for v, value in enumerate(entry):
            ax.bar(x=x_reshaped[e,v], height=value, color="black", alpha=alphas[v], yerr=error_reshaped[e,v], width=width, bottom=1,
                   linewidth=0)
else:
    for e,entry in enumerate(bars):
        for v, value in enumerate(entry):
            ax.bar(x=x_reshaped[e,v], height=value if value > 4e-2 else 4e-2, color="black", alpha=alphas[v], yerr=error_reshaped[e,v], width=width, bottom=0,
                   linewidth=0.3, error_kw=dict(lw=width/2, capsize=0, capthick=0), edgecolor="black")
ax.set_xticks(x_ticks)
ax.set_xticklabels(labels)
plt.xticks(rotation=-90, ha='center')
if scan_measure == "CV":
    ax.set_ylabel("CV")
    ax.set_ylim((0, 4))
elif scan_measure == "SD":
    ax.set_ylabel("surface conc. s.d. (pM)")
    ax.set_yscale("log")
if scan_measure == "IL-2_surf_c":
    ax.set_ylabel("surface conc. (pM)")
    ax.set_ylim((1e-1, 5e4))
    ax.set_yscale("log", base=10)
    import matplotlib as matplotlib
    locmin = matplotlib.ticker.LogLocator(base=10.0,subs=(0.2,0.4,0.6,0.8),numticks=6)
    locmaj = matplotlib.ticker.LogLocator(base=10.0,numticks=12)
    ax.yaxis.set_minor_locator(locmin)
    ax.yaxis.set_minor_formatter(matplotlib.ticker.NullFormatter())
    ax.yaxis.set_major_locator(locmaj)
    ax.yaxis.set_major_formatter(matplotlib.ticker.NullFormatter())

    plt.yticks([1e-1, 1e0, 1e1, 1e2, 1e3, 1e4], [r"10$^{-1}$", "10$^0$", r"10$^1$", r"10$^2$", r"10$^3$", r"10$^4$"])

fig.savefig(fig_path + scan_measure + ".pdf", bbox_inches='tight', transparent=True)
plt.tight_layout()
plt.show()