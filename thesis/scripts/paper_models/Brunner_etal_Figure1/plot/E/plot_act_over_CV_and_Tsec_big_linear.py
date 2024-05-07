import matplotlib.pyplot as plt
import numpy as np
import seaborn as sns
import os
import getpass
from thesis.scripts.paper_models.utilities.plotting_rc import rc_ticks
import pandas as pd

from thesis.scripts.paper_models.utilities.plot_helper import my_load_df, my_interpolation

hdd = "extra2" if os.path.exists("/extra2") else "extra"
user = getpass.getuser()

saving_string = r"/home/brunner/Documents/Current work/2024_03_15/"
if not os.path.exists(saving_string):
    os.mkdir(saving_string)

offset = 0
replicat = "replicat_index"

base_path = "/extra2/brunner/paper_models/kinetics/Figure_1C/dataframes_act_over_Tsec_large_linear_steady_state/"

spatial_cell_df, spatial_global_df =  my_load_df(base_path, offset=offset, custom_ending = "_combined")
try:
    print(np.sort(spatial_cell_df["IL-2_Tsec_fraction"].unique()))
except KeyError:
    spatial_cell_df["IL-2_Tsec_fraction"] = spatial_cell_df["fractions_Tsec"]
    print(np.sort(spatial_cell_df["IL-2_Tsec_fraction"].unique()))
try:
    print(np.sort(spatial_cell_df["IL-2_gamma"].unique()))
except KeyError:
    spatial_cell_df["IL-2_gamma"] = spatial_cell_df["misc_gamma"]
    print(np.sort(spatial_cell_df["IL-2_gamma"].unique()))
if len(spatial_cell_df["IL-2_gamma"].unique()) > 1:
    spatial_cell_df = spatial_cell_df.loc[spatial_cell_df["IL-2_gamma"] == 12.]

spatial_cell_df["IL-2_surf_c"] *= 1e3
########################################################################################################################
print("plotting")
########################################################################################################################
#%%
def EC50_calculation(E_max, E_min, k, N, R):
    return (E_max * k ** N + E_min * R ** N) / (k ** N + R ** N)

spatial_cv_act_list = []
ODE_cv_act_list = []
results = [spatial_cv_act_list, ODE_cv_act_list]
time_index = 0
for c, c_df in enumerate([spatial_cell_df]):
    c_df = c_df.loc[(c_df["time_index"] == c_df["time_index"].unique()[time_index])]
    for frac in np.sort(c_df["IL-2_Tsec_fraction"].unique()):
        frac_df = c_df.loc[(c_df["IL-2_Tsec_fraction"] == frac)]
        frac_df = frac_df.loc[frac_df["type_name"] == "Th"]
        for idx, rep in enumerate(np.sort(frac_df[replicat].unique())):
            rep_df = frac_df.loc[(frac_df[replicat] == rep)]
            try:
                assert False
                rep_df["pSTAT5"] = rep_df["IL-2_pSTAT5"]
            except:
                rep_df["pSTAT5"] = rep_df["IL-2_surf_c"] ** 3 / (
                        (EC50_calculation(E_max=125e-12, E_min=0, k=860, N=1.5, R=rep_df["IL-2_R"]) * 1e12) ** 3 +
                        rep_df["IL-2_surf_c"] ** 3).values
            act = (len(rep_df.loc[(rep_df["pSTAT5"] >= 0.5)].values)/len(rep_df.values)) * 100
            surf_c = rep_df["IL-2_surf_c"].mean()
            cv = rep_df["IL-2_surf_c"].std()/rep_df["IL-2_surf_c"].mean()
            std = rep_df["IL-2_surf_c"].std()
            if act > 0:
                results[c].append({"IL-2_Tsec_fraction" : frac,
                                 replicat: rep,
                                 "activation": act,
                                 "CV": cv,
                                 "std": std,
                                 "surface_c": surf_c})
#%%
rc_ticks['figure.figsize'] = (1.67475 * 0.85, 1.386)
sns.set_theme(context="talk", style="ticks", rc=rc_ticks)
fig, ax = plt.subplots()

for entry in results[:1]:
    df = pd.DataFrame(entry)
    # color = "#e4c137fd"
    color = "#aa0000ff"
    from matplotlib.colors import to_rgba
    from copy import deepcopy
    rgba_color = list(to_rgba(color))
    palette = []
    alphas = np.logspace(-2.5,0, len(df["IL-2_Tsec_fraction"].unique()))
    for e,entry in enumerate(np.sort(df["IL-2_Tsec_fraction"].unique())):
        alphaed_color = deepcopy(rgba_color)
        alphaed_color[-1] = alphas[e]
        palette.append(alphaed_color)
    ax = sns.scatterplot(data = df, x = "std", y = "activation", hue = "IL-2_Tsec_fraction", legend=False, palette=palette, linewidth=0.05)
plt.xlabel("surface conc. s.d. (pM)")
plt.xscale("linear")
plt.yscale("linear")
plt.xlim(0, 10)
plt.ylabel(r"pSTAT$^+$ T$_{\rm resp}$ cells (%)")
plt.xticks([0,5,10])
plt.yticks([0,50,100])
plt.title("linear")

fig.savefig(saving_string + f"Fig1_activation_over_CV_local_q_full_linear.pdf", bbox_inches='tight', transparent=True)
fig.tight_layout()
plt.show()
# #%%
# rc_ticks['figure.figsize'] = (1.67475 * 1.32, 0.81)
# sns.set_theme(context="talk", style="ticks", rc=rc_ticks)
# fig, ax = plt.subplots()
#
# for entry in results[:1]:
#     df = pd.DataFrame(entry)
#     color = "#aa0000ff"
#     from matplotlib.colors import to_rgba
#     from copy import deepcopy
#     rgba_color = list(to_rgba(color))
#     palette = []
#     alphas = np.logspace(-3,0, len(df["IL-2_Tsec_fraction"].unique()))
#     for e,entry in enumerate(np.sort(df["IL-2_Tsec_fraction"].unique())):
#         alphaed_color = deepcopy(rgba_color)
#         alphaed_color[-1] = alphas[e]
#         palette.append(alphaed_color)
#     ax = sns.scatterplot(data = df, x = "std", y = "surface_c", hue = "IL-2_Tsec_fraction", legend=False, palette=palette, linewidth=0.05)
# plt.ylabel(r"conc. avg." "\n" r"(pM)")
# plt.xlabel("surface conc. s.d. (pM)")
# plt.xscale("linear")
# plt.yscale("linear")
# plt.xlim((0,15))
# plt.xticks([0,5,10,15])
# plt.ylim(0, 30)
# plt.yticks([0, 15, 30])
#
# fig.savefig(saving_string + f"Fig1_surf_c_over_CV_local_q_full.pdf", bbox_inches='tight', transparent=True)
# fig.tight_layout()
# plt.show()

#%%
# from scipy.stats import spearmanr
# x = df["std"].values
# y = df.surface_c.values
# res = spearmanr(x,y)
# print("surface spearman:", res)

x = df["std"].values
y = df.activation.values
res = spearmanr(x,y)
print("activation spearman:", res)
