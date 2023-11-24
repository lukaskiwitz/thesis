import matplotlib.pyplot as plt
import numpy as np
import seaborn as sns
import os
from thesis.scripts.paper_models.utilities.plotting_rc import rc_ticks
from thesis.scripts.paper_models.utilities.plot_helper import my_load_df, my_interpolation
# bc = "standard"
bc = "saturated"

saving_string = r"/home/brunner/Documents/Current work/2023_11_03/boxplot_static_surf_c.pdf"
offset = 0
replicat = "replicat_index"

print("loading data")
hdd = "/extra2" if os.path.exists("/extra2") else "/extra"
path = '/extra2/brunner/paper_models/ODE/saturated/activation_q_10_R_1e4_Tsec_scan_5/'
ODE_c_df, ODE_g_df = my_load_df(path, offset = 0, custom_ending = "_combined")
ODE_c_df["IL-2_surf_c"] *= 1e3

import getpass
hdd = "extra2" if os.path.exists("/extra2") else "extra"
user = getpass.getuser()

base_path = "/extra2/brunner/paper_models/kinetics/Figure_1C/dataframes_act_over_Tsec_large_2_steady_state/"

big_c_df, big_g_df = my_load_df(base_path, offset=offset, custom_ending="")

try:
    print(np.sort(big_c_df["IL-2_Tsec_fraction"].unique()))
except KeyError:
    big_c_df["IL-2_Tsec_fraction"] = big_c_df["fractions_Tsec"]
    print(np.sort(big_c_df["IL-2_Tsec_fraction"].unique()))
try:
    print(np.sort(big_c_df["IL-2_gamma"].unique()))
except KeyError:
    big_c_df["IL-2_gamma"] = big_c_df["misc_gamma"]
    print(np.sort(big_c_df["IL-2_gamma"].unique()))
if len(big_c_df["IL-2_gamma"].unique()) > 1:
    big_c_df = big_c_df.loc[big_c_df["IL-2_gamma"] == 12.]
########################################################################################################################
print("plotting")
########################################################################################################################
#%%
def EC50_calculation(E_max, E_min, k, N, R):
    return (E_max * k ** N + E_min * R ** N) / (k ** N + R ** N)

# fracs = [0.051, 0.076, 0.126, 0.2, 0.3, 0.4]
fracs = ODE_c_df["IL-2_Tsec_fraction"].unique()
ODE_boxplot_data = []
for frac in fracs:
    ODE_frac_df = ODE_c_df.loc[(ODE_c_df["IL-2_Tsec_fraction"] == frac)]
    ODE_frac_df = ODE_frac_df.loc[ODE_frac_df["type_name"] == "Th"]

    act_list = []
    R_list = []
    for idx, rep in enumerate(np.sort(ODE_frac_df[replicat].unique())):
        tmp_df = ODE_frac_df.loc[(ODE_frac_df[replicat] == rep)]
        frac_act = tmp_df["IL-2_surf_c"].mean()
        act_list.append(frac_act)
    if len(act_list) == 0:
        print(frac)
    ODE_boxplot_data.append(np.array(act_list))

#%%
boxplot_data = []
# fracs = np.sort(big_c_df["IL-2_Tsec_fraction"].unique())
fracs = [0.02, 0.0295, 0.039, 0.058, 0.077, 0.096, 0.115, 0.1435, 0.2005, 0.248, 0.305, 0.4]

for frac in fracs:
    frac_df = big_c_df.loc[(np.abs(big_c_df["IL-2_Tsec_fraction"] - frac) < 1e-3) & (big_c_df["time_index"] == 0.)]
    print(frac_df.loc[frac_df.type_name == "Th", "IL-2_R"].mean())
    act_list_2 = []
    for rep in frac_df["replicat_index"].unique():
        rep_df = frac_df.loc[(frac_df["replicat_index"] == rep) & (frac_df.type_name == "Th")]
        rep_df["IL-2_surf_c"] *= 1e3
        act_list_2.append(rep_df["IL-2_surf_c"].mean())
    boxplot_data.append(np.array(act_list_2))
#%%
np.random.seed(1)
rc_ticks['figure.figsize'] = (1.67475 * 1.32, 1.386 * 1.27)
sns.set_theme(context = "talk", style = "ticks", rc = rc_ticks)
fig,ax = plt.subplots()

linewidth = 0.5
alpha = 0.9

color1 = "#aa0000ff"
color2 = "grey"
color3 = "blue"

flierprops = dict(marker='o', markerfacecolor=color1, markersize=0,
                  linestyle='none', markeredgecolor=color1, linewidth=linewidth)
boxprops=dict(linewidth=linewidth)
whiskerprops=dict(linewidth=linewidth)
capsprops=dict(linewidth=linewidth)

two_colours = [[color1] * len(fracs), [color2] * len(fracs), [color3] * len(fracs)]
from matplotlib.colors import to_rgba
rgba_1 = list(to_rgba(color1))
rgba_1[-1] = 0.6
rgba_2 = list(to_rgba(color2))
rgba_2[-1] = 0.6
rgba_3 = list(to_rgba(color3))
rgba_3[-1] = 0.6
rgb_colors = [[rgba_1] * len(fracs),
              [rgba_2] * len(fracs),
              [rgba_3] * len(fracs)]
four_colours = [[color1] * len(fracs) * 2, [color2] * len(fracs) * 2, [color3] * len(fracs) * 2]

for b, bdata in enumerate([boxplot_data]):
    x = [5, 10, 15, 20, 30, 40]
    x = np.array(fracs) * 100
    alphas = np.linspace(1,1, len(bdata))
    whisker_alphas = np.array([[a] * 2 for a in alphas]).flatten()
    w = 0.04 * 3
    width = lambda p, w: 10 ** (np.log10(p) + w / 2.) - 10 ** (np.log10(p) - w / 2.)
    boxp = plt.boxplot(bdata, positions=x, flierprops=flierprops, boxprops=boxprops, whiskerprops=whiskerprops, widths=width(x, w), patch_artist=True)

    for element in ["fliers", "means", "medians"]:
        for e,entry in enumerate(boxp[element]):
            plt.setp(entry, color=two_colours[b][e], alpha=alphas[e]/2)

    for element in ['boxes']:
        for e,entry in enumerate(boxp[element]):
            plt.setp(entry, facecolor=rgb_colors[b][e], edgecolor=two_colours[b][e], alpha=alphas[e]/2)

    for element in ['whiskers', "caps"]:
        for e,entry in enumerate(boxp[element]):
            plt.setp(entry, color=four_colours[b][e], linewidth=linewidth, alpha=np.clip(whisker_alphas[e] * 2,0,1))
    for e,entry in enumerate(bdata):
        factor = 1/1
        subsample = np.random.choice(entry, int(len(entry) * factor), replace=False)
        plt.scatter(x[e] + np.random.normal(0, 1, size=len(subsample))/1, subsample, color = two_colours[b][e],
                 alpha=alphas[e], s=0.75, linewidth=0)

plt.plot(ODE_c_df["IL-2_Tsec_fraction"].unique() * 100, np.mean(ODE_boxplot_data, axis=1), color=color2, linewidth=0.9)

plt.ylabel(r"surface conc. avg. (pM)")
plt.xscale("log")
plt.yscale("linear")
plt.ylim(0, 25)

plt.xlabel(r"secreting cells (%)")
plt.yticks([0,12.5,25])
plt.xlim((1.5, 49))
plt.xticks([2, 10, 40], [2, 10, 40])

fig.savefig(saving_string, bbox_inches='tight', transparent=True)
plt.tight_layout()
plt.show()