import matplotlib.pyplot as plt
# plt.rcParams['text.usetex'] = False
import numpy as np
import pandas as pd
import seaborn as sns
from plotting_rc import rc_ticks

sns.set_theme(context = "talk", style = "ticks", rc = rc_ticks)
fig, ax = plt.subplots()
# fig = plt.figure()
plt.subplots_adjust(wspace=.3)


def EC50_calculation(E_max, E_min, k, N, R):
    return (E_max * k ** N + E_min * R ** N) / (k ** N + R ** N)

path = "/extra/brunner/thesis/static/lukas/receptor_distribution/"

fig_path = "/home/brunner/Documents/Current work/2022_03_04/"

global_df = pd.read_hdf(path + 'global_df.h5', mode="r")
cell_df = pd.read_hdf(path + 'cell_df.h5', mode="r")


cell_df["IL-2_surf_c"] = cell_df["IL-2_surf_c"] * 1e3
global_df["surf_c"] *= 1e3
global_df["surf_c_std"] *= 1e3
global_df["cv"] = global_df["surf_c_std"] / global_df["surf_c"]

try:
    cell_df["pSTAT5"] = cell_df["misc_pSTAT5"]
except:
    print("pSTAT5 values not available, calculating")
    cell_df["pSTAT5"] = cell_df["IL-2_surf_c"] ** 3 / (
        (EC50_calculation(E_max=125e-12, E_min=0, k=860, N=0.55, R=cell_df["IL-2_R"]) * 1e12) ** 3 + cell_df[
                "IL-2_surf_c"] ** 3).values

frac_of_activated_cells = []
vs = cell_df["v_v"].unique()
cell_df = cell_df.loc[cell_df["type_name"] != "sec"]
SD = []
mean = []
for v in vs:
    activated_cells = len(cell_df.loc[(cell_df["v_v"] == v) & (cell_df["pSTAT5"] >= 0.5), "id_id"].unique())
    SD.append(cell_df.loc[(cell_df["v_v"] == v), "IL-2_surf_c"].std())
    mean.append(cell_df.loc[(cell_df["v_v"] == v), "IL-2_surf_c"].mean())
    frac_of_activated_cells.append(activated_cells/len(cell_df["id_id"].unique()))

vs *= 100

colour1 = "black"
ax = sns.lineplot(vs, reversed(global_df.sort_values(by = "v_v")["surf_c_std"].values), color=colour1)

ax.tick_params(axis='y', colors=colour1)
ax.yaxis.label.set_color(colour1)
ax.set(xlabel="% receptors on T$_{resp}$", ylabel="s.d. on T$_{resp}$ (pM)")


colour2 = "orange"

ax2 = ax.twinx()
ax2 = sns.lineplot(vs, reversed(np.array(frac_of_activated_cells) * 100), color=colour2)
ax2.set(ylabel="activated T$_{resp}$ (%)")
ax2.tick_params(axis='y', colors=colour2)
ax2.yaxis.label.set_color(colour2)

fig.savefig(fig_path + "Fig1_D_receptor_distribution_cv.pdf", bbox_inches='tight')

plt.tight_layout()
plt.show()


