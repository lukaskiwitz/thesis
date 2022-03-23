import matplotlib.pyplot as plt
import numpy as np
import seaborn as sns
import os
from scipy.interpolate import UnivariateSpline
from plotting_rc import rc_ticks

sns.set_theme(context = "talk", style = "ticks", rc = rc_ticks)
fig,ax = plt.subplots()

from thesis.scripts.paper_models.utilities.plot_helper import my_load_df, my_interpolation

'''
define which static sim you want to plot. 
bc defines the boundary condition, scan variable the scan.
"IL-2_Tsec_fraction" is the secreting cells scan, "IL-2_sigma" the R_lognorm scan and IL-2_KD the KD scan.
'''
bc = "standard"
# bc = "saturated"
# scan_variable = "IL-2_Tsec_fraction"
# scan_variable = "IL-2_sigma"
scan_variable = "IL-2_KD"

saving_string =r"/home/brunner/Documents/Current work/2022_03_04/" + "{m}_static_{sv}_{bc}".format(m = "yukawa", sv=scan_variable, bc = bc) + ".svg"

# plotting parameters
plot_legend = True

yscale = "linear"
xscale = "linear"
# define data start and end point
startingPoint = None
stoppingPoint = None
# maximum value
max_value = 1e5
min_value = -1
# x and ylims have to be adjusted for each plot
xlim = (0.04, 0.5)
ylim = (-1, 45)

# defines a outer layer of N cells to be ignored in plotting. Used to further limit unwanted boundary effects.
offset = 0
# which cell type to plot
cell_types = ["Th"]
# Which colors to use, dim-1 first entry is spatial, second is ODE. dim-2 first entry is concentration, second SD
y_colours = [["blue", "orange"] , ["black"]]
legend_labels = ["RD-system", "well-mixed", "s.d."]
linestyles = [["-", "-"], ["-"]]

y_variable = "surf_c"
title = "Surface concentration"

replicat = "replicat_index"

# if line smoothing is desired
c_interpolate = False
c_sf = 0.009
c_ste_sf = 0.009
std_interpolate = False
std_sf = 0.05
std_ste_sf = 0.05

########################################################################################################################
########################################################################################################################
##                                              Plotting                                                              ##
########################################################################################################################
########################################################################################################################
# enter the paths for all the scans here
hdd = "/extra2" if os.path.exists("/extra2") else "/extra"
if scan_variable == "IL-2_Tsec_fraction":
    # spatial_path = hdd + "/brunner/paper_models/statics/{bc}/Tsec_scan_no_sigma_low_R_log2/".format(bc=bc)
    # spatial_path = hdd + "/brunner/paper_models/statics/{bc}/Tsec_scan_close_to_ODE/".format(bc=bc)
    spatial_path = hdd + "/brunner/paper_models/yukawa/{bc}/Tsec_scan_low_c0_log2/".format(bc=bc)

    ODE_path = hdd + "/brunner/paper_models/ODE/{bc}/Tsec_scan_log2/".format(bc=bc)

elif scan_variable == "IL-2_sigma":
    # spatial_path = hdd + "/brunner/paper_models/statics/{bc}/R_lognorm_low_R_log2/".format(bc=bc)
    spatial_path = hdd + "/brunner/paper_models/yukawa/{bc}/R_lognorm_low_c0_log2/".format(bc=bc)
    ODE_path = hdd + "/brunner/paper_models/ODE/{bc}/R_lognorm_log2/".format(bc=bc)

elif scan_variable == "IL-2_KD":
    # spatial_path = hdd + "/brunner/paper_models/statics/saturated/KD_scan_low_R_log2/"
    spatial_path = hdd + "/brunner/paper_models/yukawa/saturated/KD_scan_low_c0_log2/"
    ODE_path = hdd + "/brunner/paper_models/ODE/saturated/KD_scan_log2/"

print("loading data from", spatial_path)
# spatial_cell_df, spatial_global_df = my_load_df(spatial_path, offset = 0, run_range = myRange, custom_ending = "_combined")
spatial_cell_df, spatial_global_df = my_load_df(spatial_path, offset = offset)

print("loading data from", ODE_path)
ODE_cell_df, ODE_global_df = my_load_df(ODE_path, offset = offset)

# adjust units
spatial_global_df["surf_c"] *= 1e3
spatial_global_df["surf_c_std"] *= 1e3

ODE_global_df["surf_c"] *= 1e3
ODE_global_df["surf_c_std"] *= 1e3

spatial_global_df["time"] = spatial_global_df["time"].div(36000)
spatial_global_df["time"] = spatial_global_df["time"].div(3600/10)
spatial_global_df = spatial_global_df.sort_values(by=replicat)

ODE_global_df["time"] = ODE_global_df["time"].div(36000)
ODE_global_df["time"] = ODE_global_df["time"].div(3600/10)
ODE_global_df = ODE_global_df.sort_values(by=replicat)

spatial_cell_df["time"] = spatial_cell_df["time"].div(3600)
spatial_cell_df["IL-2_surf_c"] = spatial_cell_df["IL-2_surf_c"].mul(1e3)

ODE_cell_df["time"] = ODE_cell_df["time"].div(3600)
ODE_cell_df["IL-2_surf_c"] = ODE_cell_df["IL-2_surf_c"].mul(1e3)

if scan_variable == "IL-2_KD":
    from thesis.scripts.patrick.ODE.driver import ODEdriver
    from thesis.scripts.patrick.ODE.parameters import p
    p["N_cells"] = len(ODE_cell_df["id"].unique())
    ODEdr = ODEdriver(p)
    try:
        if spatial_global_df.scan_name_scan_name.unique()[0] == "KD":
            spatial_global_df["IL-2_KD"] = spatial_global_df.scan_value * 1e3
            spatial_cell_df["IL-2_KD"] = spatial_cell_df.scan_value * 1e3
    except AttributeError:
        if spatial_global_df.scan_name.unique()[0] == "KD":
            spatial_global_df["IL-2_KD"] = spatial_global_df.scan_value * 1e3
            spatial_cell_df["IL-2_KD"] = spatial_cell_df.scan_value * 1e3
        elif spatial_global_df.scan_name.unique()[0] == "IL-2_KD":
            from copy import deepcopy
            p_yuk = deepcopy(p)
            p_yuk["N_cells"] = len(spatial_cell_df["id_id"].unique())
            YUKdr = ODEdriver(p_yuk)
            spatial_global_df["IL-2_KD"] = YUKdr.molecules_to_molar(spatial_global_df["IL-2_KD"].values) * 1e12
            spatial_cell_df["IL-2_KD"] = YUKdr.molecules_to_molar(spatial_cell_df["IL-2_KD"].values) * 1e12

    spatial_cell_df["IL-2_KD"] = 1 / spatial_cell_df["IL-2_KD"]
    spatial_global_df["IL-2_KD"] = 1 / spatial_global_df["IL-2_KD"]
    ODE_global_df["IL-2_KD"] = 1/(ODEdr.molecules_to_molar(ODE_global_df["IL-2_KD"].values) * 1e12)
    ODE_cell_df["IL-2_KD"] = 1/(ODEdr.molecules_to_molar(ODE_cell_df["IL-2_KD"].values) * 1e12)


########################################################################################################################
print("plotting")
########################################################################################################################
#%%


ODE_cell_df = ODE_cell_df.loc[ODE_cell_df[scan_variable] < max_value]
ODE_global_df = ODE_global_df.loc[ODE_global_df[scan_variable] < max_value]

spatial_cell_df = spatial_cell_df.loc[spatial_cell_df[scan_variable] < max_value]
spatial_global_df = spatial_global_df.loc[spatial_global_df[scan_variable] < max_value]

ODE_cell_df = ODE_cell_df.loc[ODE_cell_df[scan_variable] > min_value]
ODE_global_df = ODE_global_df.loc[ODE_global_df[scan_variable] > min_value]

spatial_cell_df = spatial_cell_df.loc[spatial_cell_df[scan_variable] > min_value]
spatial_global_df = spatial_global_df.loc[spatial_global_df[scan_variable] > min_value]

for d, dfs in enumerate([(spatial_cell_df, spatial_global_df), (ODE_cell_df, ODE_global_df)]):
    cell_df, global_df = dfs
    for c, cell_type in enumerate(cell_types):
        tmp_df = cell_df.loc[cell_df["type_name"] == cell_type]
        surf_c = []
        for s, sv in enumerate(np.sort(tmp_df[scan_variable].unique())):
            surf_c.append(tmp_df.loc[tmp_df[scan_variable] == sv, "IL-2_surf_c"].mean())

        bars_df = cell_df.loc[cell_df["type_name"] == cell_type]
        ste = []
        mean = []
        for s, scan in enumerate(np.sort(bars_df[scan_variable].unique())):
            ste.append(bars_df.loc[np.abs(bars_df[scan_variable] - scan) < 0.0001, "IL-2_" + y_variable].std() / np.sqrt(len(cell_df[replicat].unique())))
            mean.append(bars_df.loc[np.abs(bars_df[scan_variable] - scan) < 0.0001,  "IL-2_" + y_variable].mean())

        if c_interpolate == True:
            my_interpolation(np.sort(tmp_df[scan_variable].unique()), surf_c, c_sf, plot=True, label="mean",
                             color=y_colours[d][0], fill=True, std=np.array(ste), ci=0, linestyle = linestyles[d][0], legend=plot_legend)
        else:
            sns.lineplot(x = np.sort(tmp_df[scan_variable].unique()), y = surf_c,
                     color=y_colours[d][0], linestyle = linestyles[d][0], ci=0, label="mean", legend=plot_legend)
            plt.fill_between(np.sort(bars_df[scan_variable].unique()),
                         np.clip(np.array(mean) - np.array(ste), 0, None),
                         np.clip(np.array(mean) + np.array(ste), 0, None), alpha=0.15, color=y_colours[d][0])


for d, dfs in enumerate([(spatial_cell_df, spatial_global_df)]):
    cell_df, global_df = dfs
    for c, cell_type in enumerate(cell_types):
        tmp_df = cell_df.loc[cell_df["type_name"] == cell_type]
        ste = []
        std = []
        for s, scan in enumerate(np.sort(tmp_df[scan_variable].unique())):
            sliced_df = tmp_df.loc[(tmp_df[scan_variable] == scan)]
            if len(sliced_df) != 0:
                run_std = []
                for r, run in enumerate(np.sort(sliced_df[replicat].unique())):
                    run_std.append(sliced_df.loc[(sliced_df[replicat] == run), "IL-2_surf_c"].std())
                run_std = np.array(run_std)
                run_std = run_std[~np.isnan(run_std)]
                ste.append(np.std(run_std) / np.sqrt(len(sliced_df[replicat].unique())))
                std.append(np.mean(run_std))

        if std_interpolate == True:
            std = np.array(std)
            std_spl = UnivariateSpline(np.sort(tmp_df[scan_variable].unique()), std / np.max(std))
            std_spl.set_smoothing_factor(std_sf)

            ste = np.array(ste)
            ste_spl = UnivariateSpline(np.sort(tmp_df[scan_variable].unique()), ste/np.max(ste))
            ste_spl.set_smoothing_factor(std_ste_sf)

            new_x = np.linspace(np.sort(tmp_df[scan_variable].unique()).min(),
                                np.sort(tmp_df[scan_variable].unique()).max(), num=len(tmp_df[scan_variable].unique()) * 10,
                                endpoint=True)
            # new_x = np.linspace(np.sort(tmp_df[scan_variable].unique()).min(),
            #                     3, num=len(tmp_df[scan_variable].unique()) * 10,
            #                     endpoint=True)

            sns.lineplot(x = new_x, y = std_spl(new_x) * np.max(std),
                         color=y_colours[d][1], linestyle=linestyles[d][1], ci=0, label="SD", legend=plot_legend)
            plt.fill_between(new_x,
                             np.clip(std_spl(new_x)*np.max(std) - ste_spl(new_x)*np.max(ste), 0, None),
                             np.clip(std_spl(new_x)*np.max(std) + ste_spl(new_x)*np.max(ste), 0, None), alpha=0.15, color=y_colours[d][1])
        else:
            sns.lineplot(x = np.sort(tmp_df[scan_variable].unique()), y = std,
                     color=y_colours[d][1], linestyle=linestyles[d][1], ci=0, label="SD", legend=plot_legend)
            plt.fill_between(np.sort(tmp_df[scan_variable].unique()),
                         np.clip(np.array(std) - np.array(ste), 0, None),
                         np.clip(np.array(std) + np.array(ste), 0, None), alpha=0.15, color=y_colours[d][1])


if scan_variable == "IL-2_Tsec_fraction":
    xlabel = "fraction of sec. cells"
    standard = 0.1
elif scan_variable == "IL-2_sigma":
    xlabel = "Receptor heterogeneity"
    standard = 1.5
elif scan_variable == "IL-2_KD":
    xlabel = "1/K$_{D}$ (pM)"
    standard = 1/7.4

plt.axvline(standard, 0, 10)
ax.set(yscale=yscale, xscale=xscale, ylim=ylim, xlim=xlim)
ax.set_xlabel(xlabel)
ax.set_ylabel("pM")

if scan_variable == "IL-2_sigma":
    ax.set_xticks([0,1,2, 3, 4])
    ax.set_xticklabels([0, 1, 2, 3, 4])
# elif scan_variable == "IL-2_KD":
#     ax.set_xticks([0, 5, 10])
#     ax.set_xticklabels([0, 5, 10])

handles, labels = ax.get_legend_handles_labels()
labels = legend_labels
handles = [handles[0], handles[1], handles[2]]
plt.legend(handles, labels)

fig.savefig(saving_string, bbox_inches='tight')
plt.tight_layout()
plt.show()