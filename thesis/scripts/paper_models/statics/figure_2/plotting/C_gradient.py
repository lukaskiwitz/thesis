import matplotlib.pyplot as plt
import numpy as np
import seaborn as sns
import os
from plotting_rc import rc_ticks

sns.set_theme(context = "talk", style = "ticks", rc = rc_ticks)
fig,ax = plt.subplots()

from thesis.scripts.paper_models.utilities.plot_helper import my_load_df, my_interpolation

'''
define which static sim you want to plot. 
bc defines the boundary condition, scan variable the scan.
"IL-2_Tsec_fraction" is the secreting cells scan, "IL-2_sigma" the R_lognorm scan and IL-2_KD the KD scan.
'''

# bc = "standard"
bc = "saturated"

scan_variables = ["IL-2_Tsec_fraction", "IL-2_sigma", "IL-2_KD"]

saving_string = r"/home/brunner/Documents/Current work/2022_03_04/" + "gradient_{m}_static_{sv}_{bc}".format(m = "spatial_ODE_smoothed", sv=scan_variables[0], bc = bc) + ".pdf"

sv_styles = ['-','-.',':']
plot_names = ["sec. scan", "R scan", "KD scan"]
# define the standard values (fold-change = 0) for each of the scan_variables
standards = [0.05, 1.5, 1/7.4]

# which models to plot

# models = ["loc. q and R", "loc. q", "well_mixed"]
models = ["loc. q and R", "well_mixed"]
model_styles = ["-", "--"]
# Which colors to use for each scan_variable. Multiple lists are the models
model_colors = [["red", "green", "blue"], ["red", "green", "blue"], ["blue"]]
model_alphas = [1,0.5]

plot_legend = True


# plotting parameters

yscale = "linear"
yscale_base = None
xscale = "log"
xscale_base = 2
# define data start and end point
startingPoint = None
stoppingPoint = None
# maximum value
max_value = 1e5
min_value = -1

xlim = (2 ** -2.1, 2 ** 2.1)
ylim = (None, None)

# defines a outer layer of N cells to be ignored in plotting. Used to further limit unwanted boundary effects.
offset = 0
# which cell type to plot in which colour. The second colour is the SD colour
cell_types = ["Th"]

scan_measure = "IL-2_surf_c"

replicat = "replicat_index"

# if line smoothing is desired, sfs are the smoothing factors
interpolate = True
sfs = [0.002, 0.01, 0.01]

########################################################################################################################
########################################################################################################################
##                                              Plotting                                                              ##
########################################################################################################################
########################################################################################################################

# load runs, apply offset, merge dataframes
dataframes = [[] for x in models]
print("loading data")
from funcs_for_plotting import get_path, scale_dataframes
hdd = "/extra2" if os.path.exists("/extra2") else "/extra"
for m,model in enumerate(models):
    for sv, scan_variable in enumerate(scan_variables):
        m_sv_path = get_path(bc, hdd, model, scan_variable)
        c_df, g_df = my_load_df(m_sv_path, offset = 0, run_range = [0], custom_ending = "_combined")
        # scales the dataframes. Depends on your desired units. Currently everything is in pM
        c_df, g_df = scale_dataframes(c_df, g_df, model, scan_variable, sv, min_value, max_value, standards)
        dataframes[m].append([c_df, g_df])

########################################################################################################################
print("plotting")
########################################################################################################################
#%%
for m,model in enumerate(models):
    for sv,scan_variable in enumerate(scan_variables):
        c_df = dataframes[m][sv][0]
        g_df = dataframes[m][sv][1]

        ste = []
        grad = []
        for s, scan in enumerate(np.sort(c_df[scan_variable].unique())):
            sliced_df = g_df.loc[(g_df[scan_variable] == scan)]
            run_grad = []
            try:
                for r, run in enumerate(np.sort(g_df[replicat].unique())):
                    run_grad.append(sliced_df.loc[(sliced_df[replicat] == run), "Gradient"])
            except KeyError: # ODE system does not have a gradient
                run_grad = [0]

            ste.append(np.std(run_grad) / np.sqrt(len(g_df[replicat].unique())))
            grad.append(np.mean(run_grad))

        x = np.sort(c_df[scan_variable].unique())
        x /= x[np.abs(x - standards[sv]).argmin()]
        grad = np.array(grad) * 1e3

        if interpolate == True:
            new_x, spl, ax = my_interpolation(x, np.array(grad), smoothing_factor=sfs[sv], label=str(model), color=model_colors[m][sv], linestyle = model_styles[m], alpha=model_alphas[m])
        else:
            plt.plot(x, grad, label=str(model), color=model_colors[m][sv], linestyle = model_styles[m], alpha=model_alphas[m])

x_axis_label = "log fold change"

plt.ylabel("gradient")
plt.xlabel(x_axis_label)
plt.xscale(xscale, basex=xscale_base)
plt.yscale(yscale, basey=yscale_base)
plt.ylim(ylim)
plt.xlim(xlim)


import matplotlib
locmaj = matplotlib.ticker.LogLocator(base=2,numticks=1000)
ax.xaxis.set_major_locator(locmaj)
locmin = matplotlib.ticker.LogLocator(base=2,subs=np.arange(1,10)*.1 + 1,numticks=1000)
ax.xaxis.set_minor_locator(locmin)
ax.xaxis.set_minor_formatter(matplotlib.ticker.NullFormatter())

ax.set_xticks([1/4, 1/2, 1, 2, 4])
ax.set_xticklabels(["-2","", "0","", "2"])

import matplotlib.lines as mlines
handles = []
labels = ["sec. cells", "R hetero.", "1/K$_D$", "well-mixed"]
if plot_legend == True:
    for sv, scan_variable in enumerate(scan_variables):
        handles.append(mlines.Line2D([], [], color=model_colors[0][sv], alpha = model_alphas[0], linestyle = model_styles[0]))
    handles.append(mlines.Line2D([], [], color="black", alpha = model_alphas[1], linestyle = model_styles[1]))

    ax.figure.legend(handles, labels)  # (loc=(0.385,0.2))

# save the plots

fig.savefig(saving_string, bbox_inches='tight')
plt.tight_layout
plt.show()