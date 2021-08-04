from  thesis.main.MyPlotter import Plotter
import matplotlib.pyplot as plt
import numpy as np
import sys
import os
from thesis.main.my_debug import message

if len(sys.argv) > 1:
    path = sys.argv[1]
    IMGPATH = os.path.join(path,"cli_images")

else:
    from parameters import path, IMGPATH

path = path
plotter = Plotter(path)
plotter.calc_cell_activation()
message("dataset loaded")



a = 8.3
b = np.sqrt(2) * a


plotter.label_replacement.update({

    "D":"Diffusion coefficient",
    "sec_amax":"$a_{max}$",
    "sec_q":"IL-2 / (cell * s)",
    "f_sec":"% secretors",
    "f_abs":"% consumers",
    "abs_R":"IL-2R / cell",
    "Kc":"IL-2R EC50",
    "kd":"IL-2 decay",
    plotter.scan_index_key:"parameter fold-change",
    "Concentration":"Concentration (nM)",
    "run": "total",
    "run:scan_sample:SimContainer:run": "time series",
    "run:scan_sample:SimContainer:run:step": "timestep",
    "run:scan_sample:update_sim_container": "update sim_container",
    "run:write_element_tree": "write xml file",
    "run:scan_sample":"scan sample"

})

cv_lim = [0,3]
c_lim = [1e-3,1]

plotter.subplots(3,4, figsize = (a,b/2), external_legend = "axes")
plotter.global_steady_state_plot("Concentration", style = "numeric_linear",ci = "sem", hue = plotter.scan_name_key, ylim = c_lim,legend="brief", ylog=True,average=True)
plotter.global_steady_state_plot("CV", style = "numeric_linear",ci = "sem",hue=plotter.scan_name_key, legend = False, ylog=False, ylim = cv_lim,average=True)
plotter.global_steady_state_plot("SD", style = "numeric_linear",ci = "sem", hue = plotter.scan_name_key, legend=False, ylog=True, average=True)
plotter.global_steady_state_plot("SD", style = "numeric_linear",ci = "sem", hue = plotter.scan_name_key, legend=False, ylog=True, average=True)

plotter.filter = lambda df: df.loc[df["type_name"] == "abs"]

plotter.cell_steady_state_plot("IL-2_surf_c", style = "numeric_linear",ci = "sem",hue = plotter.scan_name_key, ylim = c_lim,legend=False, ylog=True)
plotter.global_steady_state_plot("surf_c_cv", style = "numeric_linear", ci = "sem",hue = plotter.scan_name_key, legend=False, ylog=False,ylim = cv_lim,average=True)
plotter.global_steady_state_plot("surf_c_std", style = "numeric_linear",ci = "sem", hue = plotter.scan_name_key, legend=False, ylog=True,average=True)
plotter.cell_steady_state_plot("activation", style = "numeric_linear",ci = "sem",hue = plotter.scan_name_key, ylim = c_lim,legend=False, ylog=False)



plotter.filter = lambda df: df.loc[df["type_name"] == "sec"]

plotter.cell_steady_state_plot("IL-2_surf_c", style = "numeric_linear",ci = "sem",hue = plotter.scan_name_key, ylim = c_lim,legend=False, ylog=True)
plotter.global_steady_state_plot("surf_c_cv", style = "numeric_linear", ci = "sem",hue = plotter.scan_name_key, legend=False, ylog=False,ylim = cv_lim,average=True)
plotter.global_steady_state_plot("surf_c_std", style = "numeric_linear",ci = "sem", hue = plotter.scan_name_key, legend=False, ylog=True,average=True)
plotter.cell_steady_state_plot("activation", style = "numeric_linear",ci = "sem",hue = plotter.scan_name_key, ylim = c_lim,legend=False, ylog=False)

plotter.make_legend()
plotter.savefig(IMGPATH + "collection_full.pdf")
plt.show()



f_1 = lambda df: df.loc[
    (df[plotter.scan_name_key] == "ratio")
    ]

f_2 = lambda df: df.loc[
    (df[plotter.scan_name_key] == "sec_q")|
    (df[plotter.scan_name_key] == "abs_R")|
    (df[plotter.scan_name_key] == "f_sec")|
    (df[plotter.scan_name_key] == "f_abs")
    ]

f_3 = lambda df: df.loc[
    (df[plotter.scan_name_key] == "D")|
    (df[plotter.scan_name_key] == "kd")|
    (df[plotter.scan_name_key] == "Koff")|
    (df[plotter.scan_name_key] == "kendo")
    ]


plotter.subplots(3,2, figsize = (2* a/3,b/2), external_legend = "axes")

plotter.global_steady_state_plot("Concentration", style = "numeric_linear",ci ="sem", hue = plotter.scan_name_key, legend="brief",xlog=True, ylog=True,ylim = c_lim,average=True, filter = f_1)
plotter.global_steady_state_plot("CV", style = "numeric_linear",ci = "sem", hue=plotter.scan_name_key, legend = False, xlog=True, ylog=False, ylim =cv_lim,average=True, filter = f_1)


plotter.global_steady_state_plot("Concentration", style = "numeric_linear",ci ="sem",hue = plotter.scan_name_key, legend="brief", ylog=True,ylim = c_lim,average=True,filter = f_2)
plotter.global_steady_state_plot("CV", style = "numeric_linear",ci ="sem",hue=plotter.scan_name_key, legend = False, ylog=False, ylim =cv_lim,average=True, filter = f_2)

plotter.global_steady_state_plot("Concentration", style = "numeric_linear",ci ="sem",hue = plotter.scan_name_key, legend="brief", ylog=True,ylim = c_lim,average=True,filter = f_3)
plotter.global_steady_state_plot("CV", style = "numeric_linear",ci ="sem",hue=plotter.scan_name_key, legend = False, ylog=False, ylim = cv_lim,average=True,filter = f_3)
plotter.make_legend()
plotter.savefig(IMGPATH + "collection_seperate.pdf")
plt.show()

select = {"task": [
    # "run",
    "run:scan_sample",
    "run:scan_sample:SimContainer:run",
    "run:scan_sample:SimContainer:run:step",
    "run:scan_sample:update_sim_container",
    "run:write_element_tree"
]}


plotter.subplots(1, 2, figsize=(a,b/4), external_legend="figure")
plotter.timing_barplot("task", select=select, legend="brief")
plotter.timing_lineplot("duration",x_name="scan_index", hue="task", select=select)
# plotter.timing_timelineplot(select=select)
plotter.make_legend()
plotter.savefig(IMGPATH + "timing.pdf")
plotter.show()

d = plotter.timing_df
T = d.loc[d["task"] =="run"]["duration"]
for t in T:
    print("total runtime {t}".format(t = t/60**2))
message("done")


metrics = {"min_distance": np.min, "mean_distance": np.mean, "max_distance": np.max}
plotter.compute_cell_distance_metric("sec", metric_dict=metrics, round = 0)

cell_df = plotter.cell_df
cell_df = cell_df.loc[(cell_df["numeric_linear"] == False) & (cell_df["scan_value"] == 1)  & (cell_df[plotter.scan_name_key] == "D")]

abs_f = lambda df: df.loc[df.type_name == "abs"]
sec_f = lambda df: df.loc[df.type_name == "sec"]

import seaborn as sns
fig,ax = plt.subplots(2,2)
ax = np.ravel(ax)

sns.boxplot(x="type_name", y = "IL-2_surf_c", data=cell_df, ax = ax[0])
# sns.stripplot(x="type_name", y = "IL-2_surf_c", data=cell_df, ax = ax[0], color="gray",size=2)

sns.boxplot(x="type_name", y = "activation", data=cell_df,ax = ax[1])
# sns.stripplot(x="type_name", y = "activation", data=cell_df,ax = ax[1], color="gray",size=2)

sns.boxplot(x="min_distance", y = "IL-2_surf_c", data=abs_f(cell_df), ax = ax[2])
# sns.stripplot(x="min_distance", y = "IL-2_surf_c", data=abs_f(cell_df), ax = ax[2], color="gray",size=2)

sns.boxplot(x="min_distance", y = "activation", data=abs_f(cell_df),ax = ax[3])
# sns.stripplot(x="min_distance", y = "activation", data=abs_f(cell_df),ax = ax[3], color="gray",size=2)

ax[1].set_xlabel("cell type")
ax[2].set_xlabel("cell type")

ax[2].set_xlabel("distance to closest secreting cell $(\mu m)$")
ax[3].set_xlabel("distance to closest secreting cell $(\mu m)$")
plt.tight_layout()
plt.savefig(os.path.join(IMGPATH,"boxplots.pdf"))
plt.show()