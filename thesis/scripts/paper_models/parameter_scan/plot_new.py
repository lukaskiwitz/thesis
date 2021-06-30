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
message("dataset loaded")



a = 8.3
b = np.sqrt(2) * a


plotter.subplots(3,3, figsize = (a,b/2), external_legend = "axes")


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

plotter.global_steady_state_plot("Concentration", style = "numeric_linear",ci = "sem", hue = plotter.scan_name_key, legend="brief", ylog=True,average=True)
plotter.global_steady_state_plot("CV", style = "numeric_linear",ci = "sem",hue=plotter.scan_name_key, legend = False, ylog=False, ylim = [0, 5],average=True)
plotter.global_steady_state_plot("SD", style = "numeric_linear",ci = "sem", hue = plotter.scan_name_key, legend=False, ylog=True, average=True)
plotter.filter = lambda df: df.loc[df["type_name"] == "abs"]

plotter.cell_steady_state_plot("IL-2_surf_c", style = "numeric_linear",ci = "sem",hue = plotter.scan_name_key, legend=False, ylog=True)
plotter.global_steady_state_plot("surf_c_cv", style = "numeric_linear", ci = "sem",hue = plotter.scan_name_key, legend=False, ylog=False,ylim = [0, 5],average=True)
plotter.global_steady_state_plot("surf_c_std", style = "numeric_linear",ci = "sem", hue = plotter.scan_name_key, legend=False, ylog=True,average=True)

plotter.filter = lambda df: df.loc[df["type_name"] == "sec"]

plotter.cell_steady_state_plot("IL-2_surf_c", style = "numeric_linear",ci = "sem",hue = plotter.scan_name_key, legend=False, ylog=True)
plotter.global_steady_state_plot("surf_c_cv", style = "numeric_linear", ci = "sem",hue = plotter.scan_name_key, legend=False, ylog=False,ylim = [0, 5],average=True)
plotter.global_steady_state_plot("surf_c_std", style = "numeric_linear",ci = "sem", hue = plotter.scan_name_key, legend=False, ylog=True,average=True)


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

plotter.global_steady_state_plot("Concentration", style = "numeric_linear",ci ="sem", hue = plotter.scan_name_key, legend="brief",xlog=True, ylog=True,average=True, filter = f_1)
plotter.global_steady_state_plot("CV", style = "numeric_linear",ci = "sem", hue=plotter.scan_name_key, legend = False, xlog=True, ylog=False, ylim = [0, 1],average=True, filter = f_1)


plotter.global_steady_state_plot("Concentration", style = "numeric_linear",ci ="sem",hue = plotter.scan_name_key, legend="brief", ylog=True,average=True,filter = f_2)
plotter.global_steady_state_plot("CV", style = "numeric_linear",ci ="sem",hue=plotter.scan_name_key, legend = False, ylog=False, ylim = [0, 1],average=True, filter = f_2)

plotter.global_steady_state_plot("Concentration", style = "numeric_linear",ci ="sem",hue = plotter.scan_name_key, legend="brief", ylog=True,average=True,filter = f_3)
plotter.global_steady_state_plot("CV", style = "numeric_linear",ci ="sem",hue=plotter.scan_name_key, legend = False, ylog=False, ylim = [0, 1],average=True,filter = f_3)
plotter.make_legend()
plotter.savefig(IMGPATH + "collection_seperate.pdf")
plt.show()


# select = {"task": [
#     # "run",
#     "run:scan_sample",
#     "run:scan_sample:SimContainer:run",
#     "run:scan_sample:SimContainer:run:step",
#     "run:scan_sample:update_sim_container",
#     "run:write_element_tree"
# ]}


# plotter.subplots(1, 2, figsize=(a,b/4), external_legend="figure")
# plotter.timing_barplot("task", select=select, legend="brief")
# plotter.timing_lineplot("duration",x_name="scan_index", hue="task", select=select)
# # plotter.timing_timelineplot(select=select)
# plotter.make_legend()
# plotter.savefig(IMGPATH + "timing.pdf")
# plotter.show()
#
# d = plotter.timing_df
# T = d.loc[d["task"] =="run"]["duration"]
# for t in T:
#     print("total runtime {t}".format(t = t/60**2))
# message("done")
