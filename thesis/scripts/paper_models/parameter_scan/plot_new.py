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

cv_lim = [0,2]
c_lim = [1e-3,10]

plotter.subplots(3,4, figsize = (a,b/2), external_legend = "axes")
plotter.global_steady_state_plot("Concentration", style = "numeric_linear",ci = "sem", hue = plotter.scan_name_key, ylim = c_lim,legend="brief", ylog=True,average=True)
plotter.global_steady_state_plot("CV", style = "numeric_linear",ci = "sem",hue=plotter.scan_name_key, legend = False, ylog=False, ylim = cv_lim,average=True)
plotter.global_steady_state_plot("SD", style = "numeric_linear",ci = "sem", hue = plotter.scan_name_key, legend=False, ylog=True, average=True)
plotter.empty_plot()
plotter.steady_state_count(hue = plotter.scan_name_key,style="type_name", subtitle = "cell fraction\n(for distance scan)", relative=True, filter= lambda df:df.loc[(df[plotter.scan_name_key] == "distance")])
plotter.steady_state_count(hue = plotter.scan_name_key,style="type_name", subtitle = "cell fraction\n(for ratio scan)", relative=True, filter= lambda df:df.loc[(df[plotter.scan_name_key] == "ratio")])
plotter.cell_steady_state_plot("IL-2_q",subtitle = "systemic secretion\n(for ratio scan)",filter= lambda df:df.loc[(df.type_name == "sec") & (df[plotter.scan_name_key] == "ratio")],cummulative=True, legend="brief")
plotter.cell_steady_state_plot("IL-2_R",subtitle = "systemic IL-2R\n(for ratio scan)", filter= lambda df:df.loc[(df.type_name == "abs") & (df[plotter.scan_name_key] == "ratio")],cummulative=True, legend="brief")
plotter.cell_histogramm("IL-2_R",subtitle = "IL-2R distribution\n(for fold-change = 1)",filter = lambda df:df.loc[(df.type_name == "abs") & (df.scan_value == 1)] , distplot_kwargs={"bins":100}, ylog=False, xlog=False, xlim=[0,5e4])

plotter.make_legend()
plotter.savefig(IMGPATH + "collection_full.pdf")
plt.show()

