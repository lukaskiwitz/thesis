import matplotlib.pyplot as plt
import numpy as np

from parameters import path, IMGPATH
from thesis.main.MyPlotter import Plotter
from thesis.main.my_debug import message

plotter = Plotter(path)
message("dataset loaded")

plotter.label_replacement.update({

    "D": "Diffusion coefficient",
    "sec_amax": "$a_{max}$",
    "sec_q": "IL-2 / (cell * s)",
    "f_sec": "% secretors",
    "f_abs": "% consumers",
    "f_new": "ratio",
    "numeric_linear": "linear",
    "abs_R": "IL-2R / cell",
    "Kc": "IL-2R EC50",
    "kd": "IL-2 decay",
    plotter.scan_index_key: "parameter fold-change",
    "Concentration": "Concentration (nM)",
    "run": "total",
    "run:scan_sample:SimContainer:run": "time series",
    "run:scan_sample:SimContainer:run:step": "timestep",
    "run:scan_sample:update_sim_container": "update sim_container",
    "run:write_element_tree": "write xml file",
    "run:scan_sample": "scan sample"

})

c_lim = [0, 0.05]
sd_lim = [0, 0.01]
cv_lim = [0, 1]

b = 8.3 / 5
a = np.sqrt(2) * b

dashes = {False: (5, 2), True: (1, 0)}
plotter.subplots(3, 3, figsize=(3 * a, 3 * b), external_legend="axes")

ci = "sem"
plotter.global_steady_state_plot("Concentration", hue=plotter.scan_name_key, ci=ci,
                                 legend="brief", ylog=False, xlog=False, average=True, ylim=c_lim, dashes=dashes)
plotter.global_steady_state_plot("CV", hue=plotter.scan_name_key, legend=False, ci=ci,
                                 ylog=False, xlog=False, ylim=cv_lim, average=True, dashes=dashes)
plotter.global_steady_state_plot("SD", hue=plotter.scan_name_key, legend=False, ci=ci,
                                 ylog=False, xlog=False, average=True, dashes=dashes, ylim=sd_lim)
plotter.filter = lambda df: df.loc[df["type_name"] == "abs"]

plotter.cell_steady_state_plot("IL-2_surf_c", hue=plotter.scan_name_key, ci=ci, legend=False,
                               ylog=False, xlog=False, dashes=dashes, ylim=c_lim)
plotter.global_steady_state_plot("surf_c_cv", hue=plotter.scan_name_key, ci=ci, legend=False,
                                 ylog=False, xlog=False, ylim=cv_lim, average=True, dashes=dashes)
plotter.global_steady_state_plot("surf_c_std", hue=plotter.scan_name_key, ci=ci, legend=False,
                                 ylog=False, xlog=False, average=True, dashes=dashes, ylim=sd_lim)

plotter.filter = lambda df: df.loc[df["type_name"] == "sec"]

plotter.cell_steady_state_plot("IL-2_surf_c", hue=plotter.scan_name_key, ci=ci, legend=False,
                               ylog=False, xlog=False, dashes=dashes, ylim=c_lim)
plotter.global_steady_state_plot("surf_c_cv", hue=plotter.scan_name_key, ci=ci, legend=False,
                                 ylog=False, xlog=False, ylim=cv_lim, average=True, dashes=dashes)
plotter.global_steady_state_plot("surf_c_std", hue=plotter.scan_name_key, ci=ci, legend=False,
                                 ylog=False, xlog=False, average=True, dashes=dashes, ylim=sd_lim)

plotter.make_legend()
plotter.savefig(IMGPATH + "collection.pdf")
plt.show()
