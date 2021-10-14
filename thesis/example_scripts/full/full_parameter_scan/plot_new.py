import matplotlib.pyplot as plt
import numpy as np

from parameters import path, IMGPATH
from thesis.main.MyPlotter import Plotter
from thesis.main.my_debug import message

path = path
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

c_lim = [1e-4, 10]
cv_lim = [0, 2]

b = 8.3 / 5
a = np.sqrt(2) * b

dashes = {False: (5, 2), True: (1, 0)}
plotter.subplots(3, 3, figsize=(3 * a, 3 * b), external_legend="axes")

f = lambda df: df.loc[~df[plotter.scan_name_key].isin(["f_new", "Kc"])]

plotter.filter = f

ci = "sem"

plotter.global_steady_state_plot("Concentration", style="numeric_linear", hue=plotter.scan_name_key, ci=ci,
                                 legend="brief", ylog=True, average=True, dashes=dashes)
plotter.global_steady_state_plot("CV", style="numeric_linear", hue=plotter.scan_name_key, legend=False, ci=ci,
                                 ylog=False, ylim=cv_lim, average=True, dashes=dashes)
plotter.global_steady_state_plot("SD", style="numeric_linear", hue=plotter.scan_name_key, legend=False, ci=ci,
                                 ylog=True, average=True, dashes=dashes)
plotter.filter = lambda df: f(df.loc[df["type_name"] == "abs"])

plotter.cell_steady_state_plot("IL-2_surf_c", style="numeric_linear", hue=plotter.scan_name_key, ci=ci, legend=False,
                               ylog=True, dashes=dashes)
plotter.global_steady_state_plot("surf_c_cv", style="numeric_linear", hue=plotter.scan_name_key, ci=ci, legend=False,
                                 ylog=False, ylim=cv_lim, average=True, dashes=dashes)
plotter.global_steady_state_plot("surf_c_std", style="numeric_linear", hue=plotter.scan_name_key, ci=ci, legend=False,
                                 ylog=True, average=True, dashes=dashes)

plotter.filter = lambda df: f(df.loc[df["type_name"] == "sec"])

plotter.cell_steady_state_plot("IL-2_surf_c", style="numeric_linear", hue=plotter.scan_name_key, ci=ci, legend=False,
                               ylog=True, dashes=dashes)
plotter.global_steady_state_plot("surf_c_cv", style="numeric_linear", hue=plotter.scan_name_key, ci=ci, legend=False,
                                 ylog=False, ylim=cv_lim, average=True, dashes=dashes)
plotter.global_steady_state_plot("surf_c_std", style="numeric_linear", hue=plotter.scan_name_key, ci=ci, legend=False,
                                 ylog=True, average=True, dashes=dashes)

plotter.make_legend()
plotter.savefig(IMGPATH + "collection.pdf")
plt.show()

# f = lambda df: df.loc[
#     (df[plotter.scan_name_key] == "sec_q")|
#     (df[plotter.scan_name_key] == "abs_R")
# ]
#
# plotter.subplots(2, 2, figsize=(2 * a ,2 * b), external_legend="axes")
# plotter.filter = f
# abs = lambda df: df.loc[df["type_name"] == "abs"]
# sec = lambda df: df.loc[df["type_name"] == "sec"]
#
# plotter.calc_cell_activation(n_il2=4)
#
# plotter.global_steady_state_plot("Concentration", style="numeric_linear", hue=plotter.scan_name_key, legend=False,
#                                  ci=ci, ylog=True, average=True)
# plotter.global_steady_state_plot("CV", style="numeric_linear", hue=plotter.scan_name_key, legend=False, ylog=False,
#                                  ci=ci, ylim=cv_lim, average=True)
# plotter.cell_steady_state_plot("activation", style="numeric_linear", hue=plotter.scan_name_key, legend=False,
#                                ylog=False, ci=ci, ylim=[0, 1], average=True, filter=abs)
#
# plotter.cell_steady_state_plot("activation", style="numeric_linear", hue=plotter.scan_name_key, legend="brief",
#                                ylog=False, ci=ci, ylim=[0, 1], average=True, filter=sec)
# plotter.make_legend()
# plotter.savefig(IMGPATH + "collection_small_qR.pdf")
# plt.show()
#
# f = lambda df: df.loc[
#     (df[plotter.scan_name_key] == "f_new")
# ]
#
# plotter.subplots(2, 1, figsize=(a,2*b), external_legend="axes")
# plotter.filter = f
#
# plotter.global_steady_state_plot("Concentration", style="numeric_linear", hue=plotter.scan_name_key, legend="brief",
#                                  ci=ci, xlog=True, ylog=True, average=True, ylim=c_lim)
# plotter.global_steady_state_plot("CV", style="numeric_linear", hue=plotter.scan_name_key, legend=False, xlog=True,
#                                  ci=ci, ylog=False, ylim=cv_lim, average=True)
# plotter.make_legend()
# plotter.savefig(IMGPATH + "collection_small_fractions.pdf")
# plt.show()
#
# f = lambda df: df.loc[
#     (df[plotter.scan_name_key] == "sec_q") |
#     (df[plotter.scan_name_key] == "abs_R") |
#     (df[plotter.scan_name_key] == "f_new")
#     ]
#
#
# plotter.subplots(2, 1, figsize=(a, 2*b), external_legend="axes")
# plotter.filter = f
#
# plotter.global_steady_state_plot("Concentration", style="numeric_linear", ci=ci, ylim=c_lim, hue=plotter.scan_name_key,
#                                  legend="brief", ylog=True, average=True)
# plotter.global_steady_state_plot("CV", style="numeric_linear", ci=ci, hue=plotter.scan_name_key, legend=False,
#                                  ylog=False, ylim=cv_lim, average=True)
# plotter.make_legend()
# plotter.savefig(IMGPATH + "collection_small_fqR.pdf")
# plt.show()
#
# f = lambda df: df.loc[
#     (df[plotter.scan_name_key] == "D") |
#     (df[plotter.scan_name_key] == "kd") |
#     (df[plotter.scan_name_key] == "Kc")
#     ]
#
#
#
# plotter.subplots(2, 1, figsize=(a, 2*b), external_legend="axes")
# plotter.filter = f
#
# plotter.global_steady_state_plot("Concentration", style="numeric_linear", hue=plotter.scan_name_key, legend="brief",
#                                  ci=ci, ylog=True, average=True, ylim=c_lim)
# plotter.global_steady_state_plot("CV", style="numeric_linear", hue=plotter.scan_name_key, legend=False, ylog=False,
#                                  ci=ci, ylim=cv_lim, average=True)
# plotter.make_legend()
# plotter.savefig(IMGPATH + "useless_parameters.pdf")
# plt.show()
#
# f = lambda df: df.loc[
#     (df[plotter.scan_name_key] == "kd")
# ]
#
#
# plotter.subplots(2, 1, figsize=(a,2* b), external_legend="axes")
# plotter.filter = f
#
# plotter.global_steady_state_plot("Concentration", style="numeric_linear", hue=plotter.scan_name_key, legend="brief",
#                                  ci=ci, ylog=True, average=True)
# plotter.global_steady_state_plot("CV", style="numeric_linear", hue=plotter.scan_name_key, legend=False, ylog=False,
#                                  ci=ci, ylim=cv_lim, average=True)
# plotter.make_legend()
# plotter.savefig(IMGPATH + "decay.pdf")
# plt.show()
#
#
#
# f = lambda df: df.loc[
#     (df[plotter.scan_name_key] == "f_new")|
#     (df[plotter.scan_name_key] == "D")
#     ]
# #
#
# plotter.scan_scale = np.linspace(0, 1, 20)
# plotter.subplots(2, 1, figsize=(a, 2*b), external_legend="axes")
# plotter.filter = f
# f = lambda df: df.loc[df.type_name == "naive"]
# # f_abs = lambda df: df.loc[df[plotter.scan_name_key] == "f_abs"]
#
# plotter.steady_state_count(style="type_name", hue=plotter.scan_name_key, legend="brief", ylim=[0, 1], ylog=False, ci=ci,
#                            xlog=True, relative=True, )
# plotter.steady_state_count(style="type_name", hue=plotter.scan_name_key, legend="brief", ylim=[0, 1], ylog=False, ci=ci,
#                            xlog=True, relative=True, )
#
# plotter.make_legend()
# plotter.savefig(IMGPATH + "fractions.pdf")
# plt.show()
#
# select = {"task": [
#     # "run",
#     "run:scan_sample",
#     "run:scan_sample:SimContainer:run",
#     "run:scan_sample:SimContainer:run:step",
#     "run:scan_sample:update_sim_container",
#     "run:write_element_tree"
# ]}
#
#
# plotter.subplots(2, 2, figsize=(2*a, 2*b), external_legend="axes")
# plotter.timing_barplot("task", select=select, legend="brief")
# plotter.timing_lineplot("duration", x_name="scan_index", hue="task", select=select)
# # plotter.timing_timelineplot(select=select)
# plotter.make_legend()
# plotter.savefig(IMGPATH + "timing.pdf")
# plotter.show()
#
# d = plotter.timing_df
# T = d.loc[d["task"] == "run"]["duration"]
# for t in T:
#     print("total runtime {t}".format(t=t / 60 ** 2))
# message("done")
