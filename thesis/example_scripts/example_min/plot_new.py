from  thesis.main.MyPlotter import Plotter
from parameters import path, IMGPATH
import matplotlib.pyplot as plt
import numpy as np
import seaborn as sns
from scipy.constants import golden

# path = [
#     '/extra/lukas/example_min/new',
#     '/extra/lukas/example_min/new_1',
# ]

plotter = Plotter(path)

plotter.scan_scale = np.linspace(0,1,5)
# plotter.scan_scale = np.logspace(-1,1,10)

plotter.label_replacement.update({

    "kd":"$k_d$",
    plotter.time_key: "time (a.u)",
    "sec_amax":"$a_{max}$",
    "sec_q":"$q$",
    plotter.scan_index_key:"parameter fold-change",
    "Concentration":"Concentration (nM)",
    "Gradient": "Gradient (nM$ / \mu m)$",
    "SD": "SD (nM)",
    "IL-2_surf_c":"IL-2 surface concentration (nM)",
})

plotter.color_dict = {
    "IL-2":"tab:red",
    "IL-2_surf_c":"tab:red",
    "sec": "tab:orange",
    "abs": "tab:green",
    "default":"tab:gray",

}

l = 8.3 * (1-0.1)#page length x
n =  4#rows
m = 3#columns



plotter.subplots(n,m, figsize = (l,l), external_legend = "axes")

plotter.global_time_series_plot("Concentration", hue = "field_name", legend="brief", style=plotter.scan_index_key)
plotter.global_time_series_plot("Gradient", hue = "field_name")
plotter.global_time_series_plot("SD", hue = "field_name")
plotter.count_plot(hue= "type_name", relative=True, legend="brief", style=plotter.scan_index_key)


plotter.cells_time_series_plot("IL-2_surf_c",hue="type_name")
plotter.cell_slice_plot("IL-2_surf_c", hue= "type_name", palette_name = "viridis")


# plotter.cells_time_series_plot("default_score_init_norm",hue= "type_name")
# plotter.cells_time_series_plot("sec_score_init_norm",hue= "type_name")
# plotter.cells_time_series_plot("abs_score_init_norm",hue= "type_name")

plotter.filter = lambda df: df.loc[df[plotter.scan_name_key] == "fractions"]

# plotter.cell_histogramm("IL-2_surf_c", hue = "type_name", distplot_kwargs = {"bins":100,"kde":False, "norm_hist":True})
plotter.global_steady_state_plot("Concentration", hue = "path_name")
plotter.steady_state_count(hue="path_name",relative=True)
plotter.cell_steady_state_plot("IL-2_surf_c",hue="path_name")
plotter.make_legend()
plotter.savefig(IMGPATH+"collection.pdf")
plotter.show()

# plotter.subplots(2,2,(8,6), external_legend = False)
# n = 4
# ec50 = 10
# plotter.cell_histogramm("IL-2_surf_c" , overlay = True, distplot_kwargs = {"bins":100})
# plotter.function_twinx_overlay(lambda x:  x**n /(ec50**n + x**n), plot_args = ("r-",), ylabel="pSTAT5 activation")
# plotter.timing_barplot("task",hue = "name")
# plotter.timing_timelineplot()
#
# plotter.make_legend()
# plotter.savefig(IMGPATH + "hist.pdf")
# plotter.show()

# plotter.subplots(2,2,(8,6), external_legend = False)
#
# plotter.timing_barplot("task",hue = "name")
# plotter.timing_lineplot("duration", hue ="task",legend="brief")
# plotter.timing_timelineplot()
#
# plotter.make_legend()
# plotter.savefig(IMGPATH + "hist.pdf")
# plotter.show()

