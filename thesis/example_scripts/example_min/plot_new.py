from  thesis.main.MyPlotter import Plotter
from parameters import path, IMGPATH
import matplotlib.pyplot as plt
import numpy as np
import seaborn as sns
from scipy.constants import golden
plotter = Plotter(path)

l = 8.3 * (1-0.1)#page length x
n =  4#rows
m = 3#columns


plotter.subplots(n,m, figsize = (l,l))
# plotter.scan_scale = np.linspace(0.01,0.2,5)
plotter.scan_scale = np.logspace(-1,1,10)

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
plotter.global_time_series_plot("Concentration", hue = "field_name", legend="brief", style=plotter.scan_index_key)
plotter.global_time_series_plot("Gradient", hue = "field_name")
plotter.global_time_series_plot("SD", hue = "field_name")
plotter.count_plot(hue= "type_name", relative=True, legend="brief", style=plotter.scan_index_key)


plotter.cells_time_series_plot("IL-2_surf_c",hue="type_name")

# plotter.cell_slice_plot("IL-2_surf_c", hue= "type_name", palette_name = "viridis")


plotter.cells_time_series_plot("default_score_init_norm",hue= "type_name")
plotter.cells_time_series_plot("sec_score_init_norm",hue= "type_name")
plotter.cells_time_series_plot("abs_score_init_norm",hue= "type_name")

plotter.cell_histogramm("IL-2_surf_c",t = [0], hue = "type_name",bins = 10, kde = False, norm_hist = True)
plotter.global_steady_state_plot("Concentration", hue = "field_name")
plotter.steady_state_count(hue="type_name",relative=True)
plotter.cell_steady_state_plot("IL-2_surf_c",hue="type_name")

plotter.make_legend()
plotter.show()
