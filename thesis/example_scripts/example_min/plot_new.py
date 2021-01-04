import numpy as np

from parameters import path, IMGPATH
from thesis.main.MyPlotter import Plotter

plotter = Plotter(path)



a = 8.3
b = np.sqrt(2) * a
a = a*0.9
b = 0.5*b

plotter.subplots(4,3, figsize=(a,b),external_legend="figure")
# plotter.scan_scale = np.linspace(0.01,0.2,5)
plotter.scan_scale = np.logspace(-1, 1, 10)
plotter.tight_layout_settings["pad"] = 1
plotter.label_replacement.update({

    "kd": "$k_d$",
    plotter.time_key: "time (a.u)",
    "sec_amax": "$a_{max}$",
    "sec_q": "$q$",
    plotter.scan_index_key: "Scan index",
    "Concentration": "Concentration (nM)",
    "Gradient": "Gradient (nM$ / \mu m)$",
    "SD": "SD (nM)",
    "IL-2_surf_c": "IL-2 surface \n concentration (nM)",
})

plotter.color_dict = {
    "IL-2": "tab:red",
    "IL-2_surf_c": "tab:red",
    "sec": "tab:orange",
    "abs": "tab:green",
    "default": "tab:gray",

}
plotter.global_time_series_plot("Concentration", hue="field_name", legend="brief", style=plotter.scan_index_key)
plotter.global_time_series_plot("Gradient", hue="field_name")
plotter.global_time_series_plot("SD", hue="field_name")

plotter.count_plot(hue="type_name", relative=True, legend="brief", style=plotter.scan_index_key)
plotter.cells_time_series_plot("IL-2_surf_c", hue="type_name")

# plotter.cell_slice_plot("IL-2_surf_c", hue= "type_name", palette_name = "viridis")


# plotter.cells_time_series_plot("default_score_init_norm",hue= "type_name")
# plotter.cells_time_series_plot("sec_score_init_norm",hue= "type_name")
# plotter.cells_time_series_plot("abs_score_init_norm",hue= "type_name")

# plotter.cell_histogramm("IL-2_surf_c",t = [0], hue = "type_name",bins = 10, kde = False, norm_hist = True)

plotter.make_legend()
plotter.savefig(IMGPATH + "collections.pdf")
plotter.show()
