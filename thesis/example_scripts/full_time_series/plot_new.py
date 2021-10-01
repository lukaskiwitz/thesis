import numpy as np

from parameters import path, IMGPATH
from thesis.main.MyPlotter import Plotter

# path = [
#     '/extra/lukas/full_parameter_scan/new',
#     '/extra/lukas/full_parameter_scan/new_1',
# ]

plotter = Plotter(path)

plotter.scan_scale = s = 10
scan_space = np.concatenate([np.logspace(-1, 0, int(s / 2)), np.logspace(0, 1, int(s / 2) + 1)[1:]])

plotter.label_replacement.update({

    "kd": "$k_d$",
    plotter.time_key: "time (a.u)",
    "sec_amax": "$a_{max}$",
    "sec_q": "$q$",
    plotter.scan_index_key: "secretion rate",
    "Concentration": "Concentration (nM)",
    "Gradient": "Gradient (nM$ / \mu m)$",
    "SD": "SD (nM)",
    "IL-2_surf_c": "IL-2 surface concentration (nM)",
})

plotter.color_dict = {
    "IL-2": "tab:red",
    "IL-2_surf_c": "tab:red",
    "sec": "tab:orange",
    "abs": "tab:green",
    "default": "tab:gray",

}

MM_PER_INCH = 2.54 * 10
margin_a = 10
margin_b = 10

a = (128 - 2 * margin_a) / MM_PER_INCH
b = (96 - 2 * margin_b) / MM_PER_INCH

plotter.subplots(2, 2, figsize=(a, b), external_legend="axes")

plotter.global_time_series_plot("Concentration", hue=plotter.scan_index_key)
plotter.count_plot(style="type_name", hue=plotter.scan_index_key, ylog=False, legend="brief")
plotter.cell_steady_state_plot("IL-2_surf_c", hue=plotter.scan_name_key, ylog=False)
plotter.steady_state_count(style="type_name", legend="brief")
plotter.make_legend()
plotter.savefig(IMGPATH + "collection.pdf")
plotter.show()

plotter.ruse_plot(IMGPATH)
