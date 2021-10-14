import numpy as np

from parameters import path, IMGPATH
from thesis.main.MyPlotter import Plotter

# path = [
#     '/extra/lukas/distance_scan/new',
#     '/extra/lukas/distance_scan/new_1',
# ]

plotter = Plotter(path)

plotter.label_replacement.update({

    "kd": "$k_d$",
    plotter.time_key: "time (a.u)",
    "sec_amax": "$a_{max}$",
    "sec_q": "$q$",
    plotter.scan_index_key: "cell-cell distance_scan ($\mu m$)",
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

plotter.scan_scale = np.linspace(5, 15, 5)
MM_PER_INCH = 2.54 * 10
margin_a = 10
margin_b = 10

a = (128 - 2 * margin_a) / MM_PER_INCH
b = (96 - 2 * margin_b) / MM_PER_INCH
c_lim = [0, 0.02]
plotter.subplots(2, 2, figsize=(a, b), external_legend="axes")

plotter.global_steady_state_plot("Concentration", hue=plotter.scan_name_key, style=plotter.model_name_key, ylog=False,
                                 xlog=False, legend="brief",
                                 ci="sd", average=True, ylim=c_lim)
plotter.cell_steady_state_plot("IL-2_surf_c", hue=plotter.scan_name_key, style=plotter.model_name_key, ylog=False,
                               xlog=False, ylim=c_lim)
plotter.steady_state_count(hue=plotter.scan_name_key, style="type_name", legend="brief", xlog=False)
plotter.make_legend()
plotter.savefig(IMGPATH + "collection.pdf")
plotter.show()
