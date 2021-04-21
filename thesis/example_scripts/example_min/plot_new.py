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

plotter.subplots(2,2,(8,6), external_legend = "axes")

plotter.global_steady_state_plot("Concentration", hue = plotter.scan_name_key, ylog=True, legend="brief")
plotter.cell_steady_state_plot("IL-2_surf_c",hue = plotter.scan_name_key, ylog=True)
plotter.steady_state_count(hue=plotter.scan_name_key, style="type_name",legend="brief")
plotter.make_legend()
plotter.savefig(IMGPATH+"collection.pdf")
plotter.show()


