import numpy as np

from parameters import path, IMGPATH
from thesis.main.MyPlotter import Plotter

# path = '/extra/kiwitz/thesis_models/receptor_distribution_box_300/f_new_test_3_non_linear/'
# IMGPATH = os.path.join(path,"images") + "/"

plotter = Plotter(path, groups=["v_v"])
s = len(plotter.global_df[plotter.scan_index_key].unique())
# plotter.scan_scale = np.linspace(0, 1, s)
# plotter.scan_ticks = [0, 0.25, 0.5, 0.75, 1]

plotter.label_replacement.update({

    "kd": "$k_d$",
    "sec_amax": "$a_{max}$",
    "sec_q": "$q$",
    plotter.scan_index_key: "fraction of R on sec",
    "Concentration": "Concentration (nM)",
    "abs": "Responder",
    "sec": "Secretor",
    "naive": "Inert",
    "pSTAT5": "Activation",
    "IL-2_surf_c": "surface IL-2 (nM)",
    "IL-2_R": "IL-2R/cell",
    "numeric_ratio": "avg. R per cell",
    "v_v": "fraction of R in sec",
    "activation": "% pSTAT5",
    "run": "total",
    "run:scan_sample:SimContainer:run": "time series",
    "run:scan_sample:SimContainer:run:step": "timestep",
    "run:scan_sample:update_sim_container": "update sim_container",
    "run:write_element_tree": "write xml file",
    "run:scan_sample": "scan sample"

})

plotter.cell_df["activation"] = 0
plotter.calc_cell_activation()

metrics = {"min_distance": np.min, "mean_distance": np.mean, "max_distance": np.max}
plotter.compute_cell_distance_metric("sec", metric_dict=metrics)

plotter.scan_scale = plotter.global_df.v_v.unique()

f = lambda df: df
stat_max = lambda df: df.loc[df.activation > 0.5]

b = 8.3 / 2
a = np.sqrt(2) * b * 1.2

plotter.subplots(2, 3, external_legend="axes", figsize=(a, b))
plotter.filter = f

plotter.global_steady_state_plot("Concentration", "v_v", xlog=False, style=plotter.model_name_key)
# plotter.global_steady_state_plot("SD","v_v",xlog=False,style=True, legend="brief")

plotter.filter = lambda df: df.loc[(df.numeric_linear == True)]

plotter.cell_steady_state_plot("activation", x_name="v_v", hue="type_name", style=plotter.model_name_key,
                               xlog=False,
                               ci=None, xlim=[0, 1], legend="brief", categorical=False)
plotter.cell_plot("v_v", "activation", filter=lambda df: df.loc[(df.type_name == "abs")],
                  style=plotter.model_name_key,
                  xlog=False, hue="min_distance", palette_name="plasma", condition=stat_max, count=True,
                  ci=None, legend=None, xlim=[0, 1], categorical=False, overlay=True)


plotter.make_legend()
plotter.savefig(IMGPATH+"cell_plot.pdf")
plotter.show()
