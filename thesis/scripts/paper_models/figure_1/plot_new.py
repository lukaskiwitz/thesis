import os

import numpy as np

from parameters import IMGPATH, path
from thesis.main.MyPlotter import Plotter

if os.path.isdir("/extra2"):
    extra = "extra2"
else:
    extra = "extra"

# path = "/extra/lukas/paper_figure1_B/test_1"
# IMGPATH = path + "/images/"

groups = [
    "numeric_linear",
]

plotter = Plotter(path, groups=groups)

R_M = 860
ec50_min = 0
ec50_max = 125 * 1e-3
n_R = 0.55
n_il2 = 2

plotter.calc_cell_activation(
    R_M=R_M, max=ec50_max, min=ec50_min, n_R=n_R, n_il2=n_il2
)
plotter.scan_scale = {0: "L", 1: "NL"}

plotter.label_replacement.update({

    "kd": "$k_d$",
    plotter.time_key: "time (a.u)",
    "sec_amax": "$a_{max}$",
    "sec_q": "$q$",
    plotter.scan_index_key: "",
    "Concentration": "[IL-2] (nM)",
    "SD": "$\sigma$",
    "CV": "$C_v$",
    "IL-2_surf_c": "surface IL-2", #$\langle$ IL-2 $\\rangle_s$(nM)",
    "surf_c": "surface IL-2",#"$\langle$ IL-2 $\\rangle_s$(nM)",
    "surf_c_std": "$\sigma(\langle$ IL-2 $\\rangle_s)$(nM)",
    "surf_c_cv": "$C_v(\langle$ IL-2 $\\rangle_s$)",
    "Gradient": "$\\nabla $[IL-2]",
    "activation": "activation probability",
    "numeric_linear":"linear BCs",
})

plotter.color_dict = {
    "IL-2": "tab:red",
    "IL-2_surf_c": "tab:red",
    "sec": "tab:orange",
    "abs": "tab:green",
    "default": "tab:gray",

}


# plotter.calc_cell_activation(R_M=1e5)


def bar_plot(a, b):
    plotter.subplots(1, 7, figsize=(a, b), external_legend="figure")

    def get_ylim(df, columns):
        if not plotter.scan_index_key in columns:
            columns.append(plotter.scan_index_key)
        return [0, df.loc[:, columns].groupby([plotter.scan_index_key]).mean().max().max()]

    ylim = get_ylim(plotter.global_df, ["SD", "Concentration", "surf_c", "surf_c_std"])
    ylim[-1] *= 1.2

    plotter.global_steady_state_barplot(["Concentration"], x_name="numeric_linear", t_mean=True, ylim=ylim,
                                        y_ticks=True)
    plotter.global_steady_state_barplot(["SD"], x_name="numeric_linear", t_mean=True, ylim=ylim, y_ticks=False)
    plotter.global_steady_state_barplot(["surf_c"], x_name="numeric_linear", t_mean=True, ylim=ylim, y_ticks=False)
    plotter.global_steady_state_barplot(["surf_c_std"], x_name="numeric_linear", t_mean=True, ylim=ylim, y_ticks=False)

    ylim = get_ylim(plotter.global_df, ["Gradient"])
    ylim[-1] *= 1.2
    plotter.global_steady_state_barplot(["Gradient"], x_name="numeric_linear", t_mean=True, ylim=ylim, y_ticks=True)

    ylim = get_ylim(plotter.global_df, ["CV", "surf_c_cv"])
    ylim[-1] *= 1.2

    plotter.global_steady_state_barplot(["CV"], x_name="numeric_linear", t_mean=True, ylim=ylim, y_ticks=True)
    plotter.global_steady_state_barplot(["surf_c_cv"], x_name="numeric_linear", t_mean=True, ylim=ylim, y_ticks=False)

    plotter.savefig(IMGPATH + "overview_barplot.pdf")


def activation_function(a, b):
    plotter.subplots(1, 1, figsize=(a, b), external_legend="figure")


    filter_abs = lambda df: df.loc[(df["type_name"] == "abs")]
    filter_scan_index = lambda df: df.loc[
        (df["IL-2_bc_type"] == "R_saturation")
    ]

    plotter.filter = lambda df: filter_scan_index(filter_abs(df))

    plotter.cell_scatter_plot(["IL-2_surf_c", "activation"], hue="type_name", marker=".", s=0.01)
    plotter.gca().set_xlim([0, 0.1])

    plotter.make_legend()
    plotter.savefig(IMGPATH + "act_function.pdf")


def activation_historamm(a, b):
    plotter.subplots(2, 2, figsize=(a, b), external_legend="figure")

    for bc, linear in [("R_saturation", False), ("linear", True)]:

        filter_abs = lambda df: df.loc[(df["type_name"] == "abs")]
        filter_scan_index = lambda df: df.loc[
            (df["IL-2_bc_type"] == bc)
        ]

        c = plotter.global_df.loc[
            (plotter.global_df["numeric_linear"] == linear)

        ]["Concentration"].mean()

        # c = plotter.cell_df.loc[
        #     (plotter.cell_df["type_name"] == "abs") &
        #     (plotter.cell_df["IL-2_bc_type"] == bc)
        #     ]["IL-2_surf_c"].mean()

        cum_q = plotter.cell_df.loc[
            (plotter.cell_df["IL-2_bc_type"] == bc)
        ]["IL-2_q"].sum()



        R_abs = plotter.cell_df.loc[
            (plotter.cell_df["type_name"] == "abs") &
            (plotter.cell_df["IL-2_bc_type"] == bc)
            ]["IL-2_R"].mean()


        mean_activation = plotter.activation(c, R_abs, R_M=R_M, min=ec50_min, max=ec50_max, n_R=n_R, n_il2=n_il2)

        print("mean c: {c}, mean activation: {a}".format(c=c, a=mean_activation))
        plotter.filter = lambda df: filter_scan_index(filter_abs(df))

        xlim = [0, np.quantile(plotter.cell_df["IL-2_surf_c"],0.95)]
        ylim = [0, 50]

        bins = 500

        plotter.cell_histogramm("IL-2_surf_c", ylim=ylim, distplot_kwargs={"bins": bins},xlim = xlim)
        plotter.cell_activation_histogramm("IL-2_surf_c", ylim=[0, 0.2], bins=bins / 10, relative=True, twinx=True,
                                           color="red", xlim=xlim)

        plotter.cell_histogramm("IL-2_surf_c", ylim=ylim, distplot_kwargs={"bins": bins}, xlim=xlim)
        plotter.cell_activation_histogramm("IL-2_surf_c", ylim=[0, 1], bins=bins / 10, cummulative=True, relative=True,
                                           twinx=True, color="green", xlim=xlim)
        plotter.gca().plot(
            np.linspace(xlim[0], xlim[1], 100),
            100 * [mean_activation], "--"
        )

    plotter.make_legend()
    plotter.savefig(IMGPATH + "activation_histogramm.pdf")


def radial_niche_plot(a, b):

    plotter.subplots(1, 2, figsize=(a, b), external_legend=False)
    c_max = np.quantile(plotter.cell_df["IL-2_surf_c"],0.95)
    plotter.cell_radial_niche_plot("IL-2_surf_c", "sec",  xlim=[1,20],ylim=[0,c_max], cell_radius = 20,ci = "sd")
    plotter.cell_radial_niche_plot("activation", "sec", xlim=[1, 20],ylim=[0,1.2], cell_radius=20, ci ="sd")

    plotter.make_legend()
    plotter.savefig(IMGPATH+"niche_plot.pdf")
    plotter.show()



a = 8.3 * 0.75
b = np.sqrt(2) * a

# bar_plot(a, b / 5)
radial_niche_plot(a,b/4)
# activation_function(a,b/2)
# activation_historamm(a,b/2)

