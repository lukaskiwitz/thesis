import glob
import os
from copy import deepcopy

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import seaborn as sns

from thesis.main.MyPlotter import Plotter
from thesis.main.my_debug import message


def do_plotting(path_list, mode=None, filter=lambda df: df, file_prefix=None, hue=None, plotter=None,
                palette_name="seismic"
                , color_dict=None, scan_ticks=None):
    b = 1.45
    a = 7 / 6 * b

    if plotter is None:
        plotter = Plotter(path_list, load_dataframes={"cell": True})

        plotter.cell_df[plotter.time_key] = plotter.cell_df[plotter.time_key] / 3600
        try:
            plotter.global_df[plotter.time_key] = plotter.global_df[plotter.time_key] / 3600
        except:
            pass
        """
        legend font 6
        label font 8
        A,B,C font 12
        title font 10 
        tick label 6
        
        """
        fs = 8
        rc_ticks = {
            "text.usetex": False,
            "font.size": fs,
            "legend.fontsize": 6,
            "legend.title_fontsize": 10,
            "xtick.labelsize": 6,
            "ytick.labelsize": 6,
            "axes.labelsize": 8,
            "axes.titlesize": 10,
            "lines.linewidth": 0.5,
            "axes.linewidth": 0.4,
            "lines.markersize": 3,
            "xtick.major.size": 2.5,
            "xtick.major.width": 0.5,
            "xtick.major.pad": 1,
            "xtick.minor.size": 1.5,
            "xtick.minor.width": 0.5,
            "xtick.minor.pad": 1,
            "ytick.major.size": 2.5,
            "ytick.major.width": 0.5,
            "ytick.major.pad": 1,
            "ytick.minor.size": 1.5,
            "ytick.minor.width": 0.5,
            "ytick.minor.pad": 1,
            "axes.labelpad": 1,
            "axes.titlepad": 1
        }
        # sns.set_context("paper")
        plt.rcParams.update(rc_ticks)

        plotter.t_max = plotter.cell_df[plotter.time_key].max()
        plotter.color_dict = {
            'naive': (0.9058823529411765, 0.1607843137254902, 0.5411764705882353),
            'sec': (0.4, 0.6509803921568628, 0.11764705882352941),
            'treg': (0.9019607843137255, 0.6705882352941176, 0.00784313725490196),
        }
        plotter.label_replacement.update({
            "D": "Diffusion coefficient",
            "sec_amax": "$a_{max}$",
            "sec_q": "IL-2 / (cell * s)",
            "f_sec": "% secretors",
            "f_abs": "% consumers",
            "abs_R": "IL-2R / cell",
            "Kc": "IL-2R EC50",
            "kd": "IL-2 decay",
            plotter.scan_index_key: "$\phi$",
            "fractions_treg": "% Treg",
            "fractions_sec": "% Tsec",
            "IL-2_surf_c": "surf. IL-2",
            "Concentration": "Concentration (nM)",
            "run": "total",
            "run:scan_sample:SimContainer:run": "time series",
            "run:scan_sample:SimContainer:run:step": "timestep",
            "run:scan_sample:update_sim_container": "update sim_container",
            "run:write_element_tree": "write xml file",
            "run:scan_sample": "scan sample"

        })
        activation_threshold = 3e3
        plotter.calc_cell_activation(
            R_M=plotter.cell_df["misc_pSTAT_k"],
            n_R=plotter.cell_df["misc_pSTAT_N"],
            n_il2=3
        )
        plotter.cell_df["IL-2_surf_c"] = plotter.cell_df["IL-2_surf_c"] * 1e3
    os.makedirs(IMGPATH, exist_ok=True)
    message("dataset loaded")

    if isinstance(color_dict, dict):
        message("color udpated applied")
        plotter.color_dict.update(color_dict)

    if scan_ticks is not None:
        message("scan ticks updated")
        plotter.scan_ticks.update(scan_ticks)

    confidence_intervall = "sem"
    if mode == "steady_state":

        marker = None
        plotter.subplots(3, 3, figsize=(3 * a, 3 * b), external_legend="axes")
        plotter.filter = filter
        xlim = [0, None]
        HUE = "fractions_treg" if hue is None else deepcopy(hue)

        f = lambda df: df.loc[df.type_name == "naive"]
        plotter.cell_steady_state_plot("IL-2_surf_c", ylog=False, xlog=False, hue=HUE, xlim=xlim,
                                       marker=marker, filter=f, ylim=[0, None], legend="brief",
                                       categorical_palette=True, palette_name=palette_name,
                                       ci=confidence_intervall)
        plotter.empty_plot()
        plotter.empty_plot()
        plotter.empty_plot()
        # plotter.steady_state_count(style="type_name", x_name="scan_value", legend=None, xlog=False,hue=HUE, xlim=xlim,
        #                            marker=marker, ylim=[0, None])

        plotter.empty_plot()
        plotter.empty_plot()

        cf = lambda df: df.loc[(df.type_name.isin(["naive"])) & (
                (df[plotter.time_key] == plotter.t_max) |
                (df[plotter.time_index_key] == 0)
        )]

        df = cf(plotter.cell_df)
        c = [
            plotter.time_index_key,
            plotter.scan_name_key,
            plotter.model_index_key,
            plotter.replicat_index_key,
            plotter.scan_index_key,
            "path_name",
            "path",
            "id_id",
            "field_name"
        ]
        c = [cc for cc in c if cc in df.columns and len(df[cc].unique()) > 1]
        ind = pd.MultiIndex.from_frame(df.loc[:, c])
        df = df.set_index(ind)
        df = df.drop(ind.names, axis=1)

        df_0 = df.loc[0]
        df_t = df.loc[plotter.time_index_max]

        df_0, df_t = df_0.loc[df_0.index.intersection(df_t.index)], df_t.loc[df_t.index.intersection(df_0.index)]

        fold_change = df_t["IL-2_R"] / df_0["IL-2_R"]

        df = df_t
        df["activated"] = (fold_change <= 1 / (0.5 * df_t["misc_gamma"])) | (fold_change >= 0.5 * df_t["misc_gamma"])

        df = df.reset_index()
        df_total = df.groupby(
            [plotter.scan_index_key, plotter.scan_name_key, plotter.replicat_index_key, HUE]).count().reset_index()

        def cf2(df):
            return len(df.loc[df.activated == True])

        df = df.groupby([plotter.scan_index_key, plotter.scan_name_key, plotter.replicat_index_key, HUE]).apply(
            cf2).reset_index()

        df["IL-2_R"] = (df[0] / df_total["IL-2_R"]) * 100
        ax, df, palette, hue = plotter.prepare_plot(df, hue=HUE)

        df, ci = plotter.compute_ci(df,
                                    [plotter.scan_index_key, plotter.replicat_index_key, plotter.model_index_key,
                                     plotter.time_index_key, hue, None],
                                    ci=confidence_intervall, estimator=None, y_names=["IL-2_R"])

        sns.lineplot(x=plotter.scan_index_key, y="IL-2_R", data=df, ax=ax, palette=palette, legend=None, marker=marker,
                     hue=hue, ci=ci)
        plotter.make_legend_entry(ax)
        plotter.finalize_steady_state_plot(ax, "IL-2_R", None, xlim, False, False, x_name=plotter.scan_index_key)
        ax.set_ylabel("% of activated Th cells")
        ax.set_ylim([0, None])

        ax, df, palette, hue = plotter.prepare_plot(plotter.cell_df, hue=HUE)
        df = df.loc[df.type_name == "naive"]
        df = df.loc[df.time_index == plotter.time_index_max]

        df_total = df.groupby([plotter.scan_index_key, plotter.time_index_key, plotter.time_key,
                               plotter.scan_name_key, plotter.replicat_index_key]).count().reset_index()

        df["activated"] = df.activation >= 0.5
        df = df.groupby([plotter.scan_index_key, plotter.time_index_key, plotter.time_key, plotter.scan_name_key,
                         plotter.replicat_index_key, HUE]).apply(cf2).reset_index()

        df["activated"] = 100 * df[0] / df_total["activation"]
        # df["activated"] = df["activation"]

        df, ci = plotter.compute_ci(df,
                                    [plotter.scan_index_key, plotter.replicat_index_key, plotter.model_index_key,
                                     plotter.time_index_key, hue, None],
                                    ci=confidence_intervall, estimator=None, y_names=["activated"])

        sns.lineplot(x=plotter.scan_index_key, y="activated", data=df, ax=ax, palette=palette, legend=None,
                     marker=marker,
                     hue=hue, ci=ci, )

        plotter.make_legend_entry(ax)
        plotter.finalize_steady_state_plot(ax, "activation", None, xlim, False, False, x_name=plotter.scan_index_key)
        ax.set_ylabel("% pSTAT5$^+$")
        ax.set_ylim([0, None])

        ax, df, palette, hue = plotter.prepare_plot(plotter.cell_df, hue=HUE)
        df = df.loc[df.type_name == "naive"]

        df = df.groupby([plotter.scan_index_key, plotter.time_index_key, plotter.time_key, plotter.scan_name_key, hue],
                        as_index=False).mean()
        df, ci = plotter.compute_ci(df,
                                    [plotter.scan_index_key, plotter.replicat_index_key, plotter.model_index_key,
                                     plotter.time_index_key, hue, None],
                                    ci=confidence_intervall, estimator=None, y_names=["IL-2_R"])
        sns.lineplot(x=plotter.scan_index_key, y="IL-2_R", data=df, ax=ax, palette=palette, legend=None,
                     marker=marker,
                     hue=hue, ci=ci)

        plotter.make_legend_entry(ax)
        plotter.finalize_steady_state_plot(ax, "IL-2_R", None, xlim, False, False, x_name=plotter.scan_index_key)
        ax.set_ylabel("IL-2R on Th-cells (molec/cell)")
        # ax.set_ylim([0, None])

        plotter.make_legend()
        plt.tight_layout()
        if file_prefix is not None:
            plt.savefig(IMGPATH + str(file_prefix) + "_global_collection.pdf")
        else:
            plt.savefig(IMGPATH + "global_collection.pdf")
        plt.show()

    if mode == "time_series":

        df = plotter.cell_df
        scan_indices = df.scan_value.unique()
        scan_indices.sort()
        df["activated"] = df["IL-2_R"] > activation_threshold
        cdf = df.groupby([plotter.time_key, "activated", "type_name", plotter.scan_name_key,
                          plotter.path_name]).count().reset_index()

        for scan_name in df[plotter.scan_name_key].unique():
            fig, ax = plt.subplots(len(scan_indices), 3, figsize=(a, b / 6 * len(scan_indices)))
            # fig,ax = plt.subplots(2,3,figsize = (a, b/6 * 2))
            if len(ax.shape) == 1:
                ax = np.expand_dims(ax, axis=1).T

            for i, scan_value in enumerate(scan_indices):
                f1 = lambda df: df.loc[
                    (df[plotter.scan_index_key] == scan_value) & (df[plotter.scan_name_key] == scan_name)]
                f2 = lambda df: df.loc[
                    (df[plotter.scan_index_key] == scan_value) & (df[plotter.scan_name_key] == scan_name)]

                sns.lineplot(x=plotter.time_key, y="IL-2_surf_c", data=f1(plotter.cell_df), units="id", estimator=None,
                             linewidth=0.25, ax=ax[i][0], hue="type_name")
                sns.lineplot(x=plotter.time_key, y="IL-2_R", data=f2(plotter.cell_df), units="id", estimator=None,
                             linewidth=0.25, ax=ax[i][1], hue="type_name")

                cf = lambda df: df.loc[
                    (df[plotter.scan_index_key] == scan_value) & (df.type_name.isin(["naive", "sec"])) & (
                                df.activated == True) & (
                            df[plotter.scan_name_key] == scan_name)]
                sns.lineplot(x=plotter.time_key, y="IL-2_R", data=cf(cdf), ax=ax[i][2], hue="type_name")
                ax[i][2].set_ylabel("Number of activated Th")
                ax[i][0].set_ylim([0, 1])

                ax[i][1].set_yscale("log")

            plt.tight_layout()
            plt.title(scan_name)
            if file_prefix is not None:
                plt.savefig(IMGPATH + str(file_prefix) + scan_name + "_kinetic_collection.pdf")
            else:
                plt.savefig(IMGPATH + scan_name + "_kinetic_collection.pdf")
            plt.show()

    return plotter



# path_list = [p + "/" for p in glob.glob("/".join(path.split("/")[:-1]) + "*/*")]

# path_list = [path]

# plotter = Plotter(path_list)
# try:
#     plotter = Plotter(path_list)
#     # print(plotter.cell_df)
#     plotter.calc_cell_activation(
#                 R_M = plotter.cell_df["misc_pSTAT_k"],
#                 n_R = plotter.cell_df["misc_pSTAT_N"],
#                 n_il2=3,
#             )
#
#     n = len(plotter.cell_df[plotter.scan_name_key].unique())
#     m = int(np.ceil(n/2))
#     a = 8.3
#     b = np.sqrt(2)*8.3
#
#     plotter.subplots(m,2, external_legend = "axes",figsize = (a/4 * m, b/6*2))
#     f = lambda df: df.loc[df.scan_name_scan_name.str.contains(r"fsec")]
#     for i, scan_name in enumerate(f(plotter.cell_df).scan_name_scan_name.unique()):
#
#
#         type_filter = lambda df: df.loc[df.type_name.isin(["naive","sec"])]
#         scan_filter = lambda df : df.loc[df.scan_name_scan_name == scan_name]
#         # scan_filter = lambda df: df.loc[df.scan_name_scan_name.str.contains(r"cs_1_ec50_1_k_860_ht_1_treg_0_sg_2.5")]
#         f = lambda df: type_filter(scan_filter(df))
#
#         legend = "brief" if i == 0 else None
#
#         plotter.cell_radial_niche_plot(
#             "activation",
#             center_type="sec",
#             ylim=[0,1],
#             xlim = [0,200],
#             hue = "scan_value",
#             style = "type_name",
#             legend =legend,
#             filter = f,
#             plot_filter = lambda df: df.loc[df.type_name.isin(["naive"])],
#             additional_categories = ["fractions_sec"],
#             ci = None
#         )
#         plt.gca().set_title(scan_name)
#     plotter.make_legend()
#     plotter.savefig(IMGPATH+"radials.pdf")
#     plotter.show()
# except:
#     pass
plotter = None

from matplotlib import cm as color_map

cm1 = color_map.get_cmap('hot', 5)
cm2 = color_map.get_cmap('summer', 5)

color_dict = {
    "fractions_treg": {
        0.025: cm1(0),
        0.05: cm1(1),
        0.2: cm1(2),
        0.0: (0, 0.0, 0.0, 0.25)
    },
    'fractions_sec': {
        0.01: (0, 0.0, 0.0, 0.25),
        0.05: cm2(0),
        0.1: cm2(1),
        0.2: cm2(2)
    }
}
scan_ticks = {
    0: 0,
    0.25: 0.25,
    0.5: 0.5,
    0.75: 0.75,
    0.9: 0.9
}

# f = lambda a,x: -a*np.power(x,2)
# x = np.linspace(0,2,100)
# fig,ax = plt.subplots(2,2)
# ax= np.ravel(ax)
#
# for i,a in enumerate([0,0.025,0.05,0.2]):
#     ax[0].plot(x,f(a,x), color = color_dict["fractions_treg"][a],linewidth = 2)
# for i,a in enumerate([0.01,0.05,0.1,0.2]):
#     ax[1].plot(x,f(a,x), color = color_dict["fractions_sec"][a],linewidth = 2)
# plt.show()

from parameters import path
import os

# path = os.path.join(path,"../dataframes_IL-7_testing_steady_state/")
path = os.path.join(path,"../dataframes_multiple_apc_dry_run_steady_state/")

IMGPATH = os.path.join(path, "images/")
os.makedirs(IMGPATH, exist_ok=True)

g = lambda df: df.loc[
    (df["fractions_treg"].isin([0, 0.025, 0.05, 0.2]))
    & (df.scan_value.isin(df.scan_value.unique()[:-1]))
    ]

f = lambda df: g(df.loc[df.scan_name_scan_name.str.contains(r"treg_cs_strength_scan")])
plotter = do_plotting(path, mode="steady_state", filter=f, hue="fractions_treg", file_prefix="treg_scan",
                      color_dict=color_dict, scan_ticks=scan_ticks)


g = lambda df: df.loc[
    (df["fractions_sec"].isin([0.01,0.05,0.1,0.2]))
    & (df.scan_value.isin(df.scan_value.unique()[:-1]))
    ]

f = lambda df: g(df.loc[df.scan_name_scan_name.str.contains(r"sec_cs_strength_scan")])
plotter = do_plotting(path, mode="steady_state", filter=f, hue="fractions_sec", file_prefix="tsec_scan",
                      color_dict=color_dict, scan_ticks=scan_ticks)

# path = "/extra/kiwitz/Treg_competition_parallel_patrick_parameters/200_3D/dataframes_20220227_run_replicats_new_parameters_box_replicats_tsec_cluster_steady_state/"
# IMGPATH = os.path.join(path,"images/")
# os.makedirs(IMGPATH,exist_ok=True)
#
# f = lambda df: df.loc[(df.scan_value.isin(df.scan_value.unique()[:-1]))]
# plotter = do_plotting(path, mode = "steady_state", filter = f, hue = "fractions_sec",color_dict=color_dict, scan_ticks =scan_ticks)
#
