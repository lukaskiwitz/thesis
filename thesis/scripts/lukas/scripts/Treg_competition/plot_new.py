import glob
import os

import matplotlib.pyplot as plt
import numpy as np

from parameters import path, scan_name
from thesis.main.MyPlotter import Plotter
from thesis.main.my_debug import message


# path = "/extra/kiwitz/Treg_competition_patricks_solver/150_3D/20211123_cluster_antigen_stimulus_response_2/"
# IMGPATH = path + "images/"


def do_plotting(path_list, mode=None, filter=lambda df: df):
    IMGPATH = "/".join(path_list[-1].split("/")[0:-1]) + "/images/"

    print(IMGPATH)
    # return None
    plotter = Plotter(path_list)
    plotter.cell_df[plotter.time_key] = plotter.cell_df[plotter.time_key] / 3600
    plotter.global_df[plotter.time_key] = plotter.global_df[plotter.time_key] / 3600
    plotter.t_max = plotter.global_df[plotter.time_key].max()
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
        plotter.scan_index_key: "q (molec/(cell*s))",
        "Concentration": "Concentration (nM)",
        "run": "total",
        "run:scan_sample:SimContainer:run": "time series",
        "run:scan_sample:SimContainer:run:step": "timestep",
        "run:scan_sample:update_sim_container": "update sim_container",
        "run:write_element_tree": "write xml file",
        "run:scan_sample": "scan sample"

    })
    activation_threshold = 3e3

    os.makedirs(IMGPATH, exist_ok=True)
    message("dataset loaded")
    a = 8.3
    b = np.sqrt(2) * a

    import seaborn as sns

    if mode == "steady_state":
        marker = None
        plotter.subplots(2, 3, figsize=(a, b / 4), external_legend="axes")
        plotter.filter = filter
        xlim = [0, None]

        plotter.global_steady_state_plot("Concentration", ylog=False, xlog=False, hue=plotter.scan_name_key, xlim=xlim,
                                         marker=marker)
        plotter.global_steady_state_plot("SD", ylog=False, xlog=False, hue=plotter.scan_name_key, xlim=xlim,
                                         marker=marker)
        plotter.global_steady_state_plot("CV", ylog=False, xlog=False, hue=plotter.scan_name_key, xlim=xlim,
                                         marker=marker)
        plotter.steady_state_count(style="type_name", legend="brief", xlog=False, hue=plotter.scan_name_key, xlim=xlim,
                                   marker=marker, ylim=[0, None])

        # plotter.cell_steady_state_plot(y_name="IL-2_surf_c", xlim = [0,80],
        #                                hue = plotter.scan_name_key, filter = lambda df:df.loc[(df[plotter.time_index_key] == 0) & (df["type_name"] == "naive")],
        #                                xlog=False, marker=marker)

        cf = lambda df: df.loc[(df.type_name.isin(["naive"])) & (df[plotter.time_key] == plotter.t_max)]
        df = plotter.cell_df

        df["activated"] = df["IL-2_R"] > activation_threshold

        df = cf(df)
        df_total = df.groupby([plotter.scan_index_key, plotter.time_index_key, plotter.time_key,
                               plotter.scan_name_key]).count().reset_index()

        def cf2(df):
            return len(df.loc[df.activated == True])

        df = df.groupby(
            [plotter.scan_index_key, plotter.time_index_key, plotter.time_key, plotter.scan_name_key]).apply(
            cf2).reset_index()

        df["IL-2_R"] = df[0] / df_total["IL-2_R"]
        ax, df, palette, hue = plotter.prepare_plot(df, hue=plotter.scan_name_key)
        sns.lineplot(x=plotter.scan_index_key, y="IL-2_R", data=df, ax=ax, palette=palette, legend=None, marker=marker,
                     hue=hue)
        plotter.make_legend_entry(ax)
        plotter.finalize_steady_state_plot(ax, "IL-2_R", None, xlim, False, False, x_name=plotter.scan_index_key)
        ax.set_ylabel("No of activated Th cells")
        ax.set_ylim([0, None])

        plotter.make_legend()
        plt.tight_layout()
        plt.savefig(IMGPATH + "global_collection.pdf")
        plt.show()

    if mode == "time_series":
        df = plotter.cell_df
        scan_indices = df.scan_value.unique()
        scan_indices.sort()
        df["activated"] = df["IL-2_R"] > activation_threshold
        cdf = df.groupby(["scan_value", plotter.time_key, "activated", "type_name", plotter.scan_name_key,
                          plotter.path_name]).count().reset_index()

        for scan_name in df[plotter.scan_name_key].unique():
            fig, ax = plt.subplots(len(scan_indices), 3, figsize=(a, b / 6 * len(scan_indices)))
            # fig,ax = plt.subplots(2,3,figsize = (a, b/6 * 2))
            for i, scan_value in enumerate(scan_indices):
                f1 = lambda df: df.loc[(df["scan_value"] == scan_value) & (df[plotter.scan_name_key] == scan_name)]
                f2 = lambda df: df.loc[(df["scan_value"] == scan_value) & (df[plotter.scan_name_key] == scan_name)]

                sns.lineplot(x=plotter.time_key, y="IL-2_surf_c", data=f1(plotter.cell_df), units="id", estimator=None,
                             linewidth=0.25, ax=ax[i][0], hue="type_name")
                sns.lineplot(x=plotter.time_key, y="IL-2_R", data=f2(plotter.cell_df), units="id", estimator=None,
                             linewidth=0.25, ax=ax[i][1], hue="type_name")

                cf = lambda df: df.loc[
                    (df["scan_value"] == scan_value) & (df.type_name.isin(["naive", "sec"])) & (
                                df.activated == True) & (
                            df[plotter.scan_name_key] == scan_name)]
                sns.lineplot(x=plotter.time_key, y="IL-2_R", data=cf(cdf), ax=ax[i][2], hue="type_name")
                ax[i][2].set_ylabel("Number of activated Th")

            plt.tight_layout()
            plt.title(scan_name)
            plt.savefig(IMGPATH + scan_name + "_kinetic_collection.pdf")
            plt.show()

    return plotter


path_list = glob.glob(os.path.join(
    "/".join(path.split("/")[:-2])
    , scan_name + "*"))

# all_paths = []
# for path in path_list:
#
#     refine_paths = glob.glob(os.path.join(path, "refine_*"))
#     refine_paths = [re.search(r"refine_(?P<n>\d)", rp) for rp in refine_paths]
#     refine_paths = [rp.string for rp in refine_paths if (int(rp.groupdict()["n"]) in [0, 1, 2, 3, 4, 5])]
#     all_paths = all_paths + refine_paths
#
#
#     plotter = do_plotting(refine_paths,  mode = "time_series")
#
#
# plotter = do_plotting(all_paths, mode = "steady_state")
#
# from parameters import path
# path = "/".join(path.split("/")[:-2])
# plotter.cell_df.to_hdf(os.path.join(path, "cell_df.h5"), key="df", mode="w")
# plotter.global_df.to_hdf(os.path.join(path, "global_df.h5"), key="df", mode="w")
#
# df = plotter.cell_df
# df = df.loc[df["time"] == plotter.t_max]
# df.to_hdf(os.path.join(path, "cell_df_ss.h5"), key="df", mode="w")


from parameters import path

path = "/".join(path.split("/")[:-2])

plotter = Plotter(path)
df = plotter.cell_df
df = df.loc[df["time"] == plotter.t_max]
df.to_hdf(os.path.join(path, "cell_df_ss.h5"), key="df", mode="w")
f = lambda df: df.loc[df.scan_name_scan_name.str.contains(r"cs_(0.0|1.0)_ec50_(1|3).0_k_860.0_ht_1.0")]

f = lambda df: df.loc[df.scan_name_scan_name.str.contains(r"cs_0.0")]
plotter = do_plotting([path], mode="steady_state", filter=f)

f = lambda df: df.loc[df.scan_name_scan_name.str.contains(r"cs_1.0")]
plotter = do_plotting([path], mode="steady_state", filter=f)
print("done")

# cdf = plotter.cell_df.groupby(["scan_name_scan_name","scan_index","scan_value","type_name","time_index"]).count().reset_index()
# cdf = f(cdf.loc[cdf.scan_value == cdf.scan_value.iloc[0]])
# sns.lineplot(x = "scan_value", y = "IL-2_D",data = cdf, hue = "type_name", style = "scan_name_scan_name")
# plt.show()

columns = [k for k, v in df.dtypes.items() if v == np.dtype('O')]
df[:, columns] = df.loc[:, columns] = df[columns].applymap(str)
