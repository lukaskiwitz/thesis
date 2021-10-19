import os
from typing import List, Dict

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import seaborn as sns
from matplotlib.patches import Rectangle
from scipy.spatial import distance_matrix

from thesis.main.my_debug import warning


class Plotter:

    def __init__(self, path, groups=[]) -> None:

        self.main_title = ""
        self.time_key: str = "time"
        self.time_index_key: str = "time_index"
        self.replicat_index_key: str = "replicat_index"

        self.n = 2
        self.m = 2
        self.fig = None
        self.gridspec = None

        self.t_max = 10
        self.model_index_key: str = "model_index"
        self.model_name_key: str = "model_name"
        self.scan_name_key: str = "scan_name_scan_name"
        self.scan_index_key: str = "scan_value"
        self.path_name: str = "path_name"
        self.scan_ticks = []
        self.legend_axes = None
        self.legend_figure = None

        self.external_legend = True
        self.legend_entries = {}
        self.max_scans = 5
        self.gridspec_index = 0

        self.scan_scale: List = []
        self.groups = groups
        self.color_dict: Dict = {}
        self.style_dict = {
        }

        self.global_df: pd.DataFrame = pd.DataFrame()
        self.cell_df: pd.DataFrame = pd.DataFrame()
        self.timing_df: pd.DataFrame = pd.DataFrame()
        self.ruse_df: pd.DataFrame = None
        self.means: pd.DataFrame = pd.DataFrame()
        self.counts: pd.DataFrame = pd.DataFrame()

        self.filter = lambda df: df

        self.label_replacement = {
            "type_name": "Cell Type",
            "field_name": "Cytokine",
            self.scan_name_key: "Scan Name",
            "time": "Time",
            "time_index": "Time",
            "scan_index": "parameter fold-change"
        }
        self.tight_layout_settings = {
            "pad": 2,
        }
        self.rc = {'font.size': 7,
                   'axes.labelsize': 7.0,
                   'axes.titlesize': 7.0,
                   'xtick.labelsize': 6.0,
                   'ytick.labelsize': 6.0,
                   'legend.fontsize': 7,
                   'axes.linewidth': 0.5,
                   'grid.linewidth': 0.5,
                   'lines.linewidth': 1.0,
                   'lines.markersize': 5.0,
                   'patch.linewidth': 0.8,
                   'xtick.major.width': 1,
                   'ytick.major.width': 1,
                   'xtick.minor.width': 0.5,
                   'ytick.minor.width': 0.5,
                   'xtick.major.size': 2,
                   'ytick.major.size': 2,
                   'xtick.minor.size': 1,
                   'ytick.minor.size': 2,
                   'axes.formatter.usemathtext': True
                   }

        sns.set_context("paper", rc=self.rc)
        self.load(path, groups=groups)

        self._prepare_color_dict()

        self.time_index_max, self.t_max = self.get_max_time_index()

    def update_rc(self, rc):

        self.rc.update(rc)

        sns.set_context("paper", rc=self.rc)

    def subplots(self, n, m, figsize=(10, 5), external_legend="axes", gridspec_args=None, reset_filter=True) -> None:

        if self.fig is not None:
            plt.close(self.fig)

        if self.legend_figure is not None:
            plt.close(self.legend_figure)

        if reset_filter:
            self.filter = lambda df: df

        if self.fig:
            plt.close()

        if gridspec_args is None:
            gridspec_args = {}

        self.legend_entries = {}

        self.n = n
        self.m = m
        self.gridspec_index = 0

        from typing import Tuple
        if external_legend == "axes":
            self.fig = plt.figure(figsize=((figsize[0] / m) * (m + 1), figsize[1]))
            self.gridspec = self.fig.add_gridspec(n, m + 1, **gridspec_args)
            self.legend_axes = self.fig.add_subplot(self.gridspec[:, -1])
        elif external_legend == "figure" or isinstance(external_legend, Tuple):
            self.fig = plt.figure(figsize=figsize)
            if isinstance(external_legend, Tuple):
                self.legend_figure = plt.figure(figsize=external_legend)
                external_legend = "figure"
            else:
                self.legend_figure = plt.figure(figsize=(figsize[0] / m, figsize[1] / n))
            self.gridspec = self.fig.add_gridspec(n, m, **gridspec_args)
        else:
            self.fig = plt.figure(figsize=figsize)
            self.legend_figure = None
            self.gridspec = self.fig.add_gridspec(n, m, **gridspec_args)

        self.external_legend = external_legend

    def gca(self):
        return self.fig.gca()

    def gcf(self):
        return self.fig

    def get_subplot_axes(self, overlay=False, gs_slice=None):

        if overlay is False and self.gridspec_index + 1 > (self.n * self.m):
            raise Exception

        o = np.mod(self.gridspec_index, self.m)
        i = np.floor_divide(self.gridspec_index, self.m)

        if not overlay:
            self.gridspec_index += 1

        if gs_slice is None:
            subplot = self.fig.add_subplot(self.gridspec[i, o])
        else:
            subplot = self.fig.add_subplot(self.gridspec[gs_slice])
        return subplot

    def make_legend(self) -> None:

        if self.external_legend == "axes":
            handles = list(self.legend_entries.keys())
            labels = list(self.legend_entries.values())
            labels = self.replace_labels(labels)

            self.legend_axes.legend(handles, labels, loc="center left")
            self.legend_axes.axis('off')

        elif self.external_legend == "figure":
            ax = self.legend_figure.gca()
            handles = list(self.legend_entries.keys())
            labels = list(self.legend_entries.values())
            labels = self.replace_labels(labels)

            ax.legend(handles, labels, loc="center left")
            ax.axis('off')

    def _prepare_color_dict(self) -> None:
        fields = self.global_df["field_name"].unique()
        cell_types = self.cell_df["type_name"].unique()
        t = self.cell_df[self.time_key].unique()

        try:
            scan_names = self.global_df[self.scan_name_key].unique()
        except:
            self.color_dict = self.get_color_dict(np.concatenate([fields, t, cell_types]))
        else:
            self.color_dict = self.get_color_dict(np.concatenate([fields, t, cell_types, scan_names]))

    def get_categorical_color(self, kvp, palette_name="Dark2") -> Dict:

        key, value = kvp

        if kvp in self.color_dict.keys():
            return self.color_dict[kvp]
        else:
            keys = list(self.color_dict.keys())
            keys.append(kvp)
            self.color_dict = self.get_color_dict(keys, palette_name=palette_name)
            return self.color_dict[kvp]

    def get_continuous_color(self, key, palette_name="viridids"):

        if key in self.color_dict.keys():
            return self.color_dict[key]
        else:
            self.color_dict[key] = palette_name
            return palette_name

    def get_palette(self, df, key, palette_name="Dark2", categorical=None) -> Dict:

        if key is None:
            return None

        values = df[key].unique()

        if len(values) > 5 or (key in self.color_dict.keys() and isinstance(self.color_dict[key], str)):
            return self.get_continuous_color(key, palette_name=palette_name)
        else:
            p = {}
            for v in values:
                p[v] = self.get_categorical_color((key, v), palette_name=palette_name)

        return p

    @staticmethod
    def get_color_dict(keys, palette_name="Dark2") -> Dict:

        palette = sns.color_palette(palette_name, len(keys))
        color_dict = {}
        for i, k in enumerate(keys):
            color_dict[k] = palette[i]
        return color_dict

    def activation(self, c, R, R_M=860, max=0.125, min=0, n_R=0.55, n_il2=4):

        def rec(R):
            ec50 = (max - min) * (1 - np.power(R, n_R) / (np.power(R_M, n_R) + np.power(R, n_R))) + min
            return ec50

        ec50 = rec(R)
        a = c ** n_il2 / (ec50 ** n_il2 + c ** n_il2)

        if isinstance(R, float) and R == 0:
            return 0
        if isinstance(R, float) and isinstance(c, float):
            return a
        else:
            a[R == 0] = 0
            return a

    def calc_cell_activation(self, R_M=860, max=0.125, min=0, n_R=0.55, n_il2=4):

        # mask = self.cell_df["scan_name_scan_name"] == name

        act = self.activation(
            self.cell_df["IL-2_surf_c"], self.cell_df["IL-2_R"],
            R_M=R_M, max=max, min=min, n_R=n_R, n_il2=n_il2)

        self.cell_df["activation"] = act

    def load(self, path, groups=[]) -> None:

        if isinstance(path, str):
            path = [path]

        acc = {}
        for p in path:

            result = self.load_single_sim(p, groups)

            for k, v in result.items():
                name_series = pd.Series([p.split("/")[-1]] * len(v))
                v[self.path_name] = name_series

                if k in acc.keys():
                    acc[k] = acc[k].append(v)
                else:
                    acc[k] = v

        for k, v in acc.items():
            v.index = pd.RangeIndex(0, len(v))
            self.__setattr__(k, v)

    def load_single_sim(self, path, groups) -> None:

        assert os.path.exists(path)
        ruse_df = pd.DataFrame()
        try:
            global_df: pd.DataFrame = pd.read_hdf(os.path.join(path, "global_df.h5"))
            global_df = self.reset_scan_index(global_df)

            cell_df: pd.DataFrame = pd.read_hdf(os.path.join(path, "cell_df.h5"), mode="r")
            if os.path.exists(os.path.join(path, "cell_constants_df.h5")):
                cell_constants: pd.DataFrame = pd.read_hdf(os.path.join(path, "cell_constants_df.h5"), mode="r")
            else:
                cell_constants = pd.DataFrame()

            self.cell_constants = cell_constants

            try:
                timing_df: pd.DataFrame = pd.read_hdf(os.path.join(path, "timing_df.h5"), mode="r")
                timing_df = self.reset_scan_index(timing_df)
            except FileNotFoundError:
                pass
            except KeyError:
                pass

            try:
                ruse_df: pd.DataFrame = pd.read_hdf(os.path.join(path, "records/ruse.h5"), mode="r")
                ruse_df = self.reset_scan_index(ruse_df)
            except FileNotFoundError:
                pass

            # self.timing_df = timing_df
        except FileNotFoundError as e:
            warning("{df} dataframe was not found".format(df=str(e)))

        if self.scan_index_key not in global_df.columns:
            self.scan_index_key = "scan_index"

        groups = groups + [
            "type_name",
            self.time_key,
            self.time_index_key,
            self.scan_index_key,
            self.model_index_key,
            self.model_name_key,
            # "IL-2_surf_c",
            # "x",
            # "y",
            # "z",
            # "id_id",
            # "id"
        ]
        if self.scan_name_key in global_df.columns:
            groups.append(self.scan_name_key)

        for g in groups:
            if g in cell_constants:
                cell_df[g] = cell_constants[g]
        cell_df = self.reset_scan_index(cell_df)

        cell_df["id"] = cell_df["id"].astype("int")
        cell_df["id_id"] = cell_df["id_id"].astype("int")

        cell_df["x"] = cell_df["x"].astype("float")
        cell_df["y"] = cell_df["y"].astype("float")
        cell_df["z"] = cell_df["z"].astype("float")

        if self.scan_name_key in list(cell_df.columns) and self.scan_name_key not in groups:
            groups.append(self.scan_name_key)

        grouped_cells = cell_df.groupby(groups, as_index=True)

        means = grouped_cells.mean()
        means.reset_index(inplace=True)
        # self.means = means

        counts = grouped_cells.count()
        counts.reset_index(inplace=True)
        counts["n"] = counts["id"]
        # counts = counts.drop(columns=counts.columns.drop(["n"]))

        # groups = groups + [
        #     "IL-2_surf_c",
        #     "x",
        #     "y",
        #     "z",
        #     "id_id",
        #     "id"
        # ]

        groups.remove("type_name")

        total = pd.Series(np.zeros(len(counts)))

        for i, g in counts.groupby(groups):
            n = g["n"].sum()
            for o in g.index:
                total.iloc[o] = n

        counts.reset_index(inplace=True)
        counts["n_rel"] = counts["n"] / total
        # self.counts = counts

        return {"global_df": global_df, "cell_df": cell_df, "means": means, "timing_df": timing_df, "counts": counts,
                "ruse": ruse_df}

    def get_max_time_index(self) -> float:

        g = self.global_df[self.time_index_key].max()
        g = [g, np.unique(self.global_df.loc[self.global_df[self.time_index_key] == g]["time"])[0]]

        c = self.cell_df[self.time_index_key].max()
        c = [c, np.unique(self.cell_df.loc[self.cell_df[self.time_index_key] == c]["time"])[0]]

        return np.max([g, c], axis=0)

    def get_label(self, key) -> str:

        if key in self.label_replacement.keys():
            return self.label_replacement[key]
        else:
            return key

    # noinspection PyUnusedLocal
    @staticmethod
    def format_scan_index_ticklabels(scan_scale, distance, round_n) -> List:
        return scan_scale

    def get_scan_ticks(self) -> (List, List):

        scan_axis_name = self.scan_index_key

        axis_df = self.global_df.loc[self.global_df[scan_axis_name].notna()]
        ticks = axis_df[self.scan_index_key].unique()

        labels = axis_df[scan_axis_name]
        return np.array(ticks), np.array(labels)

        # scale = self.scan_scale
        # ticks = range(len(self.scan_scale))
        # labels = []  # np.round(self.scan_scale,2)
        #
        # if isinstance(scale, Dict):
        #     return list(scale.keys()), list(scale.values())
        # else:
        #     # noinspection PyUnusedLocal
        #     for i, e in enumerate(scale):
        #         labels.append("")
        #
        #     if len(self.scan_ticks) > 0:
        #
        #         for i, tick in enumerate(self.scan_ticks):
        #             arg = np.argmin(np.abs(tick - scale))
        #             labels[arg] = np.round(self.scan_scale, 2)[arg]
        #
        #     else:
        #         for i, e in enumerate(scale):
        #             if i % 10 == 0:
        #                 labels[i] = np.round(self.scan_scale, 2)[i]

        # return ticks, labels

    def replace_labels(self, labels) -> Dict:

        for i, l in enumerate(labels):
            if l in list(self.label_replacement.keys()):
                labels[i] = self.label_replacement[l]
        return labels

    def format_scan_index(self, index_column) -> List:
        return 1 * np.round(self.scan_scale[index_column.apply(lambda x: int(x))], 1)

    def reduce_df(self, df, index_name) -> pd.DataFrame:

        indices = df[index_name].unique()
        if len(indices) <= 1:
            return df

        e = int(len(indices) / self.max_scans)
        e = e if e >= 1 else 1
        indices = indices[0::e]
        result = pd.DataFrame()
        for i in indices:
            result = result.append(df.loc[(df[index_name] == i)])

        return result

    def reset_scan_index(self, df) -> pd.DataFrame:

        df["raw_scan_index"] = df["scan_index"]

        if "scan_name" in list(df.columns):
            df[self.scan_name_key] = df["scan_name"]

        if not self.scan_name_key in list(df.columns):
            return df
        else:
            offsets = {}
            for scan_name in df[self.scan_name_key].unique():
                offsets[scan_name] = df.loc[df[self.scan_name_key] == scan_name]["scan_index"].min()

            for i, o in offsets.items():
                mask = (df[self.scan_name_key] == i)
                df.loc[mask, "scan_index"] = df["scan_index"] - o

            return df

    def make_legend_entry(self, ax) -> None:

        if ax.get_legend() is not None:
            handles = ax.get_legend().legendHandles
            labels = [i._text for i in ax.get_legend().texts]
            for i, l in enumerate(labels):
                try:
                    labels[i] = np.round(float(l), 2)
                except ValueError:
                    continue




        else:
            handles, labels = ax.get_legend_handles_labels()

        if self.external_legend:
            try:
                ax.get_legend().remove()
                for i, h in enumerate(handles):
                    self.legend_entries[h] = labels[i]
            except:
                pass
        elif len(handles) > 0:
            labels = self.replace_labels(labels)
            ax.legend(handles, labels)

    def savefig(self, path) -> None:

        import os
        from functools import reduce

        os.makedirs(os.path.split(path)[0], exist_ok=True)
        self.fig.tight_layout(**self.tight_layout_settings)
        self.fig.suptitle(self.main_title, fontsize=16)
        self.fig.savefig(path, dpi=300)
        if self.legend_figure and self.external_legend == "figure":
            s = path.split(".")
            legend_path = reduce(lambda x, y: x + y, s[0:-1]) + "_legend." + s[-1]
            self.legend_figure.savefig(legend_path)

    def show(self) -> None:
        self.fig.tight_layout(**self.tight_layout_settings)
        self.fig.suptitle(self.main_title, fontsize=16)
        if self.legend_figure and self.external_legend == "figure":
            self.legend_figure.show()
        self.fig.show()

    def prepare_plot(self, df, hue, reduce=False, **kwargs):

        if self.filter:
            try:
                df = self.filter(df)
            except KeyError as e:
                warning("Key: {k} was not found in dataframe. Could not apply filter.".format(k=str(e)[1:-1]))
        if "filter" in kwargs:
            try:
                df = kwargs["filter"](df)
            except KeyError as e:
                warning("Key: {k} was not found in dataframe. Could not apply filter.".format(k=str(e)[1:-1]))

        if "select" in kwargs and isinstance(kwargs["select"], Dict):
            for k, v in kwargs["select"].items():
                df = df.loc[df[k].isin(v)]

        if isinstance(hue, Dict):
            hue_values = list(hue.values())[0]
            hue = list(hue.keys())[0]
            df = df.loc[df[hue].isin(hue_values)]

        if "twinx" in kwargs and kwargs["twinx"]:
            ax = self.fig.gca().twinx()
        else:
            ax = self.get_subplot_axes(**split_kwargs(kwargs, ["overlay", "gs_slice"]))
        if reduce:
            df = self.reduce_df(df, self.scan_index_key)
        palette = self.get_palette(df, hue, **split_kwargs(kwargs, ["palette_name"]))

        if "subtitle" in kwargs:
            ax.set_title(kwargs["subtitle"])
        return ax, df, palette, hue

    def prepare_twinx_plot(self, df, y_names, **kwargs):

        ax = self.get_subplot_axes()
        ax = [ax, ax.twinx()]
        df = self.reduce_df(df, self.scan_index_key)

        palette = {}
        for y in y_names:
            palette[y] = self.get_categorical_color(y)

        return ax, df, palette

    def empty_plot(self):
        self.get_subplot_axes()
        return None

    def global_time_series_plot(self, y_name, legend=False, hue=None, style=None, ylim=None, **kwargs) -> None:

        ax, df, palette, hue = self.prepare_plot(self.global_df, hue, reduce=True, **kwargs)

        sns.lineplot(x=self.time_key, y=y_name, data=df, hue=hue, style=style,
                     ax=ax,
                     legend=legend,
                     palette=palette)

        self.make_legend_entry(ax)

        ax.set_xlabel(self.get_label(self.time_key))
        ax.set_xlim([0, self.t_max])
        ax.set_ylim(ylim)
        ax.set_ylabel(self.get_label(y_name))

    def global_steady_state_plot(self, y_name, x_name=None, legend=False, ci="sd", hue=None, style=None, ylog=False,
                                 xlog=True, ylim=None, average=False, estimator=None, dashes=True,
                                 **kwargs) -> None:

        ax, df, palette, hue = self.prepare_plot(self.global_df, hue, **kwargs)

        df, ci = self.compute_ci(
            df,
            [self.scan_index_key, self.time_index_key, hue, style],
            ci=ci, estimator=estimator, y_names=[y_name])

        if x_name is None:
            x_name = self.scan_index_key

        if not average:
            df = df.loc[df[self.time_key] == self.t_max]

        # for i,sdf in df.groupby(hue, as_index=True):
        #     mean = sdf.groupby([x_name],as_index = True).mean()
        #     mean = mean.reset_index()
        #     sd = sdf.groupby([x_name],as_index = True).std()
        #     sd = sd.reset_index()
        #     ax.errorbar(mean[x_name], mean[y_name], yerr = sd[y_name], color = palette[i])

        if style in df.columns:
            sns.lineplot(x=x_name, y=y_name, data=df, hue=hue, ax=ax, legend=legend, palette=palette,
                         style=style, ci=ci, dashes=dashes)
        else:
            sns.lineplot(x=x_name, y=y_name, data=df, hue=hue, ax=ax, legend=legend, palette=palette,
                         ci=ci, dashes=dashes)

        self.finalize_steady_state_plot(ax, y_name, ylim, ylog, xlog, x_name=x_name)

    def global_steady_state_barplot(self, y_names, x_name=None, legend=False, hue=None, ylim=None, bar_spacing=1.1,
                                    cat_spacing=1.2,
                                    barwidth=0.1, norm=False, y_ticks=True, t_mean=False, **kwargs) -> None:

        ax, df, palette, hue = self.prepare_plot(self.global_df, hue, **kwargs)
        if t_mean == False:
            df = df.loc[df[self.time_key] == self.t_max]

        if not isinstance(y_names, List):
            y_names = [y_names]

        x_ticks = []
        scan_ticks = []
        x = [0]
        y = []
        w = []

        for y_name in y_names:

            spacing = cat_spacing

            if x_name is None:
                x_values = self.scan_scale
            else:
                x_values = {}
                for v in df[x_name].unique():
                    x_values[v] = v

            for k, v in x_values.items():
                x.append(x[-1] + barwidth * spacing)
                scan_ticks.append(v)
                spacing = bar_spacing
                w.append(barwidth)
                m = df[y_name].max() if norm else 1
                if x_name is None:
                    y.append(
                        df.loc[df[self.scan_index_key] == k][y_name].iloc[0] / m
                    )
                else:
                    y.append(
                        df.loc[df[x_name] == k][y_name].iloc[0] / m
                    )

            t = np.mean(x[-len(self.scan_scale):])
            x_ticks.append(t)
        x = x[1:]

        ax.bar(x, height=y, width=w)
        ax.set_xticks(x)
        ax.set_xticklabels(scan_ticks)

        if not y_ticks:
            ax.set_yticklabels([])

        ax2 = ax.twiny()
        ax2.set_xlim(ax.get_xlim())
        ax2.set_xticks(x_ticks)
        ax2.set_xticklabels(self.replace_labels(y_names))

        ax2.xaxis.set_ticks_position("top")
        ax2.xaxis.set_label_position("top")
        ax2.spines["top"].set_position(("axes", 1))
        ax2.set_frame_on(False)
        ax2.tick_params(length=0)

        if ylim:
            ax.set_ylim(ylim)

    def finalize_steady_state_plot(self, ax, y_name, ylim, ylog, xlog, x_name=None):

        self.make_legend_entry(ax)
        ax.set_xlabel(self.get_label(self.scan_index_key))
        ax.set_ylabel(self.get_label(y_name))
        if ylog:
            ax.set_yscale("log")

        if x_name is None:
            ticks, labels = self.get_scan_ticks()
            ax.set_xticks(ticks)
            ax.set_xticklabels(labels)
        else:
            ax.set_xlabel(self.get_label(x_name))

        from matplotlib.ticker import LogLocator, AutoLocator, \
            ScalarFormatter

        if xlog:
            ax.xaxis.set_major_locator(LogLocator())
            ax.set_xscale("log")
        else:
            ax.xaxis.set_major_locator(AutoLocator())

        ax.xaxis.set_major_formatter(ScalarFormatter())

        if ylim:
            ax.set_ylim(ylim)

    def cell_steady_state_plot(self, y_name, x_name=None, legend=False, hue=None, style=None, ylog=False, xlog=True,
                               cummulative=False, estimator=None,
                               ylim=None, ci="sd", dashes=True, **kwargs) -> None:
        ax, df, palette, hue = self.prepare_plot(self.cell_df, hue, **kwargs)

        df, ci = self.compute_ci(
            df,
            [self.scan_index_key, self.replicat_index_key, self.model_index_key, self.time_index_key, hue, style],
            ci=ci, estimator=estimator, y_names=[y_name])

        if cummulative:
            if hue and style:
                df = df.groupby([self.time_index_key, self.scan_index_key, hue, style]).sum()
            elif hue:
                df = df.groupby([self.time_index_key, self.scan_index_key, hue]).sum()
            elif style:
                df = df.groupby([self.time_index_key, self.scan_index_key, style]).sum()
            else:
                df = df.groupby([self.time_index_key, self.scan_index_key]).sum()
            df = df.reset_index()
        if x_name is None:
            x_name = self.scan_index_key
        sns.lineplot(x=x_name, y=y_name, data=df, hue=hue, ax=ax, style=style, legend=legend,
                     palette=palette, ci=ci, dashes=dashes)

        self.finalize_steady_state_plot(ax, y_name, ylim, ylog, xlog, x_name=x_name)

    def cell_steady_state_barplot(self, y_name, legend=False, hue=None, style=None, ylog=False, cummulative=False,
                                  ylim=None, ci="sd", y_ticks=True, **kwargs) -> None:
        ax, df, palette, hue = self.prepare_plot(self.cell_df, hue, **kwargs)
        df = df.loc[df[self.time_key] == self.t_max]

        sns.barplot(x=self.scan_index_key, y=y_name, data=df, ax=ax, ci=ci)
        ax.set_ylabel("")
        ax.set_title(self.get_label(y_name))
        ticks, labels = self.get_scan_ticks()
        ax.set_xticklabels(labels)
        ax.set_xlabel(self.get_label(self.scan_index_key))

        if not y_ticks:
            ax.set_yticklabels([])
        if ylim:
            ax.set_ylim(ylim)

    def single_cell_steady_state_plot(self, y_name, x_name=None, legend=False, hue=None, style=None, ylog=False,
                                      xlog=True,
                                      cummulative=False,
                                      ylim=None, linewidth=0.1, units="id", **kwargs):
        ax, df, palette, hue = self.prepare_plot(self.cell_df, hue, **kwargs)

        if cummulative:
            if hue and style:
                df = df.groupby([self.time_index_key, self.scan_index_key, hue, style]).sum()
            elif hue:
                df = df.groupby([self.time_index_key, self.scan_index_key, hue]).sum()
            elif style:
                df = df.groupby([self.time_index_key, self.scan_index_key, style]).sum()
            else:
                df = df.groupby([self.time_index_key, self.scan_index_key]).sum()
            df = df.reset_index()

        if x_name is None:
            x_name = self.scan_index_key

        sns.lineplot(x=x_name, y=y_name, data=df, hue=hue, ax=ax, style=style, legend=legend,
                     palette=palette, estimator=None, units=units, linewidth=linewidth)

        self.finalize_steady_state_plot(ax, y_name, ylim, ylog, xlog, x_name=x_name)

    def steady_state_count(self, legend=None, hue=None, style=None, relative=False, ylog=False, ylim=None, xlog=True,
                           ci="sd",
                           **kwargs):

        ax, df, palette, hue = self.prepare_plot(self.counts, hue, **kwargs)
        if relative:
            y = "n_rel"
            if ylim is None:
                ylim = [0, 1]
        else:
            y = "n"
            if ylim is None:
                ylim = False

        sns.lineplot(x=self.scan_index_key, y=y, hue=hue, data=df, ax=ax, style=style, ci=ci,
                     legend=legend,
                     palette=palette)

        self.make_legend_entry(ax)

        if ylog:
            ax.set_yscale("log")

        if ylim:
            ax.set_ylim(ylim)

        if relative:
            ax.set_ylabel("Fraction of cells")
        else:
            ax.set_ylabel("Number of cells")

        self.finalize_steady_state_plot(ax, y, ylim, ylog, xlog)

        # ax.set_ylabel(self.get_label(y))
        #
        # ticks, labels = self.get_scan_ticks()
        # ax.set_xticks(ticks)
        # ax.set_xticklabels(labels)

    def count_plot(self, legend=None, hue=None, style=None, relative=False, ci="sd", **kwargs) -> None:

        ax, df, palette, hue = self.prepare_plot(self.counts, hue, reduce=True, **kwargs)

        if relative:
            y = "n_rel"
            ylim = [0, 1]
        else:
            y = "n"
            ylim = False

        ax.set_ylabel("Number of cells")

        sns.lineplot(x=self.time_key, y=y, hue=hue, data=df, ax=ax, style=style, ci=ci,
                     legend=legend,
                     palette=palette)

        self.make_legend_entry(ax)
        if ylim:
            ax.set_ylim(ylim)

        if relative:
            ax.set_ylabel("Fraction of cells")
        else:
            ax.set_ylabel("Number of cells")

        ax.set_xlabel(self.get_label(self.time_key))
        ax.set_xlim([0, self.t_max])

    def cells_time_series_plot(self, y_name, legend=False, ylim=None, hue=None, style=None, ci="sd", **kwargs) -> None:

        ax, df, palette, hue = self.prepare_plot(self.cell_df, hue, **kwargs)

        sns.lineplot(x=self.time_key, y=y_name, hue=hue, data=df, ax=ax, style=style,
                     legend=legend, palette=palette, ci=ci)

        self.make_legend_entry(ax)

        ax.set_ylabel(self.get_label(y_name))
        ax.set_xlabel(self.get_label(self.time_key))
        ax.set_xlim([0, self.t_max])
        if ylim:
            ax.set_ylim(ylim)

    def cells_time_series_twinx(self, y_names, legend=False, ylim=None, colors_axis=True, style=None, ci="sd",
                                **kwargs) -> None:

        ax, df, palette = self.prepare_twinx_plot(self.cell_df, y_names, **kwargs)

        for i, y in enumerate(y_names):

            ax_l = ax[i]

            if colors_axis:
                ax_l.tick_params(axis='y', colors=palette[y])
                ax_l.yaxis.label.set_color(palette[y])

            sns.lineplot(x=self.time_key, y=y, data=df, ax=ax_l, style=style,
                         legend=legend, ci=ci, color=palette[y])

            self.make_legend_entry(ax_l)

            ax_l.set_ylabel(self.get_label(y))
            ax_l.set_xlabel(self.get_label(self.time_key))
            ax_l.set_xlim([0, self.t_max])
            if not ylim is None:
                ax_l.set_ylim(ylim[i])

    def cell_displacement_plot(self, legend=False, hue=None, n=False, palette_name="Dark2",
                               color=None) -> None:

        def displacement(df):

            def r0(df):
                d0 = df.loc[df["time_index"] == 0]

                x0 = pd.Series(list(d0["x"])[0] * np.ones(len(df)))
                x0.index = df.index
                y0 = pd.Series(list(d0["y"])[0] * np.ones(len(df)))
                y0.index = df.index
                z0 = pd.Series(list(d0["z"])[0] * np.ones(len(df)))
                z0.index = df.index

                df["x0"] = x0
                df["y0"] = y0
                df["z0"] = z0

                return df

            df = df.groupby("id").apply(r0)

            df["dr"] = np.sqrt(
                np.power(df["x"] - df["x0"], 2) + np.power(df["y"] - df["y0"], 2) + np.power(df["z"] - df["z0"], 2))
            return df

        cell_df = self.reduce_df(self.cell_df, self.scan_index_key)
        if n:
            ids = cell_df["id"].unique()
            cell_df = cell_df.loc[cell_df["id"].isin(ids[np.random.randint(0, len(ids), n)])]

        df = cell_df.groupby([self.scan_index_key, self.scan_name_key]).apply(displacement)
        ax, df, palette, hue = self.prepare_plot(df, hue)
        ax = self.get_subplot_axes()

        palette = self.get_palette(df, hue, palette_name=palette_name)

        if color:
            sns.lineplot(x="time", y="dr", data=df, units="id", ax=ax, estimator=None, color=color, legend=legend)
        else:
            sns.lineplot(x="time", y="dr", data=df, units="id", ax=ax, estimator=None, hue=hue, legend=legend,
                         palette=palette)

        self.make_legend_entry(ax)

        ax.set_ylabel(self.get_label("dr"))
        ax.set_xlabel(self.get_label("time"))

    def cell_plot(self, x_name, y_name, legend=False, hue=None, style=None, ci="sd", time=None, ylim=None, xlim=None,
                  condition=lambda df: df, count=False, **kwargs) -> None:

        if time:
            cell_df = self.cell_df.loc[self.cell_df["time"].isin(time)]
        else:
            cell_df = self.cell_df

        ax, df, palette, hue = self.prepare_plot(cell_df, hue, **kwargs)

        df = condition(df)
        if count:
            grp = [x_name]
            if hue is not None:
                grp.append(hue)
            if style is not None:
                grp.append(style)

            df = df.groupby(grp).count()
            df = df.reset_index()
            sns.lineplot(x=x_name, y="id_id", hue=hue, data=df, ax=ax, style=style, ci=ci,
                         legend=legend, palette=palette)
        else:
            sns.lineplot(x=x_name, y=y_name, hue=hue, data=df, ax=ax, style=style, ci=ci,
                         legend=legend, palette=palette)

        self.make_legend_entry(ax)

        ax.set_ylabel(self.get_label(y_name))
        ax.set_xlabel(self.get_label(x_name))
        if ylim:
            ax.set_ylim(ylim)
        if xlim:
            ax.set_xlim(xlim)

    def cell_plot_twinx(self, x_name, y_names, legend=False, ylim=None, style=None, color_axis=True, ci="sd",
                        time=None, marker=None, **kwargs) -> None:

        ax, df, palette = self.prepare_twinx_plot(self.cell_df, y_names, **kwargs)
        if time:
            df = df.loc[df["time"].isin(time)]

        for i, y in enumerate(y_names):

            ax_l = ax[i]
            if color_axis:
                ax_l.tick_params(axis='y', colors=palette[y])
                ax_l.yaxis.label.set_color(palette[y])
            sns.lineplot(x=x_name, y=y, data=df, ax=ax_l, style=style, ci=ci, legend=legend, color=palette[y],
                         marker=marker)

            self.make_legend_entry(ax_l)

            ax_l.set_ylabel(self.get_label(y))
            ax_l.set_xlabel(self.get_label(x_name))
            if not ylim is None:
                ax_l.set_ylim(ylim[i])

    def cell_slice_plot(self, y_name, axis_name="x", legend=False, hue=None, style=None, ci="sd", **kwargs) -> None:

        ax, df, palette, hue = self.prepare_plot(self.cell_df, hue, **kwargs)

        sns.lineplot(x=axis_name, y=y_name, hue=hue, data=df, ax=ax, style=style,
                     legend=legend, palette=palette, ci=ci)

        self.make_legend_entry(ax)
        ax.set_ylabel(self.get_label(y_name))
        ax.set_xlabel(self.get_label(axis_name))

    def cell_histogramm(self, x_name, t=None, hue=None, xlim=None, ylim=None, quantiles=None, distplot_kwargs=None,
                        xlog=True, ylog=True, **kwargs):

        if quantiles is None:
            quantiles = [0, 1]
        if distplot_kwargs is None:
            distplot_kwargs = {}

        distplot_kwargs.update({"kde": False})

        ax, df, palette, hue = self.prepare_plot(self.cell_df, hue, **kwargs)

        if t:
            df = df.loc[df[self.time_key].isin(t)]
        if hue:
            for i, h in enumerate(df[hue].unique()):
                x = df.loc[df[hue] == h]
                if "bins" in distplot_kwargs:
                    distplot_kwargs["bins"] = int(len(x) / 5) if distplot_kwargs["bins"] > len(x) / 5 else \
                        distplot_kwargs["bins"]

                sns.distplot(x[x_name], color=palette[h], **distplot_kwargs, ax=ax, norm_hist=False)
        else:
            if "bins" in distplot_kwargs:
                distplot_kwargs["bins"] = int(len(df) / 5) if distplot_kwargs["bins"] > len(df) / 5 else \
                    distplot_kwargs["bins"]
            print(df[x_name].mean())
            sns.distplot(df[x_name], **distplot_kwargs, ax=ax, norm_hist=False)

        self.make_legend_entry(ax)

        if ylog:
            ax.set_yscale('log')
        if xlog:
            ax.set_xscale('log')
        if ylim is not None:
            ax.set_ylim(ylim)

        if xlim is not None:
            ax.set_xlim(xlim)
        else:
            try:
                ax.set_xlim([df[x_name].quantile(quantiles[0]), df[x_name].quantile(quantiles[1])])
            except ValueError:
                pass

        ax.set_ylabel("absolute frequency")
        ax.set_xlabel(self.get_label(x_name))

    def cell_activation_histogramm(self, x_name, cummulative=False, relative=False, color="red", showmax=None, t=None,
                                   bins=100, hue=None, xlim=None, ylim=None, xlog=True, ylog=True, **kwargs):

        ax, df, palette, hue = self.prepare_plot(self.cell_df, hue, **kwargs)

        if t:
            df = df.loc[df[self.time_key].isin(t)]

        df["bins"], b = pd.cut(df[x_name], bins, retbins=True)
        bins = df.groupby(["bins"]).count().iloc[:, 0]
        if relative:
            bins = bins / np.sum(bins)

        act = df.groupby(["bins"]).mean()["activation"] * bins
        if cummulative:
            act = np.cumsum(act)

        act = np.array(act)
        nans = np.argwhere(np.isnan(act))

        b = b[1:]
        b = np.delete(b, nans)
        act = np.delete(act, nans)

        import matplotlib.lines as lines

        ax.plot(b, act, "-", color=color)

        if showmax is not None:
            q = np.quantile(np.array(act), [showmax])

            x = b[np.argmin(np.abs(act - q))]
            print(x)
            line = lines.Line2D([x, x], [0, 1])
            ax.add_artist(line)

        if xlim is not None:
            ax.set_xlim(xlim)
        if ylim is not None:
            ax.set_ylim(ylim)
        if ylog:
            ax.set_yscale('log')
        if xlog:
            ax.set_xscale('log')

        if relative:
            ax.set_ylabel("% activation")
        else:
            ax.set_ylabel("n activated")

        ax.set_xlabel(self.get_label(x_name))

    def _compute_radial_profile(self, df, y_names, center_func=lambda df: df.loc[df["type_name"] == "sec"], n_workers=8,
                                chunksize=32):

        if not isinstance(y_names, List):
            y_names = [y_names]

        def init(_y_names, _groups, _center_func):
            global y_names
            y_names = _y_names

            global groups
            groups = _groups

            global center_func
            center_func = _center_func

        from multiprocessing import Pool
        with Pool(n_workers, initializer=init, initargs=(y_names, self.groups, center_func)) as p:
            l = list(df.groupby(["raw_scan_index", self.time_index_key]))

            c = int(len(l) / n_workers)
            c = c if c > 0 else 1

            chunksize = c if c < chunksize else chunksize
            print("running radials for l: {l} and chunksize {cs} on {n} workers".format(l=len(l), cs=chunksize,
                                                                                        n=n_workers))
            full_result = pd.concat(p.starmap(run_single_step, l, chunksize=chunksize))

        return full_result

    def compute_radial_profiles(self, y_names, center_func=lambda df: df.loc[df["type_name"] == "sec"], n_workers=8,
                                chunksize=32):

        self.radials_df = self._compute_radial_profile(self.cell_df, y_names, center_func=center_func,
                                                       n_workers=n_workers, chunksize=chunksize)

    def compute_cell_distance_metric(self, source_type_name, metric_dict={"mean_distance": np.mean}):

        for k in metric_dict.keys():
            if k in self.cell_df.columns:
                self.cell_df = self.cell_df.drop(columns=[k])

        result = []
        for i, gf in self.cell_df.groupby(["raw_scan_index", self.time_index_key]):
            print(i)
            ids, r = self.compute_distance_metric(gf, source_type_name, metrics=metric_dict)
            d = {"id_id": ids, "raw_scan_index": gf["raw_scan_index"],
                 self.time_index_key: gf[self.time_index_key]}
            d.update(r)
            result.append(pd.DataFrame(d))

        result = pd.concat(result)

        self.cell_df = self.cell_df.merge(result, on=["id_id", "raw_scan_index", self.time_index_key])

    def compute_distance_metric(self, df, type_name, metrics):

        ids, r = get_distance_matrix(df)
        source_ids = np.array(df.loc[df.type_name == type_name]["id_id"])
        folded = r[np.argwhere(np.isin(ids, source_ids)).ravel()]

        d = {}
        for k, v in metrics.items():
            d[k] = np.apply_along_axis(v, 0, folded)
        return ids, d

    def cell_radial_niche_plot(self, y_name, center_type, hue=None, style=None, xlim=None, ylim=None, ylog=False,
                               ci="sd", cell_radius=None, legend=None, estimator=None, plot_filter=lambda df: df,
                               **kwargs):

        df = self.cell_df

        # cell = df.loc[df.type_name == "sec"].iloc[5]
        # df = df.loc[(df.type_name == "abs") | (df.id_id == cell.id_id)]

        if hasattr(self, "radials_df") and hue in self.radials_df.columns and y_name in self.radials_df.columns:
            ax, df, palette, hue = self.prepare_plot(self.radials_df, hue, **kwargs)

        else:
            print("recomputing radial profiles because hue or y_name wasn't found in cache")
            ax, df, palette, hue = self.prepare_plot(df, hue, **kwargs)
            df = self._compute_radial_profile(df, [y_name, "type_name"],
                                              center_func=lambda df: df.loc[df["type_name"] == center_type])
        df = plot_filter(df)

        df, ci = self.compute_ci(
            df,
            [self.scan_index_key, self.time_index_key, "type_name", "distance", hue, style],
            ci=ci, estimator=estimator, y_names=[y_name])

        if cell_radius:
            df["distance"] = df["distance"] / cell_radius

        sns.lineplot(x="distance", y=y_name, hue=hue, style=style, data=df, ax=ax, legend=legend, ci=ci,
                     palette=palette)

        if ylog:
            ax.set_yscale("log")

        if xlim:
            ax.set_xlim(np.array(xlim))
        if ylim:
            ax.set_ylim(ylim)

        if cell_radius:
            ax.set_xlabel("distance from secreting cell (cell radius)")
        else:
            ax.set_xlabel(r"r $\mu m $")

        ax.set_ylabel(self.get_label(y_name))

        self.make_legend_entry(ax)

    # def _cell_radial_niche_plot(self, y_name, center_type, hue = None, style = None, xlim = None, ylim = None, ylog = False, ci = "sd", cell_radius = None, legend = None, estimator = None, **kwargs):
    #
    #     ax, df, palette, hue = self.prepare_plot(self.cell_df, hue, **kwargs)
    #
    #     def get_distance_matrix(df):
    #
    #         df = df.groupby("id").first()
    #         df = df.reset_index()
    #         ids = np.array(df["id"], dtype=int)
    #
    #         XX = np.array([
    #             np.array(df["x"], dtype=float),
    #             np.array(df["y"], dtype=float),
    #             np.array(df["z"], dtype=float),
    #         ]).T
    #
    #         r = distance_matrix(XX, XX, p=2)
    #         return ids, r
    #
    #
    #     ids, r = get_distance_matrix(df)
    #
    #     if cell_radius:
    #         r = r / cell_radius
    #
    #     final_result = pd.DataFrame(columns=["distance", y_name])
    #
    #     groups = [self.time_index_key,self.scan_index_key]
    #
    #     if hue:
    #         if not hue== "type_name":
    #             groups += [hue]
    #     if style:
    #         if not style == "type_name":
    #             groups += [style]
    #
    #     for o, cells in df.groupby(groups):
    #
    #         secretors_step = cells.loc[df["type_name"] == center_type]
    #
    #         for i, sec in secretors_step.iterrows():
    #             i = np.where(ids == sec["id"])[0][0]
    #
    #             result = pd.DataFrame({
    #                 "id": ids,
    #                 "distance": r[i],
    #                 "scan_index": sec["scan_index"],
    #                 y_name: cells[y_name],
    #                 "activation": cells["activation"],
    #                 "type_name": cells["type_name"]
    #             })
    #
    #             for g in groups:
    #                 result[g] = sec[g]
    #
    #             final_result = final_result.append(result)
    #
    #
    #     final_result,ci = self.compute_ci(final_result,[self.scan_index_key,self.time_index_key,"type_name","distance",hue,style], ci = ci, estimator=estimator)
    #
    #     sns.lineplot(x="distance", y=y_name, data=final_result, ci=ci, ax=ax, style=style, hue=hue ,legend=legend)
    #
    #     if ylog:
    #         ax.set_yscale("log")
    #
    #     if xlim:
    #         ax.set_xlim(np.array(xlim))
    #     if ylim:
    #         ax.set_ylim(ylim)
    #     if cell_radius:
    #         ax.set_xlabel("distance from secreting cell (cell radius)")
    #     else:
    #         ax.set_xlabel(r"r $\mu m $")
    #
    #     ax.set_ylabel(self.get_label(y_name))
    #
    #     self.make_legend_entry(ax)

    def compute_ci(self, df, group_by_columns, ci="sd", estimator=None, y_names=None):

        if ci in ["sd", None] or isinstance(ci, float) or isinstance(ci, int):
            return df, ci
        elif ci == "sem":
            for g in group_by_columns:
                if g is None:
                    group_by_columns.remove(g)

            if y_names is not None:
                df = df[list(set(group_by_columns + y_names))]

            gb = df.groupby(list(set(group_by_columns)))

            if estimator is None:
                return gb.mean().reset_index(), "sd"
            else:
                return gb.agg(estimator).reset_index(), "sd"
        else:
            return df, ci

    def cell_heatmap(self, x_name, y_name, z_name, filter={}, cmap="viridis", v_range=None, c_lines=None, levels=100,
                     xlog=False, ylog=False, accumulator=lambda groupby: groupby.mean(), **kwargs):

        hue = None

        ax, df, palette, hue = self.prepare_plot(self.cell_df, hue, **kwargs)

        gb = [x_name, y_name] + list(filter.keys())

        df = accumulator(df.groupby(gb))

        df = df.reset_index()
        for key in filter.keys():
            df = df.loc[df[key].isin(filter[key])]

        piv = df.pivot(y_name, x_name, z_name)
        # piv = piv.reindex(index=piv.index[::-1])

        x = np.array(piv.columns)
        y = np.array(piv.index)
        z = np.array(piv)

        n_ticks = 5
        ext = [min(x), max(x), min(y), max(y)]

        # cs = ax.imshow(z, extent=ext, aspect = "auto",origin="lower", cmap = cmap)
        from matplotlib.colors import Normalize
        from matplotlib.cm import ScalarMappable

        if v_range is None:
            v_range = [np.nanmin(z), np.nanmax(z)]

        levels = np.linspace(v_range[0], v_range[1], 100)

        cs = ax.contourf(x, y, z, vmin=v_range[0], vmax=v_range[1], cmap=cmap, levels=levels)
        for c in cs.collections:
            c.set_edgecolor("face")

        if c_lines:
            if not isinstance(c_lines, List):
                c_lines = np.linspace(v_range[0], v_range[1], c_lines)
            c_lines_sc = ax.contour(x, y, z, vmin=v_range[0], vmax=v_range[1], colors="black", levels=c_lines)
            ax.clabel(c_lines_sc, fmt='%2.3f', colors='black', fontsize=4)

        if ylog:
            ax.semilogy()
        if xlog:
            ax.semilogx()

        x = np.linspace(ext[0], ext[1], n_ticks)
        y = np.linspace(ext[2], ext[3], n_ticks)

        ax.set_xticks(x)
        ax.set_yticks(y)

        norm = Normalize(vmin=v_range[0], vmax=v_range[1])
        mappable = ScalarMappable(norm, cmap=cmap)

        if self.external_legend:
            plt.colorbar(mappable, label=self.get_label(z_name), ax=ax, cax=self.legend_axes)
        else:
            plt.colorbar(mappable, ax=ax)

        ax.set_ylabel(self.get_label(y_name))
        ax.set_xlabel(self.get_label(x_name))

        # ax.xaxis.set_major_formatter(StrMethodFormatter("{x:.1g}"))
        # ax.yaxis.set_major_formatter(StrMethodFormatter("{x:.0e}"))

        ax.set_xlabel(self.get_label(x_name))
        ax.set_ylabel(self.get_label(y_name))

    def global_heatmap(self, x_name, y_name, z_name, acc_filter={}, cmap="viridis", v_range=None, c_lines=None,
                       levels=100, xlog=False, ylog=False, accumulator=lambda groupby: groupby.mean(), mask=None,
                       **kwargs):

        hue = None

        ax, df, palette, hue = self.prepare_plot(self.global_df, hue, **kwargs)

        gb = [x_name, y_name] + list(acc_filter.keys())

        df = accumulator(df.groupby(gb))

        df = df.reset_index()
        for key in acc_filter.keys():
            df = df.loc[df[key].isin(acc_filter[key])]

        if mask is not None:
            df[z_name] = df.loc[df["success"] == False][z_name].replace(mask, np.nan)

        piv = df.pivot(y_name, x_name, z_name)

        # piv = piv.reindex(index=piv.index[::-1])

        x = np.array(piv.columns)
        y = np.array(piv.index)
        z = np.array(piv)

        n_ticks = 5
        ext = [min(x), max(x), min(y), max(y)]

        # cs = ax.imshow(z, extent=ext, aspect = "auto",origin="lower", cmap = cmap)
        from matplotlib.colors import Normalize
        from matplotlib.cm import ScalarMappable

        if v_range is None:
            v_range = [np.min(z), np.max(z)]
        levels = np.linspace(v_range[0], v_range[1], 100)
        cs = ax.contourf(x, y, z, vmin=v_range[0], vmax=v_range[1], cmap=cmap, levels=levels)

        for c in cs.collections:
            c.set_edgecolor("face")

        if c_lines:
            ax.contour(x, y, z, vmin=v_range[0], vmax=v_range[1], colors="red", levels=c_lines)
        if ylog:
            ax.semilogy()
        if xlog:
            ax.semilogx()

        x = np.linspace(ext[0], ext[1], n_ticks)
        y = np.linspace(ext[2], ext[3], n_ticks)

        ax.set_xticks(x)
        ax.set_yticks(y)

        norm = Normalize(vmin=v_range[0], vmax=v_range[1])
        mappable = ScalarMappable(norm, cmap=cmap)

        if self.external_legend:
            plt.colorbar(mappable, label=self.get_label(z_name), ax=ax, cax=self.legend_axes)
        else:
            plt.colorbar(mappable, ax=ax)

        ax.set_ylabel(self.get_label(y_name))
        ax.set_xlabel(self.get_label(x_name))

        # ax.xaxis.set_major_formatter(StrMethodFormatter("{x:.1g}"))
        # ax.yaxis.set_major_formatter(StrMethodFormatter("{x:.0e}"))

        ax.set_xlabel(self.get_label(x_name))
        ax.set_ylabel(self.get_label(y_name))

    def _global_heatmap(self, x_name, y_name, z_name, filter={}, accumulator=lambda groupby: groupby.mean(), **kwargs):

        hue = None

        ax, df, palette, hue = self.prepare_plot(self.global_df, hue, **kwargs)

        gb = [x_name, y_name] + list(filter.keys())

        df = accumulator(df.groupby(gb))

        df = df.reset_index()
        for key in filter.keys():
            df = df.loc[df[key].isin(filter[key])]

        piv = df.pivot(y_name, x_name, z_name)
        sns.heatmap(piv, cbar_kws={"label": self.get_label(z_name)})

        ax.set_yticklabels([round(float(i._text), 1) for i in ax.get_yticklabels()])
        ax.set_xticklabels([round(float(i._text), 1) for i in ax.get_xticklabels()])

        ax.set_xlabel(self.get_label(x_name))
        ax.set_ylabel(self.get_label(y_name))

    def function_twinx_overlay(self, f, hue=None, ylabel=None, sub_divisions=10, plot_args=(), **kwargs):

        ax, df, palette, hue = self.prepare_plot(self.cell_df, hue, **kwargs)
        xticks = ax.get_xticks()

        x = np.linspace(np.min(xticks), np.max(xticks), len(xticks) * sub_divisions)
        ax = ax.twinx()

        if ylabel:
            ax.set_ylabel(ylabel)
        else:
            ax.set_ylabel("f")

        y = np.apply_along_axis(f, 0, x)

        ax.plot(x, y, *plot_args)

    def function_plot(self, f, hue=None, xlim=None, ylim=None, plot_kwargs={}, **kwargs):

        ax, df, palette, hue = self.prepare_plot(self.cell_df, hue, **kwargs)
        ticks = ax.get_xticks()
        x = np.linspace(min(ticks), max(ticks), 100)
        y = np.apply_along_axis(f, 0, x)

        self.make_legend_entry(ax)
        ax.plot(x, y, **plot_kwargs)
        if xlim:
            ax.set_xlim(xlim)
        if ylim:
            ax.set_ylim(ylim)

    def plot_twinx_overlay(self, y, legend_name=None, hue=None, plot_args=(), y_label=None, **kwargs):

        ax, df, palette, hue = self.prepare_plot(self.cell_df, hue, **kwargs)
        xticks = ax.get_xticks()
        old_ax = ax
        ax = ax.twinx()

        line = ax.plot(xticks, y, *plot_args)

        if y_label:
            ax.set_ylabel(y_label)
        if legend_name:
            handles, labels = old_ax.get_legend_handles_labels()
            handles += [line[0]]
            labels += ["Overlay: " + legend_name]
            old_ax.legend(handles, labels)

            self.make_legend_entry(old_ax)

    def cell_scatter_plot(self, names, t=None, hue=None, m=0.05, legend=None, marker="o", s=0.1, **kwargs):

        from matplotlib.lines import Line2D
        ax, df, palette, hue = self.prepare_plot(self.cell_df, hue, **kwargs)

        if t:
            df = df.loc[df[self.time_key].isin(t)]
        if hue:
            handels = []
            labels = []

            for i, h in enumerate(df[hue].unique()):
                x = df.loc[df[hue] == h]
                ax.scatter(x[names[0]], x[names[1]], marker=marker, color=palette[(hue, h)], s=s)
                labels.append(str(h))
                handels.append(
                    Line2D([0], [0], marker=marker, color=palette[h], markersize=1,
                           markerfacecolor=palette[(hue, h)], markeredgewidth=1, linewidth=1)
                )

            ax.legend(handels, labels)
        else:
            ax.scatter(df[names[0]], df[names[1]], marker=marker, s=s)

        self.make_legend_entry(ax)

        ax.set_ylim([(1 - m) * df[names[1]].min(), (1 + m) * df[names[1]].max()])
        ax.set_xlim([(1 - m) * df[names[0]].min(), (1 + m) * df[names[0]].max()])

        ax.set_ylabel(self.get_label(names[1]))
        ax.set_xlabel(self.get_label(names[0]))

    def cell_density_plot(self, x_name, t=None, hue=None, legend=False, **kwargs):
        ax, df, palette, hue = self.prepare_plot(self.cell_df, hue, **kwargs)

        if t:
            df = df.loc[df[self.time_key].isin(t)]

        if hue:
            for i, h in enumerate(df[hue].unique()):
                x = df.loc[df[hue] == h]
                sns.kdeplot(data=x[x_name], color=palette[h], cut=0, ax=ax, legend=legend, **kwargs)
        else:
            sns.kdeplot(data=df[x_name], cut=0, ax=ax, legend=legend, **kwargs)

        self.make_legend_entry(ax)

        ax.set_xlim([df[x_name].min(), 1 * df[x_name].max()])
        ax.set_ylabel("cell density")
        ax.set_xlabel(self.get_label(x_name))

    def counts_global_twinx(self, y_names, legend=False, colors_axis=False, ci="sd"):
        pass

    def timing_barplot(self, y_name, ci="sd", legend=False, **kwargs):

        hue = y_name
        ax, df, palette, hue = self.prepare_plot(self.timing_df, hue, **kwargs)
        sns.barplot(y=y_name, x="duration", data=df, ci=ci, hue=hue, orient="h", palette=palette, dodge=False, ax=ax)
        ax.set_xlabel("time(s)")
        ax.set_yticklabels(self.replace_labels([i.get_text() for i in ax.get_yticklabels()]))
        ax.set_xlim([0, df["duration"].quantile(0.99)])

        self.make_legend_entry(ax)
        if not legend and ax.get_legend():
            ax.get_legend().remove()

    def timing_lineplot(self, y_name, x_name=None, hue=None, style=None, ci="sd", ylim=None, legend=False, **kwargs):

        ax, df, palette, hue = self.prepare_plot(self.timing_df, hue, **kwargs)

        if x_name is None:
            x_name = self.scan_index_key
        if ylim is not None:
            ax.set_ylim(ylim)

        sns.lineplot(x=x_name, y=y_name, data=df, ci=ci, hue=hue, style=style, legend=legend, ax=ax, palette=palette)

        # ax.set_xlabel("time(s)")
        # ax.set_xlabel(self.get_label(x_name))
        # ax.set_ylabel(self.get_label(y_name))
        self.make_legend_entry(ax)

    # noinspection PyUnusedLocal
    def timing_timelineplot(self, **kwargs):
        hue = None

        ax, df, palette, hue = self.prepare_plot(self.timing_df, hue, **kwargs)
        df["level"] = df["task"].map(lambda x: x.count(":"))
        df["cycle"] = pd.Series(np.zeros(len(df)))
        cycle_start = "run:scan_sample:SimContainer:run:step"

        cycles = df.loc[df["task"] == cycle_start]

        for i, c in cycles.iterrows():
            df["cycle"][(df["start"] >= c["start"]) & (df["end"] < c["end"])] = int(i - 1)

        colors = {key: self.get_categorical_color(key) for key in df["task"].unique()}

        level_size = df.groupby(["level"], as_index=False).count()

        def get_y(level, i, h=0.01):

            n = level_size.iloc[row["level"]]["task"]

            return 0.5 + level + h * i

        y_max = 0
        level_counter = {i: 0 for i in range(len(level_size))}

        for o, c in df.groupby("cycle"):

            print(len(c))
            for i, row in c.iterrows():
                width = row["duration"]

                start = row["start"] - cycles.iloc[int(row["cycle"])]["start"]
                end = row["end"] - cycles.iloc[int(row["cycle"])]["end"]

                y = get_y(row["level"], level_counter[row["level"]], h=1 / level_size.iloc[row["level"]]["task"])
                level_counter[row["level"]] = level_counter[row["level"]] + 1

                y_max = y if y > y_max else y_max

                # bar = Rectangle((start, y), width, 1/level_size.iloc[row["level"]]["task"], color=colors[row["task"]])
                bar = Rectangle((start, y), width, 0.1, color=colors[row["task"]])
                ax.add_patch(bar)

        from matplotlib.lines import Line2D
        custom_lines = []
        for k, c in colors.items():
            custom_lines.append(Line2D([0], [0], color=c, lw=4))

        ax.legend(custom_lines, colors.keys())
        self.make_legend_entry(ax)

        ax.set_xlim(0, cycles["duration"].max() * 1.2)
        ax.set_ylim([0, y_max * 1.1])

    def ruse_plot(self, IMGPATH):

        df = self.ruse
        show = ["time_index", "scan_index", "ru_utime", "ru_stime", "ru_minflt", "ru_oublock", "ru_nvcsw", "ru_nivcsw",
                "ru_maxrss", "ru_inblock"]
        for c in df.columns:
            if c not in show:
                df.drop(c, inplace=True, axis=1)

        grad = ["ru_utime", "ru_stime", "ru_minflt", "ru_oublock", "ru_nvcsw", "ru_nivcsw"]

        for k in grad:
            df[k] = np.gradient(df[k])

        a = 8.3 * 0.5
        b = np.sqrt(2) * a * 0.7

        fig, ax = plt.subplots(4, 2, figsize=(a, b), sharex=True)
        ax = np.ravel(ax)

        for i, c in enumerate(df.columns):
            if c in ["time_index", "scan_index"]:
                continue
            axl = ax[i]
            if len(df["time_index"].unique()) > 1:
                sns.lineplot(x="time_index", y=c, data=df, ax=axl, hue="scan_index", legend=False)
            else:
                sns.lineplot(x="scan_index", y=c, data=df, ax=axl, legend=False)
            # axl.set_xlabel(self.get_label(axl.get_xlabel()))
            # axl.set_ylabel(self.get_label(axl.get_label()))

        plt.tight_layout()
        plt.savefig(os.path.join(IMGPATH, "ruse.pdf"))


def split_kwargs(kwargs, keys):
    result = {}
    for k in keys:
        if k in kwargs.keys():
            result[k] = kwargs[k]

    return result


def run_single_step(i, dfg):
    rsi, ti = i

    time_index_key = "time_index"
    scan_name_key = "scan_name_scan_name"
    scan_index_key = "scan_index"

    chunk_result = []
    ids, r = get_distance_matrix(dfg)

    for i, center in center_func(dfg).iterrows():
        single_center_result = pd.DataFrame()
        center_id = center["id"]
        i = np.where(ids == center_id)[0]
        distances = r[:, i][:, 0]
        sort_index = np.argsort(distances)

        distances = np.take_along_axis(distances, sort_index, axis=0)
        single_center_result["distance"] = distances
        properties = ["id_id"] + ["raw_scan_index", time_index_key, scan_name_key,
                                  scan_index_key] + y_names + groups
        for p in properties:
            single_center_result[p] = np.take_along_axis(np.array(dfg[p]), sort_index, axis=0)

        chunk_result.append(single_center_result)

    return pd.concat(chunk_result)


def get_distance_matrix(df):  # for single replicate

    ids = np.array(df["id"], dtype=int)
    XX = np.array([
        np.array(df["x"], dtype=float),
        np.array(df["y"], dtype=float),
        np.array(df["z"], dtype=float),
    ]).T

    r = distance_matrix(XX, XX, p=2)
    return ids, r
