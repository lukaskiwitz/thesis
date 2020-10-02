from typing import List, Dict

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import seaborn as sns


class Plotter:

    def __init__(self, path):

        self.time_key: str = "time"
        self.time_index_key: str = "time_index"

        self.t_max = 10
        self.scan_name_key: str = "scan_name_scan_name"
        self.scan_index_key: str = "scan_index"
        self.legend_axes = None
        self.external_legend = True
        self.legend_entries = {}
        self.max_scans = 5
        self.gridspec_index = 0

        self.scan_scale: List = []
        self.color_dict: Dict = {}
        self.style_dict = {

        }
        self.label_replacement = {
            "type_name": "Cell Type",
            "field_name": "Cytokine",
            self.scan_name_key: "Scan Name",
            "time": "Time",
            "time_index": "Time",
            "scan_index": "parameter fold-change"
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
                   'xtick.major.width': 0.5,
                   'ytick.major.width': 0.5,
                   'xtick.minor.width': 0.8,
                   'ytick.minor.width': 0.8,
                   'xtick.major.size': 1,
                   'ytick.major.size': 1,
                   'xtick.minor.size': 0.5,
                   'ytick.minor.size': 0.5
                   }

        sns.set_context("paper", rc=self.rc)
        self.load(path)

        self._prepare_color_dict()

        self.time_index_max, self.t_max  = self.get_max_time_index()

    def subplots(self, n, m, figsize=(10, 5), external_legend=True, gridspec_args={}):

        self.external_legend = external_legend

        self.n = n
        self.m = m
        self.gridspec_index = 0

        self.fig = plt.figure(figsize=figsize)
        if external_legend:
            self.gridspec = self.fig.add_gridspec(n, m + 1, **gridspec_args)
            self.legend_axes = self.fig.add_subplot(self.gridspec[:, -1])
        else:
            self.gridspec = self.fig.add_gridspec(n, m, **gridspec_args)

    def get_subplot_axes(self, overlay=False, gs_slice=None, **kwargs):

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

    def make_legend(self):

        if self.external_legend:
            handles = list(self.legend_entries.keys())
            labels = list(self.legend_entries.values())
            labels = self.replace_labels(labels)

            self.legend_axes.legend(handles, labels, loc="center left")
            self.legend_axes.axis('off')

    def _prepare_color_dict(self):
        fields = self.global_df["field_name"].unique()
        cell_types = self.cell_df["type_name"].unique()
        t = self.cell_df[self.time_key].unique()

        try:
            scan_names = self.global_df[self.scan_name_key].unique()
        except:
            self.color_dict = self.get_color_dict(np.concatenate([fields, t, cell_types]))
        else:
            self.color_dict = self.get_color_dict(np.concatenate([fields, t, cell_types, scan_names]))

    def get_color(self, key, palette_name="Dark2"):

        if key in self.color_dict.keys():
            return self.color_dict[key]
        else:
            keys = list(self.color_dict.keys())
            keys.append(key)
            self.color_dict = self.get_color_dict(keys, palette_name=palette_name)
            return self.color_dict[key]

    def get_color_dict(self, keys, palette_name="Dark2"):

        palette = sns.color_palette(palette_name, len(keys))
        color_dict = {}
        for i, k in enumerate(keys):
            color_dict[k] = palette[i]
        return color_dict

    def get_palette(self, df, key, palette_name="Dark2"):

        if key is None:
            return None

        keys = df[key].unique()

        p = {}
        for k in keys:
            p[k] = self.get_color(k, palette_name=palette_name)

        return p

    def load(self, path) -> None:

        global_df: pd.DataFrame = pd.read_hdf(path + "global_df.h5", mode="r")
        global_df = self.reset_scan_index(global_df)
        self.global_df = global_df

        cell_df: pd.DataFrame = pd.read_hdf(path + "cell_df.h5", mode="r")
        cell_df = self.reset_scan_index(cell_df)
        self.cell_df = cell_df

        groups = [

            "type_name",
            self.time_key,
            self.time_index_key,
            self.scan_index_key,
        ]
        if self.scan_name_key in list(cell_df.columns):
            groups.append(self.scan_name_key)

        grouped_cells = self.cell_df.groupby(groups, as_index=True)

        means = grouped_cells.mean()
        means.reset_index(inplace=True)
        self.means = means

        counts = grouped_cells.count()
        counts["n"] = counts["id"]
        counts = counts.drop(columns=counts.columns.drop(["n"]))

        counts = counts.reset_index()
        groups.remove("type_name")
        total = pd.Series(np.zeros(len(counts)))
        for i, g in counts.groupby(groups):
            n = g["n"].sum()
            for o in g.index:
                total.iloc[o] = n

        counts.reset_index(inplace=True)
        counts["n_rel"] = counts["n"] / total
        self.counts = counts

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

    def format_scan_index_ticklabels(self, scan_scale, distance, round_n) -> List:
        # return scan_scale

        my_scale = np.round(scan_scale, 2)

        return my_scale

    def replace_labels(self, labels) -> Dict:

        for i, l in enumerate(labels):
            if l in list(self.label_replacement.keys()):
                labels[i] = self.label_replacement[l]
        return labels

    def format_scan_index(self, index_column):
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

        if not self.scan_name_key in list(df.columns):
            return df
        else:
            offsets = {}
            for scan_name in df[self.scan_name_key].unique():
                offsets[scan_name] = df.loc[df[self.scan_name_key] == scan_name][self.scan_index_key].min()

            for i, o in offsets.items():
                mask = (df[self.scan_name_key] == i)
                df.loc[mask, self.scan_index_key] = df[self.scan_index_key] - o

            return df

    def make_legend_entry(self, ax):

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
            plt.legend(handles, labels)

    def savefig(self, path):

        import os
        os.makedirs(os.path.dirname(path), exist_ok=True)
        self.fig.tight_layout(pad=0)
        self.fig.savefig(path)

    def show(self):
        self.fig.tight_layout(pad=0)
        self.fig.show()

    def prepare_plot(self, df, hue, **kwargs):

        ax = self.get_subplot_axes(**kwargs)
        df = self.reduce_df(df, self.scan_index_key)
        palette = self.get_palette(df, hue)
        return ax, df, palette

    def prepare_twinx_plot(self, df, y_names, **kwargs):

        ax = self.get_subplot_axes(**kwargs)
        ax = [ax, ax.twinx()]
        df = self.reduce_df(df, self.scan_index_key)

        palette = {}
        for y in y_names:
            palette[y] = self.get_color(y)

        return ax, df, palette

    def global_time_series_plot(self, y_name, legend=False, hue=None, style=None, **kwargs) -> None:

        ax, df, palette = self.prepare_plot(self.global_df, hue, **kwargs)

        sns.lineplot(x=self.time_key, y=y_name, data=df, hue=hue, style=style,
                     ax=ax,
                     legend=legend,
                     palette=palette)

        self.make_legend_entry(ax)

        ax.set_xlabel(self.get_label(self.time_key))
        ax.set_xlim([0, self.t_max])
        ax.set_ylabel(self.get_label(y_name))

    def global_steady_state_plot(self, y_name, legend=False, hue=None, style=None, ylog = False, **kwargs) -> None:

        ax, df, palette = self.prepare_plot(self.global_df, hue, **kwargs)
        df = self.global_df.loc[self.global_df[self.time_key] == self.t_max]

        if style in df.columns:
            sns.lineplot(x=self.scan_index_key, y=y_name, data=df, hue=hue, ax=ax, legend=legend, palette=palette,
                         style=style)
        else:
            sns.lineplot(x=self.scan_index_key, y=y_name, data=df, hue=hue, ax=ax, legend=legend, palette=palette)

        self.make_legend_entry(ax)

        ax.set_xlabel(self.get_label(self.scan_index_key))
        ax.set_ylabel(self.get_label(y_name))
        if ylog:
            ax.set_yscale("log")

        ax.set_xticks(range(len(self.scan_scale)))
        ax.set_xticklabels(self.format_scan_index_ticklabels(self.scan_scale, 1, 1))

    def cell_steady_state_plot(self, y_name, legend=False, hue=None, style=None, ylog = False, **kwargs):
        ax, df, palette = self.prepare_plot(self.cell_df, hue, **kwargs)

        sns.lineplot(x=self.scan_index_key, y=y_name, data=df, hue=hue, ax=ax, legend=legend, palette=palette)

        self.make_legend_entry(ax)

        ax.set_xlabel(self.get_label(self.scan_index_key))
        ax.set_ylabel(self.get_label(y_name))
        if ylog:
            ax.set_yscale("log")

        ax.set_xticks(range(len(self.scan_scale)))
        ax.set_xticklabels(self.format_scan_index_ticklabels(self.scan_scale, 1, 1))

    def steady_state_count(self, legend=None, hue=None, style=None, relative=False, ylog = False, ci="sd", **kwargs):

        ax, df, palette = self.prepare_plot(self.counts, hue, **kwargs)
        if relative:
            y = "n_rel"
            ylim = [0, 1]
        else:
            y = "n"
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

        ax.set_ylabel(self.get_label(y))

        ax.set_xticks(range(len(self.scan_scale)))
        ax.set_xticklabels(self.format_scan_index_ticklabels(self.scan_scale, 0, 1))

    def count_plot(self, legend=None, hue=None, style=None, relative=False, ci="sd", **kwargs) -> None:

        ax, df, palette = self.prepare_plot(self.counts, hue, **kwargs)

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

        ax.set_xlim([0, self.t_max])
        if ylim:
            ax.set_ylim(ylim)

        if relative:
            ax.set_ylabel("Fraction of cells")
        else:
            ax.set_ylabel("Number of cells")

        ax.set_xlabel(self.get_label(self.time_key))
        ax.set_xlim([0, self.t_max])

    def cells_time_series_plot(self, y_name, legend=False, ylim=False, hue=None, style=None, ci="sd", **kwargs) -> None:

        ax, df, palette = self.prepare_plot(self.cell_df, hue, **kwargs)

        sns.lineplot(x=self.time_key, y=y_name, hue=hue, data=df, ax=ax, style=style,
                     legend=legend, palette=palette, ci=ci)

        self.make_legend_entry(ax)

        ax.set_ylabel(self.get_label(y_name))
        ax.set_xlabel(self.get_label(self.time_key))
        ax.set_xlim([0, self.t_max])
        if ylim:
            ax.set_ylim(ylim)

    def cells_time_series_twinx(self, y_names, legend=False, ylim=False, colors_axis=True, style=None, ci="sd",
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
            if ylim:
                ax_l.set_ylim(ylim[i])

    def cell_displacement_plot(self, legend=False, ylim=False, hue=None, style=None, n=False, palette_name="Dark2",
                               color=None, **kwargs) -> None:

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
        ax, df, palette = self.prepare_plot(df, hue, **kwargs)
        ax = self.get_subplot_axes(**kwargs)

        palette = self.get_palette(df, hue, palette_name=palette_name)

        if color:
            sns.lineplot(x="time", y="dr", data=df, units="id", ax=ax, estimator=None, color=color, legend=legend)
        else:
            sns.lineplot(x="time", y="dr", data=df, units="id", ax=ax, estimator=None, hue=hue, legend=legend,
                         palette=palette)

        self.make_legend_entry(ax)

        ax.set_ylabel(self.get_label("dr"))
        ax.set_xlabel(self.get_label("time"))

    def cell_plot(self, x_name, y_name, legend=False, ylim=False, hue=None, style=None, ci="sd", palette_name="Dark2",
                  time=None, **kwargs) -> None:

        if time:
            cell_df = self.cell_df.loc[self.cell_df["time"].isin(time)]
        ax, df, palette = self.prepare_plot(cell_df, hue, **kwargs)

        sns.lineplot(x=x_name, y=y_name, hue=hue, data=df, ax=ax, style=style, ci=ci,
                     legend=legend, palette=palette)

        self.make_legend_entry(ax)

        ax.set_ylabel(self.get_label(y_name))
        ax.set_xlabel(self.get_label(x_name))
        # if ylim:
        #     ax.set_ylim(ylim)

    def cell_plot_twinx(self, x_name, y_names, legend=False, ylim=False, style=None, color_axis=True, ci="sd",
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
            if ylim:
                ax_l.set_ylim(ylim[i])

    def cell_slice_plot(self, y_name, axis_name="x", legend=False, hue=None, style=None, ci="sd", palette_name="Dark2",
                        **kwargs) -> None:

        ax, df, palette = self.prepare_plot(self.cell_df, hue, **kwargs)

        sns.lineplot(x=axis_name, y=y_name, hue=hue, data=df, ax=ax, style=style,
                     legend=legend, palette=palette, ci=ci)

        self.make_legend_entry(ax)
        ax.set_ylabel(self.get_label(y_name))
        ax.set_xlabel(self.get_label(axis_name))

    def cell_histogramm(self, x_name, t=None, hue=None, palette_name="Dark2", **kwargs):

        ax, df, palette = self.prepare_plot(self.cell_df, hue, **kwargs)

        if t:
            df = df.loc[df[self.time_key].isin(t)]
        if hue:
            for i, h in enumerate(df[hue].unique()):
                x = df.loc[df[hue] == h]
                sns.distplot(x[x_name], color=palette[h], **kwargs)
        else:
            sns.distplot(df[x_name], **kwargs)

        self.make_legend_entry(ax)
        ax.set_ylabel("H")
        ax.set_xlabel(self.get_label(x_name))

    def cell_scatter_plot(self, names, t=None, hue=None, palette_name="Dark2", m=0.05, **kwargs):

        ax, df, palette = self.prepare_plot(self.cell_df, hue, **kwargs)

        if t:
            df = df.loc[df[self.time_key].isin(t)]
        if hue:
            for i, h in enumerate(df[hue].unique()):
                x = df.loc[df[hue] == h]
                sns.scatterplot(x[names[0]], x[names[1]], **kwargs, color=palette[h], ax=ax)
        else:
            sns.scatterplot(df[names[0]], df[names[1]], **kwargs, ax=ax)

        self.make_legend_entry(ax)

        ax.set_ylim([(1 - m) * df[names[1]].min(), (1 + m) * df[names[1]].max()])
        ax.set_xlim([(1 - m) * df[names[0]].min(), (1 + m) * df[names[0]].max()])

        ax.set_ylabel(self.get_label(names[1]))
        ax.set_xlabel(self.get_label(names[0]))

    def cell_density_plot(self, x_name, t=None, hue=None, palette_name="Dark2", legend=False, **kwargs):
        ax, df, palette = self.prepare_plot(self.cell_df, hue, **kwargs)

        if t:
            df = df.loc[df[self.time_key].isin(t)]

        if hue:
            for i, h in enumerate(df[hue].unique()):
                x = df.loc[df[hue] == h]
                sns.kdeplot(data=x[x_name], color=palette[h], cut=0, ax=ax, legend=legend, **kwargs)
        else:
            sns.kdeplot(data=df[x_name], cut=0, ax=ax, legend=legend, **kwargs)

        self.make_legend_entry(ax)

        ax.set_xlim([df[x_name].min(), (1) * df[x_name].max()])
        ax.set_ylabel("cell density")
        ax.set_xlabel(self.get_label(x_name))

    def counts_global_twinx(self, y_names, legend = False, colors_axis = False, ci = "sd"):
        pass