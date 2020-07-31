from typing import List, Dict

import numpy as np
import pandas as pd
import seaborn as sns


class Plotter:

    def __init__(self,path):

        self.time_key: str = "time"
        self.time_index_key: str = "time_index"
        self.t_max = 10
        self.scan_name_key: str = "scan_name_scan_name"
        self.scan_index_key: str = "scan_index"

        self.max_scans = 5

        self.scan_scale: List = []
        self.color_dict: Dict = {
        }
        self.style_dict = {

        }
        self.label_replacement = {
            "type_name": "Cell Type",
            "field_name": "Cytokine",
            "IL-2": "IL-2",
            "time": "Time",
            "time_index": "Time"
        }


        self.load(path)

        self.t_max = self.get_max_time_index()


    def get_color(self,key):
        if key in self.color_dict.keys():
            return self.color_dict[key]
        else:
            keys = list(self.color_dict.keys())
            keys.append(key)
            self.color_dict = self.get_color_dict(keys)
            return self.color_dict[key]


    def get_color_dict(self,keys):

        palette = sns.color_palette("Dark2", len(keys))
        color_dict = {}
        for i, k in enumerate(keys):
            color_dict[k] = palette[i]
        return color_dict

    def get_palette(self,keys):

        p = {}
        for k in keys:
            p[k] = self.get_color(k)

        return p

    def load(self,path) -> None:

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

        grouped_cells = self.cell_df.groupby(groups,as_index=True)

        means = grouped_cells.mean()
        means.reset_index(inplace=True)
        self.means = means

        counts = grouped_cells.count()
        counts["n"] = counts["id"]
        counts = counts.drop(columns=counts.columns.drop(["n"]))

        counts = counts.reset_index()
        groups.remove("type_name")
        total = pd.Series(np.zeros(len(counts)))
        for i,g in counts.groupby(groups):
            n = g["n"].sum()
            for o in g.index:
                total.iloc[o] = n

        counts.reset_index(inplace=True)
        counts["n_rel"] = counts["n"]/total
        self.counts = counts

    def get_max_time_index(self) -> float:

        g = self.global_df["time_index"].max()
        c = self.cell_df["time_index"].max()

        return np.max([g,c])


    def get_label(self, key) -> str:

        if key in self.label_replacement.keys():
            return self.label_replacement[key]
        else:
            return None

    def format_x_ticklabels(self, scan_scale, distance, round_n) -> List:
        return scan_scale
        my_scale = [str(np.round(scan_scale[0], round_n))]

        for i in np.arange(1, len(scan_scale)):
            i = int(i)
            if np.abs(scan_scale[i - 1] - scan_scale[i]) > distance:
                my_scale.append(str(np.round(scan_scale[i], round_n)))
            else:
                my_scale.append("")
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
        indices = indices[0::int(len(indices) / self.max_scans)]
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


    def global_time_series_plot(self, fig, ax, y_name, y_label, legend="brief", hue_key ="field_name") -> None:
        global_df = self.reduce_df(self.global_df, self.scan_index_key)
        # global_df["scan_index"] = format_scan_index(global_df["scan_index"])
        palette = self.get_palette(global_df[hue_key].unique())

        sns.lineplot(x=self.time_index_key, y=y_name, data=global_df, hue=hue_key, style=self.scan_index_key,
                     ax=ax,
                     legend=legend,
                     palette=palette)
        handles, labels = ax.get_legend_handles_labels()
        if legend:
            ax.legend(handles, self.replace_labels(labels), loc="upper right")
        ax.set_xlabel(self.get_label(self.time_key))
        ax.set_xlim([0, self.t_max])
        ax.set_ylabel(y_label)

    def global_steady_state_plot(self, fig, ax, x_name, y_name, legend="brief", hue="type_name",
                                 leg_loc="upper right", style="scan_name_scan_name") -> None:

        df = self.global_df.loc[self.global_df[self.time_key] == self.t_max]
        if style in df.columns:
            sns.lineplot(x=x_name, y=y_name, data=df, hue=hue, ax=ax, legend=legend, palette=self.color_dict,
                         style=style)
        else:
            sns.lineplot(x=x_name, y=y_name, data=df, hue=hue, ax=ax, legend=legend, palette=self.color_dict)
        handles, labels = ax.get_legend_handles_labels()
        if legend:
            ax.legend(handles, self.replace_labels(labels), loc=leg_loc)

        # ax.set_xlabel(x_label)
        # ax.set_ylabel(y_label)
        # ax.set_yscale("log")
        # ax.set_xticks(self.scan_scale)
        # ax.set_xticklabels(self.format_x_ticklabels(self.scan_scale, 1, 1))

    def count_plot(self,fig, ax, legend="brief", ylim=False, hue_key="type_name", relative = False) -> None:


        counts = self.reduce_df(self.counts, self.scan_index_key)
        palette = self.get_palette(counts[hue_key].unique())

        y = "n_rel" if relative else "n"
        ax.set_ylabel("Number of cells")

        sns.lineplot(x=self.time_key, y=y, hue=hue_key, data=counts, ax=ax, style="scan_index", ci=None,
                     legend=legend,
                     palette=palette)

        handles, labels = ax.get_legend_handles_labels()

        if legend:
            ax.legend(handles, self.replace_labels(labels), loc="upper right")
        ax.set_xlim([0, self.t_max])
        if ylim:
            ax.set_ylim(ylim)

        if relative:
            ax.set_ylabel("Fraction of cells")
        else:
            ax.set_ylabel("Number of cells")

        ax.set_xlabel(self.get_label(self.time_index_key))
        ax.set_xlim([0, self.t_max])

    def cells_time_series_plot(self, fig, ax, score_name, y_label, legend="brief", ylim=False, hue_key ="type_name") -> None:
        cell_df = self.reduce_df(self.cell_df, "scan_index")
        palette = self.get_palette(cell_df[hue_key].unique())

        sns.lineplot(x="time", y=score_name, hue=hue_key, data=cell_df, ax=ax, style="scan_index", ci=None,
                     legend=legend, palette=palette)
        handles, labels = ax.get_legend_handles_labels()

        ax.set_ylabel(y_label)
        if legend:
            ax.legend(handles, self.replace_labels(labels), loc="upper right")

        ax.set_xlim([0, self.t_max])
        if ylim:
            ax.set_ylim(ylim)

