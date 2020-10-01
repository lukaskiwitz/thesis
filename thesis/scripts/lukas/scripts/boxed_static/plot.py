import os
import sys


import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import seaborn as sns
from parameters import path
import matplotlib.image as mpimg

def global_plot(fig, ax, global_df, y_name, y_label, legend="brief"):
    global_df = reduce_df(global_df, "scan_index")
    # global_df["scan_index"] = format_scan_index(global_df["scan_index"])

    sns.lineplot(x="time_index", y=y_name, data=global_df, hue="field_name", style="scan_index", ax=ax, legend=legend,
                 palette=color_dict)
    handles, labels = ax.get_legend_handles_labels()
    if legend:
        ax.legend(handles, labels_replace(labels), loc="upper right")
    ax.set_xlabel(time_label)
    ax.set_xlim([0, t_max])
    ax.set_ylabel(y_label)


def scan_score_plot(fig, ax, cell_df, score_name, y_label, legend="brief", ylim=False):
    cell_df = reduce_df(cell_df, "scan_index")

    # cell_df["scan_index"] = format_scan_index(cell_df["scan_index"])
    sns.lineplot(x="time", y=score_name, hue="type_name", data=cell_df, ax=ax, style="scan_index", ci=None,
                 legend=legend, palette=color_dict)
    handles, labels = ax.get_legend_handles_labels()
    if legend:
        ax.legend(handles, labels_replace(labels), loc="upper right")

    ax.set_xlim([0, t_max])
    if ylim:
        ax.set_ylim(ylim)
    ax.set_ylabel(y_label)
    ax.set_xlabel(time_label)


def count_plot(fig, ax, counts, y_name, y_label, legend="brief", ylim=False):
    counts = reduce_df(counts, "scan_index")
    sns.lineplot(x="time", y=y_name, hue="type_name", data=counts, ax=ax, style="scan_index", ci=None, legend=legend,
                 palette=color_dict)
    handles, labels = ax.get_legend_handles_labels()
    if legend:
        ax.legend(handles, labels_replace(labels), loc="upper right")
    ax.set_xlim([0, t_max])
    if ylim:
        ax.set_ylim(ylim)
    ax.set_ylabel(y_label)

    ax.set_xlabel(time_label)
    ax.set_xlim([0, t_max])
    ax.set_xlabel(time_label)


def steady_state_plot(fig, ax, df, time_name, x_name, y_name, x_label, y_label, legend="brief", hue="type_name",
                      leg_loc="upper right"):

    df = df.loc[df[time_name] == t_max]
    sns.lineplot(x=x_name, y=y_name, data=df, hue=hue, ax=ax, legend=legend, palette=color_dict)
    handles, labels = ax.get_legend_handles_labels()
    if legend:
        ax.legend(handles, labels_replace(labels), loc=leg_loc)

    ax.set_xlabel(x_label)
    ax.set_ylabel(y_label)
    ax.set_xticks(scan_scale)
    ax.set_yscale("log")
    ax.set_xticklabels(format_x_ticklabels(scan_scale, 1, 1))


def format_x_ticklabels(scan_scale, distance, round_n):
    return scan_scale
    my_scale = [str(np.round(scan_scale[0], round_n))]

    for i in np.arange(1, len(scan_scale)):
        i = int(i)
        if np.abs(scan_scale[i - 1] - scan_scale[i]) > distance:
            my_scale.append(str(np.round(scan_scale[i], round_n)))
        else:
            my_scale.append("")
    return my_scale


def labels_replace(labels,rep = None):
    if rep is None:
            rep = {
            "type_name": "Cell Type",
            "field_name": "Cytokine",
            "IL-2": "IL-2",
            "time": time_label,
            "time_index": time_label
            }

    for i, l in enumerate(labels):
        if l in list(rep.keys()):
            labels[i] = rep[l]
    return labels


def format_scan_index(index_column):
    return 1 * np.round(scan_scale[index_column.apply(lambda x: int(x))], 1)


def reduce_df(df, index_name):
    indices = df[index_name].unique()
    if len(indices) <= 2:
        return df
    indices = indices[0::int(len(indices) / 3)]
    result = pd.DataFrame()
    for i in indices:
        result = result.append(df.loc[(df[index_name] == i)])

    return result


def reset_scan_index(df):
    offsets = {}
    for scan_name in df["scan_name_scan_name"].unique():
        offsets[scan_name] = df.loc[df["scan_name_scan_name"] == scan_name]["scan_index"].min()

    for i,o in offsets.items():
        mask = (df["scan_name_scan_name"] == i)
        df.loc[mask,"scan_index"] = df["scan_index"] - o

    return df

IMGPATH = path
os.makedirs(IMGPATH, exist_ok=True)
# scan_index = 0
sns.set_context("paper", font_scale=1, rc={
    "lines.markersize": 4,
    "lines.linewidth": 2
}
                )

global_df: pd.DataFrame = pd.read_hdf(path + "global_df.h5", mode="r")
global_df = reset_scan_index(global_df)
cell_df: pd.DataFrame = pd.read_hdf(path + "cell_df.h5", mode="r")

grouped_cells = cell_df.groupby(["type_name", "time_index", "scan_index", "time"], as_index=True)

means = grouped_cells.mean()

means.reset_index(inplace=True)

std = means  # grouped_cells.std()

std.reset_index(inplace=True)

counts = grouped_cells.count()
counts["n"] = counts["id"]
counts = counts.drop(columns=counts.columns.drop(["n"]))
counts.reset_index(inplace=True)

t_max = global_df["time_index"].max()

time_label = "time (a.u.)"
scan_label = "scan index"
scan_scale = range(len(global_df["scan_index"].unique()))
c_unit = "nM"
cell_max = 1.2 * (1 + counts.max()["n"])
color_dict = {
    "naive": "blue",
    "sec": "red",
    "abs": "green",
    "IL-2": "yellow",
    "IL-6": "brown",
    "IFNg": "black",
    "D":"green",
    "kd":"red",
    "naive_amax":"green",
    "sec_q":"green",
    "sec_amax":"green",
    0:"red",
    1:"blue"
}


for i, k in enumerate(color_dict.keys()):
    color_dict[k] = sns.color_palette("muted", len(color_dict))[i]

fig, ax = plt.subplots(2, 2, sharex=False, sharey=False)
# for i,df in global_df.groupby("scan_name_scan_name"):
#     """scan plot """

steady_state_plot(fig, ax[0][0], global_df, "time_index", "scan_index", "Concentration", scan_label,
                  "t={t_max}. Cytokine ({c})".format(c=c_unit, t_max=t_max), hue="scan_name_scan_name", legend=False)
steady_state_plot(fig, ax[1][0], global_df, "time_index", "scan_index", "Gradient", scan_label,
                  "t={t_max} Gradient ({c}/$\mu m$)".format(c=c_unit, t_max=t_max), hue="scan_name_scan_name", legend=False)
steady_state_plot(fig, ax[1][1], global_df, "time_index", "scan_index", "SD", scan_label,
                  "t={t_max} SD of Nodal Values".format(c=c_unit, t_max=t_max), hue="scan_name_scan_name", legend="brief")
# steady_state_plot(fig, ax[0][1], counts, "time", "scan_index", "n", scan_label,
#                   "t={t_max} Number of cells".format(t_max=t_max), leg_loc="upper left")
#
plt.tight_layout()
plt.savefig(IMGPATH + "maximum_plot.pdf")
plt.show()

