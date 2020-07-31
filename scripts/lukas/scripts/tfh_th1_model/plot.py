import os
import sys

import matplotlib
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import seaborn as sns
from parameters import path, IMGPATH

if len(sys.argv) > 1:
    if sys.argv[1] == "hide":
        matplotlib.use('Agg')

def plot_7(fig,ax,cell_df, t=None,legend="full"):


    m = 0.99
    cell_df = reduce_df(cell_df,"scan_index")
    if t:
        cell_df = cell_df.loc[(cell_df["time"] == t)]
    cell_df = cell_df.loc[(cell_df["type_name"] == "Tn")]
    cell_df["scan_index"] = format_scan_index(cell_df["scan_index"])
    sns.scatterplot(x="IL-6_surf_c", y="IL-2_surf_c", data=cell_df,hue="scan_index",ax=ax,palette="Paired",legend=legend)

    ax.add_artist(plt.Line2D([il6_threshold,il6_threshold],[0,1],linestyle='--'))
    ax.add_artist(plt.Line2D([0,1],[il2_threshold,il2_threshold],linestyle='--'))
    # ax  = g.ax_joint
    # ax.add_artist(plt.Rectangle(
    #     [
    #     il6_threshold,
    #     il2_threshold
    #      ],
    #     (
    #         (ax.get_xlim()[1])-il6_threshold)*m,
    #     (-ax.get_ylim()[1]+il2_threshold)*m,
    #     fill=False,edgecolor="red",linewidth=5
    #     )
    # )

    handles, labels = ax.get_legend_handles_labels()
    if legend:
        ax.legend(handles, labels_replace(labels), loc="lower left")

    ax.set_xlabel("Surface IL-6 (nM)")
    ax.set_title("t = {t}".format(t=t))
    ax.set_ylabel("Surface IL-2 (nM)")

def plot_8(fig,ax):
    fig, ax = plt.subplots(2, 2, sharex=False, sharey=False)
    sns.distplot(cell_df.loc[(cell_df["type_name"] == "Tn") & (cell_df["time"] == 0)]["IL-2_surf_c"], ax=ax[0][0])
    sns.distplot(cell_df.loc[(cell_df["type_name"] == "Tn") & (cell_df["time"] == 0)]["IL-6_surf_c"], ax=ax[0][1])
    sns.distplot(cell_df.loc[(cell_df["type_name"] == "Tn") & (cell_df["time"] == 0)]["IFNg_surf_c"], ax=ax[1][0])

    # ax[0][0].add_artist(plt.Line2D([il2_threshold, il2_threshold], [0, 1], linestyle='--'))
    # ax[0][1].add_artist(plt.Line2D([il6_threshold, il6_threshold], [0, 1], linestyle='--'))
    # ax[1][0].add_artist(plt.Line2D([ifng_threshold, ifng_threshold], [0, 1], linestyle='--'))

    plt.show()

    # tfh = cell_df.loc[
    #     (cell_df["IL-6_surf_c"] > il6_threshold)&
    #     (cell_df["IL-2_surf_c"] < il2_threshold)&
    #     (cell_df["type_name"] == "Tn")
    #             ].groupby(["time","scan_index","type_name"]).count()["id"].reset_index()
    # th1 = cell_df.loc[
    #     (cell_df["IFNg_surf_c"] > ifng_threshold) &
    #     (cell_df["type_name"] == "Tn")
    #     ].groupby(["time","scan_index","type_name"]).count()["id"].reset_index()
    #
    # fig,ax = plt.subplots(2,2,sharex=True,sharey=True)
    # sns.lineplot(x="time",y="id",data=tfh,ax=ax[0][0])
    # sns.lineplot(x="time",y="id",data=th1,ax=ax[0][1])
    # plt.show()


def global_plot(fig, ax, global_df, y_name, y_label, legend="brief"):
    global_df = reduce_df(global_df, "scan_index")
    global_df["scan_index"] = format_scan_index(global_df["scan_index"])

    sns.lineplot(x="time_index", y=y_name, data=global_df, hue="field_name", style="scan_index", ax=ax, legend=legend,
                 palette=color_dict)
    handles, labels = ax.get_legend_handles_labels()
    if legend:
        ax.legend(handles, labels_replace(labels), loc="upper right")
    ax.set_xlabel(time_label)
    ax.set_xlim([0, t_max])
    ax.set_ylabel(y_label)

def scan_score_plot(fig, ax, cell_df, score_name, y_label, legend="brief",ylim=False):

    cell_df = reduce_df(cell_df,"scan_index")

    cell_df["scan_index"] = format_scan_index(cell_df["scan_index"])
    sns.lineplot(x="time", y=score_name, hue="type_name", data=cell_df,ax=ax,style ="scan_index",ci=None,legend=legend,palette=color_dict)
    handles, labels = ax.get_legend_handles_labels()
    if legend:
        ax.legend(handles, labels_replace(labels), loc="upper right")

    ax.set_xlim([0, t_max])
    if ylim:
        ax.set_ylim(ylim)
    ax.set_ylabel(y_label)
    ax.set_xlabel(time_label)

def count_plot(fig, ax, counts, y_name, y_label, legend="brief",ylim=False):

    counts = reduce_df(counts,"scan_index")
    sns.lineplot(x="time", y=y_name, hue="type_name", data=counts, ax=ax, style="scan_index", ci=None, legend=legend,palette=color_dict)
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

def maximum_global_plot(fig,ax, global_df,y_name ,x_label, y_label, legend="brief"):

    global_max = global_df.groupby(["scan_index","field_name"],as_index=False).max()
    sns.lineplot(x="scan_index",y=y_name,data=global_max, hue="field_name",palette=color_dict,ax=ax,legend=legend)
    handles, labels = ax.get_legend_handles_labels()
    if legend:
        ax.legend(handles, labels_replace(labels), loc="upper right")

    ax.set_xlabel(x_label)
    ax.set_ylabel(y_label)
    ax.set_xticks(scan_scale)
    ax.set_xticklabels(format_x_ticklabels(scan_scale, 1, 1))

def steady_state_plot(fig,ax,df, time_name,x_name, y_name, x_label, y_label, legend="brief",hue="type_name",leg_loc="upper right"):

    df = df.loc[df[time_name] == t_max]
    sns.lineplot(x=x_name, y=y_name, data=df, hue=hue, ax=ax,legend=legend,palette=color_dict)
    handles, labels = ax.get_legend_handles_labels()
    if legend:
        ax.legend(handles, labels_replace(labels), loc=leg_loc)

    ax.set_xlabel(x_label)
    ax.set_ylabel(y_label)
    ax.set_xticks(scan_scale)
    ax.set_xticklabels(format_x_ticklabels(scan_scale,1,1))

def format_x_ticklabels(scan_scale,distance, round_n):

    my_scale = [str(np.round(scan_scale[0],round_n))]

    for i in np.arange(1, len(scan_scale)):
        i = int(i)
        if np.abs(scan_scale[i - 1] - scan_scale[i]) > distance:
            my_scale.append(str(np.round(scan_scale[i], round_n)))
        else:
            my_scale.append("")
    return my_scale

def labels_replace(labels):
    rep = {
        "scan_index": r"$q_s$ (molecules/s)",
        "scan_index": r"$q_s$ (molecules/s)",
        "type_name": "Cell Type",
        "field_name": "Cytokine",
        "IL-2": "IL-2",
        "IL-6":"IL-6",
        "IFNg":r"IFN$_\gamma$",
        "Tn":r"T$_n$",
        "Th1": r"Th$_1$",
        "Tfh": r"Tf$_h$",
        "time":time_label,
        "time_index":time_label
    }
    for i,l in enumerate(labels):
        if l in list(rep.keys()):
            labels[i] = rep[l]
    return labels

def format_scan_index(index_column):

    return 1*np.round(scan_scale[index_column.apply(lambda x: int(x))],1)

def reduce_df(df,index_name):
    return df
    if "scan_index" in globals() and globals()["scan_index"]:
        df = df.loc[
            (df[index_name] == globals()["scan_index"])
        ]
    else:
        df = df.loc[(
                      (df[index_name] == 0) |
                      (df[index_name] == 4) |
                      (df[index_name] == 9)
                  )]

    return df


os.makedirs(IMGPATH, exist_ok=True)
# scan_index = 0
sns.set_context("paper", font_scale=1, rc={
    "lines.markersize": 4,
    "lines.linewidth": 2
}
                )

il2_threshold = 0.05
il6_threshold = 0.055
infg_threshold = 0.034

global_df: pd.DataFrame = pd.read_hdf(path + "global_df.h5", mode="r")

cell_df: pd.DataFrame = pd.read_hdf(path + "cell_df.h5", mode="r")

stats_df = pd.read_hdf(path + "cell_stats_df.h5", mode="r")
grouped_cells = cell_df.groupby(["type_name", "time_index", "scan_index", "time"], as_index=True)

means = grouped_cells.mean()

means.reset_index(inplace=True)

std = grouped_cells.std()

std.reset_index(inplace=True)

counts = grouped_cells.count()
counts["n"] = counts["id"]
counts = counts.drop(columns=counts.columns.drop(["n"]))
counts.reset_index(inplace=True)


t_max = 24
scan_label = r"$q_s$ (molecules/s)"
scan_label = r"$D$ $(\mu m^2/s)$"
qs = 1
time_label = "time (a.u.)"
c_unit = "nM"
cell_max = 1.2*(1+counts.max()["n"])
scan_scale = np.logspace(-1,1,10)
color_dict = {
    "Tn":"green",
    "Tfh":"red",
    "Th1":"blue",
    "IL-2":"yellow",
    "IL-6":"brown",
    "IFNg":"black"

}

for i,k in enumerate(color_dict.keys()):
    color_dict[k] = sns.color_palette("Dark2",len(color_dict))[i]

# """global plot"""
#
fig, ax = plt.subplots(2, 2, sharex=False, sharey=False)
global_plot(fig, ax[0][0], global_df, "Concentration", "avg. Cytokine ({c})".format(c=c_unit), legend=False)
global_plot(fig, ax[0][1], global_df, "Gradient", "avg. Gradient ({c}/dm)".format(c=c_unit), legend=False)
global_plot(fig, ax[1][0], global_df, "SD", r"SD of nodal values ({c})".format(c=c_unit))
count_plot(fig, ax[1][1], counts, "n", r"Number of cells", ylim=[0, cell_max])
plt.tight_layout()
plt.savefig(IMGPATH + "global_plot.pdf")
plt.show()

"""scan plot """
fig, ax = plt.subplots(2, 2, sharex=True, sharey=False)
steady_state_plot(fig, ax[0][0], global_df, "time_index", "scan_index", "Concentration", scan_label,
                  "t={t_max}. Cytokine ({c})".format(c=c_unit, t_max=t_max), hue="field_name")
maximum_global_plot(fig, ax[1][0], global_df, "Gradient", scan_label, "max. Gradient ({c}/dm)".format(c=c_unit),
                    legend=False)
maximum_global_plot(fig, ax[1][1], global_df, "SD", scan_label, "max. SD of Nodal Values".format(c=c_unit),
                    legend=False)
steady_state_plot(fig, ax[0][1], counts, "time", "scan_index", "n", scan_label,
                  "t={t_max} Number of cells".format(t_max=t_max), leg_loc="upper left")

plt.tight_layout()
plt.savefig(IMGPATH + "maximum_plot.pdf")
plt.show()

score_max = 30
"""Cluster Plot"""
fig, ax = plt.subplots(2, 2, sharex=True, sharey=False)
scan_score_plot(fig, ax[0][0], cell_df, "Tfh_score_norm", r"Tf$_h$ score (normalized)", legend=False,
                ylim=[0, score_max])
scan_score_plot(fig, ax[1][0], cell_df, "Th1_score_norm", r"Th$_1$ score (normalized)", legend=False,
                ylim=[0, score_max])
scan_score_plot(fig, ax[0][1], cell_df, "Tn_score_norm", r"T$_n$ score (normalized)", legend="brief",
                ylim=[0, score_max])
count_plot(fig, ax[1][1], counts, "n", r"Number of cells", legend=False, ylim=[0, cell_max])
plt.tight_layout()
plt.savefig(IMGPATH + "score_plot.pdf")
plt.show()

"""Cluster Plot"""
fig, ax = plt.subplots(2, 2, sharex=True, sharey=False)
scan_score_plot(fig, ax[0][0], cell_df, "Tfh_score_init_norm", r"Tf$_h$ score (normalized)", legend=False,ylim=[0,score_max])
scan_score_plot(fig, ax[1][0], cell_df, "Th1_score_init_norm", r"Th$_1$ score (normalized)", legend=False,ylim=[0,score_max])
scan_score_plot(fig, ax[0][1], cell_df, "Tn_score_init_norm", r"T$_n$ score (normalized)", legend="brief",ylim=[0,score_max])
count_plot(fig, ax[1][1], counts, "n", r"Number of cells", legend=False,ylim=[0,cell_max])
plt.tight_layout()
plt.savefig(IMGPATH+"score_plot.pdf")
plt.show()

"""joint plot"""

fig, ax = plt.subplots(2, 2, sharex=True, sharey=True)
plot_7(fig, ax[0][0], cell_df, t=1, legend=False)
plot_7(fig, ax[0][1], cell_df, t=3, legend=False)
plot_7(fig, ax[1][0], cell_df, t=6, legend=False)
plot_7(fig, ax[1][1], cell_df, t=10)
plt.tight_layout()
plt.savefig(IMGPATH + "hist_plot.pdf")
plt.show()

plot_8(None, None)
plt.savefig(IMGPATH + "hist.pdf")
plt.show()

# fig, ax = plt.subplots(2, 2, sharex=True, sharey=True)
# sns.lineplot(x = "time_index", y = "IL-2_surf_c", hue = "type_name", data=cell_df ,ax=ax[0][0],units = "id", estimator=None, linewidth = 0.1)
# sns.lineplot(x = "time_index", y = "IL-6_surf_c", hue = "type_name", data=cell_df ,ax=ax[0][1],units = "id", estimator=None, linewidth = 0.1)
# sns.lineplot(x = "time_index", y = "IFNg_surf_c", hue = "type_name", data=cell_df ,ax=ax[1][0],units = "id", estimator=None, linewidth = 0.1)
# plt.savefig(IMGPATH+"trajecories.pdf")
# plt.show()

# plt.figure()
# sns.lineplot(x = "time_index", y = "x", hue = "type_name", data=cell_df.loc[
#     # (cell_df["type_name"] == "Tfh")&
#     (
#         (cell_df["scan_index"] == 0)|
#         (cell_df["scan_index"] == 0)
#     )
# ],units="id",estimator=None)
# plt.show()
