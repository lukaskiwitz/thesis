import os

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import seaborn as sns

from parameters import path


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


IMGPATH = path
os.makedirs(IMGPATH, exist_ok=True)
# scan_index = 0
sns.set_context("paper", font_scale=1, rc={
    "lines.markersize": 4,
    "lines.linewidth": 2
}
                )

global_df: pd.DataFrame = pd.read_hdf(path + "global_df.h5", mode="r")

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
# color_dict = {
#     "default": "blue",
#     "sec": "red",
#     "abs": "green",
#     "IL-2": "yellow",
#     "IL-6": "brown",
#     "IFNg": "black",
#     0:"red",
#     1:"blue"
# }
color_dict = {
    0:"red",
    1:"blue"
}

for i, k in enumerate(color_dict.keys()):
    color_dict[k] = sns.color_palette("muted", len(color_dict))[i]

palette=color_dict

# """global plot"""
# fig, ax = plt.subplots(2, 2, sharex=False, sharey=False)
# global_plot(fig, ax[0][0], global_df, "Concentration", "avg. Cytokine ({c})".format(c=c_unit), legend=False)
# global_plot(fig, ax[0][1], global_df, "Gradient", r"avg. Gradient ({c}/$\mu m$)".format(c=c_unit), legend=False)
# # global_plot(fig, ax[0][1], global_df, "fast_grad", r"avg. Gradient ({c}/$\mu m$)".format(c=c_unit), legend=False)
# global_plot(fig, ax[1][0], global_df, "SD", r"SD of nodal values ({c})".format(c=c_unit))
# count_plot(fig, ax[1][1], counts, "n", r"Number of cells", ylim=[0, cell_max])
# plt.tight_layout()
# plt.savefig(IMGPATH + "global_plot.pdf")
# plt.show()

# if len(global_df["scan_index"].unique()) > 1:
#     """scan plot """
#     fig, ax = plt.subplots(2, 2, sharex=True, sharey=False)
#     steady_state_plot(fig, ax[0][0], global_df, "time_index", "scan_index", "Concentration", scan_label,
#                       "t={t_max}. Cytokine ({c})".format(c=c_unit, t_max=t_max), hue="field_name")
#     steady_state_plot(fig, ax[1][0], global_df, "time_index", "scan_index", "Gradient", scan_label,
#                       "t={t_max} Gradient ({c}/$\mu m$)".format(c=c_unit, t_max=t_max), hue="field_name", legend=False)
#     steady_state_plot(fig, ax[1][1], global_df, "time_index", "scan_index", "SD", scan_label,
#                       "t={t_max} SD of Nodal Values".format(c=c_unit, t_max=t_max), hue="field_name", legend=False)
#     steady_state_plot(fig, ax[0][1], counts, "time", "scan_index", "n", scan_label,
#                       "t={t_max} Number of cells".format(t_max=t_max), leg_loc="upper left")
#
#     plt.tight_layout()
#     plt.savefig(IMGPATH + "maximum_plot.pdf")
#     plt.show()
#
# score_max = 10
# """Cluster Plot"""
# fig, ax = plt.subplots(2, 2, sharex=True, sharey=False)
# scan_score_plot(fig, ax[0][0], cell_df, "abs_score_norm", r"abs score t-1(normalized)", legend=False,
#                 ylim=[0, score_max])
# scan_score_plot(fig, ax[1][0], cell_df, "sec_score_norm", r"sec score t-1 (normalized)", legend=False,
#                 ylim=[0, score_max])
#
# scan_score_plot(fig, ax[0][1], cell_df, "abs_score_init_norm", r"abs score t=0 (normalized)", legend=False,
#                 ylim=[0, score_max])
# scan_score_plot(fig, ax[1][1], cell_df, "sec_score_init_norm", r"sec score t=0 (normalized)", legend="full",
#                 ylim=[0, score_max])
#
# plt.tight_layout()
# plt.savefig(IMGPATH + "score_plot.pdf")
# plt.show()

# l = cell_df["scan_index"].unique()
# rows = int(np.ceil(len(l) / 2))
# fig, ax = plt.subplots(rows, 2, figsize=(15, rows * 4), sharey=True, sharex=True)
# ax = ax.flatten()
#
# for i, v in enumerate(l):
#     ax[i].set_title("scan_index: {v}".format(v=v))
#     legend = "brief"
#     sns.lineplot(ax=ax[i], x="time", y="IL-2_surf_c", data=cell_df.loc[cell_df["scan_index"] == v], hue="type_name",
#                  units="id", estimator=None, linewidth=0.5, legend=legend, palette=color_dict)
#
# plt.savefig(IMGPATH + "trajectory.pdf")
# plt.show()
#


linear_cell_df = cell_df.loc[cell_df["scan_index"] == 0]
sat_cell_df = cell_df.loc[cell_df["scan_index"] == 1]

fig, ax = plt.subplots(3, 2, figsize = (10,10))
ax = np.ravel(ax)

ax_l = ax[0]
sns.barplot(x="scan_index", y="Concentration", data=global_df, ax=ax_l,palette=palette
)
ax_l.set_ylabel("avg. Cytokine ({c})".format(c=c_unit))
ax_l.set_xticklabels(["linear","saturated"])


ax_l = ax[1]
sns.barplot(x="scan_index", y="Gradient", data=global_df, ax=ax_l,palette=palette)
ax_l.set_ylabel(r"avg. Gradient ({c}/$\mu m$)".format(c=c_unit))
ax_l.set_xticklabels(["linear","saturated"])

ax_l = ax[2]
sns.barplot(x="scan_index", y="SD", data=global_df, ax=ax_l,palette=palette)
ax_l.set_ylabel(r"SD of nodal values ({c})".format(c=c_unit))
ax_l.set_xticklabels(["linear","saturated"])

ax_l = ax[4]
sns.distplot(linear_cell_df["IL-2_surf_c"],ax=ax_l,color=palette[0])
sns.distplot(sat_cell_df["IL-2_surf_c"],ax=ax_l,color=palette[1])
ax_l.legend(["linear","saturated"])
ax_l.set_xlabel("cell surface concentration (nM)")

ax_l = ax[3]
sns.barplot(x = "scan_index", y = "surf_c_std", data=global_df,ax=ax_l,palette=palette)
ax_l.set_xticklabels(["linear","saturated"])
ax_l.set_ylabel("cell surface std")

# ax_l = ax[6]
# ax_l.imshow(scan_0_img)
# ax_l.axis("off")
#
# ax_l = ax[7]
# ax_l.imshow(scan_1_img)
# ax_l.axis("off")




# scan_0_img = mpimg.imread(path+"images/scan_0/IL-2/0.png")
# scan_1_img = mpimg.imread(path+"images/q_scan/IL-2/0.png")
#
# fig,ax = plt.subplots(1,2,figsize = (10,5))
# ax = np.ravel(ax)
# ax_l = ax[0]
# ax_l.set_title("linear")
# ax_l.imshow(scan_0_img)
# ax_l.axis("off")
# ax_l = ax[1]
# ax_l.set_title("saturated")
# ax_l.imshow(scan_1_img)
# ax_l.axis("off")
# plt.tight_layout(
# )
# plt.savefig(IMGPATH + "2d.pdf",dpi = 600)
# plt.show()

def f_amax(u,amax, Kc):

    return amax * u / (Kc + u)

def R_sat(u,cell):

    k_on = cell["IL-2_k_on"]/60**2
    R = cell["IL-2_R"]
    Kc = cell["IL-2_Kc"]


    return k_on * R * Kc * u / (Kc + u)

def linear(u,cell):

    k_on = cell["IL-2_k_on"]/60**2
    R = cell["IL-2_R"]

    return k_on* R * u


type_means = cell_df.groupby("type_name").mean()

x = np.linspace(cell_df["IL-2_surf_c"].min(), cell_df["IL-2_surf_c"].max(), 100)

ax_l = ax[5]
for i,t in type_means.iterrows():
    y = linear(x,t)

    ax_l.plot(x,y,color=palette[0],linewidth = 1)

for i,t in type_means.iterrows():
    y = R_sat(x,t)

    ax_l.plot(x,y,color=palette[1],linewidth = 1)


cell_df["uptake"] = -cell_df["IL-2_boundary"] + cell_df["IL-2_q"]

sns.scatterplot(x = "IL-2_surf_c", y = "uptake", data=cell_df,hue="scan_index",ax=ax_l,palette=palette,s =50,style="type_name")
ax_l.set_ylim([0,1.5*cell_df["uptake"].max()])

handles, labels = ax_l.get_legend_handles_labels()

ax_l.legend(handles, labels_replace(labels,rep={
    "0":"linear",
    "1":"saturated",
    "type_name":"Cell type",
    "scan_index":"uptake kinetic"
}
), loc="upper right")

ax_l.set_ylabel(r"uptake term $(\frac{N}{s})$")
ax_l.set_xlabel(r"cell surface concentration (nM)")

plt.show()
plt.tight_layout()
plt.savefig(IMGPATH + "comp.pdf")
