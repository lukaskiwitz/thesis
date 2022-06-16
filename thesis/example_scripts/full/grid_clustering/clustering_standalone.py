import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from matplotlib.colors import to_rgb

from thesis.cellBehaviourUtilities.bridson_sampling import bridson
from thesis.cellBehaviourUtilities.grid_clustering import make_clusters

"""
This script demonstrates the usage of the grid clustering utilities, without an actual simulation.

# # """


def get_cluster_df(fractions=[0.9, 0.1], cluster_strengths=[0, 1], s=150, cell_d=12, apc_d=120, N=3):
    assert len(fractions) == len(cluster_strengths)
    apc_margin = apc_d / 2
    if N == 3:
        cell_grid_positions = bridson(30, [0, 0, 0], [s, s, s], density_function=lambda x: cell_d)
        apcs = bridson(30, [apc_d / 2, apc_d / 2, apc_d / 2], [s - apc_d / 2, s - apc_d / 2, s - apc_d / 2],
                       density_function=lambda x: apc_d)
        cell_type = make_clusters(cell_grid_positions, apcs, fractions, cluster_strengths)
        df = pd.DataFrame(
            {"x": cell_grid_positions[:, 0], "y": cell_grid_positions[:, 1], "z": cell_grid_positions[:, 2],
             "type": cell_type})
    elif N == 2:
        cell_grid_positions = bridson(30, [0, 0], [s, s], density_function=lambda x: cell_d)
        apcs = bridson(30, [apc_d / 2, apc_d / 2], [s - apc_d / 2, s - apc_d / 2], density_function=lambda x: apc_d)
        cell_type = make_clusters(cell_grid_positions, apcs, fractions, cluster_strengths)
        df = pd.DataFrame(
            {"x": cell_grid_positions[:, 0], "y": cell_grid_positions[:, 1], "z": 0,
             "type": cell_type})
    else:
        raise ValueError("N must be 2 or 3 dimension")
    # apcs = get_apc_positions(cell_grid_positions, no_apcs=4)
    # apcs = np.array([[100, 100, 100]])

    return df, apcs


def vis(df, apcs, file_name="cluster"):
    flat = len(df["z"].unique()) == 1

    names = ["apc", "naive", "sec", "treg", "stuff"]
    colors = ["gray", to_rgb("#a7a7a732"), to_rgb("#f2643b32"), to_rgb("#417cff32")]
    fig = plt.figure(figsize=(1.3, 1.3))
    fig = plt.figure(figsize=(10, 10))

    if flat:
        ax = fig.gca()
    else:
        ax = fig.add_subplot(projection='3d')
    ax.set_xlim([df["x"].min(), df["x"].max()])
    ax.set_ylim([df["y"].min(), df["y"].max()])
    if flat:

        for i, d in df.groupby(["type"]):
            if np.isin(i, [1, 2, 3]):
                ax.scatter(d["x"], d["y"], s=10, color=colors[int(d["type"].unique())])
        ax.scatter(apcs[:, 0], apcs[:, 1], s=250, color="green")

    else:
        ax.set_zlim([df["z"].min(), df["z"].max()])
        for i, d in df.groupby(["type"]):
            if np.isin(i, [2, 3]):
                ax.scatter(d["x"], d["y"], d["z"], s=250, color=colors[int(d["type"].unique())])
        ax.scatter(apcs[:, 0], apcs[:, 1], apcs[:, 2], s=500, color="green")

    ax.tick_params(left=False, right=False, top=False, bottom=False)
    ax.set_yticklabels([])
    ax.set_xticklabels([])
    plt.tight_layout()
    # plt.legend(names)
    print(df.groupby("type").count()["x"] / len(df))


# cs = 1
# for cs in [0,0.5,1]:
#     np.random.seed(0)
#     df, apcs = get_cluster_df(s = 200, fractions = [0.9,0.1],cluster_strengths= [0,cs],cell_d=12,apc_d=120, N = 2)
#     vis(df,apcs)
#     file_name="single_clusters_{cs}".format(cs = cs)
#     plt.savefig("/home/kiwitz/Seafile/System Immunology/spatial paper/{fn}.pdf".format(fn = file_name))
#     plt.show()

cs = 1
np.random.seed(0)
df, apcs = get_cluster_df(s=200, fractions=[0.8, 0.1, 0.1], cluster_strengths=[0, 1, cs], cell_d=12, apc_d=60, N=2)
vis(df, apcs)
file_name = "two_clusters_{cs}".format(cs=cs)
# plt.savefig("/home/kiwitz/Seafile/System Immunology/spatial paper/{fn}.pdf".format(fn = file_name))
plt.show()

""""""
# s = 200
# cell_grid_positions = get_cell_grid_positions(s, s, s)

# cell_grid_positions = bridson(30,[0,0,0],[s,s,s], density_function=lambda x: 12)

# apcs = np.array([
#     np.mean(cell_grid_positions,axis = 0),
# ])

# apcs = np.array([
#     [130,130,130],
#     [260,260,260]
#
# ])
# fractions = [0.9,0.1]
# assert np.sum(fractions) <= 1
# n = 4
# x = np.linspace(0,1,n**2)
# # x = [1]
# y = []
#
# for cs in x:
#     print("running")
#     cluster_strengths = [0,cs]
#     cell_type = make_clusters(cell_grid_positions, apcs, fractions, cluster_strengths)
#
#     df = pd.DataFrame(
#         {"x": cell_grid_positions[:, 0], "y": cell_grid_positions[:, 1], "z": cell_grid_positions[:, 2],
#          "type": cell_type})
#     d = df.loc[df.type.isin([1,2])]
#     X = np.array([d["x"], d["y"], d["z"]]).T
#     d["sill"] = silhouette_samples(X, labels = d["type"])
#     d["cs"] = cs
#     # y.append(d.groupby("type",as_index=False).mean())
#     y.append(d)
#
# df = pd.concat(y)
#
# rc_ticks = {
#             "lines.linewidth": 0.5,
#             "axes.linewidth": 0.4,
#             "lines.markersize": 3,
#             "xtick.major.size": 2.5,
#             "xtick.major.width":0.5,
#             "xtick.major.pad":1,
#             "xtick.minor.size": 1.5,
#             "xtick.minor.width": 0.5,
#             "xtick.minor.pad": 1,
#             "ytick.major.size": 2.5,
#             "ytick.major.width": 0.5,
#             "ytick.major.pad": 1,
#             "ytick.minor.size": 1.5,
#             "ytick.minor.width": 0.5,
#             "ytick.minor.pad": 1,
#             "axes.labelpad":1,
#             "axes.titlepad":1
#         }
# sns.set_context("paper",font_scale = 0.7)
# plt.rcParams.update(rc_ticks)
# b = 1.3
# a = 7 / 6 * b
# plt.figure(figsize = (a,b))
# ax = plt.gca()
# sns.lineplot(x = "cs", y = "sill", data=df.loc[df.type == 2])
# # plt.legend(["background","tregs"])
# plt.xlabel("clustering strength")
# plt.ylabel("silhouette score")
# plt.tight_layout()
# plt.savefig("/home/kiwitz/Seafile/System Immunology/spatial paper/silhouette_200.pdf")
#
# plt.show()
# #
