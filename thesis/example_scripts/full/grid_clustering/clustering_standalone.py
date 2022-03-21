import time

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd

from thesis.cellBehaviourUtilities.bridson_sampling import bridson
from thesis.cellBehaviourUtilities.grid_clustering import make_clusters

"""
This script demonstrates the usage of the grid clustering utilities, without an actual simulation.

"""
s = 500
"""get cell grid. In simulation this would be generated from entity positions."""
# cell_grid_positions = get_cell_grid_positions(220, 220, 220, distance=10)
cell_grid_positions = bridson(30, [0, 0], [s, s], density_function=lambda x: 12)

"""Creates some apcs"""
# apcs = get_apc_positions(cell_grid_positions, no_apcs=4)
# apcs = np.array([[100, 100, 100]])

apcs = bridson(30, [20, 20], [s - 20, s - 20], density_function=lambda x: 80)

"""fractions for each cell type and corresponding cluster strengths"""
fractions = [0.94, 0.05, 0.01]
cluster_strengths = [0, 0.9, 1]

start = time.time()
cell_type = make_clusters(cell_grid_positions, apcs, fractions, cluster_strengths)
print(time.time() - start)

# df = pd.DataFrame(
#     {"x": cell_grid_positions[:, 0], "y": cell_grid_positions[:, 1], "z": cell_grid_positions[:, 2], "type": cell_type})
df = pd.DataFrame(
    {"x": cell_grid_positions[:, 0], "y": cell_grid_positions[:, 1], "z": 0, "type": cell_type})

"""creates 3D plot"""
names = ["apc", "naive", "sec", "treg", "stuff"]
colors = ["gray", "gray", "blue", "red"]
fig = plt.figure(figsize=(10, 10))
ax = fig.add_subplot(projection='3d')
ax.set_xlim([df["x"].min(), df["x"].max()])
ax.set_ylim([df["y"].min(), df["y"].max()])
ax.set_zlim([df["z"].min(), df["z"].max()])
ax.scatter(apcs[:, 0], apcs[:, 1], s=200, color="green")
for i, d in df.groupby(["type"]):
    if np.isin(i, [1, 2, 3]):
        ax.scatter(d["x"], d["y"], d["z"], s=50, color=colors[int(d["type"].unique())])

plt.legend(names)
plt.show()
print(df.groupby("type").count()["x"] / len(df))

""""""

# from sklearn.metrics import silhouette_samples
# import seaborn as sns
# s = 300
# # cell_grid_positions = get_cell_grid_positions(s, s, s)
#
# cell_grid_positions = bridson(20,[0,0,0],[s,s,s], density_function=lambda x: 20)
# # apcs = bridson(20,[0,0,0],[s,s,s], density_function=lambda x: 170)
#
# apcs = np.array([
#     np.mean(cell_grid_positions,axis = 0),
# ])
#
# # apcs = np.array([
# #     [130,130,130],
# #     [260,260,260]
# #
# # ])
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
#
# sns.lineplot(x = "cs", y = "sill", hue = "type", data=df)
# plt.legend(["background","tregs"])
# plt.xlabel("clustering strength")
# plt.ylabel("bg. vs. treg silhouette score")
# plt.savefig("/home/kiwitz/Desktop/silhouette_200.pdf")
# plt.show()
#
#
# fig,ax = plt.subplots(n,n,figsize = (10,10))
# ax = np.ravel(ax)
#
# for i,cs in enumerate(df.cs.unique()):
#     ax[i].set_title(str(round(cs,2)))
#     ax[i].set_xlim([0,s])
#     sns.distplot(df.loc[(df.cs == cs) & (df.type == 2)]["x"], bins = 20,ax = ax[i])
#     sns.distplot(df.loc[(df.cs == cs) & (df.type == 2)]["y"], bins=20, ax=ax[i])
#     sns.distplot(df.loc[(df.cs == cs) & (df.type == 2)]["z"], bins=20, ax=ax[i])
# plt.tight_layout()
# plt.show()
