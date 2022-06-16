import os

import matplotlib.pyplot as plt
import numpy as np

from loader import load_csv
from model import my_nn, my_nn_plot

path = "../../data/hauser/"
IMGPATH = "/home/lukas/hauser_img_kde/nearest_neighbour/"
os.makedirs(IMGPATH, exist_ok=True)

r_max = 50
n_neighbours = 50

cells = load_csv(path)
result = {}
for type_name in cells["type_name"].unique():
    fig, ax = plt.subplots(1, 2, figsize=(6, 3))
    ax = np.ravel(ax)
    X, df, kg_con, kg_dist = my_nn(cells, type_name, r_max=r_max, n_neighbours=n_neighbours)
    my_nn_plot(X, df, kg_con, kg_dist, ax_list=ax, markersize=2, edgewidth=0.1)
    ax[0].set_title(type_name + "_connectivity matrix")
    ax[1].set_title(type_name + "_connectivity graph")
    plt.tight_layout(
    )
    plt.savefig(IMGPATH + type_name + "_img.pdf")

    result[type_name] = [X, df, kg_con, kg_dist]

fig, ax = plt.subplots(1, 2, figsize=(6, 3))
ax = np.ravel(ax)

x = []
for k, v in result.items():
    kg_con = v[2]
    x.append(np.count_nonzero(kg_con, axis=0))
ax[0].boxplot(x)
ax[0].set_xticklabels(list(result.keys()))
ax[0].set_ylabel("mean number of neighbors")

x = []
for k, v in result.items():
    kg_con = v[2]
    x.append(kg_dist[kg_dist != 0])

ax[1].boxplot(x)
ax[1].set_xticklabels(list(result.keys()))
ax[1].set_ylabel("mean distance to neighbour")
plt.savefig(IMGPATH + "stats.pdf")
