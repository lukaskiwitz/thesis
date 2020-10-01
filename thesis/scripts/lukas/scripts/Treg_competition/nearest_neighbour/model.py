import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import seaborn as sns
import getpass
import KDEpy
import os
import multiprocessing as mp
# from loader import load_csv
from sklearn.neighbors import NearestNeighbors
import networkx as nx
from matplotlib.lines import Line2D


def my_nn(cells, type_name, n_neighbours = 10, r_max = 100):

    if type_name is None:
        pass
    else:
        cells = cells.loc[cells["type_name"] == type_name]

    x = cells["x"]
    y = cells["y"]

    X = np.array([x,y]).T

    nn = NearestNeighbors(radius=r_max,n_neighbors = n_neighbours, metric="euclidean")
    nn.fit(X)

    # con = nn.kneighbors_graph(mode = "connectivity").toarray()


    dist = nn.kneighbors_graph(mode = "distance").toarray()
    con = np.where(dist != 0 ,1,0)

    con = np.where(dist < r_max, 1, 0) * con
    dist= dist*con
    return X,cells,con,dist,nn

def my_nn_plot(X, cells, con, dist, ax_list = None, edgewidth = 1, edgecolor ="blue", markersize= 10):

    g = nx.Graph(con)
    if ax_list is None:
        fig,ax = plt.subplots(2,2)

        ax = np.ravel(ax)
    else:
        ax = ax_list

    g = nx.Graph(con)

    for e in g.edges:
        x1 = X[e[0]]
        x2 = X[e[1]]

        x = [x1[0], x2[0]]
        y = [x1[1], x2[1]]

        line = Line2D(x, y, linewidth=edgewidth, color=edgecolor)
        ax[1].add_artist(line)

    ax[1].axis("on")
    ax[1].tick_params(left=True, bottom=True, labelleft=True, labelbottom=True)

    ax[0].matshow(con)
    sns.scatterplot(x="x", y="y", data=cells, s=markersize, ax=ax[1])


