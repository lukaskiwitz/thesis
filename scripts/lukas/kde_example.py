import json
import multiprocessing as mp
import os
from copy import deepcopy
from math import ceil
from typing import List, Dict

import KDEpy
import fenics as fcs
import lxml.etree as et
import mpi4py.MPI as MPI
import numpy as np
import pandas as pd
from scipy.constants import N_A
import scipy
import matplotlib.pyplot as plt
import seaborn as sns

def plot(ax1, ax2, ax3 , K):

    kernel = KDEpy.FFTKDE(K, norm=2, bw=bw).fit([0])
    grid_points = 100
    grid, y = kernel.evaluate(grid_points)
    x = np.unique(grid)
    # y = points.reshape(grid_points).T
    ax1.plot(x,y)





    kernel = KDEpy.FFTKDE(K, norm = 2, bw=bw).fit(data)
    grid_points = 1000
    grid, points = kernel.evaluate(grid_points)
    x, y = np.unique(grid[:, 0]), np.unique(grid[:, 1])
    z = points.reshape(grid_points, grid_points).T

    ax2.contourf(x,y,z,100)
    # ax1.colorbar()
    ax2.scatter(data[:,0],data[:,1])
    # ax2.set_xlim([-s,s])
    # ax2.set_ylim([-s,s])

    ax3.hist(np.ravel(z),bins=10)

s = 500
N = 10
bw = 50
data = np.random.uniform(-s,s,(N,2))
IMGPATH = "/home/lukas/lab_meeting_images/KDE/"
fig ,ax = plt.subplots(3,3)

plot(ax[0][0],ax[0][1],ax[0][2],K="gaussian")
plot(ax[1][0],ax[1][1],ax[1][2],K="box")
plot(ax[2][0],ax[2][1],ax[2][2],K="tri")

ax[0][0].set_title("Kernel Function")
ax[0][1].set_title("Resulting Density \n Estimator")
ax[0][2].set_title("Histogram of \n point values")
plt.tight_layout()
os.makedirs(IMGPATH,exist_ok=True)
# plt.savefig(IMGPATH+"subplot.pdf",dpi=200)
plt.show()