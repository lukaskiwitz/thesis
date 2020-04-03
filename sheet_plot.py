import getpass

import fenics as fcs
import matplotlib.pyplot as plt
import numpy as np
import mpi4py.MPI as MPI
import pandas as pd
from matplotlib.animation import FuncAnimation, FFMpegWriter
import matplotlib.animation as animation
import matplotlib  as mpl
from my_debug import message
import os
from PostProcess import PostProcessor

user = getpass.getuser()
# path = "/extra/{u}/scan_example_small/".format(u=user)
path = "/extra/{u}/dev_testing/".format(u=user)
ext_cache="/extra/{u}/dev_testing_ext_cache/".format(u=user)




def field_update(t,scan_index,field_name):
    t = t+1
    mesh = fcs.Mesh()
    with fcs.XDMFFile(ext_cache+"mesh" + ".xdmf") as f:
        f.read(mesh)
    mesh = mesh
    V = fcs.FunctionSpace(mesh,"P",1)

    field = path + 'scan_{s}/sol/distplot/field_{f}_{t}_distPlot.h5'.format(t=t,f=field_name,s=scan_index)

    u: fcs.Function  = fcs.Function(V)
    with fcs.HDF5File(MPI.COMM_WORLD, field, "r") as f:
        f.read(u, "/" + field_name)

    u.set_allow_extrapolation(True)
    fig.clear()

    rec_mesh = fcs.RectangleMesh(fcs.Point(-1,-1),fcs.Point(1,1),100,100)
    rev_V = fcs.FunctionSpace(rec_mesh,"P",1)
    u_slice = fcs.interpolate(u,rev_V)*1e9
    p = plt.colorbar(fcs.plot(u_slice))


    # cbar = mpl.colorbar.ColorbarBase(plt.gca(), cmap=cmap, norm=mpl.colors.Normalize(vmin=0.06, vmax=0.08))
    return p

cell_df: pd.DataFrame = pd.read_hdf(path+"cell_df.h5",mode="r")
cell_df = cell_df.loc[cell_df["z"] == 0]

def update(t,scan_df,field_name):
    # cell_p = []
    # cell_t = []
    plt.gca().remove()
    t_df = scan_df.loc[scan_df["time_index"] == t]
    field_update(t, 0, field_name)
    for row_i, row in t_df.iterrows():
        n = row["type_name"]
        if True:
            p = [
                row["x"],
                row["y"]
            ]
            color = "blue"
            if n == "Tfh":
                color = "red"
            elif n == "Th1":
                color = "green"
            c = plt.Circle(p, 0.05, color=color)
            plt.gca().add_artist(c)
    plt.xlim([-1,1])
    plt.ylim([-1, 1])

for field_name in ["IL-2","IL-6","IFNg"]:
    for s in cell_df["scan_index"].unique():
        save_path = path+"render/scan_{s}/".format(s=int(s))
        os.makedirs(save_path,exist_ok=True)

        message("writing scan to {p}".format(p=save_path))
        scan_df = cell_df.loc[cell_df["scan_index"] == s]
        t_max = 24#scan_df["t"].max()
        message("t_max: {t}".format(t=t_max))
        fig = plt.figure()
        ax = plt.gca()
        ani = FuncAnimation(fig, update, frames=range(t_max),fargs=(scan_df,field_name))
        ani.save(save_path+"anim_{f}.mp4".format(f=field_name), fps = 1)

