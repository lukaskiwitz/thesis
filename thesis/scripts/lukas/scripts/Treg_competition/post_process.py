import os
import sys
import matplotlib.pyplot as plt
sys.path.append("/home/lukas/thesis/main/")
sys.path.append("/home/lukas/thesis/scenarios/")

import fenics as fcs
import numpy as np
from parameters import path

os.environ["LOG_PATH"] = path

from thesis.main.PostProcess import PostProcessor, PostProcessComputation


class my_post_process(PostProcessComputation):
    """Defines a custom computation"""

    def __init__(self):
        """this name will appear in output dataframe"""
        self.name = "autocorr"

    def __call__(self, u, grad, c_conv, grad_conv, mesh_volume, **kwargs):

        from PostProcessUtil import get_rectangle_plane_mesh
        from scipy.signal import correlate2d, fftconvolve
        path = kwargs["path"]
        ti = kwargs["time_index"]
        si = kwargs["scan_index"]
        cell_df = kwargs["cell_df"]

        dim = (50,50)
        grid = get_rectangle_plane_mesh(u,res = dim)
        u.set_allow_extrapolation(True)

        grid_V = fcs.FunctionSpace(grid, "P", 1)
        grid_u = fcs.interpolate(u, grid_V)

        values = grid_u.compute_vertex_values().reshape(dim[0]+1,dim[1]+1)
        x = grid.coordinates().T[0].reshape(dim[0]+1,dim[1]+1)[0]
        y = grid.coordinates().T[1].reshape(dim[0]+1,dim[1]+1)[:,0]

        fig = plt.figure()
        effectors = cell_df.loc[cell_df["type_name"] == "effector"]
        c_list = []
        for i,tsec in effectors.iterrows():

            r = np.array([tsec["x"],tsec["y"]])

            idx = (np.abs(x - r[0])).argmin()

            idy = (np.abs(y - r[1])).argmin()
            r_grid = [x[idx],y[idy]]
            xwindow = np.min([len(x[0:idx]),len(x[idx:-1])])
            ywindow = np.min([len(y[0:idy]), len(y[idy:-1])])


            xmin = idx-xwindow
            xmax = idx + xwindow
            ymin = idy - ywindow
            ymax = idy + ywindow

            v = values[xmin:xmax, ymin:ymax]

            v = v*c_conv

            c = fftconvolve(np.conj(np.flip(v)),v, mode="same")
            pad = [
                [xmin,len(x) - xmax]
                ,[ymin,len(y) - ymax]
            ]
            c = np.pad(c,pad,mode="constant")
            c = c.T
            c_list.append(c)

        ax = fig.gca()
        ax.contourf(x,y,np.max(c_list,axis=0),200)



        for i,c in effectors.iterrows():
            ax.add_artist(plt.Circle((c["x"],c["y"]), 5, color="red"))


        imgpath = path+"autocorr/"
        os.makedirs(imgpath,exist_ok=True)
        fig.savefig(imgpath+"s_{si}_t_{ti}".format(si=si,ti=ti))
        plt.show()



        return 0


"""number of threads can be passed as first cli argument"""
if len(sys.argv) > 1:
    threads = int(sys.argv[1])
else:
    threads = 1

"""
setting filepath to look for sim results. This is setup so that it works on the itb computers.
"""

pp = PostProcessor(path)
pp.unit_length_exponent = -6

"""appends a custom calculation.
default computations are defined in PostProcess.py"""
pp.image_settings = {
    "cell_colors": "Dark2",
    "cell_color_key": "type_name",
    "round_legend_labels":4,
    "legend_title":"",
    "dpi":350,
}
pp.computations.append(my_post_process())

"""carries out the operations in pp.computations in parallel and stores the result in xml file"""
pp.run_post_process(threads, make_images=True, kde=True)

pp.global_dataframe.to_hdf(path + 'global_df.h5', key="data", mode="w")
pp.cell_dataframe.to_hdf(path + "cell_df.h5", key="df", mode="w")
