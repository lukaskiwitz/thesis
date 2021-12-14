import os
import sys

import fenics as fcs
import matplotlib.pyplot as plt
import numpy as np

from parameters import path
from thesis.main.PostProcess import PostProcessor, PostProcessComputation

"""number of threads can be passed as first cli argument"""
if len(sys.argv) == 2:
    n_processes = int(sys.argv[1])
else:
    n_processes = 4

if len(sys.argv) == 3:
    n_processes = int(sys.argv[1])
    path = str(sys.argv[2])
else:
    pass




class my_post_process(PostProcessComputation):
    """Defines a custom computation"""

    def __init__(self):
        """this name will appear in output dataframe"""
        self.name = "autocorr"

    def __call__(self, u, grad, c_conv, grad_conv, mesh_volume, **kwargs):

        from thesis.main.PostProcessUtil import get_rectangle_plane_mesh
        from scipy.signal import fftconvolve

        path = kwargs["path"]
        ti = kwargs["time_index"]
        si = kwargs["scan_index"]
        cell_df = kwargs["cell_df"]

        dim = (50, 50)
        grid, lim = get_rectangle_plane_mesh(u, res=dim)
        u.set_allow_extrapolation(True)

        grid_V = fcs.FunctionSpace(grid, "P", 1)
        grid_u = fcs.interpolate(u, grid_V)

        values = grid_u.compute_vertex_values().reshape(dim[0]+1,dim[1]+1)
        x = grid.coordinates().T[0].reshape(dim[0]+1,dim[1]+1)[0]
        y = grid.coordinates().T[1].reshape(dim[0]+1,dim[1]+1)[:,0]

        fig = plt.figure()
        effectors = cell_df.loc[cell_df["type_name"] == "sec"]
        c_list = []
        imgpath = path + "autocorr/"
        for i,tsec in effectors.iterrows():

            r = np.array([tsec["x"],tsec["y"]])

            idx = (np.abs(x - r[0])).argmin()

            idy = (np.abs(y - r[1])).argmin()
            r_grid = [x[idx],y[idy]]
            xwindow = np.min([len(x[0:idx]),len(x[idx:-1])])
            ywindow = np.min([len(y[0:idy]), len(y[idy:-1])])

            xmin = idx - xwindow
            xmax = idx + xwindow
            ymin = idy - ywindow
            ymax = idy + ywindow

            v = values[xmin:xmax, ymin:ymax]

            v = v * c_conv

            c = fftconvolve(np.conj(np.flip(v)), v, mode="same")
            # saves in cell centered coordiante frame; might be more convenient
            # np.save(imgpath + "./autocorr_{s}_{t}_{id}".format(s=si, t=ti, id=tsec["id"]), c)

            pad = [
                [xmin, len(x) - xmax]
                , [ymin, len(y) - ymax]
            ]
            c = np.pad(c, pad, mode="constant")
            c = c.T
            # saves in global coordiante frame
            np.save(imgpath + "./autocorr_{s}_{t}_{id}".format(s=si, t=ti, id=tsec["id"]), c)
            c_list.append(c)

        ax = fig.gca()
        ax.contourf(x,y,np.max(c_list,axis=0),200)


        for i,c in effectors.iterrows():
            ax.add_artist(plt.Circle((c["x"],c["y"]), 5, color="red"))



        os.makedirs(imgpath,exist_ok=True)
        fig.savefig(imgpath+"s_{si}_t_{ti}".format(si=si,ti=ti))
        plt.show()



        return 0


"""number of threads can be passed as first cli argument"""
if len(sys.argv) > 1:
    n_processes = int(sys.argv[1])
else:
    n_processes = 4

"""
setting filepath to look for sim results. This is setup so that it works on the itb computers.
"""


def to_rgb(h):
    return [int(h[2 * i:2 * i + 2], 16) / 255 for i in range(3)]


pp = PostProcessor(path)
pp.unit_length_exponent = -6

"""appends a custom calculation.
naive computations are defined in PostProcess.py"""
pp.image_settings = {
    "cell_colors": "Dark2",
    "cell_color_key": "type_name",
    "round_legend_labels": 4,
    "legend_title": "",
    "dpi": 350,
    "paraview_settings":
        {
            "orientation_axes": False,
            "render_view_size": [400, 400],  # images size
            "cell_type_title_font_size": 20,  # legend font
            "cell_type_label_font_size": 20,  # legend font
            "slice_zoom": 2,  # orthographic "zoom" from slice
            "slice_origin": [0, 0, 0],  # origin for slice from center of bounding box
            "slice_normal": [0, 0, 1],  # normal direction, also defines camera perspective
            "layer_distance": 20,  # distance of cell layer in simulation
            "axis_title_font_size": 0,  # slice label
            "axis_label_font_size": 0,  # slice label
            "number_format": '%2.2f',  # color bar number format
            "axis_ticks": False,  # slice view
            "axis_edges": False,  # slice view
            "field_color_preset": "Blue Orange (divergent)",  # concentration color scheme; paraview preset names
            "color_bar_range": [0, 0.1],
            # concentration range for color bar (nM), can be commented out.single value sets lower bound.
            # "color_bar_range": [0.01,0.04],  # concentration range for color bar (nM), can be commented out.
            "opacity_range": [0.005, 0.1],  # opacity range for volume rendering (nM), can be commented out
            "volume_camera_pos": [900, 70, 45],
            # camera position in spherical coords. around center of data bounding box
            "volume_raytracing": True,  # Raytracing on off
            "volume_raytracing_progressive_passes": 0,  # better image quality; use with caution
            "volume_raytracing_samples": 4,  # better image quality
            "volume_raytracing_ambient_samples": 4,
            "volume_raytracing_light_scale": 1,
            "marker_view_uniform_opacity": False,  # draw cell type marker with uniform opacity in marker image
            "lookup": {
                "1": ["naive", to_rgb("cbcbcbff"), 0.5],
                # mesh_function_value:[name,(r,g,b),opacity] opacity only works with raytracing
                "2": ["sec", to_rgb("f2643bff"), 1],
                "3": ["abs", to_rgb("417cffff"), 0.5],
                "4": ["tregs", to_rgb("0f00b0"), 0.5]
            },
        }
}

# pp.computations = (ParaviewRender)
"""carries out the operations in pp.computations in parallel and stores the result in xml file"""
pp.run_post_process(n_processes)
