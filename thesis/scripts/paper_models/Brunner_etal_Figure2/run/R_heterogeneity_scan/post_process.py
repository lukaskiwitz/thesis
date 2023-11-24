import os
import sys

import fenics as fcs
import matplotlib.pyplot as plt
import numpy as np

from parameters import path
from thesis.main.PostProcess import PostProcessor, PostProcessComputation


if len(sys.argv) == 3:
    path = str(sys.argv[1])
    n_processes = int(sys.argv[2])
elif len(sys.argv) == 4:
    path = str(sys.argv[1])
    n_processes = int(sys.argv[2])
    time_index = int(sys.argv[3])
else:
    assert False

class my_post_process(PostProcessComputation):
    """Defines a custom computation"""

    def __init__(self):
        """this name will appear in output dataframe"""
        self.name = "Gradient"

    def __call__(self, u, grad, c_conv, grad_conv, mesh_volume, **kwargs):
        V = u.function_space()
        mesh = V.mesh()
        degree = V.ufl_element().degree()
        W = fcs.VectorFunctionSpace(mesh, 'P', degree)
        grad_u = fcs.project(grad(u), W).vector().array().reshape(-1,3)
        grad_u = grad
        # h1: float = fcs.sqrt(fcs.assemble(fcs.dot(grad*grad_conv, grad*grad_conv)  * fcs.dX)) / mesh_volume

        g = np.reshape(grad.vector().vec().array, (u.vector().vec().size, 3))
        g = np.transpose(g)

        my_grad = np.sqrt(np.power(g[0], 2) + np.power(g[1], 2) + np.power(g[2], 2))
        my_grad = np.mean(my_grad) * grad_conv
        return my_grad


"""
setting filepath to look for sim results. This is setup so that it works on the itb computers.
"""
def to_rgb(h):
    return [int(h[2*i:2*i+2],16)/255 for i in range(3)]


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
            "orientation_axes":False,
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
            "color_bar_range": [0,0.1],  # concentration range for color bar (nM), can be commented out.single value sets lower bound.
            # "color_bar_range": [0.01,0.04],  # concentration range for color bar (nM), can be commented out.
            "opacity_range": [0.005, 0.1],  # opacity range for volume rendering (nM), can be commented out
            "volume_camera_pos": [900, 70, 45],
            # camera position in spherical coords. around center of data bounding box
            "volume_raytracing": True,  # Raytracing on off
            "volume_raytracing_progressive_passes": 0,  # better image quality; use with caution
            "volume_raytracing_samples": 4,  # better image quality
            "volume_raytracing_ambient_samples":4,
            "volume_raytracing_light_scale":1,
            "marker_view_uniform_opacity": False,  # draw cell type marker with uniform opacity in marker image
            "lookup": {
                "1": ["Th", to_rgb("cbcbcbff"), 0.5], # mesh_function_value:[name,(r,g,b),opacity] opacity only works with raytracing
                "2": ["Tsec", to_rgb("f2643bff"), 1],
                "3": ["abs", to_rgb("417cffff"), 0.5],
                "4": ["tregs", to_rgb("0f00b0"), 0.5]
            },
        }
}

# pp.computations = (ParaviewRender)
# pp.computations = []

"""carries out the operations in pp.computations in parallel and stores the result in xml file"""
if "time_index" in list(globals().keys()):
    pp.run_post_process(n_processes, time_indices=[time_index])
else:
    pp.run_post_process(n_processes)
