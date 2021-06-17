import os
import sys
import fenics as fcs

sys.path.append("/home/lukas/thesis/main/")
sys.path.append("/home/lukas/thesis/scenarios/")

import numpy as np
from parameters import path

os.environ["LOG_PATH"] = path

from thesis.main.PostProcess import PostProcessor, PostProcessComputation


def to_rgb(h):
    return [int(h[2*i:2*i+2],16)/255 for i in range(3)]

class my_post_process(PostProcessComputation):
    """Defines a custom computation"""

    def __init__(self):
        """this name will appear in output dataframe"""
        self.name = "h1_norm"

    def __call__(self, u, grad, c_conv, grad_conv, mesh_volume, **kwargs):
        """Calculates the h1 norm, another measure for inhomogeneity"""
        h1: float = fcs.sqrt(fcs.assemble(fcs.dot(grad*grad_conv, grad*grad_conv)  * fcs.dX)) / mesh_volume
        return h1


"""number of threads can be passed as first cli argument"""
if len(sys.argv) > 1:
    threads = int(sys.argv[1])
else:
    threads = 4

"""
setting filepath to look for sim results. This is setup so that it works on the itb computers.
"""

pp = PostProcessor(path)
pp.unit_length_exponent = -6


"""appends a custom calculation.
default computations are defined in PostProcess.py"""
pp.computations.append(my_post_process())

"""Imaging settings"""
pp.visual_conversion  = 1e3# converts concentrations to pM for paraview visualisation


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
            "axis_title_font_size": 15,  # slice label
            "axis_label_font_size": 15,  # slice label
            "number_format": '%2.2f',  # color bar number format
            "axis_ticks": True,  # slice view
            "axis_edges": True,  # slice view
            "field_color_preset": "Blue Orange (divergent)",  # concentration color scheme; paraview preset names
            # "color_bar_range": [0,0.05],  # concentration range for color bar (nM), can be commented out.single value sets lower bound.
            # "color_bar_range": [0.01,0.04],  # concentration range for color bar (nM), can be commented out.
            # "opacity_range": [0.001, 0.05],  # opacity range for volume rendering (nM), can be commented out
            "volume_camera_pos": [800, 70, -110],# camera position in spherical coords. around center of data bounding box
            "volume_raytracing": True,  # Raytracing on off
            "volume_raytracing_progressive_passes": 0,  # better image quality; use with caution
            "volume_raytracing_samples": 8,  # better image quality
            "volume_raytracing_ambient_samples":2,
            "volume_raytracing_light_scale":1,
            "marker_view_uniform_opacity": True,  # draw cell type marker with uniform opacity in marker image
            "lookup": {
                "1": ["default", to_rgb("cbcbcbff"), 0.5], # mesh_function_value:[name,(r,g,b),opacity] opacity only works with raytracing
                "2": ["sec", to_rgb("f2643bff"), 1],
                "3": ["abs", to_rgb("417cffff"), 0.5],
            },
        }
}

"""overwrites settings for individual cytokines (redundant for a single cytokine)"""
pp.image_settings_fields = {
    "IL-2":
        {
            "paraview_settings":
                {
                "field_color_preset": "Blue Orange (divergent)",
                # "color_bar_range": [0,30],
                # "opacity_range": [1, 5],
                }
        }

}

"""carries out the operations in pp.computations in parallel and stores the result in xml file"""
pp.run_post_process(threads, render_paraview=False, make_images=False, kde=False)