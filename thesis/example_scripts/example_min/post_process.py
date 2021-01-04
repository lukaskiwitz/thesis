import os
import sys

sys.path.append("/home/lukas/thesis/main/")
sys.path.append("/home/lukas/thesis/scenarios/")

import numpy as np
from parameters import path

os.environ["LOG_PATH"] = path

from thesis.main.PostProcess import PostProcessor, PostProcessComputation


class my_post_process(PostProcessComputation):
    """Defines a custom computation"""

    def __init__(self):
        """this name will appear in output dataframe"""
        self.name = "fast_grad"

    def __call__(self, u, grad, c_conv, grad_conv, mesh_volume, **kwargs):
        # gradient: float = fcs.assemble(fcs.sqrt(fcs.dot(grad, grad)) * fcs.dX) * grad_conv / mesh_volume

        g = np.reshape(grad.vector().vec().array, (u.vector().vec().size, 3))
        g = np.transpose(g)

        my_grad = np.sqrt(np.power(g[0], 2) + np.power(g[1], 2) + np.power(g[2], 2))
        my_grad = np.mean(my_grad) * grad_conv

        return my_grad


"""number of n_processes can be passed as first cli argument"""
if len(sys.argv) > 1:
    threads = int(sys.argv[1])
else:
    threads = os.cpu_count()

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
    "round_legend_labels": 4,
    "legend_title": "",
    "dpi": 350,
    "paraview_settings":
        {
            "orientation_axes":False,
            "render_view_size": [400, 400],  # images size
            "cell_type_title_font_size": 20,  # legend font
            "cell_type_label_font_size": 20,  # legend font
            "slice_zoom": 1.3,  # orthographic "zoom" from slice
            "slice_origin": [0, 0, 0],  # origin for slice from center of bounding box
            "slice_normal": [0, 0, 1],  # normal direction for slice plane; also defines camera perspective
            "layer_distance": 20,  # distance of cell layer in simulation
            "axis_title_font_size": 15,  # slice label
            "axis_label_font_size": 15,  # slice label
            "number_format": '%2.2f',  # color bar number format
            "axis_ticks": True,  # slice view
            "axis_edges": True,  # slice view
            "field_color_preset": "Cool to Warm",  # concentration color scheme; paraview preset names
            # "color_bar_range": [0.01],  # concentration range for color bar (nM), can be commented out.single value sets lower bound.
            # "color_bar_range": [0.01,0.04],  # concentration range for color bar (nM), can be commented out.
            # "opacity_range": [0, 0.1],  # opacity range for volume rendering (nM), can be commented out
            "volume_camera_pos": [600, 70, -110], # camera position in spherical coords.  around center of data bounding box
            "volume_raytracing": True,  # Raytracing on off
            "volume_raytracing_progressive_passes": 0,  # better image quality; use with caution
            "volume_raytracing_samples": 2,  # better image quality
            "marker_view_uniform_opacity": True,  # draw cell type marker with uniform opacity in marker image
            "lookup": {
                "1": ["default", (0.8, 0.8, 0.8), 0.25], # mesh_function_value:[name,(r,g,b),opacity]; opacity only works with raytracing
                "2": ["sec", (0.8, 0.2, 0.2), 1],
                "3": ["abs", (0.2, 0.2, 0.8), 1]
            },
        }
}
# pp.computations.append(my_post_process())

"""carries out the operations in pp.computations in parallel and stores the result in xml file"""
pp.run_post_process(threads, render_paraview=True)
