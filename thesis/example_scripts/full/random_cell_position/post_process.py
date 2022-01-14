import os
import sys
import fenics as fcs
import shutil

sys.path.append("/home/brunner/thesis/thesis/main/")
sys.path.append("/home/brunner/thesis/thesis/scenarios/")

import numpy as np
from copy import deepcopy
from parameters import path, ext_cache

print(path)

os.environ["LOG_PATH"] = path

from thesis.main.PostProcess import PostProcessor, PostProcessComputation


class my_post_process(PostProcessComputation):
    """Defines a custom computation"""

    def __init__(self):
        """this name will appear in output dataframe"""
        self.name = "h1_norm"

    def __call__(self, u, grad, c_conv, grad_conv, mesh_volume, **kwargs):
        V = u.function_space()
        mesh = V.mesh()
        degree = V.ufl_element().degree()
        W = fcs.VectorFunctionSpace(mesh, 'P', degree)

        h1: float = fcs.sqrt(fcs.assemble(fcs.dot(grad*grad_conv, grad*grad_conv)  * fcs.dX)) / mesh_volume

        return h1


"""number of threads can be passed as first cli argument"""
if len(sys.argv) > 1:
    threads = int(sys.argv[1])
else:
    threads = 6


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
}

pp.computations.append(my_post_process())

"""carries out the operations in pp.computations in parallel and stores the result in xml file"""
pp.run_post_process(threads, kde=False)
