import os
import sys
import fenics as fcs
import shutil

sys.path.append("/home/brunner/thesis/thesis/main/")
sys.path.append("/home/brunner/thesis/thesis/scenarios/")

import numpy as np
from parameters import path, ext_cache
print("loading from", path)

os.environ["LOG_PATH"] = path

from thesis.main.PostProcess import PostProcessor, PostProcessComputation, FenicsScalarFieldComputation
from thesis.main.PostProcessUtil import get_mesh_volume


class my_post_process(FenicsScalarFieldComputation):
    """Defines a custom computation"""
    name = "h1_norm"

    def __call__(self):

        mesh_volume = get_mesh_volume(self.u.function_space().mesh())
        h1: float = fcs.sqrt(fcs.assemble(fcs.dot(self.grad*self.grad_conv, self.grad*self.grad_conv)  * fcs.dX)) / mesh_volume

        return h1


"""number of threads can be passed as first cli argument"""
if len(sys.argv) > 1:
    threads = int(sys.argv[1])
else:
    threads = 6


# user = getpass.getuser()
# path = "/extra/brunner/thesis/kinetic/q_fraction_k_factor/"

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

pp.computations.append(my_post_process)

"""carries out the operations in pp.computations in parallel and stores the result in xml file"""
pp.run_post_process(threads, kde=False)

# if run >= 1:
#     for i in range(40): #N is the number of fractions scanned over
#         shutil.rmtree(path + "scan_" + str(i), ignore_errors=True)
#     shutil.rmtree(path + "solver_tmpil2", ignore_errors=True)
#     shutil.rmtree(path + "tmp", ignore_errors=True)
#     shutil.rmtree(path + "records", ignore_errors=True)
#     shutil.rmtree(path + "images", ignore_errors=True)
#     shutil.rmtree(path + "cache", ignore_errors=True)
#     try:
#         os.remove(path + "log.scan")
#     except FileNotFoundError:
#         pass
#     try:
#         os.remove(path + "timing.pdf")
#     except FileNotFoundError:
#         pass
#     try:
#         os.remove(path + "timing.csv")
#     except FileNotFoundError:
#         pass