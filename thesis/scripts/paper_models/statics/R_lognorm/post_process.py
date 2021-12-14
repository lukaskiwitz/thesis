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
        # grad_u = fcs.project(grad(u), W).vector().array().reshape(-1,3)
        # grad_u = grad
        h1: float = fcs.sqrt(fcs.assemble(fcs.dot(grad*grad_conv, grad*grad_conv)  * fcs.dX)) / mesh_volume
        # conc: float = fcs.sqrt(fcs.assemble(fcs.dot(u, u) * fcs.dX)) * c_conv / mesh_volume
        # print(gradient, conc)

        # g = np.reshape(grad.vector().vec().array, (u.vector().vec().size, 3))
        # g = np.transpose(g)
        #
        # my_grad = np.sqrt(np.power(g[0], 2) + np.power(g[1], 2) + np.power(g[2], 2))
        # my_grad = np.mean(my_grad) * grad_conv
        # h1 = 0
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

pp.computations.append(my_post_process())

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