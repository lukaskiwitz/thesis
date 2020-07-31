import os
import sys

sys.path.append("/home/brunner/thesis/thesis/main/")
sys.path.append("/home/brunner/thesis/thesis/scenarios/")

import numpy as np
from parameters_q_fraction import path_kinetic2
path=path_kinetic2
# path = "/extra/brunner/thesis/kinetic/q_fraction_test/"

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

"""carries out the operations in pp.computations in parallel and stores the result in xml file"""
pp.run_post_process(threads, make_images=False)

pp.global_dataframe.to_hdf(path + 'global_df.h5', key="data", mode="w")
pp.cell_dataframe.to_hdf(path+"cell_df.h5", key="df", mode="w")
