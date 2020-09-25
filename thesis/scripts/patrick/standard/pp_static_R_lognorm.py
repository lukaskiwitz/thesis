import getpass
import sys

import fenics as fcs
import numpy as np
import os
import shutil

from PostProcess import PostProcessor


class my_post_process:
    """Defines a custom computation"""

    def __init__(self):
        """this name will appear in output dataframe"""
        self.name = "maxGradient"

    def __call__(self, u, grad, c_conv, grad_conv, mesh_volume, V=None, V_vec=None) -> float:

        """__call__ must have this call signature.
        returns the maximum gradient
        """

        gradient = fcs.project(fcs.sqrt(fcs.dot(grad, grad)),V)
        return np.max(np.array(gradient.vector()))*grad_conv


"""number of threads can be passed as first cli argument"""
if len(sys.argv) > 1:
    threads = int(sys.argv[1])
else:
    threads = 6

for j in range(15):
    """
    setting filepath to look for sim results. This is setup so that it works on the itb computers.
    """
    # user = getpass.getuser()
    path = "/extra/brunner/thesis/static/R_lognorm_multi/run" + str(j) + "/"
    os.environ["LOG_PATH"] = path

    """
    setting filepath to look for sim results. This is setup so that it works on the itb computers.
    """
    pp = PostProcessor(path)
    pp.unit_length_exponent = -6

    """appends a custom calculation.
    default computations are defined in PostProcess.py"""

    pp.cell_colors = "Dark2"
    pp.cell_color_key = "type_name"
    pp.legend_title = "Cell Type"

    """carries out the operations in pp.computations in parallel and stores the result in xml file"""
    pp.run_post_process(threads, make_images=False)

    pp.global_dataframe.to_hdf(path + 'global_df.h5', key="data", mode="w")
    pp.cell_dataframe.to_hdf(path + "cell_df.h5", key="df", mode="w")

    if j > 0:
        for i in range(20): #20 is the number of fractions scanned over
            shutil.rmtree(path + "scan_" + str(i), ignore_errors=True)
        shutil.rmtree(path + "solver_tmp", ignore_errors=True)
        try:
            os.remove(path + "log.scan")
        except FileNotFoundError:
            pass
        try:
            os.remove(path + "timing.pdf")
        except FileNotFoundError:
            pass
        try:
            os.remove(path + "timing.csv")
        except FileNotFoundError:
            pass