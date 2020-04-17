import sys

import fenics as fcs
import numpy as np
from parameters import path

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

        gradient = fcs.project(fcs.sqrt(fcs.dot(grad, grad)), V)
        return np.max(np.array(gradient.vector())) * grad_conv


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

# pp.computations.append(my_post_process())

"""carries out the operations in pp.computations in parallel and stores the result in xml file"""
pp.write_post_process_xml(threads)

"""collects the post processing result from xml file in dataframes"""
pp.make_dataframes(kde=True)
pp.global_dataframe.to_hdf(path + 'global_df.h5', key="data", mode="w")
pp.cell_dataframe.to_hdf(path + "cell_df.h5", key="df", mode="w")
