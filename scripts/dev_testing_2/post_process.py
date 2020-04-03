import getpass
from PostProcess import PostProcessor
import sys
import fenics as fcs
import numpy as np

class my_post_process:

    def __init__(self):
        self.name = "maxGradient"

    def __call__(self, u, grad, c_conv, grad_conv, mesh_volume, V = None, V_vec = None) -> float:

        gradient = fcs.project(fcs.sqrt(fcs.dot(grad, grad)),V)
        return np.max(np.array(gradient.vector()))*grad_conv

if len(sys.argv) > 1:
    threads = int(sys.argv[1])
else:
    threads = 4


user = getpass.getuser()
path = "/extra/{u}/dev_testing/".format(u=user)

pp = PostProcessor(path)
pp.unit_length_exponent = -6
pp.computations.append(my_post_process())

pp.write_post_process_xml(threads)

pp.make_dataframes(kde=True)
pp.global_dataframe.to_hdf(path + 'global_df.h5', key="data", mode="w")
pp.cell_dataframe.to_hdf(path+"cell_df.h5", key="df", mode="w")
pp.cell_stats.to_hdf(path+"cell_stats_df.h5", key="df", mode="w")