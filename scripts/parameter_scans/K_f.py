#!/home/kiwitz/anaconda3/envs/fenicsproject/bin/python

import numpy as np

from parameter_scan_test import get_parameter_space, path_prefix

name = "K_f"
# fraction of receptors on secreting cells
pList = get_parameter_space("fraction", 1)
def get_fraction_from_K(K):
        return  (K - 2*np.sqrt(10)*np.sqrt(K-K**2))/(-40+41*K)

f_list = get_fraction_from_K(pList)

# plt.plot(pList,f_list)
# plt.show()
path = path_prefix+name+"/"
scan = [{"fraction": i[0]} for i in f_list]
# execute_scan(scan, p_c, T, path, ext_cache)
