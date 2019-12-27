#!/home/kiwitz/anaconda3/envs/fenicsproject/bin/python

from parameter_scan_setup import execute_scan
from parameter_scan_test import get_parameter_space,T,n, ext_cache,path_prefix,p_c

import numpy as np
import matplotlib.pyplot as plt
name = "K_f"
# fraction of receptors on secreting cells
pList = get_parameter_space("fraction", 1)
def get_fraction_from_p(p):
        a = p_c["fraction"]
        n = 1500
        # st = a*p_c["R_h"]*n + (1- a)*p_c["R_l"]*n
        st = a * 400 * n + (1 - a) * 10 * n
        print("total {st}".format(st=st))
        return  {
                "R_il2_s":(-p*st)/((-1+a)),
                "R_ils_f":-((-1+p)*st)/(a)
        }

f_list = [get_fraction_from_p(p[0]) for p in pList]
for i in f_list:
        print(i)
        # vs = list(i.values())
        # print(vs[0]+vs[1])
# print(f_list)

# plt.plot(pList,f_list)
# plt.show()
# path = path_prefix+name+"/"
# scan = [{"fraction": i[0]} for i in f_list]
# execute_scan(scan, p_c, T, path, ext_cache)
