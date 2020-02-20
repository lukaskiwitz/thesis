#!/home/kiwitz/anaconda3/envs/fenicsproject/bin/python

from parameter_scan_setup import execute_scan
from parameter_scan_test import get_parameter_space,T, ext_cache,path_prefix,p_c

name = "Diffusion"
# Diffusion
pList = get_parameter_space("D", 1)
path = path_prefix+name+"/"
scan = [{"D": i[0]} for i in pList]
execute_scan(scan, p_c, T, path, ext_cache)
