#!/home/kiwitz/anaconda3/envs/fenicsproject/bin/python

from parameter_scan_setup import execute_scan
from parameter_scan_test import get_parameter_space,T,n, ext_cache,path_prefix,p_c

name = "Ril2F"
# Receptor number on il2 non secretors
pList = get_parameter_space("R_h", 1)
path = path_prefix+name+"/"
scan = [{"R_il2_f": i[0]} for i in pList]
execute_scan(scan, p_c, T, path, ext_cache)
