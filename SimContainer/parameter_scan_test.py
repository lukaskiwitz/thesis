# -*- coding: utf-8 -*-

from bcFunctions import cellBC_il2, cellBC_il6, outerBC_il2, outerBC_il6, absorbing_il2
import mpi4py.MPI as MPI
from run_scan import run
from scipy.constants import N_A
import BC as bc
import numpy as np
import fenics as fcs
import matplotlib as mpl
from parameter_scan_setup import execute_scan
mpl.use('Agg')

def get_parameter_space(key,f):
    return np.transpose([
        np.linspace(p_c[key] * f, p_c[key] * 1/f, n)
    ])

T = list(range(1))
n = 10
ext_cache: str = "/extra/kiwitz/extCache/"
path_prefix: str = "/extra/kiwitz/sensitivity/"
p_c = {
        "k_on": 10e9*111.6/60**2,  # 111.6 per hour
        "rho": 0.05,  # mu
        "D": (10**0.5*0.01)**2,  # mu² per s
        "R_h": 400*N_A**-1*10e9,
        "R_l": 10*N_A**-1*10e9,
        "kd": 0,  # 0.1/(60*2),
        "q_h": 10*N_A**-1*10e9,
        "q_l": 1*N_A**-1*10e9,
        "radius": 2,
        "N": 6,
        "fraction": 0.25
    }
# Diffusion
pList = get_parameter_space("D",0.1)
path = path_prefix+"Diffusion"+"/"
scan = [{"D": i[0]} for i in pList]
execute_scan(scan, p_c, T, path, ext_cache)



# fraction
pList = get_parameter_space("fraction",0.1)
path = path_prefix+"fraction"+"/"
scan = [{"fraction": i[0]} for i in pList]
execute_scan(scan, p_c, T, path, ext_cache)

# k_d
pList = get_parameter_space("kd",0.1)
path = path_prefix+"kd"+"/"
scan = [{"kd": i[0]} for i in pList]
execute_scan(scan, p_c, T, path, ext_cache)

# k_on
pList = get_parameter_space("k_on",0.1)
path = path_prefix+"kON"+"/"
scan = [{"k_on": i[0]} for i in pList]
execute_scan(scan, p_c, T, path, ext_cache)

# fraction
pList = get_parameter_space("q_h",0.1)
path = path_prefix+"q_il_s"+"/"
scan = [{"q_il2_s": i[0]} for i in pList]
execute_scan(scan, p_c, T, path, ext_cache)

# Receptor number on il2 secretors
pList = get_parameter_space("R_h",0.1)
path = path_prefix+"R_il2_s"+"/"
scan = [{"R_il2_s": i[0]} for i in pList]
execute_scan(scan, p_c, T, path, ext_cache)

# Receptor number on il2 secretors
pList = get_parameter_space("R_h",0.1)
path = path_prefix+"R_il2_s"+"/"
scan = [{"R_il2_s": i[0]} for i in pList]
execute_scan(scan, p_c, T, path, ext_cache)

# Receptor number on il2 non secretors
pList = get_parameter_space("R_h",0.1)
path = path_prefix+"R_il2_f"+"/"
scan = [{"R_il2_f": i[0]} for i in pList]
execute_scan(scan, p_c, T, path, ext_cache)
