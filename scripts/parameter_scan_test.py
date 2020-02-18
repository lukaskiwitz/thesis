# -*- coding: utf-8 -*-

from scipy.constants import N_A
import numpy as np
import matplotlib as mpl
from parameter_scan_setup import execute_scan
# mpl.use('Agg')


def get_parameter_space(key: str, e: int) -> np.ndarray:
    return np.transpose([
        np.logspace(-e, e, n, base=10)*p_c[key]
    ])


T = list(range(5))
n = 25
ext_cache: str = "/extra/kiwitz/extCache/"
path_prefix: str = "/extra/kiwitz/sensitivity_new_test/"
p_c = {
        "k_on": 1e9*111.6/60**2,  # 111.6 per hour
        "rho": 0.05,  # mu
        "D": (10**0.5*0.01)**2,  # mu² per s
        "R_h": 400*N_A**-1*1e9,
        "R_l": 10*N_A**-1*1e9,
        "kd":  0.1/(60*2),
        "q_h": 2*N_A**-1*1e9,
        "q_l": 0.2*N_A**-1*1e9,
        "radius": 2,
        "N": 6,
        "fraction": 0.1
    }
