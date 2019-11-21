# -*- coding: utf-8 -*-

from scipy.constants import N_A
import numpy as np
import matplotlib as mpl
from parameter_scan_setup import execute_scan
mpl.use('Agg')


# def get_parameter_space(key: str, e: int) -> np.ndarray:
#     return np.transpose([
#         np.logspace(-e, e, n, base=10)*p_c[key]
#     ])


def get_parameter_space(key: str, e: int) -> np.ndarray:
    return np.transpose([
        np.logspace(-e, e, n, base=10)*p_c[key]
    ])


T = list(range(5))
n = 50
ext_cache: str = "/extra/kiwitz/extCache/"
path_prefix: str = "/extra/kiwitz/sensitivity_fullReceptor/"
p_c = {
        "k_on": 10e9*111.6/60**2,  # 111.6 per hour
        "rho": 0.05,  # mu
        "D": (10**0.5*0.01)**2,  # mu² per s
        "R_h": 4000*N_A**-1*10e9,
        "R_l": 100*N_A**-1*10e9,
        "kd":  0.1/(60*2),
        "q_h": 10*N_A**-1*10e9,
        "q_l": 1*N_A**-1*10e9,
        "radius": 2,
        "N": 6,
        "fraction": 0.1
    }
