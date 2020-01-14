#from bcFunctions import  outerBC_il2, outerBC_il6
import mpi4py.MPI as MPI
from run_scan import run
import BC as bc
import numpy as np
import fenics as fcs
import matplotlib as mpl
from typing import List, Dict
import sys
from PostProcess import  PostProcessor
# mpl.use('Agg')

comm = MPI.COMM_WORLD
rank = comm.Get_rank()
size = comm.Get_size()


def execute_scan(scan: List[Dict],p_c: Dict, T:List[float], path:str, ext_cache: str) -> None:
    if len(sys.argv) == 1:
        print("started with no parameters")
        return
    def bD(n):
        return n/10

    x = np.linspace(-bD(15), bD(15), int(15))
    y = np.linspace(-bD(10), bD(10), int(10))
    z = np.linspace(-bD(10), bD(10), int(10))

    d_x = x[-1]+0.2
    d_y = y[-1]+0.2
    d_z = z[-1]+0.2

    p_d = {"x": x,
           "y": y,
           "z": z,
           "d_x": d_x,
           "d_y": d_y,
           "d_z": d_z
           }
    p_sim = {  # default values
        "R_il2": 0,
        "q_il2": 0,
        "R_il6": 0,
        "q_il6": 0
    }

    p_il2 = {
        "R_il2_s": p_c["R_l"],  # secretors
        "R_il2_f": p_c["R_h"],  # fraction
        "R_il2_n": p_c["R_l"],  # normal
        "R_il2_b": p_c["R_h"],  # boundary
        "q_il2_s": p_c["q_h"],
        "q_il2_f": p_c["q_l"],
        "q_il2_n": p_c["q_l"],
        "q_il2_b": p_c["q_l"]
    }
    p_il6 = {
        "R_il6_s": p_c["R_l"],
        "R_il6_f": p_c["R_l"],
        "R_il6_n": 0,
        "R_il6_b": 0,
        "q_il6_s": p_c["q_h"]*0.5,
        "q_il6_f": p_c["q_l"],
        "q_il6_n": p_c["q_l"],
        "q_il6_b": p_c["q_h"],
    }

    domainBC = [
        bc.outerIntegral(lambda u, p: fcs.Constant(0),
                         "!near(x[0],-{d_x})".format(d_x=d_x), field_quantity="il2"),
        bc.outerIntegral(lambda u, p: fcs.Constant(0),
                         "near(x[0],-{d_x})".format(d_x=d_x), field_quantity="il2"),
        # bc.outerIntegral(outerBC_il2,
        #                  "near(x[0],-{d_x})".format(d_x=d_x), field_quantity="il2"),
        bc.outerIntegral(lambda u, p: fcs.Constant(0),
                         "near(x[0],-{d_x})".format(d_x=d_x), field_quantity="il6"),
        bc.outerIntegral(lambda u, p: fcs.Constant(0),
                         "!near(x[0],-{d_x})".format(d_x=d_x), field_quantity="il6")
    ]

    run(
        {**p_d, **p_c, **p_il2, **p_il6, **p_sim},
        T,
        domainBC,
        path,
        extCache=ext_cache,
        scan=scan
    )

