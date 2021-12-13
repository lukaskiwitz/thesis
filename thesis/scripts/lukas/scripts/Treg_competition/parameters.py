import getpass

import numpy as np
import os
"""defines cytokines. Currently their parameters can be changed using the scan sample interface.
"field_quantity" determines which boundary contitions interact with which fields, 
the name is only used in IO/post processing."""

cytokines = [
    {
        "name": "IL-2",
        "field_quantity": "il2",
        "k_on": 111.6,  # receptor binding constant 1/(nM*h),
        "D": 10,  # Diffusion constant mu^2
        "kd": 0.1,  # cytokine decay in medium 1/h,
        "Kc": 0.01
    }
]

"""Sets up cells types. 
The first entry is the default cell type. The "fraction" entry is meaningless.
"""

cell_types_dict = [
    {"name": "naive",
     "fraction": 1,
     "il2": {"amax": 0, "q": 0, "bc_type":"amax_saturation"},  # [Receptor number per cell, secretion in molecules/s]
     "internal_solver": "ResponseTimeSolver"
     },
    {"name": "sec",
     "fraction": 0.02,
     "il2": {"amax": 1, "q": 100,"bc_type":"amax_saturation"},
     "internal_solver": "ResponseTimeSolver",
     "bw":25,
     "cluster_strength":0.8
     },
    {"name": "Treg",
     "fraction": 0.2,
     "il2": {"amax": 100, "q": 0,"bc_type":"amax_saturation"},
     "internal_solver": "ResponseTimeSolver"
     }
]

"""defines the variable aspects of the geometry. Unit is micro meters"""
geometry = {
    "margin": 20,  # margin around the cell grid
    "distance": 20,  # distance between cell centers
    "rho": 5,  # cell radius
    "x_grid": 500,  # dimensions of the cell grid
    "y_grid": 500,
    # "z_grid":0,# comment out for single cell layer
    "norm_area": 4 * np.pi * 5 **2
}

"""
parameters regarding meshing and fenics. unit_length_exponent is necessary for calculation concentrations. 
-6  means the simulation is in micro meters.
"""
numeric = {
    "linear_solver": "gmres",
    "preconditioner": "amg",
    "linear": False,
    "krylov_atol": 1e-35,
    "krylov_rtol": 1e-5,
    "newton_atol": 1e-35,
    "newton_rtol": 1e-3,# newton method relative tolerance
    "dofs_per_node": 30000,
    "max_mpi_nodes": os.cpu_count(),
    "cells_per_worker": 50,
    "max_pool_size": os.cpu_count(),
    "min_char_length": 0.1,  # mesh
    "max_char_length": 3,  # mesh
    "unit_length_exponent": -6  # for concentration conversion
}

user = getpass.getuser()
model_name = "Treg_competition"
name = "500_2D"
scan_name = "response_test"

hdd = "extra2" if os.path.exists("/extra2") else "extra"

path = "/{hdd}/{u}/{mn}/{n}/{sn}/".format(u=user, n=name, mn=model_name, sn = scan_name,hdd = hdd)
ext_cache = r"../../{mn}_{n}_ext_cache_/".format(mn=model_name,n = name)
IMGPATH = path + "images/"
