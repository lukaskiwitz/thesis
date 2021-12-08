import getpass
import os

import numpy as np

"""defines cytokines. Currently their parameters can be changed using the scan sample interface.
"field_quantity" determines which boundary contitions interact with which fields, 
the name is only used in IO/post processing."""

cytokines = [
    {
        "name": "IL-2",
        "field_quantity": "il2",
        "k_on": 111.6,  # receptor binding constant 1/(nM*h),
        "D": 10,  # Diffusion constant mu^2
        "kd": 0.1,  # cytokine decay in medium 1/h
        "k_endo": 1.1e-3,
        "k_off": 0.83,
    }
]
"""Sets up cells types. 
The first entry is the naive cell type. The "fraction" entry is meaningless.
"""



cell_types_dict = [
    {"name": "naive",
     "fraction": 1,
     "il2": {"R": 100, "q": 0, "bc_type": "patrick_saturation", "ths": 0.01},
     "EC50": 0.1,
     "R_start": 100,
     "misc": {"sigma": 2.5,
              "gamma": 40,
              "states": [0, 0, 0, 0, 0, 0],
              "pos_half_time": 1,
              "neg_half_time": 1,
              "hill_factor": 3,
              "Km_pos": 0.5,
              "Km_neg": 0.75,
              "R_start_neg": 100,
              "R_start_pos": 100,
              "pSTAT5_signal": True,
              "KD": 7.437e-3,
              "nu": 1e-3,
              "pSTAT5": 0,
              "name": "Th"},
     "internal_solver": "kineticSolver"
     },
    {"name": "sec",
     "fraction": 0.25,
     "il2": {"R": 100, "q": 10, "bc_type": "patrick_saturation"},
     "EC50": 0.1,
     "R_start": 100,
     "misc": {
         "sigma": 2.5,
         "gamma": 40,
         "states": [0, 0, 0, 0, 0, 0],
         "pos_half_time": 1,
         "neg_half_time": 1,
         "hill_factor": 3,
         "Km_pos": 0.5,
         "Km_neg": 0.75,
         "R_start_neg": 100,
         "R_start_pos": 100,
         "pSTAT5_signal": True,
         "KD": 7.437e-3,
         "nu": 1e-3,
         "pSTAT5": 0,
         "name": "Th"},
     "internal_solver": "",
     "clustering": {"strength": 1}
     },
    {"name": "treg",
     "fraction": 0.25,
     "il2": {"R": 1000, "q": 0, "bc_type": "patrick_saturation"},
     "clustering": {"bw": 100, "strength": 0},
     "EC50": 0.1,
     "R_start": 1000,
     "misc": {"sigma": 2.5,
              "gamma": 10,
              "states": [0, 0, 0, 0, 0, 0],
              "pos_half_time": 1,
              "neg_half_time": 1,
              "hill_factor": 3,
              "Km_pos": 0.5,
              "Km_neg": 0.75,
              "R_start_neg": 1000,
              "R_start_pos": 1000,
              "pSTAT5_signal": True,
              "KD": 7.437e-3,
              "nu": 1e-3,
              "pSTAT5": 0,
              "name": "treg"},
     "internal_solver": "kineticSolver",
     "clustering": {"strength": 1}
     }
]

"""defines the variable aspects of the geometry. Unit is micro meters"""
geometry = {
    "margin": 20,  # margin around the cell grid
    "distance": 15,  # distance between cell centers
    "rho": 5,  # cell radius
    "x_grid": 110,  # dimensions of the cell grid
    "y_grid": 110,
    "z_grid": 110,  # comment out for single cell layer
    "norm_area": 4 * np.pi * 5 ** 2
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
    "newton_rtol": 1e-3,  # newton method relative tolerance
    "dofs_per_node": int(30000 / 2),
    "max_mpi_nodes": os.cpu_count(),
    "cells_per_worker": 50,
    "max_pool_size": os.cpu_count(),
    "min_char_length": 1,  # mesh
    "max_char_length": 5,  # mesh
    "unit_length_exponent": -6  # for concentration conversion
}

user = getpass.getuser()
model_name = "Treg_competition_patricks_solver_with_kevin_parameters_parallel"
name = "110_3D_fsec_scan"
scan_name = "20211208_reproducing_thurley2015_antigen_response"

hdd = "extra2" if os.path.exists("/extra2") else "extra"

path = "/{hdd}/{u}/{mn}/{n}/{sn}/".format(u=user, n=name, mn=model_name, sn=scan_name, hdd=hdd)
ext_cache = r"../../{mn}_{n}_ext_cache_/".format(mn=model_name, n=name)
IMGPATH = path + "images/"
