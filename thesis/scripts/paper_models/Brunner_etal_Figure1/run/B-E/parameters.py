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
        "KD": 7.437e-3,
        # "KD": 40e-3,
        "hill_factor": 3,
        # "k_endo": 0.00046, #np.log(2)/(25*60)
        "k_endo": 0.0011, #1/(15*60)
    }
]

"""Sets up cells types. 
The first entry is the default cell type. The "fraction" entry is meaningless.
"""

cell_types_dict = [
    {

     "name": "Tnaive",
     "fraction": 0,
     "il2": {"R": 1e2, "q": 0, "bc_type": "linear", "global_q": False},  # [Receptor number per cell, secretion in molecules/s]
     "clustering": {"strength": 0},
     "misc": {"sigma":1,
              "gamma": 1,
              "states":[],
              "hill_factor": 3,
              "Km_pos": 0.5,
              "Km_neg": 0.5,
              "pSTAT5_signal": False,
              # "KD": 7.437e-3,
              # "KD": 40e-3,
              "nu": 1e-3,
              "name": "Tnaive"},
     "internal_solver": ""
     },
    {
     "name": "Tsec",
     "fraction": 0.1,
     "il2": {"R": 0, "q": 10, "bc_type": "linear", "global_q": False}, #"R_saturation"
     "clustering": {"strength": 0},
     "misc": {"sigma": 0,
              "gamma": 1,
              "states":[],
              "hill_factor": 3,
              "Km_pos": 0.5,
              "Km_neg": 0.5,
              "pSTAT5_signal": True,
              "KD": 7.437e-3,
              "k_endo": 0.00046,
              "nu": 1e-3,
              "name": "Tsec"},
     "internal_solver": ""
     },
    {
     "name": "Th",
     "fraction": 0.9,
     "il2": {"R": 1e4, "q": 0, "bc_type": "linear", "global_q": False, "pSTAT5": 0}, # "R_saturation"
     "clustering": {"strength": 0},
     "EC50": {"EC50": 0},
     "R_start": {"R_start": 1e4},
     "misc": {"sigma": 1,
              "gamma": 1,
              "states":[0,0,0,0,0,0],
              "pos_half_time": 1,
              "neg_half_time": 1,
              "hill_factor": 3,
              "Km_pos": 0.5,
              "Km_neg": 0.5,
              "R_start_neg": 1e4,
              "R_start_pos": 1e4,
              "pSTAT5_signal": True,
              "KD": 7.437e-3,
              "k_endo": 0.00046,
              "nu": 1e-3,
              "pSTAT5": 0,
              "name": "Th",
              "EC50_N": 1.5},
     "internal_solver": ""
     },
    {
        "name": "Treg",
        "fraction": 0.0,
        "il2": {"R": 1e5, "q": 0, "bc_type": "linear", "global_q": False, "pSTAT5": 0},
        "clustering": {"strength": 0},
        "EC50": {"EC50": 0},
        "R_start": {"R_start": 1e5},
        "misc": {"sigma": 0.5,
                 "gamma": 100,
                 "states": [0, 0, 0, 0, 0, 0],
                 "pos_half_time": 1,
                 "neg_half_time": 1,
                 "hill_factor": 3,
                 "Km_pos": 0.5,
                 "Km_neg": 0.5,
                 "R_start_neg": 1e5,
                 "R_start_pos": 1e4,
                 "pSTAT5_signal": True,
                 "KD": 7.437e-3,
                 "nu": 1e-3,
                 "pSTAT5": 0,
                 "name": "Treg",
                 "EC50_N": 1.5},
        "internal_solver": "kineticSolver"
    },
]
# mode = "sheet"
mode = "grid"
# mode = ""

if mode == "sheet":

    geometry = {
        "randomize": "bridson",
        "margin": 40,  # margin around the cell grid
        "distance": 12,  # distance between cell centers
        "rho": 5,  # cell radius
        "x_grid": 180,  # dimensions of the cell grid
        "y_grid": 180,
        # "z_grid": 180,# comment out for single cell layer
        "norm_area": 4 * np.pi * 5 ** 2
    }

    numeric = {
        "linear_solver": "gmres",
        "preconditioner": "amg",
        "linear": False,
        "krylov_atol": 1e-35,
        "krylov_rtol": 1e-5,
        "newton_atol": 1e-35,
        "newton_rtol": 1e-2,  # newton method relative tolerance
        "dofs_per_node": int(100000),
        "max_mpi_nodes": os.cpu_count(),
        "cells_per_worker": 50,
        "max_pool_size": 1,
        "min_char_length": 0.2,  # mesh, smaller = finer
        "max_char_length": 3,  # mesh, smaller = finer
        "unit_length_exponent": -6  # for concentration conversion
    }
    name = "300_2D"

elif mode == "grid":
    geometry = {
        "randomize": False,
        "margin": 40,  # margin around the cell grid
        "distance": 20,  # distance between cell centers
        "rho": 5,  # cell radius
        "x_grid": 240,  # dimensions of the cell grid
        "y_grid": 240,
        "z_grid": 240,  # comment out for single cell layer
        "norm_area": 4 * np.pi * 5 ** 2,
    }

    numeric = {
        "linear_solver": "gmres",
        "preconditioner": "amg",
        "linear": False,
        "krylov_atol": 1e-35,
        "krylov_rtol": 1e-5,
        "newton_atol": 1e-35,
        "newton_rtol": 1e-2,  # newton method relative tolerance
        "dofs_per_node": int(30000),
        "max_mpi_nodes": os.cpu_count(),
        "cells_per_worker": 50,
        "max_pool_size": 1,
        "min_char_length": 0.07,  # mesh, smaller = finer
        "max_char_length": 7,  # mesh, smaller = finer
        "unit_length_exponent": -6  # for concentration conversion
    }
    name = "test"

else:
    geometry = {
        "randomize": "bridson",
        "margin": 40,  # margin around the cell grid
        "distance": 12,  # distance between cell centers
        "rho": 5,  # cell radius
        "x_grid": 180,  # dimensions of the cell grid
        "y_grid": 180,
        "z_grid": 180,  # comment out for single cell layer
        "norm_area": 4 * np.pi * 5 ** 2
    }

    numeric = {
        "linear_solver": "gmres",
        "preconditioner": "amg",
        "linear": True,
        "krylov_atol": 1e-35,
        "krylov_rtol": 1e-5,
        "newton_atol": 1e-35,
        "newton_rtol": 1e-2,  # newton method relative tolerance
        "dofs_per_node": int(150000),
        "max_mpi_nodes": os.cpu_count(),
        "cells_per_worker": 50,
        "max_pool_size": 1,
        "min_char_length": 0.07,  # mesh, smaller = finer
        "max_char_length": 7,  # mesh, smaller = finer
        "unit_length_exponent": -6  # for concentration conversion
    }
    name = "300_3D"

user = getpass.getuser()
model_name = "Figure_1C"
# scan_name = "pos_treg_0"
scan_name = "act_over_Tsec_large_linear"

hdd = "extra2" if os.path.exists("/extra2") else "extra"
path = f"/{hdd}/{user}/paper_models/kinetics/{model_name}/{scan_name}/"
ext_cache = r"../../{mn}_{n}_ext_cache/".format(mn=model_name, n=name)
IMGPATH = path + "images/"
