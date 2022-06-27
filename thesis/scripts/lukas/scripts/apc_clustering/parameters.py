import getpass
import os

import numpy as np

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

cell_types_dict = [
    {"name": "naive",
     "fraction": 1,
     "il2": {"R": 1500, "q": 0, "bc_type": "patrick_saturation", "ths": 0.01},
     "EC50": 0.1,
     "R_start": 1500,
     "misc": {
         "pSTAT_k": 860,
         "pSTAT_N": 1,
         "sigma": 1,
         "gamma": 30,
         "states": [0, 0, 0, 0, 0, 0],
         "pos_half_time": 1,
         "neg_half_time": 1,
         "hill_factor": 3,
         "Km_pos": 0.5,
         "Km_neg": 0.75,
         "R_start_neg": 1500,
         "R_start_pos": 1500,
         "pSTAT5_signal": True,
         "KD": 7.437e-3,
         "nu": 1e-3,
         "pSTAT5": 0,
         "name": "Th"},
     "internal_solver": "kineticSolver"
     },
    {"name": "treg",
     "fraction": 0,
     "il2": {"R": 1e4, "q": 0, "bc_type": "patrick_saturation"},
     "clustering": {"strength": 0},
     "EC50": 0.1,
     "R_start": 3000,
     "misc": {"sigma": 1,
              "gamma": 30,  # Fuhrmann2016 supplement
              "pSTAT_k": 860,
              "pSTAT_N": 1,
              "states": [0, 0, 0, 0, 0, 0],
              "pos_half_time": 1,
              "neg_half_time": 1,
              "hill_factor": 3,
              "Km_pos": 0.5,
              "Km_neg": 0.75,
              "R_start_neg": 3000,
              "R_start_pos": 3000,
              "pSTAT5_signal": True,
              "KD": 7.437e-3,
              "nu": 1e-3,
              "pSTAT5": 0,
              "name": "treg"},
     "internal_solver": "kineticSolver",
     },
    {"name": "sec",
     "fraction": 0.25,
     "il2": {"R": 100, "q": 10, "bc_type": "patrick_saturation"},
     "misc": {
         "pSTAT_k": 860,
         "pSTAT_N": 1},
     "internal_solver": "",
     "clustering": {"strength": 1}
     },
]
mode = "box"

if mode == "sheet":

    geometry = {
        "randomize": "bridson",
        "margin": 12,  # margin around the cell grid
        "distance": 12,  # distance between cell centers
        "rho": 5,  # cell radius
        "x_grid": 500,  # dimensions of the cell grid
        "y_grid": 500,
        # "z_grid": 200,# comment out for single cell layer
        "norm_area": 4 * np.pi * 5 ** 2
    }

    numeric = {
        "linear_solver": "gmres",
        "preconditioner": "amg",
        "linear": False,
        "krylov_atol": 1e-35,
        "krylov_rtol": 1e-5,
        "newton_atol": 1e-35,
        "newton_rtol": 1e-3,  # newton method relative tolerance
        "dofs_per_node": int(60000),
        "max_mpi_nodes": os.cpu_count(),
        "cells_per_worker": 50,
        "max_pool_size": 1,
        "min_char_length": 1,  # mesh
        "max_char_length": 5,  # mesh
        "unit_length_exponent": -6  # for concentration conversion
    }
    name = "500_2D"

else:
    geometry = {
        "randomize": "bridson",
        "margin": 12,  # margin around the cell grid
        "distance": 12,  # distance between cell centers
        "rho": 5,  # cell radius
        "x_grid": 300,  # dimensions of the cell grid
        "y_grid": 300,
        "z_grid": 300,  # comment out for single cell layer
        "norm_area": 4 * np.pi * 5 ** 2
    }

    numeric = {
        "linear_solver": "gmres",
        "preconditioner": "amg",
        "linear": False,
        "krylov_atol": 1e-35,
        "krylov_rtol": 1e-5,
        "newton_atol": 1e-35,
        "newton_rtol": 1e-3,  # newton method relative tolerance
        "dofs_per_node": int(30000),
        "max_mpi_nodes": 8,#os.cpu_count(),
        "cells_per_worker": 50,
        "max_pool_size": 1,
        "min_char_length": 1,  # mesh
        "max_char_length": 5,  # mesh
        "unit_length_exponent": -6  # for concentration conversion
    }
    name = "300_3D"

user = getpass.getuser()
model_name = "Treg_competition_parallel_patrick_parameters"
scan_name = "multiple_apc_dry_run"

hdd = "extra2" if os.path.exists("/extra2") else "extra"
path = "/{hdd}/{u}/{mn}/{n}/{sn}/".format(u=user, n=name, mn=model_name, sn=scan_name, hdd=hdd)
ext_cache = r"../../{mn}_{n}_ext_cache_/".format(mn=model_name, n=name)
IMGPATH = path + "images/"
