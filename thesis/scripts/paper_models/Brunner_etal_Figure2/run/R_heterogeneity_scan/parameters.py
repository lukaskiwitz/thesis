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
        "k_endo": 0.00046, #np.log(2)/(25*60)
        # "k_endo": 0.0011, #1/(15*60)
    }
]

"""Sets up cells types. 
The first entry is the default cell type. The "fraction" entry is meaningless.
"""

cell_types_dict = [
    {

     "name": "Tnaive",
     "fraction": 0,
     "il2": {"R": 1e2, "q": 0, "bc_type": "patrick_saturation", "global_q": False},  # [Receptor number per cell, secretion in molecules/s]
     "misc": {"sigma":1,
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
     "il2": {"R": 0, "q": 10, "bc_type": "patrick_saturation", "global_q": True},
     "misc": {"sigma":0,
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
     "il2": {"R": 1e4, "q": 0, "bc_type": "patrick_saturation", "global_q": True}, # "R_saturation"
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
        "il2": {"R": 1e5, "q": 0, "bc_type": "patrick_saturation", "global_q": False, "pSTAT5": 0},
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

"""defines the variable aspects of the geometry. Unit is micro meters"""
geometry = {
    "margin": 40,  # margin around the cell grid
    "distance": 20,  # distance between cell centers
    "rho": 5,  # cell radius
    "x_grid": 240,  # dimensions of the cell grid
    "y_grid": 240,
    "z_grid": 240,# comment out for single cell layer
    "norm_area": 4 * np.pi * 5 **2
}

boundary = [
    {"name": "box",
     "expr":"true",
     "il2":{"q":0, "R":0, "bc_type": "linear"},
     },

]

"""
parameters regarding meshing and fenics. unit_length_exponent is necessary for calculation concentrations. 
-6  means the simulation is in micro meters.
"""
numeric = {
    "linear_solver": "gmres",
    "preconditioner": "hypre_amg",
    "linear": False,
    "krylov_atol": 1e-35,
    "krylov_rtol": 1e-5,
    "newton_atol": 1e-35,
    "newton_rtol": 1e-5,
    "dofs_per_node": 40000,
    "max_mpi_nodes": int(os.cpu_count() - 1),
    "cells_per_worker": 50,
    "max_pool_size": 1,
    "min_char_length": 0.05,  # mesh, smaller = finer
    "max_char_length": 5,  # mesh, smaller = finer
    "unit_length_exponent": -6  # for concentration conversion
}

hdd = "extra2" if os.path.exists("/extra2") else "extra"

user = getpass.getuser()

model_name = "Figure_2C"
# scan_name = "Tsec_scan_constant_R"
scan_name = "R_Tsec_0.05"
path = "/{extra}/{u}/paper_models/statics/saturated/{mn}/{n}/".format(u=user, mn=model_name, n=scan_name, extra = hdd)
IMGPATH = path + "images/"

ext_cache = fr"../../{scan_name}_3D_ext_cache/"