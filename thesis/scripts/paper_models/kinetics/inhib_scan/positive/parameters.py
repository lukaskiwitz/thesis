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
        "kd": 0.1,  # cytokine decay in medium 1/s
        "hill_factor": 3,
        "k_endo": 1.1e-3,
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
     "misc": {"sigma":1,
              "states":[],
              "hill_factor": 3,
              "Km_pos": 7e-1,
              "Km_neg": 1e-7,
              "pSTAT5_signal": True,
              "KD": 7.437e-3,
              "nu": 1e-3,
              "name": "Tnaive"},
     "internal_solver": "kineticSolver"
     },
    {
     "name": "Tsec",
     "fraction": 0.35,
     "il2": {"R": 1e2, "q": 10, "bc_type": "linear", "global_q": False}, #"R_saturation"
     "misc": {"sigma": 1,
              "states":[],
              "hill_factor": 3,
              "Km_pos": 0.65,
              "Km_neg": 0.75,
              "pSTAT5_signal": True,
              "KD": 7.437e-3,
              "nu": 1e-3,
              "name": "Tsec"},
     "internal_solver": "kineticSolver"
     },
    {
     "name": "Th",
     "fraction": 0.65,
     "il2": {"R": 1e3, "q": 0, "bc_type": "linear", "global_q": False},
     "misc": {"sigma": 2.5,
              "states":[0,0,0,0,0,0],
              "pos_half_time": 1,
              "neg_half_time": 1,
              "hill_factor": 3,
              "Km_pos": 0.65,
              "Km_neg": 0.75,
              "R_start_neg": 1e3,
              "R_start_pos": 1e3,
              "pSTAT5_signal": True,
              "KD": 7.437e-3,
              "nu": 1e-3,
              "pSTAT5": 0,
              "name": "Th"},
     "internal_solver": "kineticSolver"
     },
    {
        "name": "Treg",
        "fraction": 0.1,
        "il2": {"R": 1e3, "q": 0, "bc_type": "linear", "global_q": False},
        "misc": {"sigma": 0.5,
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
                 "name": "Treg"},
        "internal_solver": "kineticSolver"
    },
]

"""defines the variable aspects of the geometry. Unit is micro meters"""
geometry = {
    "margin": 40,  # margin around the cell grid
    "distance": 20,  # distance between cell centers
    "rho": 5,  # cell radius
    "x_grid": 280,  # dimensions of the cell grid
    "y_grid": 280,
    "z_grid": 280,# comment out for single cell layer
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
    "linear": True,
    "krylov_atol": 1e-35,
    "krylov_rtol": 1e-5,
    "newton_atol": 1e-35,
    "newton_rtol": 1e-5,
    "dofs_per_node": 15000,
    "max_mpi_nodes": int(os.cpu_count()/4),
    "cells_per_worker": 50,
    "max_pool_size": int(os.cpu_count()/4),
    "min_char_length": 0.06,  # mesh, smaller = finer
    "max_char_length": 6,  # mesh, smaller = finer
    "unit_length_exponent": -6  # for concentration conversion
}

hdd = "extra2" if os.path.exists("/extra2") else "extra"

user = getpass.getuser()

model_name = "inhib_scan"
name = "positive"
path = "/{extra}/{u}/paper_models/kinetics/{mn}/{n}/".format(u=user, n=name, mn=model_name, extra = hdd)
IMGPATH = path + "images/"

ext_cache = r"../../{mn}_ext_cache/".format(mn=model_name)