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
        "k_on": 100,  # receptor binding constant 1/(nM*h),
        "D": 10,  # Diffusion constant mu^2
        "kd": 0.0  # cytokine decay in medium 1/h
    }
]

"""Sets up cells types. 
The first entry is the default cell type. The "fraction" entry is meaningless.
"""

cell_types_dict = [
    {"name": "default",
     "fraction": 1,
     "il2": {"R": 0, "q": 0, "bc_type": "linear"},  # [Receptor number per cell, secretion in molecules/s]
     "internal_solver": "kineticSolver"
     },
    {"name": "changed",
     "fraction": 0.25,
     "il2": {"R": 5000, "bc_type": "linear"},
     "internal_solver": "kineticSolver"
     },
    {"name": "Treg",
     "fraction": 0.25,
     "il2": {"R": 27000, "q": 0, "bc_type": "linear"},
     "internal_solver": "kineticSolver"
     },
]

"""defines the variable aspects of the geometry. Unit is micro meters"""
geometry = {
    "margin": 20,  # margin around the cell grid
    "distance": 20,  # distance between cell centers
    "rho": 5,  # cell radius
    "x_grid": 240,  # dimensions of the cell grid
    "y_grid": 240,
    "z_grid": 240,# comment out for single cell layer
    "norm_area": 4 * np.pi * 5 **2
}

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
    "dofs_per_node": 30000,
    "max_mpi_nodes": int(os.cpu_count()/4),
    "cells_per_worker": 50,
    "max_pool_size": int(os.cpu_count()/4),
    "min_char_length": 0.01,  # mesh
    "max_char_length": 3,  # mesh
    "unit_length_exponent": -6  # for concentration conversion
}

user = getpass.getuser()
model_name = "q_fraction"
name = "scan_name"

path = "/extra2/brunner/thesis/static/q_fraction_new_paras_multi/"
ext_cache = r"../q_fraction_large_ext_cache/"
path_kinetic = "/extra2/brunner/thesis/kinetic/q_fraction_medium_qfrac_scan_new_paras/"
IMGPATH = path + "images/"