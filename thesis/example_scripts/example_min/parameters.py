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
        "k_on": 116,  # receptor binding constant 1/(nM*h),
        "D": 10,  # Diffusion constant mu^2
        "kd": 0.1  # cytokine decay in medium 1/h
    }
]

"""Sets up cells types. T
The first entry is the default cell type. There the "fraction" entry is meaningless.
"""

cell_types_dict = [
    {"name": "default",
     "fraction": 0,
     "il2": {"R": 4000, "q": 0, "Kc": 0.01, "amax": 0, "bc_type": "linear"},
     "internal_solver": "RuleBasedSolver"
     },
    {"name": "sec",
     "fraction": 0.5,
     "il2": {"R": 4000, "q": 20, "Kc": 0.01, "amax": 1, "bc_type": "linear"},
     "internal_solver": "RuleBasedSolver"
     },
    {"name": "abs",
     "fraction": 0.5,
     "il2": {"R": 20000, "q": 0, "Kc": 0.01, "amax": 100, "bc_type": "linear"},
     "internal_solver": "RuleBasedSolver"
     }
]

"""defines the variable aspects of the geometry. Unit is micro meters"""
geometry = {
    "margin": 20,  # margin around the cell grid
    "distance": 20,  # distance between cell centers
    "rho": 5,  # cell radius
    "x_grid": 200,  # dimensions of the cell grid
    "y_grid": 200,
    # "z_grid":100,# comment out for single cell layer
    "norm_area": 4 * np.pi * 5 ** 2
}

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
    "dofs_per_node": 30000,
    "max_mpi_nodes": os.cpu_count(),
    "cells_per_worker": 50,
    "max_pool_size": os.cpu_count(),
    "min_char_length": 0.05,  # mesh
    "max_char_length": 2,  # mesh
    "unit_length_exponent": -6  # for concentration conversion
}

user = getpass.getuser()
model_name = "example_min"
name = "test"

path = "/extra/{u}/{mn}/{n}/".format(u=user, n=name, mn=model_name)
ext_cache = r"../{mn}_ext_cache/".format(mn=model_name)
IMGPATH = path + "images/"
