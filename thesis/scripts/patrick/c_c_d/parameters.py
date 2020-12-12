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
        "kd": 0.1  # cytokine decay in medium 1/h
    }
]

"""Sets up cells types. 
The first entry is the default cell type. The "fraction" entry is meaningless.
"""

cell_types_dict = [
    {"name": "Tnaive",
     "fraction": 0,
     "il2": {"R": 1e2, "q": 0, "bc_type": "R_saturation"},  # [Receptor number per cell, secretion in molecules/s]
     "internal_solver": "kineticSolver"
     },
    {"name": "Tsec",
     "fraction": 0.25,
     "il2": {"R": 1e4, "bc_type": "R_saturation"}, #"R_saturation"
     "internal_solver": "kineticSolver"
     },
    {"name": "Treg",
     "fraction": 0.75,
     "il2": {"R": 1e4, "q": 0, "bc_type": "R_saturation"},
     "internal_solver": "kineticSolver"
     },
    {"name": "blank",
     "fraction": 0.75,
     "il2": {"R": 0, "q": 0, "bc_type": "R_saturation"},
     "internal_solver": "kineticSolver"
     },
]

"""defines the variable aspects of the geometry. Unit is micro meters"""
geometry = {
    "margin": 20,  # margin around the cell grid
    "distance": 30,  # distance between cell centers
    "rho": 5,  # cell radius
    "x_grid": 300,  # dimensions of the cell grid
    "y_grid": 300,
    #"z_grid": 200,# comment out for single cell layer
    "norm_area": 4 * np.pi * 5 **2
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
    "max_mpi_nodes": int(os.cpu_count()/4),
    "cells_per_worker": 50,
    "max_pool_size": int(os.cpu_count()/4),
    "min_char_length": 0.05,  # mesh
    "max_char_length": 5,  # mesh
    "unit_length_exponent": -6  # for concentration conversion
}

user = getpass.getuser()
model_name = "q_fraction"
name = "scan_name"

path = "/extra2/brunner/thesis/static/q_fraction_new_paras_multi/"
ext_cache = r"../c-c-d_" + str(geometry["distance"]) + "/"
hdd = "/extra2" if os.path.exists("/extra2") else "/extra"
path_kinetic = hdd + "/brunner/thesis/kinetic/c_c_d/kinetic_c_c_d_" + str(geometry["distance"]) + "/"
IMGPATH = path + "images/"