import getpass
import os

import numpy as np

"""defines cytokines. Currently their parameters can be changed using the scan sample interface.
"field_quantity" determines which boundary conditions interact with which fields, 
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
The first entry is the default cell type. There the "fraction" entry is meaningless.
The second entry defines parameters attached to a cytokine field, here il2. bc_type references bcFunctions.py in main.
Misc lets you attach any paremeters to entitys/cell types. List, float and bool should all be possible.
The internal_solver attaches the user defined solver by name.
"""
R_h = 1e4
R_l = 1e2
q = 100

bc_type = "patrick_saturation"
cell_types_dict = [
    {"name": "default",
     "fraction": 0,
     "il2": {"R": R_h, "q": 5, "Kc": 0.01, "amax": 0, "ths": 0.01, "bc_type": bc_type},
     "internal_solver": "RuleBasedSolver"
     },
    {"name": "sec",
     "fraction": 0.005,
     # "fraction": 0.05,
     "il2": {"R": R_l, "q": q, "Kc": 0.01, "amax": 1, "ths": 0.01, "bc_type": bc_type},
     "internal_solver": ""
     },
    {"name": "abs",
     "fraction": 0,
     "il2": {"R": R_h, "q": 0, "Kc": 0.01, "amax": 1, "ths": 0.01, "bc_type": bc_type},
     "internal_solver": ""
     }
]

boundary = [
    {"name": "box",
     "expr": "!near(x[0],{d})",
     "il6": {"q": 0, "R": 0, "bc_type": "R_saturation"},
     "ifng": {"q": 0, "R": 0, "bc_type": "R_saturation"},
     },
    {"name": "left_boundary",
     "expr": "near(x[0],{d})",
     "il6": {"q": 0, "R": 0, "amax": 0, "bc_type": "R_saturation"},
     "ifng": {"q": 0, "R": 0, "bc_type": "R_saturation"},
     }
]

"""defines the variable aspects of the geometry. Unit is micro meters"""
geometry = {
    "margin": 20,  # margin around the cell grid
    "distance": 20,  # distance between cell centers
    "rho": 5,  # cell radius
    "x_grid": 100,  # dimensions of the cell grid (edge length in um)
    "y_grid": 100,
    "z_grid": 100,  # comment out for single cell layer
    "norm_area": 4 * np.pi * 5 ** 2  # area passed to bc function of outside boundary
}

"""
parameters regarding meshing and fenics. unit_length_exponent is necessary for calculation concentrations. 
-6  means the simulation is in micro meters.
"""
numeric = {
    "linear_solver": "gmres",
    "preconditioner": "hypre_amg",
    "linear": False,  # use linear fenics solver. If your chosen bc_type is not linear, flip it
    "krylov_atol": 1e-35,
    "krylov_rtol": 1e-5,  # linear solver relative tolerance
    "newton_atol": 1e-35,
    "newton_rtol": 1e-5,  # newton method relative tolerance
    "dofs_per_node": 7500,  # calc_boundary_values degrees of freedom per mpi node for pde solving
    "max_mpi_nodes": os.cpu_count() / 2,  # max nodes for fenics solver
    "cells_per_worker": 50,
    "max_pool_size": 4,  # max number of worker to extract boundary conditions at runtime
    "min_char_length": 1,  # mesh resolution, smaller = finer
    "max_char_length": 5,  # mesh resolution, smaller = finer
    "unit_length_exponent": -6  # for concentration conversion
}

"""
absolute paths
"""
if os.path.isdir("/extra2"):
    extra = "extra2"
else:
    extra = "extra"

user = getpass.getuser()
model_name = "full_time_series"
name = "example_1"
path = "/{extra}/{u}/{mn}/{n}/".format(u=user, n=name, mn=model_name, extra=extra)
IMGPATH = path + "images/"

"""
relative to path; best to leave this alone
"""
ext_cache = r"../{mn}_ext_cache/".format(mn=model_name)
