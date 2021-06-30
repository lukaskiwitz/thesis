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
        "k_off":0.83,
    }
]

"""Sets up cells types. T
The first entry is the default cell type. There the "fraction" entry is meaningless.
"""

R_l = 1e2
R_h = 5000
q = 10

rat = 1/24

f_sec =   lambda  v: (rat*v)/((rat*v)+1)
f_abs =  lambda v: 1/((rat*v)+1)

cell_types_dict = [
    {"name": "naive",
     "fraction": 0,
     "il2": {"q": 0, "R":R_l, "bc_type": "patrick_saturation"},
     "internal_solver": ""
     },
    {"name": "sec",
     "fraction": f_sec(1),
     "il2": {"q": q, "R": R_l, "bc_type": "patrick_saturation"},
     "internal_solver": ""
     },
    {"name": "abs",
     "fraction": f_abs(1),
     "il2": {"q": 0, "R": R_h, "bc_type": "patrick_saturation"},
     "internal_solver": ""
     }
]

"""defines the variable aspects of the geometry. Unit is micro meters"""
geometry = {
    "margin": 20,  # margin around the cell grid
    "distance": 20,  # distance between cell centers
    "rho": 5,  # cell radius
    "x_grid": 100,  # dimensions of the cell grid (edge length in um)
    "y_grid": 100,
    "z_grid":100,# comment out for single cell layer
    "norm_area": 4 * np.pi * 5 ** 2# area passed to bc function of outside boundary
}

"""
parameters regarding meshing and fenics. unit_length_exponent is necessary for calculation concentrations. 
-6  means the simulation is in micro meters.
"""
numeric = {
    "linear_solver": "gmres",
    "preconditioner": "amg",
    "linear": False,# use linear fenics solver
    "krylov_atol": 1e-35,
    "krylov_rtol": 1e-5,# linear solver relative tolerance
    "newton_atol": 1e-35,
    "newton_rtol": 1e-5,# newton method relative tolerance
    "dofs_per_node": 2*15000,#target degrees of freedom per mpi node for pde solving
    "max_mpi_nodes": int(os.cpu_count()),# max nodes for fenics solver
    "cells_per_worker": 50,
    "max_pool_size": 16,#max number of worker to extract boundary conditions at runtime
    "min_char_length": 1,  # mesh resolution
    "max_char_length": 5,  # mesh resolution
    "unit_length_exponent": -6  # for concentration conversion
}

"""
absolute paths
"""

if os.path.isdir("/extra2"):
    extra = "extra2"
else:
    extra ="extra"


user = getpass.getuser()
model_name = "boxed_static_small"
name = "test_1"
path = "/{extra}/{u}/paper_models/{mn}/{n}/".format(u=user, n=name, mn=model_name, extra = extra)
IMGPATH = path + "images/"

"""
relative to path; best to leave this alone
"""
ext_cache = r"../{mn}_ext_cache/".format(mn=model_name)

