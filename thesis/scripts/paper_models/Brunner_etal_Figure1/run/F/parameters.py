import getpass
import logging
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
        "k_endo": 0.00046,
        "k_off":0.83,
    }
]

"""Sets up cells types. T
The first entry is the default cell type. There the "fraction" entry is meaningless.
"""

rat = 1/19

f_sec =   lambda  v: (rat*v)/((rat*v)+1)
f_abs =  lambda v: 1/((rat*v)+1)

cell_types_dict = [
    {"name": "naive",
     "fraction": 0,
     "il2": {"q": 0, "R":1e2, "bc_type": "patrick_saturation"},
     "misc": {"sigma": 0,
              "gamma": 1,},
     "internal_solver": ""
     },
    {"name": "Tsec",
     "fraction": 0.05,
     "il2": {"q": 10, "R": 0, "bc_type": "patrick_saturation", "global_q": False},
     "misc": {"sigma": 0,
              "gamma": 1,
              "name": "Tsec"},
     "internal_solver": ""
     },
    {"name": "Th",
     "fraction": 0.95,
     "il2": {"q": 0, "R": 1e4, "bc_type": "patrick_saturation"},
     "misc": {"sigma": 1,
              "gamma": 1,
              "name": "Th",
              "k_endo": 0.00046,
              },
     "internal_solver": ""
     }
]

d = lambda x, v: (x-10) * v + 10
f = lambda n,d : np.ceil(n**(1/3) * d + d)
grid = f(1000,d(20,1))


"""defines the variable aspects of the geometry. Unit is micro meters"""
geometry = {
    "margin": 20,  # margin around the cell grid
    "distance": 20,  # distance between cell centers
    "rho": 5,  # cell radius
    "x_grid": grid,  # dimensions of the cell grid (edge length in um)
    "y_grid": grid,
    "z_grid":grid, # comment out for single cell layer
    "norm_area": 4 * np.pi * 5 ** 2# area passed to bc function of outside boundary
}

"""
parameters regarding meshing and fenics. unit_length_exponent is necessary for calculation concentrations. 
-6  means the simulation is in micro meters.
"""
numeric = {
    "linear_solver": "gmres",
    "preconditioner": "amg",
    "linear": False,  # use linear fenics solver
    "krylov_atol": 1e-35,
    "krylov_rtol": 1e-5,  # linear solver relative tolerance
    "newton_atol": 1e-35,
    "newton_rtol": 1e-5,  # newton method relative tolerance
    "dofs_per_node": 15000,  # calc_boundary_values degrees of freedom per mpi node for pde solving
    "max_mpi_nodes": int(os.cpu_count() / 2),  # max nodes for fenics solver
    "cells_per_worker": 50,
    "max_pool_size": 1,  # max number of worker to extract boundary conditions at runtime
    "min_char_length": 0.06,  # mesh, smaller = finer
    "max_char_length": 6,  # mesh, smaller = finer
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
# get_model_name = "boxed_static_300"

model_name = "boxed_static"
name = "parameter_scan_redo"
path = "/{extra}/{u}/paper_models/{mn}/{n}/".format(u=user, n=name, mn=model_name, extra = extra)
IMGPATH = path + "images/"

"""
relative to path; best to leave this alone
"""
ext_cache = r"../{n}_ext_cache/".format(n=name)

# logging.basicConfig(
#     filename=os.path.join(path, "debug.log"),
#     level=logging.INFO,
#     filemode="w",
#     format='%(levelname)s::%(asctime)s %(message)s',
#     datefmt='%I:%M:%S')
