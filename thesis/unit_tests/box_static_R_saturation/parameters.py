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
        "kd": 0.1  # cytokine decay in medium 1/h
    }
]

"""Sets up cells types. T
The first entry is the default cell type. There the "fraction" entry is meaningless.
"""
R_h = 1e4
R_l = 1e2
q = 10


misc = {"list":[1,2.58,1e100],"bool":True,"float":np.pi,"int":2,"str":"HELLO WORLD"}

cell_types_dict = [
    {"name": "default",
     "mc":0,
     "fraction": 0,
     "il2": {"R": 0, "q": 0, "Kc": 0.01, "amax": 0,"ths":0.01, "bc_type": "R_saturation"},
     "misc":misc,
     "internal_solver": ""
     },
    {"name": "sec",
     "mc":0,
     "fraction": 0.04,
     "il2": {"R": R_l, "q": q, "Kc": 0.01, "amax": 1, "ths":0.01, "bc_type": "R_saturation"},
      "misc":misc,
     "internal_solver": ""
     },
    {"name": "abs",
     "mc":0,
     "fraction": 0.5,
     "il2": {"R": R_h, "q": 0, "Kc": 0.01, "amax": 1,"ths":0.01, "bc_type": "R_saturation"},
     "misc": misc,
     "internal_solver": ""
     }
]

boundary = [
    {"name": "box",
     "expr":"!near(x[0],{d})",
     "il6":{"q":0, "R":0, "bc_type": "R_saturation"},
     "ifng":{"q":0, "R":0, "bc_type": "R_saturation"},
     },
    {"name": "left_boundary",
     "expr":"near(x[0],{d})",
     "il6":{"q":0, "R":0, "amax":0, "bc_type": "R_saturation"},
     "ifng":{"q":0, "R":0, "bc_type": "R_saturation"},
     }
]

"""defines the variable aspects of the geometry. Unit is micro meters"""
geometry = {
    "margin": 20,  # margin around the cell grid
    "distance": 20,  # distance between cell centers
    "rho": 5,  # cell radius
    "x_grid": 100,  # dimensions of the cell grid (edge length in um)
    "y_grid": 100,
    "z_grid": 100,# comment out for single cell layer
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
    "dofs_per_node": 15000,#target degrees of freedom per mpi node for pde solving
    "max_mpi_nodes":2,# max nodes for fenics solver
    "cells_per_worker": 50,
    "max_pool_size": 2,#max number of worker to extract boundary conditions at runtime
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
model_name = "box_static_R_saturation"
name = "example_1"
path = "/{extra}/{u}/unit_test/{mn}/{n}/".format(u=user, n=name, mn=model_name, extra = extra)
IMGPATH = path + "images/"

"""
relative to path; best to leave this alone
"""
ext_cache = r"../{mn}_ext_cache/".format(mn=model_name)

