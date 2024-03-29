import getpass

import numpy as np

"""defines cytokines. Currently their parameters can be changed using the scan sample interface.
"field_quantity" determines which boundary contitions interact with which fields, 
the name is only used in IO/post processing."""

cytokines = [
    {
        "name": "IL-2",
        "field_quantity": "il2"
    }
]

"""Sets up cells types. 
The first entry is the default cell type. The "fraction" entry is meaningless.
"""

cell_types_dict = [
    {"name": "default",
     "fraction": 1,
     "il2": [500, 1],  # [Receptor number per cell, secretion in molecules/s]
     "internal_solver": "RuleBasedSolver"
     },
    {"name": "sec",
     "fraction": 0.1,
     "il2": [100, 10],
     "internal_solver": "RuleBasedSolver"
     },
    {"name": "abs",
     "fraction": 0.1,
     "il2": [4000, 1],
     "internal_solver": "RuleBasedSolver"
     }
]

"""defines the variable aspects of the geometry. Unit is micro meters"""
geometry = {
    "margin": 20,  # margin around the cell grid
    "distance": 20,  # distance between cell centers
    "rho": 5,  # cell radius
    "x_grid": 500,  # dimensions of the cell grid
    "y_grid": 500,
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
    "krylov_atol": 1e-35,
    "krylov_rtol": 1e-5,
    "min_char_length": 0.1,  # mesh
    "max_char_length": 3,  # mesh
    "unit_length_exponent": -6  # for concentration conversion
}

user = getpass.getuser()

model_name = "mut_inhib_toy_model"
name = "scan_name"

path = "/extra/{u}/{mn}/{n}/".format(u=user, n=name, mn=model_name)
ext_cache = "../{mn}_ext_cache/".format(mn=model_name)
IMGPATH = path + "images/""../{n}_ext_cache/".format(n=name)
