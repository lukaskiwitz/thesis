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
        "kd": 0.1  # cytokine decay in medium 1/h
    }
]

"""Sets up cells types. 
The first entry is the default cell type. There the "fraction" entry is meaningless.
The second entry defines parameters attached to a cytokine field, here il2. bc_type references bcFunctions.py in main.
Misc lets you attach any paremeters to entitys/cell types. List, float and bool should all be possible.
The internal_solver attaches the user defined solver by name.
"""
R_h = 1e2
R_l = 1e2
q = 10

myList = [1,2,3]

cell_types_dict = [
    {"name": "default",
     "fraction": 0,
     "il2": {"R": 0, "q": 0, "Kc": 0.01, "amax": 0,"ths":0.01, "bc_type": "linear"},
     "misc":{"custom_para":myList},
     "internal_solver": "RuleBasedSolver"
     },
    {"name": "sec",
     "fraction": 0.04,
     "il2": {"R": R_l, "q": q, "Kc": 0.01, "amax": 1, "ths":0.01, "bc_type": "linear"},
     "misc":{"custom_para":myList},
     "internal_solver": "RuleBasedSolver"
     },
    {"name": "abs",
     "fraction": 0.96,
     "il2": {"R": R_h, "q": 0, "Kc": 0.01, "amax": 1,"ths":0.01, "bc_type": "linear"},
     "misc":{"custom_para":myList},
     "internal_solver": "RuleBasedSolver"
     }
]

"""
Sets up boundary conditions and their respective bc_type.
Here the fields "il6" and "ifng" are explicitly defined. The expression or condition is passed to fenics and should be 
written in C. near(x[0],{d}) means if value x[0], which is the x-axis, is close to d. d is calculated from the margin, 
compare line 186 in box_grid.py. Alternatively one could set a value or tolerance here directly.
If this condition returns true the defined boundary condition is applied. If none is given for the cytokine field the
standard box is applied.
"""
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
    "x_grid": 150,  # dimensions of the cell grid (edge length in um)
    "y_grid": 150,
    "z_grid": 150, # comment out for single cell layer
    "norm_area": 4 * np.pi * 5 ** 2 # area passed to bc function of outside boundary
}

"""
parameters regarding meshing and fenics. unit_length_exponent is necessary for calculation concentrations. 
-6  means the simulation is in micro meters.
"""
numeric = {
    "linear_solver": "gmres",
    "preconditioner": "hypre_amg",
    "linear": True,# use linear fenics solver. If your chosen bc_type is not linear, flip it
    "krylov_atol": 1e-35,
    "krylov_rtol": 1e-5,# linear solver relative tolerance
    "newton_atol": 1e-35,
    "newton_rtol": 1e-5,# newton method relative tolerance
    "dofs_per_node": 15000,#calc_boundary_values degrees of freedom per mpi node for pde solving
    "max_mpi_nodes": os.cpu_count()/2,# max nodes for fenics solver
    "cells_per_worker": 50,
    "max_pool_size": os.cpu_count()/2,#max number of worker to extract boundary conditions at runtime
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
    extra ="extra"


user = getpass.getuser()
model_name = "example_hendrik"
name = "example_1"
path = "/{extra}/{u}/{mn}/{n}/".format(u=user, n=name, mn=model_name, extra = extra)
IMGPATH = path + "images/"

"""
relative to path; best to leave this alone
"""
ext_cache = r"../{mn}_ext_cache/".format(mn=model_name)