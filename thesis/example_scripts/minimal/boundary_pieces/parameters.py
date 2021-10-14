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
        "Kc": 0.01,
        "k_endo": 1.1e-3,
        "k_off": 0.83,
    },
    {
        "name": "IL-6",
        "field_quantity": "il6",
        "k_on": 111.6,  # receptor binding constant 1/(nM*h),
        "D": 10,  # Diffusion constant mu^2
        "kd": 0.1,  # cytokine decay in medium 1/h
        "Kc": 0.01,
        "k_endo": 1.1e-3,
        "k_off": 0.83,
    }
]

"""Sets up cells types.
The first entry is the default cell type. There the "fraction" entry is meaningless.
"""
R_h = 1e4
R_l = 1e2
q = 1

cell_types_dict = [
    {"name": "default",
     "mc": 0,
     "fraction": 0,
     "il2": {"R": 1e2, "q": 0, "bc_type": "patrick_saturation"},
     "il6": {"R": 1e2, "q": 0, "bc_type": "patrick_saturation"},
     "internal_solver": ""
     },
    {"name": "sec",
     "mc": 0,
     "fraction": 0.1,
     "il2": {"R": 1e2, "q": 10, "bc_type": "patrick_saturation"},
     "il6": {"R": 1e2, "q": 10, "bc_type": "patrick_saturation"},
     "internal_solver": ""
     },
    {"name": "abs",
     "mc": 0,
     "fraction": 0.5,
     "il2": {"R": 1e4, "q": 0, "bc_type": "patrick_saturation"},
     "il6": {"R": 1e4, "q": 0, "bc_type": "patrick_saturation"},
     "internal_solver": ""
     }
]

boundary = [
    {"name": "left_boundary",
     "expr": "near(x[0],{p0x})",
     "il6": {"q": 5, "R": 1e2, "bc_type": "patrick_saturation"},
     "il2": {"q": 0, "R": 0, "bc_type": "patrick_saturation"},
     },
    {"name": "top_boundary",
     "expr": "near(x[1],{p1y})",
     "il6": {"q": 0, "R": 0, "bc_type": "patrick_saturation"},
     "il2": {"q": 5, "R": 1e2, "bc_type": "patrick_saturation"},
     },
    {"name": "right_boundary",
     "expr": "near(x[0],{p1x})",
     "il6": {"q": 0, "R": 1e2, "bc_type": "patrick_saturation"},
     "il2": {"q": 0, "R": 0, "bc_type": "patrick_saturation"},
     },
    {"name": "bottom_boundary",
     "expr": "near(x[1],{p0y})",
     "il6": {"q": 0, "R": 0, "bc_type": "patrick_saturation"},
     "il2": {"q": 5, "R": 1e2, "bc_type": "patrick_saturation"},
     },
    {"name": "ventral_boundary",
     "expr": "near(x[2],{p0z})",
     "il6": {"q": 0, "R": 0, "bc_type": "patrick_saturation"},
     "il2": {"q": 0, "R": 0, "bc_type": "patrick_saturation"},
     },
    {"name": "dorsal_boundary",
     "expr": "near(x[2],{p1z})",
     "il6": {"q": 0, "R": 0, "bc_type": "patrick_saturation"},
     "il2": {"q": 0, "R": 0, "bc_type": "patrick_saturation"},
     },
]

"""defines the variable aspects of the geometry. Unit is micro meters"""
geometry = {
    "margin": 20,  # margin around the cell grid
    "distance": 20,  # distance between cell centers
    "rho": 5,  # cell radius
    "x_grid": 150,  # dimensions of the cell grid (edge length in um)
    "y_grid": 150,
    "z_grid": 150,  # comment out for single cell layer
    "norm_area": 4 * np.pi * 5 ** 2  # area passed to bc function of outside boundary
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
    "dofs_per_node": 15000,  # target degrees of freedom per mpi node for pde solving
    "max_mpi_nodes": os.cpu_count() / 2,  # max nodes for fenics solver
    "cells_per_worker": 50,
    "max_pool_size": os.cpu_count() / 2,  # max number of worker to extract boundary conditions at runtime
    "min_char_length": 1,  # mesh resolution
    "max_char_length": 5,  # mesh resolution
    "unit_length_exponent": -6  # for concentration conversion
}
