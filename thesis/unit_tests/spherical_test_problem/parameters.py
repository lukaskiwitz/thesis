import getpass
import os
import numpy as np

cytokine = {
    "name": "IL-2",
    "field_quantity": "il2",
    "k_on": 116,  # receptor binding constant 1/(nM*h),
    "D": 10,  # Diffusion constant mu^2
    "kd": 0.1  # cytokine decay in medium 1/h
}

cell = {
    "name": "cell",
    "il2": {"R": 1e2, "q": 10, "bc_type": "linear"},
    "internal_solver": ""
}

boundary = {
    "name": "box",
    "expr": "true",
    "il2": {"R": 1e2, "N": 12, "scenario": "high_density"},
}
geometry = {
    "rho": 5,  # cell radius
    "radius": 20,
    "norm_area": 4 * np.pi * 20 ** 2  # area passed to bc function of outside boundary
}

numeric = {
    "linear_solver": "gmres",
    "preconditioner": "amg",
    "linear": True,  # use linear fenics solver
    "krylov_atol": 1e-35,
    "krylov_rtol": 1e-10,  # linear solver relative tolerance
    "newton_atol": 1e-35,
    "newton_rtol": 1e-10,  # newton method relative tolerance
    "dofs_per_node": 15000,  # target degrees of freedom per mpi node for pde solving
    "max_mpi_nodes": 2,  # max nodes for fenics solver
    "cells_per_worker": 50,
    "max_pool_size": os.cpu_count(),  # max number of worker to extract boundary conditions at runtime
    "min_char_length": 0.1,  # mesh resolution
    "max_char_length": 3,  # mesh resolution
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
model_name = "spherical_test_problem"
name = "test_1"
path = "/{extra}/{u}/unit_test/{mn}/{n}/".format(u=user, n=name, mn=model_name, extra=extra)
IMGPATH = path + "images/"

"""
relative to path; best to leave this alone
"""
ext_cache = r"../{mn}_ext_cache/".format(mn=model_name)
