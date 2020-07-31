import getpass

import numpy as np

cytokines = [
    {
        "name": "IL-2",
        "field_quantity": "il2",
        "k_on": 111.6,
        "D": 10,
        "kd": 0.1
    },
    {
        "name": "IL-6",
        "field_quantity": "il6",
        "k_on": 111.6,
        "D": 10,
        "kd": 0.1
    },
    {
        "name": "IFNg",
        "field_quantity": "ifng",
        "k_on": 111.6,
        "D": 10,
        "kd": 0.1
    }
]


R_h = 4000
R_l = 100
q_h = 10
q_l = 1

cell_types_dict = [
    # {"name": "Tn",
    #  "fraction": 1,
    #  "il2": [R_h, q_h, 0.06],
    #  "il6": [R_l, 0, 0.04],
    #  "ifng": [R_l, 0, 0.034],
    #  "internal_solver": "RuleBasedSolver"
    #  },
    {"name": "Tn",
     "fraction": 1,
     "il2": [R_h, q_h, 0.035],
     "il6": [R_l, 0, 0.07],
     "ifng": [R_l, 0, 0.035],
     "internal_solver": "RuleBasedSolver"
     },
    {"name": "Tfh",
     "fraction": 0.1,
     "il2": [R_h, q_l],
     "il6": [R_h, q_h],
     "ifng": [0, 0],
     "internal_solver": "RuleBasedSolver"
     },
    {"name": "Th1",
     "fraction": 0.1,
     "il2": [0, 0],
     "il6": [0, 0],
     "ifng": [R_h, q_h],
     "internal_solver": "RuleBasedSolver"
     }
]

geometry = {
    "margin": 20,
    "distance": 20,
    "rho": 5,
    "x_grid": 200,
    "y_grid": 200,
    "norm_area": 4 * np.pi * 5 ** 2
}

numeric = {
    "linear_solver": "gmres",
    "preconditioner": "hypre_amg",
    "krylov_atol": 1e-35,
    "krylov_rtol": 1e-5,
    "min_char_length": 0.1,
    "max_char_length": 3,
    "unit_length_exponent": -6
}

user = getpass.getuser()
model_name = "tfh_th1_model_medium"
name = "boundary_pdiff_scan"

path = "/extra/{u}/{mn}/{n}/".format(u=user, n=name, mn=model_name)

ext_cache = "../{mn}_ext_cache/".format(mn=model_name)
IMGPATH = path + "images/"
