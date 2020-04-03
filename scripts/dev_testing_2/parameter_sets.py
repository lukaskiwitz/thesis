from copy import deepcopy
import numpy as np
from scipy.constants import N_A
from ParameterSet import PhysicalParameter, PhysicalParameterTemplate

R = PhysicalParameterTemplate(PhysicalParameter("R", 0, to_sim=N_A ** -1 * 1e9))

k_on = PhysicalParameterTemplate(PhysicalParameter("k_on", 111.6, to_sim=1e9 / 60 ** 2, is_global=True))

q = PhysicalParameterTemplate(PhysicalParameter("q", 0, to_sim=N_A ** -1 * 1e9))


D = PhysicalParameterTemplate(
    PhysicalParameter("D", 10, to_sim=1,
                      is_global=True))

kd = PhysicalParameterTemplate(PhysicalParameter("kd", 0.1, to_sim=1 / (60 * 2), is_global=True))

cytokines = [
    {
        "name": "IL-2",
        "field_quantity": "il2"
    },
    {
        "name": "IL-6",
        "field_quantity": "il6"
    },
    {
        "name": "IFNg",
        "field_quantity": "ifng"
    }
]


R_h = 4000
R_l = 100
q_h = 10
q_l = 1

cell_types_dict = [
    {"name": "Tn",
     "fraction": 1,
     "il2": [R_h, q_h],
     "il6": [R_l, 0],
     "ifng": [R_l, 0],
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
    "margin":20,
    "distance":20,
    "rho":5,
    "x_grid":100,
    "y_grid":100
}

numeric = {
    "linear_solver":"gmres",
    "preconditioner":"hypre_amg",
    "krylov_atol": 1e-35,
    "krylov_rtol": 1e-5,
    "min_char_length": 0.1,
    "max_char_length": 3,
    "unit_length_exponent":-6
}