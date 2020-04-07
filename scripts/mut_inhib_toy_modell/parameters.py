from scipy.constants import N_A

from ParameterSet import PhysicalParameter, PhysicalParameterTemplate

"""Sets up parameter templates. This are callable object, which return a full copy of themselves 
with a new value (set in post units). This is so that conversion information has to be specified only one."""
R = PhysicalParameterTemplate(PhysicalParameter("R", 0, to_sim=N_A ** -1 * 1e9))

k_on = PhysicalParameterTemplate(PhysicalParameter("k_on", 111.6, to_sim=1e15 / 60 ** 2, is_global=True))
q = PhysicalParameterTemplate(PhysicalParameter("q", 0, to_sim=N_A ** -1 * 1e9))
D = PhysicalParameterTemplate(PhysicalParameter("D", 10, to_sim=1, is_global=True))
kd = PhysicalParameterTemplate(PhysicalParameter("kd", 0.1, to_sim=1 / (60 * 2), is_global=True))

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
    # "z_grid":100# comment out for single cell layer
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

# R = R(1000).get_in_sim_unit()
# k_on  = k_on(111.6).get_in_sim_unit()
# q = q(1).get_in_sim_unit()
# u = 1e-16
