try:
    import fenics as fcs
except RuntimeError:
    import os
    os.environ['PATH'] = '/home/brunner/anaconda3/envs/Lukas2/bin:/home/brunner/.local/bin:/home/brunner/anaconda3/condabin:/usr/local/bin:/usr/bin:/bin:/usr/local/games:/usr/games:/opt/puppetlabs/bin'
    import fenics as fcs

import getpass
import random
import sys
import os

sys.path.append("/home/brunner/thesis/thesis/main/")
sys.path.append("/home/brunner/thesis/thesis/scenarios/")

import numpy as np
from scipy.constants import N_A
# from sympy import symbols, solve
from scipy.integrate import solve_ivp

from static_parameters_R_lognorm import cytokines, cell_types_dict, geometry, numeric, path_kinetic, ext_cache


import thesis.main.StateManager as StateManager
from thesis.main.InternalSolver import InternalSolver
from thesis.main.ParameterSet import MiscParameter, ParameterCollection, ScannablePhysicalParameter, ParameterSet
from thesis.main.ScanContainer import ScanContainer, ScanSample
from thesis.main.SimContainer import SimContainer
from thesis.scripts.patrick.driver.solvers import kineticSolver
from thesis.scripts.patrick.driver.states import updateState
from box_grid_q_fraction import setup
import mpi4py.MPI as MPI

comm = MPI.COMM_WORLD
rank = comm.Get_rank()
size = comm.Get_size()

import logging
os.environ["LOG_PATH"] = path_kinetic
# logging.basicConfig(level=logging.DEBUG,filemode="w", format='%(levelname)s::%(asctime)s %(message)s', datefmt='%I:%M:%S')
logging.getLogger('FFC').setLevel(logging.ERROR)




"""Setup/Simulation"""

"""
setting filepath for simulation results. This is setup so that it works on the itb computers.
If ext_cache = "" the mesh will be cached for each field separately
"""
# ext_cache = "/extra/brunner/para_handling/static/R_lognorm/ext_cache/"
# path = "/extra/brunner/thesis/kinetic/q_fraction_k_factor/"
path = path_kinetic


user = getpass.getuser()


"""Setting up a parameters scan now has a object oriented interface. This is the container class"""
scan_container = ScanContainer()

"""the setup function is defined in an external file. It builds the SimContainer for this simulation.
The variable aspects are passed as list "cytokines, cell_types, geometry and numeric. 
These are imported as modules and can be modified in parameters.py """

sc: SimContainer = setup(cytokines, cell_types_dict, geometry, numeric, path, ext_cache)

"""Imports the parameter Templates"""
from box_grid_q_fraction import get_parameter_templates
from static_parameters_R_lognorm import numeric

templates = get_parameter_templates(numeric["unit_length_exponent"])
t_gamma = templates["gamma"]
t_hill_factor = templates["hill_factor"]
t_Tsec_fraction = templates["Tsec_fraction"]
t_R = templates["R"]
t_R_start = templates["R_start"]
t_D = templates["D"]
t_Th_fraction = templates["Th_fraction"]
t_c0 = templates["c0"]
t_pSTAT5_signal = templates["pSTAT5_signal"]
t_pSTAT5 = templates["pSTAT5"]
t_EC50 = templates["EC50"]
t_sigma = templates["sigma"]

"""Sets up Scannable parameters from parameters templates"""

Tsec_fraction = ScannablePhysicalParameter(t_Tsec_fraction(0.0), lambda x, v: v)
Th_fraction = ScannablePhysicalParameter(t_Th_fraction(0.0), lambda x, v: v)
gamma = ScannablePhysicalParameter(t_gamma(0.0), lambda x, v: v)
hill_factor = ScannablePhysicalParameter(t_hill_factor(2), lambda x, v: v)
c0 = ScannablePhysicalParameter(t_c0(0.0), lambda x, v: v)
pSTAT5_signal = ScannablePhysicalParameter(t_pSTAT5_signal(False), lambda x, v: v)
sigma = ScannablePhysicalParameter(t_sigma(1e4), lambda x, v: v)

Tsec_distribution_array = np.array([])
Th_distribution_array = np.array([])

# a = [1/x for x in reversed(np.arange(10,110,10))]
# b = [x for x in np.arange(10,110,10)]
# c = np.concatenate([a,b])
# a = [10]

c = [1]
for v in c:#np.linspace(0, 20000.0, 1): #np.logspace(-1,1,3):
    for sig in np.arange(0,1.05e4,500):
        for hill_fac in [3]:
            """Scans over parameters that are associated with a field"""
            sim_parameters = [
                ParameterCollection("IL-2", [gamma(v)], field_quantity="il2"),
                ParameterCollection("IL-2", [Th_fraction(0.75)], field_quantity="il2"),
                ParameterCollection("IL-2", [c0(8.5e-12)], field_quantity="il2"), #8.637363
                ParameterCollection("IL-2", [Tsec_fraction(0.25)], field_quantity="il2"),
                ParameterCollection("IL-2", [hill_factor(hill_fac)], field_quantity="il2"),
                ParameterCollection("IL-2", [pSTAT5_signal(False)], field_quantity="il2"),
                ParameterCollection("IL-2", [sigma(sig)], field_quantity="il2"),
                # ParameterCollection("IL-2", [kd(v)], field_quantity="il2")
            ]

            """Scans over parameters that are associated with an entity_type"""
            entity_types = [
                # default.get_updated([ParameterCollection("IL-2", [bc_type(v)])]),
                # abs.get_updated([ParameterCollection("IL-2", [bc_type(v)])]),
                # sec.get_updated([ParameterCollection("IL-2", [bc_type(v)])]),
            ]
            """Scans over parameters that are associated with the outer domain
            This is a dictionary. If boundary pieces where defined in the setup function, they can be referenced by name. 
            Here the pieces "left_boundary" and "box" are defined."""
            outer_domain_dict = {
                # "left_boundary": [ParameterCollection("IL-2",[R(v)])],
                # "box": [ParameterCollection("IL-2",[R(v)])]
            }
            """Creates container object for one sample of a the parameter scan. 
            The Lists/Dicts can be empty for default parameters."""
            sample = ScanSample(sim_parameters, entity_types, outer_domain_dict)
            scan_container.add_sample(sample)

"""signs up the internal solver with the sim container. 
It can be referenced in a cell_type definition by its name field
"""
sc.add_internal_solver(kineticSolver)

"""State Manager updates the parameters of simulation objects in accordance with scan samples defined above and 
manages the orderly IO of simulation results and metadata for post processing."""

stMan = StateManager.StateManager(path)
stMan.sim_container = sc
stMan.scan_container = scan_container

"""sets up time range"""

stMan.T = [0,1,2]

"""defines a function which is called by StateManager before a parameter scan. 
Here it is used to assign cell types
"""

uS = updateState(0,geometry,templates, Tsec_distribution_array, Th_distribution_array, offset=0)

def pre_scan(state_manager, scan_index):
    uS.step_R_lognorm(state_manager.sim_container)



def pre_step(sc, time_index, t, T):
    # calculate the average surface concentration for the first solution
    if time_index == 2:
        list = []
        for e in sc.entity_list:
            # print(e.p.get_as_dictionary())
            list.append(e.p.get_physical_parameter("surf_c", "IL-2").get_in_post_unit()*1e-9)
        if len(list) != 0:
            for e in sc.entity_list:
                e.p.update(ParameterSet("dummy", [ParameterCollection("IL-2", [c0(np.mean(list))], field_quantity="il2")] ))
    return None

stMan.sim_container.pre_step = pre_step
stMan.pre_scan = pre_scan

"""Runs the ParameterScan"""

if len(sys.argv) > 1:
    if not sys.argv[1] == "mesh":
        stMan.run()
else:

    stMan.run()
