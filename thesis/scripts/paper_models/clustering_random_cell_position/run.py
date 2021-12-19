try:
    import fenics as fcs
except RuntimeError:
    import os
    os.environ['PATH'] = '/home/brunner/anaconda3/envs/Lukas2/bin:/home/brunner/.local/bin:/home/brunner/anaconda3/condabin:/usr/local/bin:/usr/bin:/bin:/usr/local/games:/usr/games:/opt/puppetlabs/bin'
    import fenics as fcs

import logging
from copy import deepcopy


module_logger = logging.getLogger(__name__)
import numpy as np
# from sympy import symbols, solve

from parameters import cytokines, cell_types_dict, geometry, numeric, path, ext_cache, boundary
import thesis.main.StateManager as StateManager
from thesis.main.ParameterSet import ScannableParameter, PhysicalParameter
from thesis.main.ScanContainer import ScanContainer, ScanDefintion, ScanType
from thesis.scripts.paper_models.utilities.converging_clustering import conv_clustering
from thesis.scenarios.box_grid import setup
from thesis.main.my_debug import message

"""Setup/Simulation"""

"""
setting filepath for simulation results. This is setup so that it works on the itb computers.
If ext_cache = "" the mesh will be cached for each field separately
"""
# ext_cache = "/extra/brunner/para_handling/static/R_lognorm/ext_cache/"
# path = "/extra/brunner/thesis/kinetic/q_fraction_k_factor/"

# path = path

"""Setting up a parameters scan now has a object oriented interface. This is the container class"""
scan_container = ScanContainer()

"""the setup function is defined in an external file. It builds the SimContainer for this simulation.
The variable aspects are passed as list "cytokines, cell_types, geometry and numeric. 
These are imported as modules and can be modified in parameters.py """

scenario = setup(cytokines, deepcopy(cell_types_dict), boundary, geometry, numeric)

"""Imports the parameter Templates"""
parameter_pool = scenario.parameter_pool

"""Retrieves and entity type from sim container for scanning"""
Tnaive = scenario.get_entity_type_by_name("Tnaive")
Tsec = scenario.get_entity_type_by_name("Tsec")
Th = scenario.get_entity_type_by_name("Th")

"""
Sets up a parameter scan. ScannableParameter takes a function with two arguments, 
here conveniently a lambda function, to set the given PhysicalParameter according to input. 
A simpler definition would be 'lambda x, v: v'. 
scan_space sets up the values which are in the end inserted into the lambda function. ScanDefinition then combines
both while defining on which entity this scan should be applied to with ScanType.[...].
To actually run this scan setup we attach it to the scan_container with add_single_parameter_scan.
"""

t_gamma = parameter_pool.get_template("gamma")
t_global_q = parameter_pool.get_template("global_q")

gamma = ScannableParameter(t_gamma(0.0), lambda x, v: v)
global_q = ScannableParameter(t_global_q(False), lambda x, v: v)

Tsec_distribution_array = np.array([])
Th_distribution_array = np.array([])
Treg_distribution_array = np.array([])

s = 2
scan_space = np.linspace(1, 0, s)
# scan_space = [0]

cs = ScannableParameter(PhysicalParameter("strength", 0.5), lambda x, v: v)

sec_cs_def = ScanDefintion(cs, "clustering", scan_space, ScanType.ENTITY, entity_type=Tsec)
abs_cs_def = ScanDefintion(cs, "clustering", [1], ScanType.ENTITY, entity_type=Th)

scan_gamma = 40
gamma_def = ScanDefintion(gamma, "IL-2", [scan_gamma], ScanType.GLOBAL, field_quantity="il2")
Tsec_fraction = ScannableParameter(PhysicalParameter("Tsec_fraction", 1, is_global=True), lambda x, v: v)

# scan_space = np.around(np.linspace(8,40,10), 3)
# scan_space = np.around(np.linspace(2000,10000,10), 3)
# scan_space = [0.25 for x in range(3)]

# Tsec_fraction = ScannableParameter(PhysicalParameter("Tsec_fraction", 1, is_global=True), lambda x, v: v)
# q_def = ScanDefintion(q, "IL-2", scan_space, ScanType.ENTITY, field_quantity="il2", entity_type=Tsec)
# R_def = ScanDefintion(R, "IL-2", scan_space, ScanType.ENTITY, field_quantity="il2", entity_type=Tsec)
# Tsec_def = ScanDefintion(Tsec_fraction, "IL-2", scan_space, ScanType.GLOBAL, field_quantity="il2")

# scan_container.add_single_parameter_scan([q_def], scan_name="q_scan")
# scan_container.add_single_parameter_scan([R_def], scan_name="R_scan")
# scan_container.add_single_parameter_scan([Tsec_def, gamma_def], scan_name="Tsec_fraction")
scan_container.add_single_parameter_scan([sec_cs_def, abs_cs_def, gamma_def], scan_name = "clustering") #, remesh_scan_sample=True)

"""signs up the internal solver with the sim container. 
It can be referenced in a cell_type definition by its name field
"""

from thesis.cellBehaviourUtilities.cell_solver import kineticSolver
scenario.internal_solvers = [kineticSolver]

"""State Manager updates the parameters of simulation objects in accordance with scan samples defined above and 
manages the orderly IO of simulation results and metadata for post processing."""

stMan = StateManager.StateManager(path)
stMan.scenario = scenario
scenario.marker_lookup = {"Tnaive": 1, "Tsec": 2, "Th": 3, "Treg": 4}
scenario.markers = ["type_name", "IL-2_surf_c", "IL-2_R"]
stMan.scan_container = scan_container
stMan.compress_log_file = True

"""sets up time range"""

dt = 3600  # 1h
length = 2
max_T = dt * 2

myRange = np.arange(0, length)
def exp_func(x, a, b, c):
    return np.round(a * np.exp(b * x) + c, 4)
a = 2 * dt
c = -a
b = np.log((max_T - c) / a) / (length - 1)
# T = exp_func(myRange, a, b, c)
T = np.linspace(0, max_T, length)

stMan.T = T
kineticSolver.T = stMan.T


"""defines a function which is called by StateManager before a parameter scan. 
Here it is used to assign cell types
"""

# cc = conv_clustering(replicat_index = 0, geometry = geometry, parameter_pool = parameter_pool)

from thesis.scripts.paper_models.utilities.apc_clustering import update_state as uS

# def pre_scan(state_manager, scan_index):
#     cc.step(state_manager.sim_container, scale = 0.35 * scan_index)

def pre_replicat(sc, time_index, replicat_index, t, T):
    # uS.reset_draw_arrays()
    # cc.step(sc)
    uS(sc, replicat_index, parameter_pool, apcs = np.array([[80,80,80]]))

def post_replicat(sc, time_index, replicat_index, t, T):
    list = []
    for e in sc.entity_list:
        celltype = e.p.get_misc_parameter("name", "misc").get_in_post_unit()
        if celltype == "Th":
            list.append(e.p.get_physical_parameter("surf_c", "IL-2").get_in_post_unit() * 1e-9)
    message("#####", module_logger)
    message("mean concentration in M: " + str(np.mean(list)), module_logger)
    message("#####", module_logger)

    R_list = []
    for e in sc.entity_list:
        celltype = e.p.get_misc_parameter("name", "misc").get_in_post_unit()
        if celltype == "Th":
            R_list.append(e.p.get_physical_parameter("R", "IL-2").get_in_post_unit())
    message("#####", module_logger)
    message("mean R: " + str(np.mean(R_list)), module_logger)
    message("#####", module_logger)

def pre_step(sc, time_index, replicat_index, t, T):
    # calculate the average surface concentration for the first solution
    # if time_index == 2:
    #     list = []
    #     for e in sc.entity_list:
    #         # print(e.p.get_as_dictionary())
    #         list.append(e.p.get_physical_parameter("surf_c", "IL-2").get_in_post_unit()*1e-9)
    #     if len(list) != 0:
    #         for e in sc.entity_list:
    #             e.p.update(ParameterSet("dummy", [ParameterCollection("IL-2", [c0(np.median(list))], field_quantity="il2")] ))
    # return None
    pass

# stMan.pre_scan = pre_scan
# stMan.pre_step = pre_step
stMan.pre_replicat = pre_replicat
stMan.post_replicat = post_replicat

"""Runs the ParameterScan"""
stMan.run(model_names=["pde_model"], ext_cache=ext_cache, number_of_replicats=1)
