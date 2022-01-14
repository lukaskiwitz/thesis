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
import logging
from copy import deepcopy

sys.path.append("/home/brunner/thesis/thesis/main/")
sys.path.append("/home/brunner/thesis/thesis/scenarios/")

import numpy as np
from copy import deepcopy

from parameters import cytokines, cell_types_dict, geometry, numeric, path, ext_cache, boundary

os.environ["LOG_PATH"] = path
LOG_PATH = os.environ.get("LOG_PATH") if os.environ.get("LOG_PATH") else "./"
os.makedirs(LOG_PATH, exist_ok=True)
logging.basicConfig(filename=LOG_PATH + "debug.log", level=logging.INFO, filemode="w",
                    format='%(levelname)s::%(asctime)s %(message)s', datefmt='%I:%M:%S')

os.environ["LOG_PATH"] = path

import thesis.main.StateManager as StateManager
from thesis.main.ParameterSet import ScannableParameter, PhysicalParameter, PhysicalParameterTemplate, \
    MiscParameter
from thesis.main.ScanContainer import ScanContainer, ScanDefintion, ScanType
from thesis.scripts.paper_models.utilities.states import updateState
from thesis.scenarios.box_grid import setup
from thesis.main.my_debug import message
import mpi4py.MPI as MPI


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
Treg = scenario.get_entity_type_by_name("Treg")

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

scan_space = [0, 10, 100]

steps = ScannableParameter(MiscParameter("steps", 100), lambda x, v: v)
steps_def = ScanDefintion(steps, "geometry", scan_space, ScanType.GLOBAL)

scan_container.add_single_parameter_scan([steps_def], scan_name="steps", remesh_scan_sample=True)

"""State Manager updates the parameters of simulation objects in accordance with scan samples defined above and 
manages the orderly IO of simulation results and metadata for post processing."""

stMan = StateManager.StateManager(path)
stMan.scenario = scenario
scenario.marker_lookup = {"Tnaive": 1, "Tsec": 2, "Th": 3, "Treg": 4}
scenario.markers = ["type_name", "IL-2_surf_c", "IL-2_R"]
stMan.scan_container = scan_container
stMan.compress_log_file = True

"""sets up time range"""

stMan.T = [0,1]

"""defines a function which is called by StateManager before a parameter scan. 
Here it is used to assign cell types
"""

uS = updateState(0, 0, geometry,parameter_pool, Tsec_distribution_array, Th_distribution_array, Treg_distribution_array, offset=0)

def pre_scan(state_manager, scan_index):
    uS.step(state_manager.sim_container)

def pre_replicat(sc, time_index, replicat_index, t, T):
    uS.reset_draw_arrays()
    uS.step(sc)

def post_replicat(sc, time_index, replicat_index, t, T):
    list = []
    for e in sc.entity_list:
        celltype = e.p.get_misc_parameter("name", "misc").get_in_post_unit()
        if celltype == "Th":
            list.append(e.p.get_physical_parameter("surf_c", "IL-2").get_in_post_unit() * 1e-9)
    message("#####")
    message("mean concentration in M: " + str(np.mean(list)))
    message("#####")

    R_list = []
    for e in sc.entity_list:
        celltype = e.p.get_misc_parameter("name", "misc").get_in_post_unit()
        if celltype == "Th":
            R_list.append(e.p.get_physical_parameter("R", "IL-2").get_in_post_unit())
    message("#####")
    message("mean R: " + str(np.mean(R_list)))
    message("#####")



stMan.pre_scan = pre_scan

stMan.pre_replicat = pre_replicat
# stMan.post_replicat = post_replicat

"""Runs the ParameterScan"""
stMan.run(model_names=["pde_model"], ext_cache=ext_cache, number_of_replicats=1)
