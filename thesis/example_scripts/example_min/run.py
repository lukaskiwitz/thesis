import os
import sys

sys.path.append("/home/lukas/thesis/main/")
sys.path.append("/home/lukas/thesis/scenarios/")

import random

import numpy as np
from parameters import cytokines, cell_types_dict, geometry, numeric, path, ext_cache, boundary

import logging
os.environ["LOG_PATH"] = path
LOG_PATH = os.environ.get("LOG_PATH") if os.environ.get("LOG_PATH") else "./"
os.makedirs(LOG_PATH,exist_ok=True)
logging.basicConfig(filename=LOG_PATH+"debug.log",level=logging.INFO,filemode="w", format='%(levelname)s::%(asctime)s %(message)s', datefmt='%I:%M:%S')

os.environ["LOG_PATH"] = path

import thesis.main.StateManager as StateManager
from thesis.main.InternalSolver import InternalSolver
from thesis.main.ParameterSet import MiscParameter, ParameterCollection, ScannablePhysicalParameter, PhysicalParameter
from thesis.main.ScanContainer import ScanContainer, ScanSample, ScanDefintion, ScanType
from thesis.main.SimContainer import SimContainer
from thesis.scenarios.box_grid import setup
from thesis.main.assign_fractions import assign_fractions


from scipy.stats import poisson

class RuleBasedSolver(InternalSolver):
    name = "RuleBasedSolver"

    def step(self, t1, t2, dt, p, entity=None):

        il2_threshold = p.get_physical_parameter("ths", "IL-2").get_in_post_unit()
        il2 = p.get_physical_parameter("surf_c", "IL-2").get_in_post_unit()

        rand = np.random.uniform(0,1)


        def hill(c,ths,invert = False):

            l_max = 1
            n = 4
            if invert:
                return l_max * (1-(c ** n / (c ** n + ths**n))) * dt
            else:
                return l_max * ((c ** n / (c ** n + ths**n))) * dt

        il2_cdf = 1-poisson.cdf(k = 1, mu = hill(il2, il2_threshold,invert = True))

        if entity.type_name == "default":
            if il2_cdf > rand:
                entity.change_type = "sec"
        return p


def updateState(sc, t):

    """sets cell types according to the values given in fractions.
    The pseudo random seed depends on t, so that cell placement is repeatable. """

    assign_fractions(sc, t)

"""Setup/Simulation"""

"""
setting filepath for simulation results. This is setup so that it works on the itb computers.
If ext_cache = "" the mesh will be cached for each field separately
"""

"""Setting up a parameters scan now has a object oriented interface. This is the container class"""
scan_container = ScanContainer()

"""the setup function is defined in an external file. It builds the SimContainer for this simulation.
The variable aspects are passed as list "cytokines, cell_types, geometry and numeric. 
These are imported as modules and can be modified in parameters.py """

sc: SimContainer = setup(cytokines, cell_types_dict,boundary, geometry, numeric, path, ext_cache)

"""Imports the parameter Templates"""
from thesis.scenarios.box_grid import get_parameter_templates
templates = get_parameter_templates(numeric["unit_length_exponent"])

t_D = templates["D"]
t_R = templates["R"]
t_q = templates["q"]
t_kd = templates["kd"]
t_amax = templates["amax"]
t_Kc = templates["Kc"]

"""Sets up Scannable parameters from parameters templates"""
R = ScannablePhysicalParameter(t_R(1e4), lambda x, v: x * v)
q = ScannablePhysicalParameter(t_q(1), lambda x, v: x * v)
D = ScannablePhysicalParameter(t_D(10), lambda x, v: x * v)
amax = ScannablePhysicalParameter(t_amax(1), lambda x, v: x * v)
bc_type = ScannablePhysicalParameter(MiscParameter("bc_type", "linear"), lambda x, v: v)
is_linear = ScannablePhysicalParameter(MiscParameter("linear", True,is_global = True), lambda x, v: v)
kd = ScannablePhysicalParameter(t_kd(0.1), lambda x, v: x * v)
f_sec = ScannablePhysicalParameter(PhysicalParameter("sec", 0.1, is_global=True), lambda x, v: 1/(v+1))
f_abs = ScannablePhysicalParameter(PhysicalParameter("abs", 0.1, is_global=True), lambda x, v: v/(v+1))
Kc = ScannablePhysicalParameter(t_Kc(0.01), lambda x, v: x * v)

"""Retrieves and entity type from sim container for scanning"""
default = sc.get_entity_type_by_name("default")
abs = sc.get_entity_type_by_name("abs")
sec = sc.get_entity_type_by_name("sec")

s = 2

scan_space = np.concatenate([np.logspace(-1,0,int(s/2)), np.logspace(0,1,int(s/2)+1)[1:]])

f_sec_def = ScanDefintion(f_sec, "fractions", scan_space, ScanType.GLOBAL)
f_abs_def = ScanDefintion(f_abs, "fractions", scan_space, ScanType.GLOBAL)
# scan_container.add_single_parameter_scan([f_abs_def, f_sec_def], scan_name = "f")

D_def = ScanDefintion(D,"IL-2",  scan_space, ScanType.GLOBAL, field_quantity = "il2")
kd_def = ScanDefintion(kd,"IL-2",  scan_space, ScanType.GLOBAL, field_quantity = "il2")
# sec_q_def = ScanDefintion(q,"IL-2",  scan_space, ScanType.ENTITY, field_quantity = "il2", entity_type = sec)
# abs_R_def = ScanDefintion(R,"IL-2",  scan_space, ScanType.ENTITY, field_quantity = "il2", entity_type = abs)

# scan_container.add_2d_parameter_scan(
#     ([D_def],"D"),([kd_def],"kd")
# )
scan_container.add_single_parameter_scan([D_def], scan_name = "D")


"""signs up the internal solver with the sim container. 
It can be referenced in a cell_type definition by its name field
"""
sc.add_internal_solver(RuleBasedSolver)

"""State Manager updates the parameters of simulation objects in accordance with scan samples defined above and 
manages the orderly IO of simulation results and metadata for post processing."""

stMan = StateManager.StateManager(path)
stMan.sim_container = sc
sc.marker_lookup = {"default":1, "sec":2, "abs":3}#labels to apper in marker function; 0 denotes background
stMan.scan_container = scan_container
stMan.compress_log_file = True

# stMan.dt = 1
# stMan.N = 5
"""sets up time range"""

stMan.T = [0,1]
# stMan.T = np.linspace(0,10,10)


"""defines a function which is called by StateManager before a parameter scan. 
Here it is used to assign cell types
"""

def pre_scan(state_manager, scan_index):
    updateState(state_manager.sim_container, 0)


stMan.pre_scan = pre_scan

"""Runs the ParameterScan"""
stMan.run()







