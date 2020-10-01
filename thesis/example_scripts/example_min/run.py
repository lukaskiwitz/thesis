import os
import sys

sys.path.append("/home/lukas/thesis/main/")
sys.path.append("/home/lukas/thesis/scenarios/")

# if "LD_LIBRARY_PATH" in os.environ:
#     os.environ['LD_LIBRARY_PATH'] = "/home/lukas/anaconda3/envs/fenics/lib:"+os.environ['LD_LIBRARY_PATH']
# else:
#     os.environ['LD_LIBRARY_PATH'] = "/home/lukas/anaconda3/envs/fenics/lib"

import random

import numpy as np
from parameters import cytokines, cell_types_dict, geometry, numeric, path, ext_cache, boundary

import logging
os.environ["LOG_PATH"] = path
LOG_PATH = os.environ.get("LOG_PATH") if os.environ.get("LOG_PATH") else "./"
os.makedirs(LOG_PATH,exist_ok=True)
logging.basicConfig(filename=LOG_PATH+"debug.log",level=logging.DEBUG,filemode="w", format='%(levelname)s::%(asctime)s %(message)s', datefmt='%I:%M:%S')

os.environ["LOG_PATH"] = path

import thesis.main.StateManager as StateManager
from thesis.main.InternalSolver import InternalSolver
from thesis.main.ParameterSet import MiscParameter, ParameterCollection, ScannablePhysicalParameter, PhysicalParameter
from thesis.main.ScanContainer import ScanContainer, ScanSample
from thesis.main.SimContainer import SimContainer
from thesis.scenarios.box_grid import setup
import mpi4py.MPI as MPI

comm = MPI.COMM_WORLD
rank = comm.Get_rank()
size = comm.Get_size()


class RuleBasedSolver(InternalSolver):
    """Controls the internal logic of cells. An instance is created for each cell.
    Persistent values should be stored in the cell ParameterSet, values attached to the InternalSolver
    object may get lost."""

    name = "RuleBasedSolver"

    def __init__(self):
        pass

    def step(self, t1 ,t2, dt, p, entity=None):


        """
        Is called for each cell after the PDEs are solved.
        Implements cell behaviour. Here cells of type "default" have a 50% transition
        probability if il2 surface concentration is greater than the threshold.
        """



        if entity.type_name == "default":

            il2 = p.get_physical_parameter("surf_c", "IL-2").get_in_post_unit()
            il2_threshold = p.get_physical_parameter("ths", "IL-2").get_in_post_unit()

            if np.random.uniform(0, 1) > 0.5:
                if il2 > il2_threshold:
                    entity.change_type = "sec"
                else:
                    entity.change_type = "abs"
        return p


def updateState(sc, t):

    """sets cell types according to the values given in fractions.
    The pseudo random seed depends on t, so that cell placement is repeatable. """

    ran = random.Random()
    ran.seed(t)

    for i, e in enumerate(sc.entity_list):

        fractions = sc.p.get_collection("fractions")
        e.change_type = fractions.parameters[0].name

        for f in fractions.parameters[1:]:
            draw = ran.random()
            if draw > 1 - f.get_in_sim_unit():
                e.change_type = f.name
                break

        e.p.add_parameter_with_collection(MiscParameter("id", int(i)))

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
# from parameters import numeric

templates = get_parameter_templates(numeric["unit_length_exponent"])

t_D = templates["D"]
t_R = templates["R"]
t_q = templates["q"]
t_kd = templates["kd"]

"""Sets up Scannable parameters from parameters templates"""
R = ScannablePhysicalParameter(t_R(40000), lambda x, v: x * v)
q = ScannablePhysicalParameter(t_q(1), lambda x, v: x * v)
D = ScannablePhysicalParameter(t_D(10), lambda x, v: x * v)
t_amax = templates["amax"]
amax = ScannablePhysicalParameter(t_amax(100), lambda x, v: x * v)
bc_type = ScannablePhysicalParameter(MiscParameter("bc_type", "linear"), lambda x, v: v)
kd = ScannablePhysicalParameter(t_kd(0.1), lambda x, v: x * v)
f = ScannablePhysicalParameter(PhysicalParameter("sec", 0.1, is_global=True), lambda x, v: v)

"""Retrieves and entity type from sim container for scanning"""
default = sc.get_entity_type_by_name("default")
abs = sc.get_entity_type_by_name("abs")
sec = sc.get_entity_type_by_name("sec")

# scan_container.add_single_sim_parameter_scan(f,"fractions", "", [0.01,0.1,0.5], scan_name = "fractions")

for v in np.logspace(-1,1,5):

    """Scans over parameters that are associated with a field"""
    sim_parameters = [
        # ParameterCollection("IL-2", [D(v)], field_quantity="il2"),
        # ParameterCollection("IL-2", [kd(v)], field_quantity="il2"),
        # ParameterCollection("fractions", [f(v)]),
    ]

    """Scans over parameters that are associated with an entity_type"""
    entity_types = [
        default.get_updated([
            ParameterCollection("IL-2", [R(v),q(v)])
        ]),
        abs.get_updated([ParameterCollection("IL-2", [])]),
        sec.get_updated([ParameterCollection("IL-2", [])]),
    ]
    """Scans over parameters that are associated with the outer domain
    This is a dictionary. If boundary pieces where defined in the setup function, they can be referenced by name.
    Here the pieces "left_boundary" and "box" are defined."""
    outer_domain_dict = {
    }
    """Creates container object for one sample of a the parameter scan.
    The Lists/Dicts can be empty for default parameters."""
    sample = ScanSample(sim_parameters, entity_types, outer_domain_dict, scan_name="entity_scan")
    scan_container.add_sample(sample)

"""signs up the internal solver with the sim container. 
It can be referenced in a cell_type definition by its name field
"""
sc.add_internal_solver(RuleBasedSolver)

"""State Manager updates the parameters of simulation objects in accordance with scan samples defined above and 
manages the orderly IO of simulation results and metadata for post processing."""

stMan = StateManager.StateManager(path)
stMan.sim_container = sc
sc.lookup = {"default":1, "sec":2, "abs":3}#labels to apper in marker function; 0 denotes background
stMan.scan_container = scan_container
stMan.compress_log_file = True
# stMan.dt = 1
# stMan.N = 5
"""sets up time range"""

# stMan.T = [0,0.2,0.4,0.6,0.8,1,2,3,10]
stMan.T = np.linspace(0,10,10)

"""defines a function which is called by StateManager before a parameter scan. 
Here it is used to assign cell types
"""

def pre_scan(state_manager, scan_index):
    updateState(state_manager.sim_container, 0)


stMan.pre_scan = pre_scan

"""Runs the ParameterScan"""

if len(sys.argv) > 1:
    if not sys.argv[1] == "mesh":
        stMan.run()
else:

    stMan.run()
