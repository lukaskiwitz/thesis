import getpass
import random
import sys

import numpy as np

import StateManager
from InternalSolver import InternalSolver
from ParameterSet import MiscParameter, ParameterCollection, ScannablePhysicalParameter
from ScanContainer import ScanContainer, ScanSample
from SimContainer import SimContainer


class RuleBasedSolver(InternalSolver):

    """Controls the internal logic of cells. An instance is created for each cell.
    Persistent values should be stored in the cell ParameterSet, values attached to the InternalSolver
    object may get lost."""

    name = "RuleBasedSolver"

    def __init__(self):
        self.il2_threshold = 0.06

    def step(self, t, dt, p, entity=None):

        """
        Is called for each cell after the PDEs are solved.
        Implements cell behaviour. Here cells of type "default" have a 10% transition
        probability if il2 surface concentration is greater than the threshold.
        """

        if entity.type_name == "default":
            il2 = p.get_physical_parameter("surf_c", "IL-2").get_in_post_unit()
            if np.random.uniform(0, 1) > 0.5:
                if il2 > self.il2_threshold:
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
user = getpass.getuser()
path = "/extra/{u}/example_min/".format(u=user)
ext_cache = "../example_min_ext_cache/"

"""Setting up a parameters scan now has a object oriented interface. This is the container class"""
scan_container = ScanContainer()

"""the setup function is defined in an external file. It builds the SimContainer for this simulation.
The variable aspects are passed as list "cytokines, cell_types, geometry and numeric. 
These are imported as modules and can be modified in parameters.py """
from setup import setup
from parameters import cytokines, cell_types_dict, geometry, numeric

sc: SimContainer = setup(cytokines, cell_types_dict, geometry, numeric, path, ext_cache)

"""Imports the parameter Templates"""
from parameters import R, q, D, kd

"""Sets up Scannable parameters from parameters templates"""
R = ScannablePhysicalParameter(R(1000), lambda x, v: x * v)
q = ScannablePhysicalParameter(q(1), lambda x, v: x * v)
D = ScannablePhysicalParameter(D(10), lambda x, v: x * v)
kd = ScannablePhysicalParameter(kd(0.1), lambda x, v: x * v)

f = ScannablePhysicalParameter(MiscParameter("sec", 0.1, is_global=True), lambda x, v: v)

"""Retrieves and entity type from sim container for scanning"""
default = sc.get_entity_type_by_name("default")

for v in np.linspace(0.05, 0.3, 4):
    """Scans over parameters that are associated with a field"""
    sim_parameters = [
        # ParameterCollection("IL-2", [D(v)], field_quantity="il2"),
        # ParameterCollection("IL-2", [kd(v)], field_quantity="il2"),
        ParameterCollection("fractions", [f(v)]),
    ]

    """Scans over parameters that are associated with an entity_type"""
    entity_types = [
        # (default.get_updated([ParameterCollection("IL-2",[R(v)])])),
        # (default.get_updated([ParameterCollection("IL-2", [q(v)])]))
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
sc.add_internal_solver(RuleBasedSolver)

"""State Manager updates the parameters of simulation objects in accordance with scan samples defined above and 
manages the orderly IO of simulation results and metadata for post processing."""

stMan = StateManager.StateManager(path)
stMan.sim_container = sc
stMan.scan_container = scan_container

"""sets up time range"""
stMan.T = np.arange(0, 10, 1)

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
