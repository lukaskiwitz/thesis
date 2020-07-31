try:
    import fenics as fcs
except RuntimeError:
    import os
    os.environ['PATH'] = '/home/brunner/anaconda3/envs/Lukas2/bin:/home/brunner/.local/bin:/home/brunner/anaconda3/condabin:/usr/local/bin:/usr/bin:/bin:/usr/local/games:/usr/games:/opt/puppetlabs/bin'
    import fenics as fcs

import getpass
import random
import sys

import numpy as np
from parameters_R_lognorm import cytokines, cell_types_dict, geometry, numeric, path, ext_cache

import StateManager
from InternalSolver import InternalSolver
from ParameterSet import MiscParameter, ParameterCollection, ScannablePhysicalParameter
from ScanContainer import ScanContainer, ScanSample
from SimContainer import SimContainer
from box_grid_R_lognorm import setup


class kineticSolver(InternalSolver):

    """Controls the internal logic of cells. An instance is created for each cell.
    Persistent values should be stored in the cell ParameterSet, values attached to the InternalSolver
    object may get lost."""

    name = "kineticSolver"

    def __init__(self):
        self.il2_threshold = 0.0048

    def step(self, t, dt, p, entity=None):

        """
        Is called for each cell after the PDEs are solved.
        Implements cell behaviour. Here cells of type "default" have a 10% transition
        probability if il2 surface concentration is greater than the threshold.
        """

        # if entity.type_name == "default":
        #     il2 = p.get_physical_parameter("surf_c", "IL-2").get_in_post_unit()
        #     if np.random.uniform(0,1) > 0.9 and il2 > self.il2_threshold:
        #         entity.change_type = "sec"
        return p


def updateState(sc, t):

    """sets cell types according to the values given in fractions.
    The pseudo random seed depends on t, so that cell placement is repeatable. """

    ran = random.Random()
    ran.seed(t)
    E = t_R.p.get_in_sim_unit()
    var = sc.p.get_physical_parameter("sigma", "IL-2").get_in_sim_unit()
    tmp_sigma = np.sqrt(np.log(var ** 2 / E ** 2 + 1))
    mean = np.log(E) - 1 / 2 * tmp_sigma ** 2

    for i, e in enumerate(sc.entity_list):
        e.p.add_parameter_with_collection(MiscParameter("id", int(i)))
        # e.change_type = "changed"
        e.change_type = "changed"
        e.p.get_physical_parameter("R", "IL-2").set_in_sim_unit(np.random.lognormal(mean, tmp_sigma))


"""Setup/Simulation"""

"""
setting filepath for simulation results. This is setup so that it works on the itb computers.
If ext_cache = "" the mesh will be cached for each field separately
"""
# ext_cache = "/extra/brunner/para_handling/static/R_lognorm/ext_cache/"
tmp_path = path
for j in range(30):
    user = getpass.getuser()
    path = tmp_path + "run" + str(j) + "/"

    """Setting up a parameters scan now has a object oriented interface. This is the container class"""
    scan_container = ScanContainer()

    """the setup function is defined in an external file. It builds the SimContainer for this simulation.
    The variable aspects are passed as list "cytokines, cell_types, geometry and numeric. 
    These are imported as modules and can be modified in parameters.py """

    sc: SimContainer = setup(cytokines, cell_types_dict, geometry, numeric, path, ext_cache)

    """Imports the parameter Templates"""
    from box_grid_R_lognorm import t_sigma, t_R #, t_fraction

    """Sets up Scannable parameters from parameters templates"""
    # R = ScannablePhysicalParameter(R(20000), lambda x, v: x * v)
    # q = ScannablePhysicalParameter(q(60), lambda x, v: x * v)
    # D = ScannablePhysicalParameter(D(10), lambda x, v: x * v)
    # kd = ScannablePhysicalParameter(D(0), lambda x, v: x * v)
    # fraction = ScannablePhysicalParameter(t_fraction(0.0), lambda x, v: v)
    sigma = ScannablePhysicalParameter(t_sigma(0.0), lambda x, v: v)

    """Retrieves and entity type from sim container for scanning"""
    # default = sc.get_entity_type_by_name("default")

    for v in np.linspace(0, 20000.0, 20): #np.logspace(-1,1,3):

        """Scans over parameters that are associated with a field"""
        sim_parameters = [
            ParameterCollection("IL-2", [sigma(v)], field_quantity="il2"),
            # ParameterCollection("IL-2", [kd(v)], field_quantity="il2")
        ]

        """Scans over parameters that are associated with an entity_type"""
        entity_types = [
        #     (default.get_updated([ParameterCollection("IL-2",[R(v)])])),
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
    sc.add_internal_solver(kineticSolver)

    """State Manager updates the parameters of simulation objects in accordance with scan samples defined above and 
    manages the orderly IO of simulation results and metadata for post processing."""

    stMan = StateManager.StateManager(path)
    stMan.sim_container = sc
    stMan.scan_container = scan_container

    """sets up time range"""
    stMan.T = np.arange(0, 1, 1)

    """defines a function which is called by StateManager before a parameter scan. 
    Here it is used to assign cell types
    """

    def pre_scan(state_manager, scan_index):
        updateState(state_manager.sim_container, 0)

    stMan.pre_scan = pre_scan

    """Runs the ParameterScan"""
    stMan.run()