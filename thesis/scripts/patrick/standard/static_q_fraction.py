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
from parameters_q_fraction import cytokines, cell_types_dict, geometry, numeric, path, ext_cache

os.environ["LOG_PATH"] = path

import thesis.main.StateManager as StateManager
from thesis.main.InternalSolver import InternalSolver
from thesis.main.ParameterSet import MiscParameter, ParameterCollection, ScannablePhysicalParameter
from thesis.main.ScanContainer import ScanContainer, ScanSample
from thesis.main.SimContainer import SimContainer
from box_grid_q_fraction import setup
import mpi4py.MPI as MPI

comm = MPI.COMM_WORLD
rank = comm.Get_rank()
size = comm.Get_size()

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
    # ran.seed(t)

    number_of_Tn_changed = int(round(len(sc.entity_list) * sc.p.get_physical_parameter("fraction", "IL-2").get_in_sim_unit(), 0))
    draws = np.random.choice(range(len(sc.entity_list)), number_of_Tn_changed, replace=False)
    no_of_secs = 0
    no_of_Treg = 0

    from scipy.constants import N_A
    q_il2_sum = len(sc.entity_list) * 60 * N_A ** -1 * 1e9

    for i, e in enumerate(sc.entity_list):
        e.change_type = "default"
        if i in draws:
            e.change_type = "changed"
            no_of_secs += 1
        e.p.get_physical_parameter("R", "IL-2").set_in_post_unit(20000.0)


    # Treg_draws = np.random.choice(np.setdiff1d(range(len(sc.entity_list)),draws), int(round(len(sc.entity_list) * 0.25, 0)), replace=False)
    # for i, e in enumerate(sc.entity_list):
    #     if i in Treg_draws:
    #         e.change_type = "Treg"
    #         no_of_Treg += 1
    #     e.p.add_parameter_with_collection(MiscParameter("id", int(i)))

    print("Number of secreting cells: ", no_of_secs)
    # print("Number of Tregs: ", no_of_Treg)
    if no_of_secs != 0:
        print("Cells q:", (q_il2_sum*N_A * 1e-9)/no_of_secs)
    else:
        print("cells q = 0")
    for i,e in enumerate(sc.entity_list):
        if no_of_secs != 0 and e.change_type == "changed":
            e.p.get_physical_parameter("q", "IL-2").set_in_sim_unit(q_il2_sum/no_of_secs)
    # sum = 0
    # for i,e in enumerate(sc.entity_list):
    #     if e.change_type == "sec":
    #         sum += e.p.get_physical_parameter("q", "IL-2").get_in_post_unit()
    # print(sum)


"""Setup/Simulation"""

"""
setting filepath for simulation results. This is setup so that it works on the itb computers.
If ext_cache = "" the mesh will be cached for each field separately
"""

# ext_cache = "/extra/brunner/para_handling/static/q_fraction/ext_cache/"

tmp_path = path
for j in range(15):
    user = getpass.getuser()
    path = tmp_path + "run" + str(j) + "/"
    # path = tmp_path + "run" + str(j) + "/"

    """Setting up a parameters scan now has a object oriented interface. This is the container class"""
    scan_container = ScanContainer()

    """the setup function is defined in an external file. It builds the SimContainer for this simulation.
    The variable aspects are passed as list "cytokines, cell_types, geometry and numeric. 
    These are imported as modules and can be modified in parameters.py """

    sc: SimContainer = setup(cytokines, cell_types_dict, geometry, numeric, path, ext_cache)

    """Imports the parameter Templates"""
    from box_grid_q_fraction import get_parameter_templates
    from parameters_q_fraction import numeric

    templates = get_parameter_templates(numeric["unit_length_exponent"])
    t_D = templates["D"]
    t_fraction = templates["fraction"]
    # t_run = templates["run"]

    """Sets up Scannable parameters from parameters templates"""
    # R = ScannablePhysicalParameter(R(20000), lambda x, v: x * v)
    # q = ScannablePhysicalParameter(q(60), lambda x, v: x * v)
    D = ScannablePhysicalParameter(t_D(10), lambda x, v: v)
    # kd = ScannablePhysicalParameter(D(0), lambda x, v: x * v)
    fraction = ScannablePhysicalParameter(t_fraction(0.0), lambda x, v: v)
    # run = ScannablePhysicalParameter(t_run(0), lambda x, v: v)

    """Retrieves and entity type from sim container for scanning"""
    default = sc.get_entity_type_by_name("default")
    changed = sc.get_entity_type_by_name("changed")
    # Treg = sc.get_entity_type_by_name("Treg")

    for v in np.linspace(0.1,1.0,10): #np.logspace(-1,1,3):
        sim_parameters = [
            ParameterCollection("IL-2", [D(10)], field_quantity="il2"),
            ParameterCollection("IL-2", [fraction(v)], field_quantity="il2"),
        ]
        """Scans over parameters that are associated with a field"""

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
    stMan.dt = 1
    stMan.T = np.arange(0, 1, 1)

    """defines a function which is called by StateManager before a parameter scan. 
    Here it is used to assign cell types
    """


    def pre_scan(state_manager, scan_index):
        updateState(state_manager.sim_container, 0)

    stMan.pre_scan = pre_scan

    """Runs the ParameterScan"""
    stMan.run()

