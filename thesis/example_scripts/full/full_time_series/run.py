import logging
import os

import numpy as np
from scipy.stats import poisson

import thesis.main.StateManager as StateManager
from parameters import cytokines, cell_types_dict, geometry, numeric, path, boundary, ext_cache
from thesis.main.InternalSolver import InternalSolver
from thesis.main.ScanContainer import ScanContainer
from thesis.scenarios.box_grid import setup, assign_fractions

os.makedirs(path, exist_ok=True)
logging.basicConfig(
    filename=os.path.join(path, "sim.log"),
    level=logging.INFO,
    filemode="w",
    format='%(levelname)s::%(asctime)s %(message)s',
    datefmt='%I:%M:%S')

"""
A solver is an object which is usually applied on a cell/entity. This solver must have a step function which acts upon 
and returns p. p is the parameter object associated with each cell and unique for each cell.
"""


class RuleBasedSolver(InternalSolver):
    name = "RuleBasedSolver"

    def on_type_change(self, p, replicat_index, entity=None):
        pass

    def step(self, t1, t2, dt, p, entity=None, **kwargs):

        def activation(c, R, R_M=860, max=0.125, min=0, n_R=0.55, n_il2=4):

            def rec(R):
                ec50 = (max - min) * (1 - np.power(R, n_R) / (np.power(R_M, n_R) + np.power(R, n_R))) + min
                return ec50

            ec50 = rec(R)
            a = c ** n_il2 / (ec50 ** n_il2 + c ** n_il2)

            if isinstance(R, float) and R == 0:
                return 0
            if isinstance(R, float) and isinstance(c, float):
                return a
            else:
                a[R == 0] = 0
                return a

        # il2_threshold = p.get_physical_parameter("ths", "IL-2").get_in_post_unit()
        il2 = p.get_physical_parameter("surf_c", "IL-2").get_in_post_unit()
        il2R = p.get_physical_parameter("R", "IL-2").get_in_post_unit()

        def hill(c, R, invert=False, l_max=1 / 3600, R_M=860):

            n = 10
            if invert:
                return l_max * (1 - activation(c, R, R_M=R_M)) * dt
            else:
                return l_max * (activation(c, R, R_M=R_M)) * dt

        il2_cdf = 1 - poisson.cdf(k=1, mu=hill(il2, il2R, invert=False, R_M=860))
        il2_cdf_r = 1 - poisson.cdf(k=1, mu=hill(il2, il2R, invert=False, R_M=860, l_max=3 / 3600))

        if entity.type_name == "default":
            if il2_cdf > np.random.uniform(0, 1):
                entity.change_type = "sec"
            elif il2_cdf_r > np.random.uniform(0, 1):
                entity.change_type = "abs"
        return p


def update_state(sc, replicat_index):
    """sets cell types according to the values given in fractions.
    The pseudo random seed depends on t, so that cell placement is repeatable. """

    assign_fractions(sc, replicat_index)
    sc.apply_type_changes(
        replicat_index)  # normally type change are applied after this method, but this should be fine as well


"""Setup/Simulation"""

"""Setting up a parameters scan now has a object oriented interface. This is the container class"""
scan_container = ScanContainer()

"""the setup function is defined in an external file. It builds the SimContainer for this simulation.
The variable aspects are passed as list "cytokines, cell_types, geometry and numeric. 
These are imported as modules and can be modified in parameters.py """

scenario = setup(cytokines, cell_types_dict, boundary, geometry, numeric)

"""Imports the parameter Templates"""
parameter_pool = scenario.parameter_pool

"""Retrieves an entity type from 03_scenario for scanning"""
default = scenario.get_entity_type_by_name("default")
abs = scenario.get_entity_type_by_name("abs")
sec = scenario.get_entity_type_by_name("sec")

"""
signs up the internal solver with the sim container. 
It can be referenced in a cell_type definition by its name field (see parameter file).
"""
scenario.internal_solvers = [RuleBasedSolver]

"""
State Manager updates the parameters of simulation objects in accordance with scan samples defined above and 
manages the orderly IO of simulation results and metadata for post processing.
"""

stMan = StateManager.StateManager(path)
stMan.scenario = scenario
scenario.marker_lookup = {"default": 1, "sec": 2, "abs": 3}  # labels to apper in marker function; 0 denotes background
stMan.scan_container = scan_container
stMan.compress_log_file = True

"""sets up time range"""
t_unit = 3600
stMan.T = np.arange(0, t_unit * 10, t_unit / 5)

"""
Defines a function which is called by StateManager before a parameter scan. 
Here it is used to assign cell types.
With access to state_manager this can be used for any manipulation of the system pre replicate.
"""


def pre_replicat(sim_container, time_index, replicat_index, t, T):
    update_state(sim_container, replicat_index)


stMan.pre_replicat = pre_replicat

"""Runs the ParameterScan"""
stMan.run(model_names=["pde_model", "ode_model"], ext_cache=ext_cache, number_of_replicats=3)
