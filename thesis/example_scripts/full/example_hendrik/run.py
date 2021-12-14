import logging
import os

import numpy as np
from scipy.stats import poisson

import thesis.main.StateManager as StateManager
from parameters import cytokines, cell_types_dict, geometry, numeric, path, boundary
from thesis.main.InternalSolver import InternalSolver
from thesis.main.ParameterSet import ScannableParameter, PhysicalParameter
from thesis.main.ScanContainer import ScanContainer, ScanDefintion, ScanType
from thesis.scenarios.box_grid import setup, assign_fractions

os.makedirs(path, exist_ok=True)

"""
A solver is an object which is usually applied on a cell/entity. This solver must have a step function which acts upon 
and returns p. p is the parameter object associated with each cell and unique for each cell.
"""


class RuleBasedSolver(InternalSolver):
    name = "RuleBasedSolver"

    def on_type_change(self, p, replicat_index, entity=None):
        pass

    def step(self, t1, t2, dt, p, entity=None, **kwargs):

        il2_threshold = p.get_physical_parameter("ths", "IL-2").get_in_post_unit()
        il2 = p.get_physical_parameter("surf_c", "IL-2").get_in_post_unit()

        rand = np.random.uniform(0, 1)

        def hill(c, ths, invert=False):

            l_max = 1
            n = 4
            if invert:
                return l_max * (1 - (c ** n / (c ** n + ths ** n))) * dt
            else:
                return l_max * ((c ** n / (c ** n + ths ** n))) * dt

        il2_cdf = 1 - poisson.cdf(k=1, mu=hill(il2, il2_threshold, invert=True))

        if entity.type_name == "default":
            if il2_cdf > rand:
                entity.change_type = "sec"
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

"""Retrieves an entity type from scenario for scanning"""
default = scenario.get_entity_type_by_name("default")
abs = scenario.get_entity_type_by_name("abs")
sec = scenario.get_entity_type_by_name("sec")

"""
Sets up a parameter scan. ScannableParameter takes a function with two arguments, 
here conveniently a lambda function, to set the given PhysicalParameter according to input. 
A simpler definition would be 'lambda x, v: v'. 
scan_space sets up the values which are in the end inserted into the lambda function. ScanDefinition then combines
both while defining on which entity this scan should be applied to with ScanType.[...].
To actually run this scan setup we attach it to the scan_container with add_single_parameter_scan.
"""

"""log scan space centered around 1"""
s = 10
scan_space = np.concatenate([np.logspace(-1, 0, int(s / 2)), np.logspace(0, 1, int(s / 2) + 1)[1:]])

# scan over sec/abs ratio
f_sec = ScannableParameter(PhysicalParameter("sec", 1, is_global=True), lambda x, v: 1 / (v + 1))
f_abs = ScannableParameter(PhysicalParameter("abs", 1, is_global=True), lambda x, v: v / (v + 1))

f_sec_def = ScanDefintion(f_sec, "fractions", scan_space, ScanType.GLOBAL)
f_abs_def = ScanDefintion(f_abs, "fractions", scan_space, ScanType.GLOBAL)
scan_container.add_single_parameter_scan([f_abs_def, f_sec_def], scan_name="f")

"""
Templates, as defined in box_grid.py, allow for easy scanning and unit conversion.
"""
# scan over diffusion constant
t_D = parameter_pool.get_template("D")
D = ScannableParameter(t_D(10), lambda x, v: x * v)
D_def = ScanDefintion(D, "IL-2", scan_space, ScanType.GLOBAL, field_quantity="il2")
# scan_container.add_single_parameter_scan([D_def], scan_name="D")

# scan over secretion rate for sec-cells
t_q = parameter_pool.get_template("q")
q = ScannableParameter(t_q(1), lambda x, v: x * v)
sec_q_def = ScanDefintion(q, "IL-2", scan_space, ScanType.ENTITY, field_quantity="il2", entity_type=sec)
# scan_container.add_single_parameter_scan([naive_q], scan_name="q")

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
stMan.T = [0, 3600]

"""
Defines a function which is called by StateManager before a parameter scan. 
Here it is used to assign cell types.
With access to state_manager this can be used for any manipulation of the system pre scan.
"""


def pre_scan(state_manager, scan_index):
    update_state(state_manager.sim_container, 0)


stMan.pre_scan = pre_scan

"""Runs the ParameterScan"""
stMan.run()
