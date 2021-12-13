import logging
import os

import numpy as np
from scipy.stats import poisson

import thesis.main.StateManager as StateManager
from parameters import cytokines, cell_types_dict, geometry, numeric, path, boundary
from thesis.main.InternalSolver import InternalSolver
from thesis.main.ParameterSet import ScannableParameter, MiscParameter
from thesis.main.ScanContainer import ScanContainer, ScanDefintion, ScanType
from thesis.scenarios.box_grid import setup, assign_fractions, distribute_receptors


os.makedirs(path, exist_ok=True)


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

scenario = setup(cytokines, cell_types_dict, boundary, geometry, numeric)

"""scan of cell-cell-distance"""
distance = ScannableParameter(MiscParameter("distance", 20), lambda x, v: v + 10)
margin = ScannableParameter(MiscParameter("margin", 20), lambda x, v: v)
scan_space = np.linspace(5, 15, 5)
distance_def = ScanDefintion(distance, "geometry", scan_space, ScanType.GLOBAL)
margin_def = ScanDefintion(margin, "geometry", scan_space, ScanType.GLOBAL)

"""add parameter scans to scan_container"""
scan_container.add_single_parameter_scan([distance_def, margin_def], scan_name="distance", remesh_scan_sample=True)


def update_state(sc, replicat_index):
    assign_fractions(sc, replicat_index)
    sc.apply_type_changes(replicat_index)
    entity_list = sc.entity_list

    for type_name in sc.get_number_of_entities().keys():
        distribute_receptors(entity_list, replicat_index, type_name)


def pre_replicat(sim_container, time_index, replicat_index, t, T):
    update_state(sim_container, replicat_index)


stMan = StateManager.StateManager(path)
stMan.scenario = scenario

stMan.marker_lookup = {"default": 1, "sec": 2, "abs": 3}  # labels to apper in marker function; 0 denotes background

stMan.compress_log_file = True
stMan.pre_replicat = pre_replicat

stMan.scan_container = scan_container
stMan.T = [0, 1]
stMan.run()
