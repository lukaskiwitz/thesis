import logging
import os

import numpy as np

import thesis.main.StateManager as StateManager
from parameters import cytokines, cell_types_dict, geometry, path, ext_cache, numeric
from thesis.main.ParameterSet import MiscParameter, ParameterCollection, ScannableParameter, PhysicalParameter
from thesis.main.ScanContainer import ScanContainer, ScanSample
from thesis.scenarios.box_grid import setup, assign_fractions




# noinspection PyPep8Naming
def update_state(sim_container, replicat_index):
    assign_fractions(sim_container, replicat_index)
    sim_container.apply_type_changes(replicat_index)

    v = sim_container.p.get_physical_parameter("v", "v").get_in_sim_unit()
    R_per_cell = sim_container.p.get_physical_parameter("R_per_cell", "v").get_in_sim_unit()

    n_entities = sim_container.get_number_of_entities()
    sec_n = n_entities["sec"] if "sec" in n_entities.keys() else 0
    abs_n = n_entities["abs"] if "abs" in n_entities.keys() else 0

    n = sec_n + abs_n
    total_R = R_per_cell * n

    sec_R = (v * total_R / sec_n) if sec_n > 0 else 0
    abs_R = (total_R - sec_R * sec_n) / abs_n if abs_n > 0 else 0

    abs = sim_container.get_entity_type_by_name("abs")
    sec = sim_container.get_entity_type_by_name("sec")

    sim_container.add_entity_type(abs.get_updated(ParameterCollection("IL-2", [t_R(abs_R)], field_quantity="il2")))
    sim_container.add_entity_type(sec.get_updated(ParameterCollection("IL-2", [t_R(sec_R)], field_quantity="il2")))
    sim_container.reapply_entity_types(replicat_index)


scan_container = ScanContainer()
scenario = setup(cytokines, cell_types_dict, [], geometry, numeric, path)
pool = scenario.parameter_pool

t_R = pool.get_template("R")
bc_type = ScannableParameter(MiscParameter("bc_type", "linear"), lambda x, v: v)
linear_solver = ScannableParameter(MiscParameter("linear", True, is_global=True), lambda x, v: v)

"""Retrieves entity types from sim container"""
naive = scenario.get_entity_type_by_name("naive")
abs = scenario.get_entity_type_by_name("abs")
sec = scenario.get_entity_type_by_name("sec")

s = 30
R_space = np.linspace(0, 1, s)
for bc, linear in [("linear", True)]:
    RpC = 5e3
    for v in R_space:
        sample = ScanSample(
            [
                ParameterCollection("numeric", [linear_solver(linear)]),
                ParameterCollection("v", [
                        PhysicalParameter("v", v, is_global=True),
                        PhysicalParameter("R_per_cell", RpC, is_global=True),
                    ])
                ],
                [
                    abs.get_updated([ParameterCollection("IL-2", [bc_type(bc)], field_quantity="il2")]),
                    naive.get_updated([ParameterCollection("IL-2", [bc_type(bc)], field_quantity="il2")]),
                    sec.get_updated(ParameterCollection("IL-2", [bc_type(bc)], field_quantity="il2"))
                ],
            {},
            scan_name="scan_RpC_{rpc}_{l}".format(rpc=RpC, l=linear)
        )
        scan_container.add_sample(sample)

state_manager = StateManager.StateManager(path)
state_manager.scenario = scenario
state_manager.scan_container = scan_container
state_manager.dt = 1
state_manager.T = [0, 1]


def pre_replicat(sc, time_index, replicat_index, t, T):
    update_state(sc, replicat_index)


state_manager.pre_replicat = pre_replicat

"""Runs the ParameterScan"""
state_manager.run(ext_cache=ext_cache, model_names=["pde_model", "ode_model"], number_of_replicats=1)
