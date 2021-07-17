import os
import sys


sys.path.append("/home/lukas/thesis/main/")
sys.path.append("/home/lukas/thesis/scenarios/")

import random

import numpy as np
from parameters import cytokines, cell_types_dict, geometry, path, ext_cache, numeric
import logging

os.environ["LOG_PATH"] = path
LOG_PATH = os.environ.get("LOG_PATH") if os.environ.get("LOG_PATH") else "./"
os.makedirs(LOG_PATH, exist_ok=True)
logging.basicConfig(filename=LOG_PATH + "debug.log", level=logging.INFO, filemode="w",
                    format='%(levelname)s::%(asctime)s %(message)s', datefmt='%I:%M:%S')

import thesis.main.StateManager as StateManager
from thesis.main.ParameterSet import MiscParameter, ParameterCollection, ScannablePhysicalParameter, ParameterSet, \
    PhysicalParameter
from thesis.main.ScanContainer import ScanContainer, ScanSample
from thesis.main.SimContainer import SimContainer
from thesis.scenarios.box_grid import setup
from thesis.main.assign_fractions import assign_fractions
import mpi4py.MPI as MPI

from thesis.scenarios.box_grid import get_parameter_templates

comm = MPI.COMM_WORLD
rank = comm.Get_rank()
size = comm.Get_size()


# noinspection PyPep8Naming
def update_state(sim_container, t):

    assign_fractions(sim_container,t)
    sim_container.apply_type_changes()

    v = sim_container.p.get_physical_parameter("v", "v").get_in_sim_unit()
    R_per_cell = sim_container.p.get_physical_parameter("R_per_cell", "v").get_in_sim_unit()

    n_entities = sim_container.get_number_of_entites()
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
    sim_container.apply_type_changes()


scan_container = ScanContainer()

sc: SimContainer = setup(cytokines, cell_types_dict, [], geometry, numeric, path, ext_cache)

# from thesis.scenarios.box_grid import get_parameter_templates
# from parameters import numeric


templates = get_parameter_templates(numeric["unit_length_exponent"])
t_R = templates["R"]

bc_type = ScannablePhysicalParameter(MiscParameter("bc_type", "linear"), lambda x, v: v)
linear_solver = ScannablePhysicalParameter(MiscParameter("linear", True, is_global = True), lambda x, v: v)

"""Retrieves entity types from sim container"""
naive = sc.get_entity_type_by_name("naive")
abs = sc.get_entity_type_by_name("abs")
sec = sc.get_entity_type_by_name("sec")

s = 10

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
                    abs.get_updated([ParameterCollection("IL-2", [bc_type(bc)],field_quantity="il2")]),
                    naive.get_updated([ParameterCollection("IL-2", [bc_type(bc)],field_quantity="il2")]),
                    sec.get_updated(ParameterCollection("IL-2", [bc_type(bc)], field_quantity="il2"))
                ],
                {},
                scan_name="scan_RpC_{rpc}_{l}".format(rpc=RpC, l=linear)
            )
            scan_container.add_sample(sample)

print(len(scan_container.scan_samples))
stMan = StateManager.StateManager(path)
stMan.sim_container = sc
stMan.scan_container = scan_container
stMan.dt = 1
stMan.T = [0,1]#list(range(11))


# noinspection PyUnusedLocal
def pre_scan(state_manager, scan_index):
    update_state(state_manager.sim_container, 0)


stMan.pre_scan = pre_scan

"""Runs the ParameterScan"""
stMan.sim_container.save_domain()
stMan.run()