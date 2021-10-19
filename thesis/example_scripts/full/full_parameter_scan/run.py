import logging
import os

import numpy as np

import thesis.main.StateManager as StateManager
from parameters import cytokines, cell_types_dict, geometry, numeric, path, ext_cache
from parameters import f_sec as fs, f_abs as fr, f_n as fn, R_h, q
from thesis.main.ParameterSet import MiscParameter, ScannableParameter, PhysicalParameter
from thesis.main.ScanContainer import ScanContainer, ScanDefintion, ScanType
from thesis.scenarios.box_grid import setup, assign_fractions

os.makedirs(path, exist_ok=True)
logging.basicConfig(
    filename=os.path.join(path, "sim.log"),
    level=logging.INFO,
    filemode="w",
    format='%(levelname)s::%(asctime)s %(message)s',
    datefmt='%I:%M:%S')


def updateState(sc, t):
    assign_fractions(sc, t)


scan_container = ScanContainer()

scenario = setup(cytokines, cell_types_dict, [], geometry, numeric, path)
pool = scenario.parameter_pool
R = ScannableParameter(pool.get_template("R")(R_h), lambda x, v: x * v)
q = ScannableParameter(pool.get_template("q")(q), lambda x, v: x * v)
D = ScannableParameter(pool.get_template("D")(10), lambda x, v: x * v)
kd = ScannableParameter(pool.get_template("kd")(0.1), lambda x, v: x * v)
Kc = ScannableParameter(pool.get_template("Kc")(0.01), lambda x, v: x * v)
f_sec = ScannableParameter(PhysicalParameter("sec", 1, is_global=True), lambda x, v: fs(v))
f_abs = ScannableParameter(PhysicalParameter("abs", 1, is_global=True), lambda x, v: fr(v))
f_naive = ScannableParameter(PhysicalParameter("naive", 1, is_global=True), lambda x, v: fn(v))

"""Retrieves entity types from sim container"""
naive = scenario.get_entity_type_by_name("naive")
abs = scenario.get_entity_type_by_name("abs")
sec = scenario.get_entity_type_by_name("sec")

s = 10
scan_space = np.concatenate([np.logspace(-1, 0, int(s / 2)), np.logspace(0, 1, int(s / 2) + 1)[1:]])

abs_R_def = ScanDefintion(R, "IL-2", scan_space, ScanType.ENTITY, field_quantity="il2", entity_type=abs)
sec_q_def = ScanDefintion(q, "IL-2", scan_space, ScanType.ENTITY, field_quantity="il2", entity_type=sec)
D_def = ScanDefintion(D, "IL-2", scan_space, ScanType.GLOBAL, field_quantity="il2")
kd_def = ScanDefintion(kd, "IL-2", scan_space, ScanType.GLOBAL, field_quantity="il2")
Kc_def = lambda t: ScanDefintion(Kc, "IL-2", scan_space, ScanType.ENTITY, field_quantity="il2", entity_type=t)
f_sec_def = ScanDefintion(f_sec, "fractions", scan_space, ScanType.GLOBAL)
f_abs_def = ScanDefintion(f_abs, "fractions", scan_space, ScanType.GLOBAL)
f_naive_def = ScanDefintion(f_naive, "fractions", scan_space, ScanType.GLOBAL)

for bc, linear in [("linear", True), ("R_saturation", False)]:

    bc_def = lambda t: ScanDefintion(
        ScannableParameter(MiscParameter("bc_type", "linear"), lambda x, v: bc), "IL-2", scan_space, ScanType.ENTITY,
        field_quantity="il2", entity_type=t
    )
    linear_def = ScanDefintion(ScannableParameter(MiscParameter("linear", True, is_global=True), lambda x, v: linear),
                               "numeric", scan_space, ScanType.GLOBAL)

    if linear == False:
        scan_container.add_single_parameter_scan(
            [Kc_def(naive), Kc_def(abs), Kc_def(sec), bc_def(naive), bc_def(sec), bc_def(abs), linear_def],
            scan_name="Kc")

    scan_container.add_single_parameter_scan([D_def, bc_def(naive), bc_def(sec), bc_def(abs), linear_def],
                                             scan_name="D")
    scan_container.add_single_parameter_scan([kd_def, bc_def(naive), bc_def(sec), bc_def(abs), linear_def],
                                             scan_name="kd")

    scan_container.add_single_parameter_scan([sec_q_def, bc_def(naive), bc_def(sec), bc_def(abs), linear_def],
                                             scan_name="sec_q")
    scan_container.add_single_parameter_scan([abs_R_def, bc_def(naive), bc_def(sec), bc_def(abs), linear_def],
                                             scan_name="abs_R")
    scan_container.add_single_parameter_scan(
        [f_abs_def, f_sec_def, bc_def(naive), bc_def(sec), bc_def(abs), linear_def], scan_name="f_new")

state_manager = StateManager.StateManager(path)
state_manager.scenario = scenario
state_manager.scan_container = scan_container
state_manager.dt = 1
state_manager.T = [0, 1]


def pre_replicat(sc, time_index, replicat_index, t, T):
    updateState(sc, replicat_index)


state_manager.pre_replicat = pre_replicat

"""Runs the ParameterScan"""
state_manager.run(ext_cache=ext_cache, model_names=["pde_model"], number_of_replicats=3)
