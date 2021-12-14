import numpy as np

import numpy as np

import thesis.main.StateManager as StateManager
from parameters import cytokines, cell_types_dict, geometry, numeric, path, ext_cache
from thesis.main.ParameterSet import ScannableParameter, PhysicalParameter, PhysicalParameterTemplate
from thesis.main.ScanContainer import ScanContainer, ScanDefintion, ScanType
from thesis.scenarios.box_grid import setup, assign_fractions


def updateState(sc, replicat_index):
    assign_fractions(sc, replicat_index)

    sc.apply_type_changes(replicat_index)

    R = np.unique(
        [e.p.get_physical_parameter("R", "IL-2").get_in_post_unit() for e in sc.entity_list if e.type_name == "abs"])
    assert len(R) == 1
    E = R[0]
    if E == 0:
        return None
    var = 1
    np.random.seed(replicat_index)
    tmp_sigma = np.sqrt(np.log((var * E) ** 2 / E ** 2 + 1))
    mean = np.log(E) - 1 / 2 * tmp_sigma ** 2
    for e in sc.entity_list:
        if e.type_name == "abs":
            R_draw = np.random.lognormal(mean, tmp_sigma)
            e.p.get_physical_parameter("R", "IL-2").set_in_post_unit(R_draw)


scan_container = ScanContainer()

from thesis.main.MyParameterPool import MyParameterPool

custom_pool = MyParameterPool()
custom_pool.add_template(PhysicalParameterTemplate(PhysicalParameter("my_p", 1, to_sim=1e-4)))
scenario = setup(cytokines, cell_types_dict, [], geometry, numeric, custom_pool=custom_pool)
pool = scenario.parameter_pool

t_D = pool.get_template("D")
t_R = pool.get_template("R")
t_q = pool.get_template("q")
t_kd = pool.get_template("kd")
t_amax = pool.get_template("amax")
t_Kc = pool.get_template("Kc")
t_kendo = pool.get_template("k_endo")
t_koff = pool.get_template("k_off")

from parameters import R_h, q
from parameters import rat as ratio, f_sec as fs, f_abs as fr

R_constant = ScannableParameter(t_R(R_h), lambda x, v: (5e3 * 0.9) / fr(v))
R = ScannableParameter(t_R(R_h), lambda x, v: x * v)
q_constant = ScannableParameter(t_q(q), lambda x, v: (30 * 0.1) / fs(v))

q = ScannableParameter(t_q(q), lambda x, v: x * v)
D = ScannableParameter(t_D(10), lambda x, v: x * v)
kd = ScannableParameter(t_kd(0.1), lambda x, v: x * v)

kendo = ScannableParameter(t_kendo(1.1e-3), lambda x, v: x * v)
koff = ScannableParameter(t_koff(0.83), lambda x, v: x * v)

f_sec = ScannableParameter(PhysicalParameter("sec", ratio, is_global=True), lambda x, v: fs(v))
f_abs = ScannableParameter(PhysicalParameter("abs", ratio, is_global=True), lambda x, v: fr(v))

"""Retrieves entity types from sim container"""
default = scenario.get_entity_type_by_name("default")
abs = scenario.get_entity_type_by_name("abs")
sec = scenario.get_entity_type_by_name("sec")

# s = 4
# fc = 10
# e = np.log10(fc) / np.log10(10)
# scan_space = np.concatenate([np.logspace(-e, 0, int(s / 2)), np.logspace(0, e, int(s / 2))[1:]])
#
# fc = 10
# e = np.log10(fc) / np.log10(10)
# scan_space_2 = np.concatenate([np.logspace(-e, 0, int(s / 2)), np.logspace(0, e, int(s / 2))[1:]])

scan_space = [0.1, 1, 10]
scan_space_2 = [0.1, 1, 10]

f_sec_def = ScanDefintion(f_sec, "fractions", scan_space, ScanType.GLOBAL)
q_constant_def = ScanDefintion(q_constant, "IL-2", scan_space, ScanType.ENTITY, field_quantity="il2", entity_type=sec)
R_constant_df = ScanDefintion(R_constant, "IL-2", scan_space, ScanType.ENTITY, field_quantity="il2", entity_type=abs)

f_abs_def = ScanDefintion(f_abs, "fractions", scan_space, ScanType.GLOBAL)

D_def = ScanDefintion(D, "IL-2", scan_space, ScanType.GLOBAL, field_quantity="il2")
kd_def = ScanDefintion(kd, "IL-2", scan_space, ScanType.GLOBAL, field_quantity="il2")

sec_q_def = ScanDefintion(q, "IL-2", scan_space, ScanType.ENTITY, field_quantity="il2", entity_type=sec)
abs_R_def = ScanDefintion(R, "IL-2", scan_space, ScanType.ENTITY, field_quantity="il2", entity_type=abs)

kendo_def = ScanDefintion(kendo, "IL-2", scan_space, ScanType.GLOBAL, field_quantity="il2")
koff_def = ScanDefintion(koff, "IL-2", scan_space, ScanType.GLOBAL, field_quantity="il2")

from thesis.main.ParameterSet import MiscParameter

d = lambda x, v: (x - 10) * v + 10
f = lambda n, d: np.ceil(n ** (1 / 3) * d + d)

distance = ScannableParameter(MiscParameter("distance", 20), d)
margin = ScannableParameter(MiscParameter("margin", 20), d)
n = 200

x_grid = ScannableParameter(MiscParameter("x_grid", 100), lambda x, v: f(n, d(20, v)))
y_grid = ScannableParameter(MiscParameter("y_grid", 100), lambda x, v: f(n, d(20, v)))
z_grid = ScannableParameter(MiscParameter("z_grid", 100), lambda x, v: f(n, d(20, v)))

distance_def = ScanDefintion(distance, "geometry", scan_space_2, ScanType.GLOBAL)
margin_def = ScanDefintion(margin, "geometry", scan_space_2, ScanType.GLOBAL)
x_def = ScanDefintion(x_grid, "geometry", scan_space_2, ScanType.GLOBAL)
y_def = ScanDefintion(y_grid, "geometry", scan_space_2, ScanType.GLOBAL)
z_def = ScanDefintion(z_grid, "geometry", scan_space_2, ScanType.GLOBAL)

for bc, linear in [("linear", True), ("patrick_saturation", False)]:
    bc_def = lambda t: ScanDefintion(
        ScannableParameter(MiscParameter("bc_type", "linear"), lambda x, v: bc), "IL-2", scan_space,
        ScanType.ENTITY,
        field_quantity="il2", entity_type=t
    )
    linear_def = ScanDefintion(
        ScannableParameter(MiscParameter("linear", True, is_global=True), lambda x, v: linear),
        "numeric", scan_space, ScanType.GLOBAL)

    scan_container.add_single_parameter_scan(
        [margin_def, x_def, y_def, z_def, bc_def(sec), bc_def(abs), linear_def, distance_def], scan_name="distance",
        remesh_scan_sample=True)

    if linear == False:
        scan_container.add_single_parameter_scan([kendo_def,bc_def(sec),bc_def(abs),linear_def], scan_name = "kendo")
        scan_container.add_single_parameter_scan([koff_def,bc_def(sec),bc_def(abs),linear_def], scan_name = "Koff")


    scan_container.add_single_parameter_scan([D_def,bc_def(sec),bc_def(abs),linear_def],scan_name = "D")
    scan_container.add_single_parameter_scan([kd_def, bc_def(sec), bc_def(abs), linear_def], scan_name="kd")

    scan_container.add_single_parameter_scan([sec_q_def,bc_def(sec),bc_def(abs),linear_def], scan_name = "sec_q")
    scan_container.add_single_parameter_scan([abs_R_def,bc_def(sec),bc_def(abs),linear_def], scan_name = "abs_R")

    scan_container.add_single_parameter_scan(
        [f_abs_def, f_sec_def, bc_def(sec), q_constant_def, R_constant_df, bc_def(abs), linear_def], scan_name="ratio")


def pre_replicat(sc, time_index, replicat_index, t, T):
    updateState(sc, replicat_index)


state_manager = StateManager.StateManager(path)
state_manager.scenario = scenario

state_manager.marker_lookup = {"default": 1, "sec": 2, "abs": 3}

state_manager.compress_log_file = True
state_manager.pre_replicat = pre_replicat

state_manager.scan_container = scan_container
state_manager.T = [0, 1]

state_manager.run(ext_cache=ext_cache)
