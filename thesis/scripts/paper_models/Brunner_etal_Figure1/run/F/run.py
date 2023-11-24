import numpy as np

import numpy as np

import thesis.main.StateManager as StateManager
from parameters import cytokines, cell_types_dict, geometry, numeric, path, ext_cache
from thesis.main.ParameterSet import ScannableParameter, PhysicalParameter, PhysicalParameterTemplate
from thesis.main.ScanContainer import ScanContainer, ScanDefintion, ScanType
from thesis.scenarios.box_grid import setup, assign_fractions
from thesis.scripts.paper_models.utilities.states import updateState
from sympy import Integer


scan_container = ScanContainer()

from thesis.main.MyParameterPool import MyParameterPool

scenario = setup(cytokines, cell_types_dict, [], geometry, numeric)
pool = scenario.parameter_pool

t_D = pool.get_template("D")
t_R = pool.get_template("R")
t_q = pool.get_template("q")
t_kd = pool.get_template("kd")
t_amax = pool.get_template("amax")
t_Kc = pool.get_template("Kc")
t_kendo = pool.get_template("k_endo")
t_koff = pool.get_template("k_off")
from scipy.constants import N_A
t_sigma = PhysicalParameterTemplate(PhysicalParameter("sigma", 1e3, to_sim=N_A ** -1 * 1e9, is_global=True))
t_KD = pool.get_template("KD")


R = ScannableParameter(t_R(1e4), lambda x, v: x * v)
q = ScannableParameter(t_q(10), lambda x, v: x * v)
Tsec_fraction = ScannableParameter(PhysicalParameter("Tsec_fraction", 0.05, is_global=True), lambda x, v: x * v)
sigma = ScannableParameter(t_sigma(1), lambda x, v: x * v)
D = ScannableParameter(t_D(10), lambda x, v: x * v)
kd = ScannableParameter(t_kd(0.1), lambda x, v: x * v)
KD = ScannableParameter(t_KD(7.437 * 1e-3), lambda x, v: x * v)

kendo = ScannableParameter(t_kendo(0.00046), lambda x, v: x * v)
koff = ScannableParameter(t_koff(0.83), lambda x, v: x * v)


"""Retrieves entity types from sim container"""
default = scenario.get_entity_type_by_name("default")
Th = scenario.get_entity_type_by_name("Th")
Tsec = scenario.get_entity_type_by_name("Tsec")


no_of_scan_points = 9
standard_space = np.logspace(-1,1,no_of_scan_points)


R_scan_space = standard_space
q_scan_space = standard_space
Tsec_scan_space = standard_space
sigma_scan_space = standard_space
D_scan_space = standard_space
distance_scan_space = np.linspace(0.1, 5, no_of_scan_points)
kd_scan_space = standard_space
KD_scan_space = standard_space

R_def = ScanDefintion(R, "IL-2", R_scan_space, ScanType.ENTITY, field_quantity="il2", entity_type=Th)
q_def = ScanDefintion(q, "IL-2", q_scan_space, ScanType.ENTITY, field_quantity="il2", entity_type=Tsec)
Tsec_def = ScanDefintion(Tsec_fraction, "IL-2", Tsec_scan_space, ScanType.GLOBAL, field_quantity="il2")
sigma_def = ScanDefintion(sigma, "IL-2", sigma_scan_space, ScanType.ENTITY, field_quantity="il2", entity_type=Th)
D_def = ScanDefintion(D, "IL-2", D_scan_space, ScanType.GLOBAL, field_quantity="il2")

KD_def = ScanDefintion(KD, "IL-2", KD_scan_space, ScanType.ENTITY, field_quantity="il2", entity_type=Th)
# gamma_def = ScanDefintion(gamma, "IL-2", [1], ScanType.GLOBAL, field_quantity="il2")
# sigma_def = ScanDefintion(sigma, "IL-2", [1], ScanType.GLOBAL, field_quantity="il2")

kd_def = ScanDefintion(kd, "IL-2", kd_scan_space, ScanType.GLOBAL, field_quantity="il2")

# sec_q_def = ScanDefintion(q, "IL-2", scan_space, ScanType.ENTITY, field_quantity="il2", entity_type=sec)
# abs_R_def = ScanDefintion(R, "IL-2", scan_space, ScanType.ENTITY, field_quantity="il2", entity_type=Th)
#
# kendo_def = ScanDefintion(kendo, "IL-2", scan_space, ScanType.GLOBAL, field_quantity="il2")
# koff_def = ScanDefintion(koff, "IL-2", scan_space, ScanType.GLOBAL, field_quantity="il2")

from thesis.main.ParameterSet import MiscParameter

d = lambda x, v: (x - 10) * v + 10

def f(n, d):
    cells = Integer(n)**(Integer(1)/Integer(3))
    assert cells.is_integer, "only works for cubes"
    return cells * d + d

distance = ScannableParameter(MiscParameter("distance", 20), d)
margin = ScannableParameter(MiscParameter("margin", 20), d)
n = 1000

x_grid = ScannableParameter(MiscParameter("x_grid", 100), lambda x, v: f(n, d(20, v)))
y_grid = ScannableParameter(MiscParameter("y_grid", 100), lambda x, v: f(n, d(20, v)))
z_grid = ScannableParameter(MiscParameter("z_grid", 100), lambda x, v: f(n, d(20, v)))

distance_def = ScanDefintion(distance, "geometry", distance_scan_space, ScanType.GLOBAL)
margin_def = ScanDefintion(margin, "geometry", distance_scan_space, ScanType.GLOBAL)
x_def = ScanDefintion(x_grid, "geometry", distance_scan_space, ScanType.GLOBAL)
y_def = ScanDefintion(y_grid, "geometry", distance_scan_space, ScanType.GLOBAL)
z_def = ScanDefintion(z_grid, "geometry", distance_scan_space, ScanType.GLOBAL)

for bc, linear in [("patrick_saturation", False)]:
    bc_def = lambda t: ScanDefintion(
        ScannableParameter(MiscParameter("bc_type", "linear"), lambda x, v: bc), "IL-2", [1],
        ScanType.ENTITY,
        field_quantity="il2", entity_type=t
    )
    linear_def = ScanDefintion(
        ScannableParameter(MiscParameter("linear", True, is_global=True), lambda x, v: linear),
        "numeric", [1], ScanType.GLOBAL)

    scan_container.add_single_parameter_scan([R_def,bc_def(Tsec),bc_def(Th),linear_def], scan_name = "R")
    scan_container.add_single_parameter_scan([sigma_def,bc_def(Tsec),bc_def(Th),linear_def], scan_name = "sigma")
    scan_container.add_single_parameter_scan([q_def,bc_def(Tsec),bc_def(Th),linear_def], scan_name = "q")
    scan_container.add_single_parameter_scan([Tsec_def, bc_def(Tsec), bc_def(Th), linear_def], scan_name="T_sec")
    scan_container.add_single_parameter_scan([D_def,bc_def(Tsec),bc_def(Th),linear_def],scan_name = "D")
    scan_container.add_single_parameter_scan(
        [margin_def, x_def, y_def, z_def, bc_def(Tsec), bc_def(Th), linear_def, distance_def], scan_name = "distance",
        remesh_scan_sample=True)
    scan_container.add_single_parameter_scan([kd_def, bc_def(Tsec), bc_def(Th), linear_def], scan_name = "kd")
    scan_container.add_single_parameter_scan([KD_def, bc_def(Tsec), bc_def(Th), linear_def], scan_name="KD")


uS = updateState(0, 0, geometry, pool, offset=0)

def pre_replicat(sc, time_index, replicat_index, t, T):
    uS.step(sc, replicat_index)


state_manager = StateManager.StateManager(path)
state_manager.scenario = scenario

state_manager.marker_lookup = {"naive": 1, "Tsec": 2, "Th": 3}
scenario.markers = ["type_name", "IL-2_surf_c", "IL-2_R", "pSTAT5"]

state_manager.compress_log_file = True
state_manager.pre_replicat = pre_replicat

state_manager.scan_container = scan_container
state_manager.T = [0, 1]

state_manager.run(model_names=["pde_model"], ext_cache=ext_cache, number_of_replicats=10)
