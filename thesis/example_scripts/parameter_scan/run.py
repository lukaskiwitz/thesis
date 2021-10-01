import numpy as np

import thesis.main.StateManager as StateManager
from parameters import cytokines, cell_types_dict, geometry, numeric, boundary
from thesis.main.ParameterSet import ScannablePhysicalParameter
from thesis.main.ScanContainer import ScanContainer, ScanDefintion, ScanType
from thesis.scenarios.box_grid import setup, assign_fractions


def updateState(sc, t):
    assign_fractions(sc, t)


path = "/extra/kiwitz/parameter_scan_example/test_1/"
ext_cache = r"../{mn}_ext_cache/".format(mn="parameter_scan_example")

scan_container = ScanContainer()
scenario = setup(cytokines, cell_types_dict, boundary, geometry, numeric)
parameter_pool = scenario.parameter_pool

s = 10
scan_space = np.concatenate([np.logspace(-1, 0, int(s / 2)), np.logspace(0, 1, int(s / 2) + 1)[1:]])

t_D = parameter_pool.get_template("D")
D = ScannablePhysicalParameter(t_D(10), lambda x, v: x * v)
D_def = ScanDefintion(D, "IL-2", scan_space, ScanType.GLOBAL, field_quantity="il2")
scan_container.add_single_parameter_scan([D_def], scan_name="D")

t_R = parameter_pool.get_template("R")
R = ScannablePhysicalParameter(t_R(1e4), lambda x, v: x * v)
abs = scenario.get_entity_type_by_name("abs")
abs_R_def = ScanDefintion(R, "IL-2", scan_space, ScanType.ENTITY, field_quantity="il2", entity_type=abs)
scan_container.add_single_parameter_scan([abs_R_def], scan_name = "R_abs")

stMan = StateManager.StateManager(path)
stMan.scenario = scenario
stMan.scan_container = scan_container

stMan.T = [0,1]

def pre_scan(state_manager, scan_index):
    updateState(state_manager.sim_container, 0)

stMan.pre_scan = pre_scan
stMan.run()







