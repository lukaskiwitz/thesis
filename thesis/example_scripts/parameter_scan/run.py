import numpy as np

import thesis.main.StateManager as StateManager
from parameters import cytokines, cell_types_dict, geometry, numeric, boundary
from thesis.main.ParameterSet import ScannablePhysicalParameter
from thesis.main.ScanContainer import ScanContainer, ScanDefintion, ScanType
from thesis.main.SimContainer import SimContainer
from thesis.main.assign_fractions import assign_fractions
from thesis.scenarios.box_grid import setup


def updateState(sc, t):

    assign_fractions(sc, t)

path = "/extra/kiwitz/parameter_scan_example/test_1/"
ext_cache = r"../{mn}_ext_cache/".format(mn="parameter_scan_example")


scan_container = ScanContainer()
sc: SimContainer = setup(cytokines, cell_types_dict,boundary, geometry, numeric, path, ext_cache)

from thesis.scenarios.box_grid import get_parameter_templates
templates = get_parameter_templates(numeric["unit_length_exponent"])

s = 10
scan_space = np.concatenate([np.logspace(-1,0,int(s/2)), np.logspace(0,1,int(s/2)+1)[1:]])

t_D = templates["D"]
D = ScannablePhysicalParameter(t_D(10), lambda x, v: x * v)
D_def = ScanDefintion(D,"IL-2",  scan_space, ScanType.GLOBAL, field_quantity = "il2")
scan_container.add_single_parameter_scan([D_def], scan_name = "D")

t_R = templates["R"]
R = ScannablePhysicalParameter(t_R(1e4), lambda x, v: x * v)
abs = sc.get_entity_type_by_name("abs")
abs_R_def = ScanDefintion(R,"IL-2",  scan_space, ScanType.ENTITY, field_quantity = "il2", entity_type = abs)
scan_container.add_single_parameter_scan([abs_R_def], scan_name = "R_abs")

stMan = StateManager.StateManager(path)
stMan.sim_container = sc
stMan.scan_container = scan_container

stMan.T = [0,1]

def pre_scan(state_manager, scan_index):
    updateState(state_manager.sim_container, 0)

stMan.pre_scan = pre_scan
stMan.run()







