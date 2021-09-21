from parameters import cytokines, cell_types_dict, geometry, numeric, boundary
from thesis.scenarios.box_grid import setup
from thesis.main.assign_fractions import assign_fractions
from thesis.main.StateManager import StateManager

path = "/extra/kiwitz/statemanager_example/test_1/"
sc = setup(
    cytokines,
    cell_types_dict,
    boundary,
    geometry,
    numeric,
    path,
    "")

stMan = StateManager(path)
stMan.sim_container = sc
stMan.T = [0,1]

def pre_scan(state_manager, scan_index):
    assign_fractions(state_manager.sim_container, 0)

stMan.pre_scan = pre_scan
stMan.run()