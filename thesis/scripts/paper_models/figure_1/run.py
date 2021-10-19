import logging
import os

from parameters import cytokines, cell_types_dict, geometry, numeric, path, ext_cache, boundary
from thesis.main.StateManager import StateManager
from thesis.scenarios.box_grid import assign_fractions, setup

os.makedirs(path, exist_ok=True)
logging.basicConfig(
    filename=os.path.join(path, "sim.log"),
    level=logging.INFO,
    filemode="w",
    format='%(levelname)s::%(asctime)s %(message)s',
    datefmt='%I:%M:%S')


def update_state(sc, t):
    assign_fractions(sc, t)


scenario = setup(cytokines, cell_types_dict, boundary, geometry, numeric)
stMan = StateManager(path)
stMan.scenario = scenario

stMan.marker_lookup = {"default": 1, "sec": 2, "abs": 3}  # labels to apper in marker function; 0 denotes background
stMan.compress_log_file = True
stMan.T = [0, 1]


def pre_replicat(sc, time_index, replicat_index, t, T):
    update_state(sc, replicat_index)


stMan.pre_replicat = pre_replicat
stMan.run(ext_cache=ext_cache, number_of_replicats=1)
