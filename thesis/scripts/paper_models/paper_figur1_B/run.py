import os
import sys

sys.path.append("/home/lukas/thesis/main/")
sys.path.append("/home/lukas/thesis/scenarios/")

import random

import numpy as np
from parameters import cytokines, cell_types_dict, geometry, numeric, path, ext_cache, boundary
from thesis.main.assign_fractions import assign_fractions

import logging
os.environ["LOG_PATH"] = path
LOG_PATH = os.environ.get("LOG_PATH") if os.environ.get("LOG_PATH") else "./"
os.makedirs(LOG_PATH,exist_ok=True)
logging.basicConfig(filename=LOG_PATH+"debug.log",level=logging.INFO,filemode="w", format='%(levelname)s::%(asctime)s %(message)s', datefmt='%I:%M:%S')

os.environ["LOG_PATH"] = path

import thesis.main.StateManager as StateManager
from thesis.main.InternalSolver import InternalSolver
from thesis.main.ParameterSet import MiscParameter, ParameterCollection, ScannablePhysicalParameter, PhysicalParameter
from thesis.main.ScanContainer import ScanContainer, ScanSample
from thesis.main.SimContainer import SimContainer
from thesis.scenarios.box_grid import setup

def updateState(sc, t):

    assign_fractions(sc, t)

scan_container = ScanContainer()


sc: SimContainer = setup(cytokines, cell_types_dict,boundary, geometry, numeric, path, ext_cache)
stMan = StateManager.StateManager(path)
stMan.sim_container = sc
sc.marker_lookup = {"default":1, "sec":2, "abs":3}#labels to apper in marker function; 0 denotes background
stMan.scan_container = scan_container
stMan.compress_log_file = True

stMan.T = [0,1]
# stMan.T = list(range(11))

def pre_scan(state_manager, scan_index):
    updateState(state_manager.sim_container, 0)

def pre_step(sc, time_index, t, T):
    updateState(sc, time_index)

stMan.pre_scan = pre_scan
sc.pre_step = pre_step
stMan.run()




