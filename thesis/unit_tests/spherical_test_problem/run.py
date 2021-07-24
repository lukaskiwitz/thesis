import os
import numpy as np
from parameters import cytokine, geometry, numeric, path, boundary
import logging

os.environ["LOG_PATH"] = path
LOG_PATH = os.environ.get("LOG_PATH") if os.environ.get("LOG_PATH") else "./"
os.makedirs(LOG_PATH,exist_ok=True)
logging.basicConfig(filename=LOG_PATH+"debug.log",level=logging.INFO,filemode="w", format='%(levelname)s::%(asctime)s %(message)s', datefmt='%I:%M:%S')
os.environ["LOG_PATH"] = path

import thesis.main.StateManager as StateManager
from thesis.main.ParameterSet import MiscParameter, ScannablePhysicalParameter
from thesis.main.SimContainer import SimContainer
from thesis.scenarios.spherical_test_problem import setup, get_parameter_templates
from thesis.main.ScanContainer import ScanContainer, ScanDefintion, ScanType

scan_container = ScanContainer()
sc: SimContainer = setup(cytokine, boundary, geometry, numeric, path, "")
templates = get_parameter_templates(numeric["unit_length_exponent"])

max_char_length = ScannablePhysicalParameter(MiscParameter("max_char_length",5),lambda x,v : v)
scan_def_max_char = ScanDefintion(max_char_length, "numeric", [5,4,3,2,1],ScanType.GLOBAL)
# scan_def_max_char = ScanDefintion(max_char_length, "numeric", [5]*5,ScanType.GLOBAL)
scan_container.add_single_parameter_scan([scan_def_max_char],"max_char_length")

stMan = StateManager.StateManager(path)
stMan.sim_container = sc
sc.marker_lookup = {"default":1, "sec":2, "abs":3}#labels to apper in marker function; 0 denotes background
stMan.scan_container = scan_container
stMan.compress_log_file = True
sc.save_domain()

for f in stMan.sim_container.fields:
    f.remesh = True

stMan.T = [0,1]
stMan.run()




