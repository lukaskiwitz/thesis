import os
import sys

sys.path.append("/home/lukas/thesis/main/")
sys.path.append("/home/lukas/thesis/scenarios/")

import random

import numpy as np
from parameters import cytokines, cell_types_dict, geometry, numeric, path, ext_cache
import logging

os.environ["LOG_PATH"] = path
LOG_PATH = os.environ.get("LOG_PATH") if os.environ.get("LOG_PATH") else "./"
os.makedirs(LOG_PATH,exist_ok=True)
logging.basicConfig(filename=LOG_PATH+"debug.log",level=logging.INFO,filemode="w", format='%(levelname)s::%(asctime)s %(message)s', datefmt='%I:%M:%S')

import thesis.main.StateManager as StateManager
from thesis.main.ParameterSet import MiscParameter, ParameterCollection, ScannablePhysicalParameter, PhysicalParameter
from thesis.main.ScanContainer import ScanContainer, ScanSample, ScanDefintion, ScanType
from thesis.main.SimContainer import SimContainer
from thesis.scenarios.box_grid import setup
import mpi4py.MPI as MPI
from thesis.main.assign_fractions import assign_fractions
comm = MPI.COMM_WORLD
rank = comm.Get_rank()
size = comm.Get_size()

def updateState(sc, t):

    assign_fractions(sc, t)



scan_container = ScanContainer()

sc: SimContainer = setup(cytokines, cell_types_dict, [], geometry, numeric, path, ext_cache)


from thesis.scenarios.box_grid import get_parameter_templates
from parameters import numeric


templates = get_parameter_templates(numeric["unit_length_exponent"])

t_D = templates["D"]
t_R = templates["R"]
t_q = templates["q"]
t_kd = templates["kd"]
t_amax = templates["amax"]
t_Kc = templates["Kc"]
t_kendo = templates["k_endo"]
t_koff = templates["k_off"]

from parameters import R_h, q
from parameters import rat as ratio, f_sec as fs, f_abs as fr

R = ScannablePhysicalParameter(t_R(R_h), lambda x, v: x * v)
q = ScannablePhysicalParameter(t_q(q), lambda x, v: x * v)
D = ScannablePhysicalParameter(t_D(10), lambda x, v: x * v)
kd = ScannablePhysicalParameter(t_kd(0.1), lambda x, v: x * v)

kendo = ScannablePhysicalParameter(t_kendo(1.1e-3), lambda x, v: x * v)
koff = ScannablePhysicalParameter(t_koff(0.83), lambda x, v: x * v)

f_sec = ScannablePhysicalParameter(PhysicalParameter("sec", ratio, is_global=True), lambda x,v: fs(v))
f_abs = ScannablePhysicalParameter(PhysicalParameter("abs", ratio, is_global=True), lambda x,v: fr(v))


"""Retrieves entity types from sim container"""
naive = sc.get_entity_type_by_name("naive")
abs = sc.get_entity_type_by_name("abs")
sec = sc.get_entity_type_by_name("sec")

s = 10

scan_space = np.concatenate([np.logspace(-1,0,int(s/2)), np.logspace(0,1,int(s/2)+1)[1:]])


f_sec_def = ScanDefintion(f_sec, "fractions", scan_space, ScanType.GLOBAL)
f_abs_def = ScanDefintion(f_abs, "fractions", scan_space, ScanType.GLOBAL)

D_def = ScanDefintion(D,"IL-2",  scan_space, ScanType.GLOBAL, field_quantity = "il2")
kd_def = ScanDefintion(kd,"IL-2",  scan_space, ScanType.GLOBAL, field_quantity = "il2")

sec_q_def = ScanDefintion(q,"IL-2",  scan_space, ScanType.ENTITY, field_quantity = "il2", entity_type = sec)
abs_R_def = ScanDefintion(R,"IL-2",  scan_space, ScanType.ENTITY, field_quantity = "il2", entity_type = abs)


kendo_def = ScanDefintion(kendo, "IL-2", scan_space, ScanType.GLOBAL, field_quantity = "il2")
koff_def = ScanDefintion(koff, "IL-2", scan_space, ScanType.GLOBAL, field_quantity = "il2")

# bc_type = ScannablePhysicalParameter(MiscParameter("bc_type", "linear"), lambda x, v: v)
# linear_solver = ScannablePhysicalParameter(MiscParameter("linear", True, is_global = True), lambda x, v: v)

for bc, linear in [("linear",True),("patrick_saturation",False)]:

    bc_def = lambda t: ScanDefintion(
        ScannablePhysicalParameter(MiscParameter("bc_type", "linear"), lambda x, v: bc),"IL-2",scan_space,ScanType.ENTITY,
        field_quantity="il2", entity_type=t
    )
    linear_def = ScanDefintion(ScannablePhysicalParameter(MiscParameter("linear", True, is_global = True), lambda x, v: linear),
                               "numeric", scan_space, ScanType.GLOBAL)

    if linear == False:
        scan_container.add_single_parameter_scan([kendo_def,bc_def(sec),bc_def(abs),linear_def], scan_name = "kendo")
        scan_container.add_single_parameter_scan([koff_def,bc_def(sec),bc_def(abs),linear_def], scan_name = "Koff")

    scan_container.add_single_parameter_scan([D_def,bc_def(sec),bc_def(abs),linear_def],scan_name = "D")
    scan_container.add_single_parameter_scan([kd_def,bc_def(sec),bc_def(abs),linear_def], scan_name = "kd")

    scan_container.add_single_parameter_scan([sec_q_def,bc_def(sec),bc_def(abs),linear_def], scan_name = "sec_q")
    scan_container.add_single_parameter_scan([abs_R_def,bc_def(sec),bc_def(abs),linear_def], scan_name = "abs_R")
    scan_container.add_single_parameter_scan([f_abs_def, f_sec_def,bc_def(sec),bc_def(abs),linear_def], scan_name = "ratio")


stMan = StateManager.StateManager(path)
stMan.sim_container = sc
stMan.scan_container = scan_container
stMan.dt = 1

stMan.T = [0,1,3,4]
# stMan.T = list(range(11))


def pre_scan(state_manager, scan_index):
    updateState(state_manager.sim_container, 0)


def pre_step(sc, time_index, t, T):
    updateState(sc, time_index)


stMan.pre_scan = pre_scan
sc.pre_step = pre_step


"""Runs the ParameterScan"""
stMan.sim_container.save_domain()
stMan.run()

import json
records = stMan.record.gather_records()
records_path = os.path.join(stMan.path, "records")
os.makedirs(records_path,exist_ok=True)

with open(os.path.join(records_path,"dump.json"),"w") as f:
    json.dump(records,f)