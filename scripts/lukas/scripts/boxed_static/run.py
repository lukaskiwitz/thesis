import os
import sys

sys.path.append("/home/lukas/thesis/main/")
sys.path.append("/home/lukas/thesis/scenarios/")

import random

import numpy as np
from parameters import cytokines, cell_types_dict, geometry, numeric, path, ext_cache

os.environ["LOG_PATH"] = path

import thesis.main.StateManager as StateManager
from thesis.main.InternalSolver import InternalSolver
from thesis.main.ParameterSet import MiscParameter, ParameterCollection, ScannablePhysicalParameter
from thesis.main.ScanContainer import ScanContainer, ScanSample
from thesis.main.SimContainer import SimContainer
from thesis.scenarios.box_grid import setup
import mpi4py.MPI as MPI

comm = MPI.COMM_WORLD
rank = comm.Get_rank()
size = comm.Get_size()

def updateState(sc, t):


    ran = random.Random()
    ran.seed(t)

    for i, e in enumerate(sc.entity_list):

        fractions = sc.p.get_collection("fractions")
        e.change_type = fractions.parameters[0].name

        for f in fractions.parameters[1:]:
            draw = ran.random()
            if draw > 1 - f.get_in_sim_unit():
                e.change_type = f.name
                break

        e.p.add_parameter_with_collection(MiscParameter("id", int(i)))


scan_container = ScanContainer()

sc: SimContainer = setup(cytokines, cell_types_dict, geometry, numeric, path, ext_cache)


from thesis.scenarios.box_grid import get_parameter_templates
from parameters import numeric


templates = get_parameter_templates(numeric["unit_length_exponent"])

t_D = templates["D"]
t_R = templates["R"]
t_q = templates["q"]
t_kd = templates["kd"]
t_amax = templates["amax"]

R = ScannablePhysicalParameter(t_R(40000), lambda x, v: x * v)
q = ScannablePhysicalParameter(t_q(10), lambda x, v: x * v)
D = ScannablePhysicalParameter(t_D(10), lambda x, v: x * v)
amax = ScannablePhysicalParameter(t_amax(1), lambda x, v: x * v)
bc_type = ScannablePhysicalParameter(MiscParameter("bc_type", "linear"), lambda x, v: v)
is_linear = ScannablePhysicalParameter(MiscParameter("linear", True,is_global = True), lambda x, v: v)
kd = ScannablePhysicalParameter(t_kd(0.1), lambda x, v: x * v)
f = ScannablePhysicalParameter(MiscParameter("sec", 0.1, is_global=True), lambda x, v: v)

"""Retrieves entity types from sim container"""
naive = sc.get_entity_type_by_name("naive")
abs = sc.get_entity_type_by_name("abs")
sec = sc.get_entity_type_by_name("sec")


def sim_parameter_scan(scanable,collection_name, field_quantity, scan_space, scan_name = None):

    result = []
    for v in scan_space:

        sim_parameters = [
            ParameterCollection(collection_name, [scanable(v)], field_quantity=field_quantity),
        ]

        sample = ScanSample(sim_parameters, [], {},scan_name = scan_name)
        result.append(sample)
    return result

def entity_scan(entities,scanable,collection_name,field_quantity,scan_space, scan_name = None):

    result = []
    for v in scan_space:
        entity_types = []
        for e in entities:
            if not field_quantity is None:
                e = e.get_updated([ParameterCollection(collection_name, [scanable(v)], field_quantity=field_quantity)])

            e = e.get_updated([ParameterCollection(collection_name, [scanable(v)])])
            entity_types.append(e)

        sample = ScanSample([], entity_types, {}, scan_name=scan_name)
        result.append(sample)
    return result

s = 10
scan_space = np.logspace(-1,1,s)

for sample in sim_parameter_scan(D,"IL-2","il2",scan_space,scan_name = "D"):
    scan_container.add_sample(sample)
#
for sample in sim_parameter_scan(kd,"IL-2","il2",scan_space,scan_name = "kd"):
    scan_container.add_sample(sample)

for sample in entity_scan([sec],amax,"IL-2","il2",scan_space,scan_name = "sec_amax"):
    scan_container.add_sample(sample)

for sample in entity_scan([sec],q,"IL-2","il2",scan_space,scan_name = "sec_q"):
    scan_container.add_sample(sample)



# sc.add_internal_solver(RuleBasedSolver)


stMan = StateManager.StateManager(path)
stMan.sim_container = sc
stMan.scan_container = scan_container
stMan.dt = 1

stMan.T = np.arange(0, 1, 1)



def pre_scan(state_manager, scan_index):
    updateState(state_manager.sim_container, 0)


stMan.pre_scan = pre_scan

"""Runs the ParameterScan"""
stMan.sim_container.save_domain()

if len(sys.argv) > 1:
    if not sys.argv[1] == "mesh":
        stMan.run()
else:

    stMan.run()
