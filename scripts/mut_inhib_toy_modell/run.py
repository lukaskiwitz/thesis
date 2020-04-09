import getpass
import random
import sys

import numpy as np
from parameters import cytokines, cell_types_dict, geometry, numeric
from setup import setup

import StateManager
from InternalSolver import InternalSolver
from ParameterSet import MiscParameter, ParameterCollection, ScannablePhysicalParameter
from ScanContainer import ScanContainer, ScanSample
from SimContainer import SimContainer


class RuleBasedSolver(InternalSolver):
    name = "RuleBasedSolver"

    def __init__(self):
        pass

    def step(self, t, dt, p, entity=None):

        if entity.type_name == "default":
            il2 = p.get_physical_parameter("surf_c", "IL-2").get_in_post_unit()
            il2_t = p.get_misc_parameter("threshold", "IL-2").get_in_post_unit()

            if np.random.uniform(0, 1) > 0.9:
                if il2 > il2_t:
                    entity.change_type = "sec"
                else:
                    entity.change_type = "abs"
        return p


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


"""Setup/Simulation"""

user = getpass.getuser()
path = "/extra/{u}/mut_inhib_toy_model/".format(u=user)
ext_cache = "../mut_inhib_toy_model_ext_cache/"

scan_container = ScanContainer()

sc: SimContainer = setup(cytokines, cell_types_dict, geometry, numeric, path, ext_cache)

# R = ScannablePhysicalParameter(R(1000), lambda x, v: x * v)
# q = ScannablePhysicalParameter(q(1), lambda x, v: x * v)
# D = ScannablePhysicalParameter(D(10), lambda x, v: x * v)
# kd = ScannablePhysicalParameter(kd(0.1), lambda x, v: x * v)
th = ScannablePhysicalParameter(MiscParameter("threshold", 0.06), lambda x, v: x * v)
f = ScannablePhysicalParameter(MiscParameter("sec", 0.1, is_global=True), lambda x, v: v)

default = sc.get_entity_type_by_name("default")

for v in np.linspace(0.5, 1.5, 10):
    sim_parameters = [
        # ParameterCollection("IL-2", [D(100)], field_quantity="il2"),
        # ParameterCollection("IL-2", [kd(v)], field_quantity="il2")
    ]

    entity_types = [
        (default.get_updated([ParameterCollection("IL-2", [th(v)])])),
        # (default.get_updated([ParameterCollection("IL-2", [q(v)])]))
    ]

    outer_domain_dict = {
        # "left_boundary": [ParameterCollection("IL-2",[R(v)])],
        # "box": [ParameterCollection("IL-2",[R(v)])]
    }

    sample = ScanSample(sim_parameters, entity_types, outer_domain_dict)
    scan_container.add_sample(sample)

sc.add_internal_solver(RuleBasedSolver)

stMan = StateManager.StateManager(path)
stMan.sim_container = sc
stMan.scan_container = scan_container

stMan.T = np.arange(0, 30, 1)


def pre_scan(state_manager, scan_index):
    updateState(state_manager.sim_container, 0)


stMan.pre_scan = pre_scan

if len(sys.argv) > 1:
    if not sys.argv[1] == "mesh":
        stMan.run()
else:

    stMan.run()
