import getpass
import random

import numpy as np
from parameter_sets import cytokines, cell_types_dict, geometry, numeric
from setup import setup

import StateManager
from InternalSolver import InternalSolver
from ParameterSet import MiscParameter, ParameterCollection, ScannablePhysicalParameter
from ScanContainer import ScanContainer, ScanSample
from SimContainer import SimContainer


class RuleBasedSolver(InternalSolver):
    name = "RuleBasedSolver"

    def __init__(self):
        self.transition = (False, 0, "Tn")
        self.il2_threshold = 0.052
        self.il6_threshold = 0.015
        self.infg_threshold = 0.0092

    def step(self, t, dt, p, entity=None):

        if entity.type_name == "Tn":
            if self.transition[0]:
                if self.transition[1] <= 0:
                    entity.change_type = self.transition[2]
                else:
                    self.transition = (True, self.transition[1] - dt, self.transition[2])
            elif np.random.rand(1) > 0.9:  # chance for type change; uniform distribution
                draw = np.random.normal(1, 0.2)  # time to type change; normal distribution

                il2 = p.get_physical_parameter("surf_c", "IL-2").get_in_post_unit()
                il6 = p.get_physical_parameter("surf_c", "IL-6").get_in_post_unit()
                ifng = p.get_physical_parameter("surf_c", "IFNg").get_in_post_unit()

                if il2 < self.il2_threshold and il6 > self.il6_threshold:  # il2neg and il6pos
                    self.transition = (True, draw, "Tfh")
                elif ifng > self.infg_threshold:  # infg pos
                    self.transition = (True, draw, "Th1")
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


user = getpass.getuser()
path = "/extra/{u}/dev_testing/".format(u=user)
ext_cache = "../dev_testing_ext_cache/"

scan_container = ScanContainer()
sc: SimContainer = setup(cytokines, cell_types_dict, geometry, numeric, path, ext_cache)

from parameter_sets import R, q, D

R = ScannablePhysicalParameter(R(1000), lambda x, v: x * v)
q = ScannablePhysicalParameter(q(1), lambda x, v: x * v)
D = ScannablePhysicalParameter(D(10), lambda x, v: x * v)

Tn = sc.get_entity_type_by_name("Tn")

for v in [1]:  # np.logspace(-1,1,1):

    sim_parameters = [
        ParameterCollection("IL-2", [D(10)], field_quantity="il2"),
        ParameterCollection("IL-6", [D(10)], field_quantity="il6"),
        ParameterCollection("IFNg", [D(10)], field_quantity="ifng")
    ]

    entity_types = [
        # (Tn.get_updated([ParameterCollection("IL-2",[R(v)])]))
    ]

    outer_domain_dict = {
        # "left_boundary": [ParameterCollection("IL-2",[R(v)])]
    }

    sample = ScanSample(sim_parameters, entity_types, outer_domain_dict)
    scan_container.add_sample(sample)

sc.add_internal_solver(RuleBasedSolver)

stMan = StateManager.StateManager(path)
stMan.sim_container = sc
stMan.scan_container = scan_container
stMan.T = np.arange(0, 25, 1)


def pre_scan(state_manager, scan_index):
    updateState(state_manager.sim_container, 0)


stMan.pre_scan = pre_scan
stMan.run()
