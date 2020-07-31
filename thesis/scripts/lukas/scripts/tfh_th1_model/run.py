import random
import sys

import numpy as np
from parameters import cytokines, cell_types_dict, geometry, numeric
from parameters import path, ext_cache

import thesis.main.StateManager
from thesis.main.InternalSolver import InternalSolver
from thesis.main.ParameterSet import MiscParameter, ParameterCollection, ScannablePhysicalParameter, PhysicalParameter, \
    PhysicalParameterTemplate
from thesis.main.ScanContainer import ScanContainer, ScanSample
from thesis.main.SimContainer import SimContainer
from thesis.scenarios.box_grid import setup


class RuleBasedSolver(InternalSolver):
    name = "RuleBasedSolver"

    def __init__(self):

        self.transition = (False, 0, "Tn")

        # self.il2_threshold = ths["il2"]#0.06
        # self.il6_threshold = ths["il2"]#0.05
        # self.infg_threshold = ths["il2"]#0.034

    def step(self, t, dt, p, entity=None):

        il2_threshold = p.get_physical_parameter("ths", "IL-2").get_in_post_unit()  # 0.06
        il6_threshold = p.get_physical_parameter("ths", "IL-6").get_in_post_unit()  # 0.05
        ifng_threshold = p.get_physical_parameter("ths", "IFNg").get_in_post_unit()  # 0.034

        p_diff = p.get_physical_parameter("p", "probs").get_in_post_unit()

        if entity.type_name == "Tn":
            if self.transition[0]:
                if self.transition[1] <= 0:
                    entity.change_type = self.transition[2]
                else:
                    self.transition = (True, self.transition[1] - dt, self.transition[2])
            elif np.random.rand(1) > p_diff:  # chance for type change; uniform distribution
                draw = np.random.normal(1, 0.2)  # time to type change; normal distribution

                il2 = p.get_physical_parameter("surf_c", "IL-2").get_in_post_unit()
                il6 = p.get_physical_parameter("surf_c", "IL-6").get_in_post_unit()
                ifng = p.get_physical_parameter("surf_c", "IFNg").get_in_post_unit()

                if il2 < il2_threshold and il6 > il6_threshold:  # il2neg and il6pos
                    self.transition = (True, draw, "Tfh")
                elif ifng > ifng_threshold:  # infg pos
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


scan_container = ScanContainer()
sc: SimContainer = setup(cytokines, cell_types_dict, geometry, numeric, path, ext_cache)

from thesis.scenarios.box_grid import t_R, t_q

R = ScannablePhysicalParameter(t_R(4000), lambda x, v: x * v)
q = ScannablePhysicalParameter(t_q(1), lambda x, v: x * v)
# D = ScannablePhysicalParameter(t_D(10), lambda x, v: x * v)

p_diff = PhysicalParameterTemplate(PhysicalParameter("p", 0.1, to_sim=1, is_global=True))
p = ScannablePhysicalParameter(p_diff(0.1), lambda x, v: v)

Tn = sc.get_entity_type_by_name("Tn")

# for v in np.logspace(-1,1,4):
for v in [0.9]:#np.linspace(0.1, 0.9, 2):
    sim_parameters = [
        ParameterCollection("probs", [p(v)])
    ]

    entity_types = [

    ]

    outer_domain_dict = {
        "left_boundary": [
            ParameterCollection("IL-2",[R(4000)]),
            ParameterCollection("IL-6",[q(1)])
        ]
    }

    sample = ScanSample(sim_parameters, entity_types, outer_domain_dict)
    scan_container.add_sample(sample)

sc.add_internal_solver(RuleBasedSolver)

stMan = StateManager.StateManager(path)
stMan.sim_container = sc
stMan.scan_container = scan_container
stMan.T = np.arange(0, 25, 1)
stMan.dt = 1


def pre_scan(state_manager, scan_index):
    updateState(state_manager.sim_container, 0)


stMan.pre_scan = pre_scan

"""Runs the ParameterScan"""
if len(sys.argv) > 1:
    if not sys.argv[1] == "mesh":
        stMan.run()
else:

    stMan.run()
