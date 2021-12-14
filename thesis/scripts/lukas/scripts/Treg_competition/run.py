import random
import sys

import matplotlib.pyplot as plt
import numpy as np

from parameters import cytokines, cell_types_dict, geometry, numeric, path, ext_cache
from thesis.main.Entity import Cell
from thesis.main.InternalSolver import InternalSolver
from thesis.main.ParameterSet import MiscParameter, ParameterCollection, ScannableParameter, PhysicalParameter
from thesis.main.ScanContainer import ScanContainer, ScanSample
from thesis.main.SimContainer import SimContainer
from thesis.main.StateManager import StateManager
from thesis.scenarios.box_grid import setup


# if "LD_LIBRARY_PATH" in os.environ:
#     os.environ['LD_LIBRARY_PATH'] = "/home/lukas/anaconda3/envs/fenics/lib:"+os.environ['LD_LIBRARY_PATH']
# else:
#     os.environ['LD_LIBRARY_PATH'] = "/home/lukas/anaconda3/envs/fenics/lib"


class ResponseTimeSolver(InternalSolver):
    name = "ResponseTimeSolver"

    def __init__(self):
        self.il2_threshold = 0.012
        self.T = 0

    def step(self, t, dt, p, entity=None):
        import pandas as pd
        import scipy.stats as stats

        types = ["naive", "Treg", "sec"]
        Kc = p.get_physical_parameter("Kc", "IL-2").get_in_post_unit()

        def get_A():

            def f(T, il2):
                x = stats.gamma.cdf(T, 2, scale=1) * il2 / (Kc + il2)
                return x

            f0 = lambda x, c: 0

            A = [
                [f0, f0, f0],
                [f0, f0, f0],
                [f, f0, f0],
            ]

            A = pd.DataFrame(A, columns=types, index=types)
            return A

        A = get_A()

        il2 = p.get_physical_parameter("surf_c", "IL-2").get_in_post_unit()
        type_name = entity.type_name

        trans = False

        for k, g in A[type_name].items():

            if g(self.T, il2) > np.random.uniform(0, 1):
                trans = True
                entity.change_type = k
                self.T = 0
                break

        if not trans:
            self.T = self.T + dt

        return p


def updateState(sc, t):
    from thesis.main.MyKDE import get_kde_from_df, evalutate_kernel_on_grid, get_cell_df

    for i, e in enumerate(sc.entity_list):
        e.p.add_parameter_with_collection(MiscParameter("id", int(i)))

    ran = random.Random()
    ran.seed(t)
    np.random.seed(t)

    Treg_frac = sc.p.get_physical_parameter("Treg", "fractions").get_in_sim_unit()
    e_frac = sc.p.get_physical_parameter("sec", "fractions").get_in_sim_unit()
    clustering_strength = sc.get_entity_type_by_name("sec").p.get_physical_parameter("strength",
                                                                                     "clustering").get_in_sim_unit()
    bw = sc.get_entity_type_by_name("sec").p.get_physical_parameter("bw", "clustering").get_in_sim_unit()

    effectors = []
    cells = []

    def is_effector(cell):
        center = cell.center
        if ran.uniform(0, 1) < e_frac:
            return True
        else:
            return False

    for i, e in enumerate(sc.entity_list):

        if isinstance(e, Cell):

            fractions = sc.p.get_collection("fractions")
            e.change_type = fractions.parameters[0].name

            if is_effector(e):
                e.change_type = "sec"
                effectors.append(e)
            else:
                cells.append(e)

    draw = np.random.uniform(0, 1, len(sc.entity_list))

    effector_df = get_cell_df(effectors)
    cell_df = get_cell_df(cells)

    kernel, kernel_vis = get_kde_from_df(effector_df, "gaussian", bw)

    grid_points = 100
    x, y, v = evalutate_kernel_on_grid(kernel_vis, grid_points)

    plt.contourf(x, y, v[:, :], 100)
    plt.xlim([0, 300])
    plt.ylim([0, 300])
    plt.colorbar()

    for e_c in np.array([effector_df["x"], effector_df["y"]]).T:
        plt.gca().add_artist(plt.Circle(e_c, 5, color="blue"))

    treg_density = kernel.evaluate(np.transpose([cell_df["x"], cell_df["y"], cell_df["z"]]))

    n_cells = len(sc.entity_list)
    n_Treg = Treg_frac * n_cells

    sort_indices = np.flip(np.argsort(treg_density))
    cells = np.take_along_axis(np.array(cells), sort_indices, axis=0)
    treg_density = np.take_along_axis(np.array(treg_density), sort_indices, axis=0)

    from numpy.random import normal

    while n_Treg > 0 and len(cells) > 0:
        n_Treg -= 1

        high = int((1 - clustering_strength) * len(treg_density))

        i = np.random.randint(0, high) if high > 0 else 0

        cell = cells[i]
        cells = np.delete(cells, i)
        treg_density = np.delete(treg_density, i)

        cell.change_type = "Treg"
        plt.gca().add_artist(plt.Circle([cell.center[0], cell.center[1]], 5, color="red"))

    plt.show()
    print("")


"""Setup/Simulation"""

scan_container = ScanContainer()

sc: SimContainer = setup(cytokines, cell_types_dict, geometry, numeric, path, ext_cache)

from thesis.scenarios.box_grid import get_parameter_templates

templates = get_parameter_templates(numeric["unit_length_exponent"])

t_D = templates["D"]
t_R = templates["R"]
t_q = templates["q"]
t_kd = templates["kd"]
t_amax = templates["amax"]

# R = ScannableParameter(t_R(40000), lambda x, v: x * v)
q = ScannableParameter(t_q(100), lambda x, v: x * v)
D = ScannableParameter(t_D(10), lambda x, v: x * v)
kd = ScannableParameter(t_kd(0.1), lambda x, v: x * v)
amax = ScannableParameter(t_amax(100), lambda x, v: x * v)

c_s = ScannableParameter(PhysicalParameter("strength", 0.1, to_sim=1), lambda x, v: v)
f = ScannableParameter(PhysicalParameter("Treg", 0.1, is_global=True, to_sim=1), lambda x, v: v)
# bw = ScannableParameter(PhysicalParameter("bw", 10 ,to_sim = 1), lambda x,v: x*v)


default = sc.get_entity_type_by_name("naive")
effector = sc.get_entity_type_by_name("sec")
treg = sc.get_entity_type_by_name("Treg")


# for v in np.logspace(-1,1,10):
#
#     sim_parameters = [
#         ParameterCollection("IL-2", [D(v)], field_quantity="il2"),
#         # ParameterCollection("IL-2", [kd(v)], field_quantity="il2"),
#         # ParameterCollection("fractions", [f(v)]),
#     ]
#
#
#     entity_types = [
#         # (sec.get_updated([ParameterCollection("clustering",[c_s(v)])])),
#         # (sec.get_updated([ParameterCollection("clustering", [bw(v)])])),
#         # (sec.get_updated([ParameterCollection("IL-2", [q(v)])])),
#         # (treg.get_updated([ParameterCollection("IL-2",[R(v)])])),
#         # (naive.get_updated([ParameterCollection("IL-2", [q(v)])]))
#     ]
#
#     outer_domain_dict = {
#         # "left_boundary": [ParameterCollection("IL-2",[R(v)])],
#         # "box": [ParameterCollection("IL-2",[R(v)])]
#     }
#
#     sample = ScanSample(sim_parameters, entity_types, outer_domain_dict)
#     scan_container.add_sample(sample)

def sim_parameter_scan(scanable, collection_name, field_quantity, scan_space, scan_name=None):
    result = []
    assert isinstance(scanable, ScannableParameter)
    for v in scan_space:
        sim_parameters = [
            ParameterCollection(collection_name, [scanable(v)], field_quantity=field_quantity),
        ]

        sample = ScanSample(sim_parameters, [], {}, scan_name=scan_name)
        result.append(sample)
    return result


def entity_scan(entities, scanable, collection_name, field_quantity, scan_space, scan_name=None):
    result = []
    assert isinstance(scanable, ScannableParameter)
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
# scan_space = np.logspace(-1,1,s)
scan_space = np.linspace(0, 1, s)

# for sample in entity_scan([treg],amax,"IL-2","il2",scan_space,scan_name = "treg_amax"):
#     scan_container.add_sample(sample)
#
# for sample in entity_scan([sec],q,"IL-2","il2",scan_space,scan_name = "effector_q"):
#     scan_container.add_sample(sample)
#
# for sample in sim_parameter_scan(D,"IL-2","il2",scan_space,scan_name = "D"):
#     scan_container.add_sample(sample)
#
# for sample in sim_parameter_scan(kd, "IL-2", "il2", scan_space, scan_name="kd"):
#     scan_container.add_sample(sample)

# for sample in sim_parameter_scan(f,"fractions","Treg",scan_space,scan_name = "f"):
#     scan_container.add_sample(sample)
#
# for sample in entity_scan([sec],c_s,"clustering","",scan_space,scan_name = "cluster_strength"):
#     scan_container.add_sample(sample)

sc.add_internal_solver(ResponseTimeSolver)

stMan = StateManager(path)
stMan.sim_container = sc
stMan.scan_container = scan_container
stMan.dt = 1

stMan.T = np.arange(0, 30, 1)


def pre_scan(state_manager, scan_index):
    updateState(state_manager.sim_container, 0)


stMan.pre_scan = pre_scan

"""Runs the ParameterScan"""
if len(sys.argv) > 1:
    if not sys.argv[1] == "mesh":
        stMan.run()
else:

    stMan.run()
