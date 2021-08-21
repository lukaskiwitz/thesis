import os
import sys

sys.path.append("/home/lukas/thesis/main/")
sys.path.append("/home/lukas/thesis/scenarios/")

import numpy as np
from parameters import cytokines, cell_types_dict, geometry, numeric, path, ext_cache, boundary

import logging
os.environ["LOG_PATH"] = path
LOG_PATH = os.environ.get("LOG_PATH") if os.environ.get("LOG_PATH") else "./"
os.makedirs(LOG_PATH,exist_ok=True)
logging.basicConfig(filename=LOG_PATH+"debug.log",level=logging.INFO,filemode="w", format='%(levelname)s::%(asctime)s %(message)s', datefmt='%I:%M:%S')

os.environ["LOG_PATH"] = path

import thesis.main.StateManager as StateManager
from thesis.main.InternalSolver import InternalSolver
from thesis.main.ParameterSet import ScannablePhysicalParameter, PhysicalParameter,PhysicalParameterTemplate
from thesis.main.ScanContainer import ScanContainer, ScanDefintion, ScanType
from thesis.main.SimContainer import SimContainer
from thesis.scenarios.scenario_box_grid import setup
from thesis.main.assign_fractions import assign_fractions
from scipy.stats import poisson

class RuleBasedSolver(InternalSolver):
    name = "RuleBasedSolver"

    def step(self, t1, t2, dt, p, entity=None):

        il2_threshold = p.get_physical_parameter("ths", "IL-2").get_in_post_unit()
        il2 = p.get_physical_parameter("surf_c", "IL-2").get_in_post_unit()

        rand = np.random.uniform(0,1)


        def hill(c,ths,invert = False):

            l_max = 1
            n = 4
            if invert:
                return l_max * (1-(c ** n / (c ** n + ths**n))) * dt
            else:
                return l_max * ((c ** n / (c ** n + ths**n))) * dt

        il2_cdf = 1-poisson.cdf(k = 1, mu = hill(il2, il2_threshold,invert = True))

        if entity.type_name == "default":
            if il2_cdf > rand:
                entity.change_type = "sec"
        return p



"""Setup/Simulation"""

"""
setting filepath for simulation results. This is setup so that it works on the itb computers.
If ext_cache = "" the mesh will be cached for each field separately
"""

"""Setting up a parameters scan now has a object oriented interface. This is the container class"""
scan_container = ScanContainer()

"""the setup function is defined in an external file. It builds the SimContainer for this simulation.
The variable aspects are passed as list "cytokines, cell_types, geometry and numeric. 
These are imported as modules and can be modified in parameters.py """

from thesis.main.MyParameterPool import MyParameterPool
custom_pool = MyParameterPool()
custom_pool.add_template(PhysicalParameterTemplate(PhysicalParameter("my_p",1,to_sim=1e-4)))
scenario = setup(cytokines, cell_types_dict, boundary, geometry, numeric, custom_pool=custom_pool)




"""Retrieves and entity type from sim container for scanning"""
default = scenario.get_entity_type_by_name("default")
abs = scenario.get_entity_type_by_name("abs")
sec = scenario.get_entity_type_by_name("sec")


"""log scan space centered around 1"""
# s = 10
# scan_space = np.concatenate([np.logspace(-1,0,int(s/2)), np.logspace(0,1,int(s/2)+1)[1:]])
scan_space = [0.1,1,10]

#scan over sec/abs ratio

# f_sec = ScannablePhysicalParameter(PhysicalParameter("sec", 1, is_global=True), lambda x, v: 1/(v+1))
# f_abs = ScannablePhysicalParameter(PhysicalParameter("abs", 1, is_global=True), lambda x, v: v/(v+1))
#
# f_sec_def = ScanDefintion(f_sec, "fractions", scan_space, ScanType.GLOBAL)
# f_abs_def = ScanDefintion(f_abs, "fractions", scan_space, ScanType.GLOBAL)


pool = scenario.parameter_pool

#scan over diffusion constant
# t_D = pool.get_template("D")

# D = ScannablePhysicalParameter(t_D(10), lambda x, v: x * v)
# D_def = ScanDefintion(D,"IL-2",  scan_space, ScanType.GLOBAL, field_quantity = "il2")


#scan over secretion rate for sec-cells
# t_q = pool.get_template("q")
# q = ScannablePhysicalParameter(t_q(1), lambda x, v: x * v)
# sec_q_def = ScanDefintion(q,"IL-2",  scan_space, ScanType.ENTITY, field_quantity = "il2", entity_type = sec)
#
# R_boundary = ScannablePhysicalParameter(pool.get_template("R")(4e4), lambda x,v: x * v)
# domain_R_def = ScanDefintion(R_boundary,"IL-2",scan_space,ScanType.BOUNDARY, boundary_pieces_name="left_boundary", field_quantity="il2")

from thesis.main.ParameterSet import MiscParameter
distance = ScannablePhysicalParameter(MiscParameter("distance",20), lambda x, v: (x-10) * v + 10)
margin = ScannablePhysicalParameter(MiscParameter("margin",20), lambda x, v: (x-10) * v + 10)

# s = 10
# scan_space = np.linspace(10,30,s)
# scan_space = [10,20,30]


distance_def = ScanDefintion(distance,"geometry",scan_space,ScanType.GLOBAL)
margin_def = ScanDefintion(margin,"geometry", scan_space ,ScanType.GLOBAL)


# scan_container.add_single_parameter_scan([f_abs_def, f_sec_def], scan_name = "f")
# scan_container.add_single_parameter_scan([D_def], scan_name = "D")
# scan_container.add_single_parameter_scan([sec_q_def], scan_name = "q")
# scan_container.add_single_parameter_scan([domain_R_def], scan_name = "boundary_R")

for bc, linear in [("linear",True),("patrick_saturation",False)]:

    bc_def = lambda t: ScanDefintion(
        ScannablePhysicalParameter(MiscParameter("bc_type", "linear"), lambda x, v: bc),"IL-2",scan_space,ScanType.ENTITY,
        field_quantity="il2", entity_type=t
    )
    linear_def = ScanDefintion(ScannablePhysicalParameter(MiscParameter("linear", True, is_global = True), lambda x, v: linear),
                               "numeric", scan_space, ScanType.GLOBAL)

    scan_container.add_single_parameter_scan([distance_def,margin_def], scan_name = "distance", dynamic_mesh = True)

def updateState(sc, t):
    assign_fractions(sc, t)

    sc.apply_type_changes()

    R = np.unique(
        [e.p.get_physical_parameter("R", "IL-2").get_in_post_unit() for e in sc.entity_list if e.type_name == "abs"])
    assert len(R) == 1
    E = R[0]
    var = 1
    np.random.seed(t)
    tmp_sigma = np.sqrt(np.log((var * E) ** 2 / E ** 2 + 1))
    mean = np.log(E) - 1 / 2 * tmp_sigma ** 2
    for e in sc.entity_list:
        if e.type_name == "abs":
            R_draw = np.random.lognormal(mean, tmp_sigma)
            e.p.get_physical_parameter("R", "IL-2").set_in_post_unit(R_draw)


def pre_scan(state_manager, scan_index):

    updateState(state_manager.sim_container, 0)

def pre_step(sim_container, time_index, t, T):

    updateState(sim_container,time_index)

stMan = StateManager.StateManager(path)
stMan.scenario = scenario

stMan.marker_lookup = {"default":1, "sec":2, "abs":3}#labels to apper in marker function; 0 denotes background

stMan.compress_log_file = True
stMan.pre_scan = pre_scan
stMan.pre_step = pre_step

stMan.scan_container = scan_container
stMan.T = [0,1]
# stMan.T = range(10)

stMan.run()








