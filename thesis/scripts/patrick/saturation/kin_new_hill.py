try:
    import fenics as fcs
except RuntimeError:
    import os
    os.environ['PATH'] = '/home/brunner/anaconda3/envs/Lukas2/bin:/home/brunner/.local/bin:/home/brunner/anaconda3/condabin:/usr/local/bin:/usr/bin:/bin:/usr/local/games:/usr/games:/opt/puppetlabs/bin'
    import fenics as fcs

import getpass
import random
import sys
import os

sys.path.append("/home/brunner/thesis/thesis/main/")
sys.path.append("/home/brunner/thesis/thesis/scenarios/")

import numpy as np
from scipy.constants import N_A
from sympy import symbols, solve

from parameters_q_fraction import cytokines, cell_types_dict, geometry, numeric, path_kinetic, ext_cache

os.environ["LOG_PATH"] = path_kinetic

import thesis.main.StateManager as StateManager
from thesis.main.InternalSolver import InternalSolver
from thesis.main.ParameterSet import MiscParameter, ParameterCollection, ScannablePhysicalParameter, ParameterSet
from thesis.main.ScanContainer import ScanContainer, ScanSample
from thesis.main.SimContainer import SimContainer
from box_grid_q_fraction import setup
import mpi4py.MPI as MPI

comm = MPI.COMM_WORLD
rank = comm.Get_rank()
size = comm.Get_size()

class kineticSolver(InternalSolver):

    """Controls the internal logic of cells. An instance is created for each cell.
    Persistent values should be stored in the cell ParameterSet, values attached to the InternalSolver
    object may get lost."""

    name = "kineticSolver"
    #
    # def __init__(self):
    #     self.il2_threshold = 0.0048
    def step(self,t,t2,dt,p,entity=None):
        if t > 0: # since pre_step calculates c_0 in the first step we can only start at the second step
            # print(p.get_as_dictionary())
            # exit
            N = p.get_physical_parameter("hill_factor", "IL-2").get_in_post_unit()
            gamma = p.get_physical_parameter("gamma", "IL-2").get_in_sim_unit()
            # gamma = p["gamma"] #10

            eta = 1/72000 # 1/s
            # a = 10
            c_0 = p.get_physical_parameter("c0", "IL-2").get_in_post_unit() #1.90e-12 # M
            #
            # k_factor = 2# p.get_physical_parameter("k_factor", "IL-2").get_in_sim_unit()
            # k = k_factor * c_0 # M
            R_start = p.get_physical_parameter("R_start", "R_start").get_in_sim_unit()
            #
            # if c_0 == 0:
            #     nenner = 1
            # else:
            #     nenner = (gamma * c_0 ** N + k ** N) / (c_0 ** N + k ** N)
            # alpha = a * R_start * eta / nenner # R/t


            il2 = p.get_physical_parameter("surf_c", "IL-2").get_in_post_unit()*1e-9 #post is nM, neeeded in M
            R_il2 = p.get_physical_parameter("R", "IL-2").get_in_sim_unit()

            R_k_max = 4e-11 # max pSTAT5 EC50 in M
            R_k_min = 2e-12 # min pSTAT55 EC50 in M
            R_k = 1e4*N_A ** -1 * 1e9 # Receptors at half pSTAT5 maximum
            R_N = 1.2 # hill coefficient
            EC50 = (R_k_max * R_k ** R_N + R_k_min * R_il2 ** R_N) / (R_k ** R_N + R_il2 ** R_N)

            R_mean = 1e4 * N_A ** -1 * 1e9
            kmin = R_mean / gamma * eta
            kmax = R_mean * gamma * eta

            if p.get_physical_parameter("pSTAT5_signal", "IL-2").get_in_post_unit() == True:
                k_x = symbols("k_x")
                typical_pSTAT5 = c_0 ** 3 / (c_0 ** 3 + EC50 ** 3)
                k = solve((kmin * k_x ** N + kmax * typical_pSTAT5 ** N) / ((k_x ** N + typical_pSTAT5 ** N) * eta) - R_mean)[0]

                pSTAT5 = il2**3/(il2**3 + EC50**3)
                # print(pSTAT5)
                new_R_il2 = R_il2 + dt * ((kmin * k ** N + kmax * pSTAT5 ** N) / (k ** N + pSTAT5 ** N) - eta * R_il2)
                # print(new_R_il2 / (N_A ** -1 * 1e9))

                new_R_il2 = float(new_R_il2)
                p.get_physical_parameter("R", "IL-2").set_in_sim_unit(new_R_il2)
            else:
                k_x = symbols("k_x")
                try:
                    k = solve((kmin * k_x ** N + kmax * c_0 ** N) / ((k_x ** N + c_0 ** N) * eta) - R_mean)[0]
                except:  # if gamma ~ 1
                    k = 1

                new_R_il2 = R_il2 + dt * ((kmin * k ** N + kmax * il2 ** N) / (k ** N + il2 ** N) - eta * R_il2)

                new_R_il2 = float(new_R_il2)
                p.get_physical_parameter("R", "IL-2").set_in_sim_unit(new_R_il2)

        return p


def updateState(sc, t):
    """sets cell types according to the values given in fractions.
    The pseudo random seed depends on t, so that cell placement is repeatable. """
    global Tsec_distribution_array
    global Treg_distribution_array
    if len(Tsec_distribution_array) != 0:
        draws = Tsec_distribution_array
    else:
        tmp_fraction = sc.entity_list[0].p.get_physical_parameter("Tsec_fraction", "IL-2").get_in_post_unit()
        number_of_Tsec = int(round(len(sc.entity_list) * tmp_fraction))
        print("fraction = ", tmp_fraction)
        draws = np.random.choice(range(len(sc.entity_list)), number_of_Tsec, replace=False)
        Tsec_distribution_array = draws
    no_of_secs = 0
    no_of_Treg = 0

    from scipy.constants import N_A
    q_il2_sum = len(sc.entity_list) * 0.25 * 10 * N_A ** -1 * 1e9

    for i, e in enumerate(sc.entity_list):
        e.change_type = "Tnaive"
        if i in draws:
            e.change_type = "Tsec"
            no_of_secs += 1
            # e.p.get_physical_parameter("R", "IL-2").set_in_post_unit(5000)
            # print(e.p.get_physical_parameter("R", "IL-2").get_in_post_unit())
    if len(Treg_distribution_array) != 0 and 1 == 0:
        Treg_draws = Treg_distribution_array
    else:
        Treg_draws = np.random.choice(np.setdiff1d(range(len(sc.entity_list)), draws),
                                  int(round(len(sc.entity_list) * e.p.get_physical_parameter("Treg_fraction", "IL-2").get_in_sim_unit(), 0)), replace=False)
        Treg_distribution_array = Treg_draws
    for i, e in enumerate(sc.entity_list):
        if i in Treg_draws:
            e.change_type = "Treg"
            no_of_Treg += 1
        e.p.add_parameter_with_collection(MiscParameter("id", int(i)))
    print("Number of secreting cells: ", no_of_secs)
    print("Number of Tregs: ", no_of_Treg)
    if no_of_secs != 0:
        print("Cells q:", q_il2_sum / no_of_secs)
    else:
        print("cells q = 0")
    for i, e in enumerate(sc.entity_list):
        if no_of_secs != 0 and e.change_type == "Tsec":
            e.p.get_physical_parameter("q", "IL-2").set_in_sim_unit(q_il2_sum / no_of_secs)
            e.p.add_parameter_with_collection(t_R_start(1e4, in_sim=False))
        elif e.change_type == "Treg":
            e.p.add_parameter_with_collection(t_R_start(1e4, in_sim=False))
        elif e.change_type == "Tnaive":
            e.p.add_parameter_with_collection(t_R_start(1e2, in_sim=False))
    # sum = 0
    # for i,e in enumerate(sc.entity_list):
    #     if e.change_type == "sec":
    #         sum += e.p.get_physical_parameter("q", "IL-2").get_in_post_unit()
    # print(sum)
    #
    # for i, e in enumerate(sc.entity_list):
    #     e.change_type = "default"
    #     if i in draws:
    #         e.change_type = "changed"
    #         no_of_secs += 1
    #     e.p.add_parameter_with_collection(MiscParameter("id", int(i)))
    #     e.p.get_physical_parameter("R", "IL-2").set_in_post_unit(20000)
    # print("Number of secreting cells: ", no_of_secs)
    # if no_of_secs != 0:
    #     print("Cells q:", q_il2_sum/no_of_secs)
    # else:
    #     print("cells q = 0")
    # for i,e in enumerate(sc.entity_list):
    #     if no_of_secs != 0 and e.change_type == "changed":
    #         e.p.get_physical_parameter("q", "IL-2").set_in_sim_unit(q_il2_sum/no_of_secs)
    #     e.p.add_parameter_with_collection(t_R_start(20000, in_sim=False))


"""Setup/Simulation"""

"""
setting filepath for simulation results. This is setup so that it works on the itb computers.
If ext_cache = "" the mesh will be cached for each field separately
"""
# ext_cache = "/extra/brunner/para_handling/static/R_lognorm/ext_cache/"
# path = "/extra/brunner/thesis/kinetic/q_fraction_k_factor/"
path = path_kinetic


user = getpass.getuser()


"""Setting up a parameters scan now has a object oriented interface. This is the container class"""
scan_container = ScanContainer()

"""the setup function is defined in an external file. It builds the SimContainer for this simulation.
The variable aspects are passed as list "cytokines, cell_types, geometry and numeric. 
These are imported as modules and can be modified in parameters.py """

sc: SimContainer = setup(cytokines, cell_types_dict, geometry, numeric, path, ext_cache)

"""Imports the parameter Templates"""
from box_grid_q_fraction import get_parameter_templates
from parameters_q_fraction import numeric

templates = get_parameter_templates(numeric["unit_length_exponent"])
t_gamma = templates["gamma"]
t_hill_factor = templates["hill_factor"]
t_Tsec_fraction = templates["Tsec_fraction"]
t_R = templates["R"]
t_R_start = templates["R_start"]
t_D = templates["D"]
t_Treg_fraction = templates["Treg_fraction"]
t_c0 = templates["c0"]
t_pSTAT5_signal = templates["pSTAT5_signal"]
# t_amax = templates["amax"]

"""Sets up Scannable parameters from parameters templates"""
# R = ScannablePhysicalParameter(R(20000), lambda x, v: x * v)
# q = ScannablePhysicalParameter(q(60), lambda x, v: x * v)
# D = ScannablePhysicalParameter(t_D(10), lambda x, v: x * v)
# kd = ScannablePhysicalParameter(D(0), lambda x, v: x * v)
Tsec_fraction = ScannablePhysicalParameter(t_Tsec_fraction(0.0), lambda x, v: v)
Treg_fraction = ScannablePhysicalParameter(t_Treg_fraction(0.0), lambda x, v: v)
gamma = ScannablePhysicalParameter(t_gamma(0.0), lambda x, v: v)
hill_factor = ScannablePhysicalParameter(t_hill_factor(2), lambda x, v: v)
# amax = ScannablePhysicalParameter(t_amax(100), lambda x, v: x * v)
# bc_type = ScannablePhysicalParameter(MiscParameter("bc_type", "linear"), lambda x, v: v)
c0 = ScannablePhysicalParameter(t_c0(0.0), lambda x, v: v)
pSTAT5_signal = ScannablePhysicalParameter(t_pSTAT5_signal(False), lambda x, v: v)

Tsec_distribution_array = np.array([])
Treg_distribution_array = np.array([])

# a = [1/x for x in reversed(np.arange(10,110,10))]
a = [x for x in np.arange(10,110,10)]


for v in a: #np.linspace(0, 20000.0, 1): #np.logspace(-1,1,3):
    for frac in [0.25]:
        for hill_fac in [3]:
            """Scans over parameters that are associated with a field"""
            sim_parameters = [
                ParameterCollection("IL-2", [gamma(v)], field_quantity="il2"),
                ParameterCollection("IL-2", [Treg_fraction(0.75)], field_quantity="il2"),
                ParameterCollection("IL-2", [c0(8.5e-12)], field_quantity="il2"), #8.637363
                ParameterCollection("IL-2", [Tsec_fraction(frac)], field_quantity="il2"),
                ParameterCollection("IL-2", [hill_factor(hill_fac)], field_quantity="il2"),
                ParameterCollection("IL-2", [pSTAT5_signal(True)], field_quantity="il2"),
                # ParameterCollection("IL-2", [kd(v)], field_quantity="il2")
            ]

            """Scans over parameters that are associated with an entity_type"""
            entity_types = [
                # default.get_updated([ParameterCollection("IL-2", [bc_type(v)])]),
                # abs.get_updated([ParameterCollection("IL-2", [bc_type(v)])]),
                # sec.get_updated([ParameterCollection("IL-2", [bc_type(v)])]),
            ]
            """Scans over parameters that are associated with the outer domain
            This is a dictionary. If boundary pieces where defined in the setup function, they can be referenced by name. 
            Here the pieces "left_boundary" and "box" are defined."""
            outer_domain_dict = {
                # "left_boundary": [ParameterCollection("IL-2",[R(v)])],
                # "box": [ParameterCollection("IL-2",[R(v)])]
            }
            """Creates container object for one sample of a the parameter scan. 
            The Lists/Dicts can be empty for default parameters."""
            sample = ScanSample(sim_parameters, entity_types, outer_domain_dict)
            scan_container.add_sample(sample)

"""signs up the internal solver with the sim container. 
It can be referenced in a cell_type definition by its name field
"""
sc.add_internal_solver(kineticSolver)

"""State Manager updates the parameters of simulation objects in accordance with scan samples defined above and 
manages the orderly IO of simulation results and metadata for post processing."""

stMan = StateManager.StateManager(path)
stMan.sim_container = sc
stMan.scan_container = scan_container

"""sets up time range"""

# range_1 = np.arange(0,5*stMan.dt, stMan.dt*0.25)
# range_2 = np.arange(5*stMan.dt,15*stMan.dt, stMan.dt*0.5)
# range_3 = np.arange(16*stMan.dt, 30*stMan.dt, stMan.dt*2)
# range_4 = np.arange(32*stMan.dt, 120*stMan.dt, stMan.dt*5)
dt = 3600 #1h
length = 30
max_T = dt * 120

myRange = np.arange(0,length)
def exp_func(x,a,b,c):
    return a*np.exp(b*x) + c
a = 2*dt
c = -a
b = np.log((max_T-c)/a)/(length-1)

stMan.T = exp_func(myRange,a,b,c)

"""defines a function which is called by StateManager before a parameter scan. 
Here it is used to assign cell types
"""

def pre_scan(state_manager, scan_index):
    updateState(state_manager.sim_container, 0)

def pre_step(sc, time_index, t, T):
    # calculate the average surface concentration for the first solution
    if time_index == 1:
        list = []
        for e in sc.entity_list:
            # print(e.p.get_as_dictionary())
            list.append(e.p.get_physical_parameter("surf_c", "IL-2").get_in_post_unit()*1e-9)
        if len(list) != 0:
            for e in sc.entity_list:
                e.p.update(ParameterSet("dummy", [ParameterCollection("IL-2", [c0(np.mean(list))], field_quantity="il2")] ))
    return None

stMan.sim_container.pre_step = pre_step
stMan.pre_scan = pre_scan

"""Runs the ParameterScan"""

if len(sys.argv) > 1:
    if not sys.argv[1] == "mesh":
        stMan.run()
else:

    stMan.run()
