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
from parameters_R_lognorm import cytokines, cell_types_dict, geometry, numeric, path_kinetic, ext_cache

os.environ["LOG_PATH"] = path_kinetic

import thesis.main.StateManager as StateManager
from thesis.main.InternalSolver import InternalSolver
from thesis.main.ParameterSet import MiscParameter, ParameterCollection, ScannablePhysicalParameter
from thesis.main.ScanContainer import ScanContainer, ScanSample
from thesis.main.SimContainer import SimContainer
from box_grid_R_lognorm import setup
import mpi4py.MPI as MPI

comm = MPI.COMM_WORLD
rank = comm.Get_rank()
size = comm.Get_size()

class kineticSolver(InternalSolver):

    """Controls the internal logic of cells. An instance is created for each cell.
    Persistent values should be stored in the cell ParameterSet, values attached to the InternalSolver
    object may get lost."""

    name = "kineticSolver"

    # def __init__(self):
        # self.il2_threshold = 0.0048

    def step(self,t,dt,p,entity=None):
        # print(p.get_as_dictionary())
        # exit()
        N=3
        gamma = p.get_physical_parameter("gamma", "IL-2").get_in_sim_unit()


        eta = 1/72000 # 1/s
        c_0 = 17.55e-12 # M

        k_factor = p.get_physical_parameter("k_factor", "IL-2").get_in_sim_unit()
        k = k_factor * c_0 # M
        R_start = p.get_physical_parameter("R_start", "R_start").get_in_sim_unit()

        nenner = (gamma * c_0 ** N + k ** N) / (c_0 ** N + k ** N)
        alpha = R_start * eta / nenner # R/t

        il2 = p.get_physical_parameter("surf_c", "IL-2").get_in_post_unit()*1e-9

        R_il2 = p.get_physical_parameter("R", "IL-2").get_in_sim_unit()

        # surf is molar, 1e-9
        # "R_l" is Mol, 1e-13
        new_R_il2 = R_il2 + dt * (alpha * (gamma * il2 ** N + k ** N) / (il2 ** N + k ** N) - eta * R_il2)
        p.get_physical_parameter("R", "IL-2").set_in_sim_unit(new_R_il2)

        return p


def updateState(sc, t):

    """sets cell types according to the values given in fractions.
    The pseudo random seed depends on t, so that cell placement is repeatable. """

    # ran = random.Random()
    # ran.seed(t)
    if len(R_distribution_dict) != 0:
        for i, e in enumerate(sc.entity_list):
            e.p.add_parameter_with_collection(MiscParameter("id", int(i)))
            # e.change_type = "changed"
            e.change_type = "changed"
            e.p.get_physical_parameter("R", "IL-2").set_in_sim_unit(R_distribution_dict[i])
            e.p.add_parameter_with_collection(t_R_start(R_distribution_dict[i], in_sim=True))
    else:
        E = t_R.p.get_in_sim_unit()
        var = t_sigma.p.get_in_sim_unit()
        tmp_sigma = np.sqrt(np.log(var ** 2 / E ** 2 + 1))
        mean = np.log(E) - 1 / 2 * tmp_sigma ** 2

        for i, e in enumerate(sc.entity_list):
            e.p.add_parameter_with_collection(MiscParameter("id", int(i)))
            # e.change_type = "changed"
            e.change_type = "changed"
            draw = np.random.lognormal(mean, tmp_sigma)
            e.p.get_physical_parameter("R", "IL-2").set_in_sim_unit(draw)
            e.p.add_parameter_with_collection(t_R_start(draw, in_sim=True))
            R_distribution_dict[i] = draw




"""Setup/Simulation"""

"""
setting filepath for simulation results. This is setup so that it works on the itb computers.
If ext_cache = "" the mesh will be cached for each field separately
"""
# ext_cache = "/extra/brunner/para_handling/static/R_lognorm/ext_cache/"
path = path_kinetic
# path = "/extra/brunner/thesis/kinetic/R_lognorm_k_factor/"

user = getpass.getuser()


"""Setting up a parameters scan now has a object oriented interface. This is the container class"""
scan_container = ScanContainer()

"""the setup function is defined in an external file. It builds the SimContainer for this simulation.
The variable aspects are passed as list "cytokines, cell_types, geometry and numeric. 
These are imported as modules and can be modified in parameters_old.py """

sc: SimContainer = setup(cytokines, cell_types_dict, geometry, numeric, path, ext_cache)

"""Imports the parameter Templates"""
from box_grid_R_lognorm import get_parameter_templates
from parameters_R_lognorm import numeric

templates = get_parameter_templates(numeric["unit_length_exponent"])
t_gamma = templates["gamma"]
t_k_factor = templates["k_factor"]
t_fraction = templates["fraction"]
t_R = templates["R"]
t_R_start = templates["R_start"]
t_kd = templates["kd"]
t_sigma = templates["sigma"]
# t_amax = templates["amax"]

"""Sets up Scannable parameters from parameters templates"""
# R = ScannablePhysicalParameter(R(20000), lambda x, v: x * v)
# q = ScannablePhysicalParameter(q(60), lambda x, v: x * v)
# D = ScannablePhysicalParameter(D(10), lambda x, v: x * v)
kd = ScannablePhysicalParameter(t_kd(0), lambda x, v: v)
# fraction = ScannablePhysicalParameter(t_fraction(0.0), lambda x, v: v)
gamma = ScannablePhysicalParameter(t_gamma(0.0), lambda x, v: v)
k_factor = ScannablePhysicalParameter(t_k_factor(2), lambda x, v: v)
# amax = ScannablePhysicalParameter(t_amax(100), lambda x, v: x * v)
# bc_type = ScannablePhysicalParameter(MiscParameter("bc_type", "linear"), lambda x, v: v)

R_distribution_dict = {}

for v in [0.5,2]: #np.linspace(0, 20000.0, 1): #np.logspace(-1,1,3):
    for kd_value in [0.0001, 0.001, 0.01, 0.1]:
        """Scans over parameters that are associated with a field"""
        sim_parameters = [
            ParameterCollection("IL-2", [gamma(v)], field_quantity="il2"),
            ParameterCollection("IL-2", [kd(kd_value)], field_quantity="il2"),
            ParameterCollection("IL-2", [k_factor(2)], field_quantity="il2")
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
stMan.dt = 3600*20

"""sets up time range"""
stMan.T = np.arange(0, 10, 1)

"""defines a function which is called by StateManager before a parameter scan. 
Here it is used to assign cell types
"""

def pre_scan(state_manager, scan_index):
    updateState(state_manager.sim_container, 0)

stMan.pre_scan = pre_scan

"""Runs the ParameterScan"""

if len(sys.argv) > 1:
    if not sys.argv[1] == "mesh":
        stMan.run()
else:

    stMan.run()
