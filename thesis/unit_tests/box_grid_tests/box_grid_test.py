import os
import time
import unittest
from copy import deepcopy
from shutil import rmtree

import numpy as np

from thesis.main.MyPlotter import Plotter
from thesis.main.ParameterSet import ScannableParameter
from thesis.main.PostProcess import PostProcessor
from thesis.main.ScanContainer import ScanDefintion, ScanType, ScanContainer
from thesis.main.StateManager import StateManager
from thesis.scenarios.box_grid import setup as scenario_setup

cytokines = [
    {
        "name": "IL-2",
        "field_quantity": "il2",
        "k_on": 111.6,  # receptor binding constant 1/(nM*h),
        "D": 10,  # Diffusion constant mu^2
        "kd": 0.1,  # cytokine decay in medium 1/h
        "k_endo": 1.1e-3,
        "k_off": 0.83,
    }
]

R_l = 1e2
R_h = 5e3
q = 10

cell_types_dict = [
    {"name": "naive",
     "fraction": 0,
     "il2": {"q": 0, "R": R_l, "bc_type": "patrick_saturation"},
     "internal_solver": ""
     },
    {"name": "sec",
     "fraction": 0.1,
     "il2": {"q": q, "R": R_l, "bc_type": "patrick_saturation"},
     "internal_solver": ""
     },
    {"name": "abs",
     "fraction": 0.8,
     "il2": {"q": 0, "R": R_h, "bc_type": "patrick_saturation"},
     "internal_solver": ""
     }
]

geometry = {
    "margin": 20,  # margin around the cell grid
    "distance": 20,  # distance between cell centers
    "rho": 5,  # cell radius
    "x_grid": 100,  # dimensions of the cell grid (edge length in um)
    "y_grid": 100,
    "z_grid": 100,  # comment out for single cell layer
    "norm_area": 4 * np.pi * 5 ** 2  # area passed to bc function of outside boundary
}

numeric = {

    "linear_solver": "gmres",
    "preconditioner": "amg",
    "linear": False,  # use linear fenics solver
    "krylov_atol": 1e-35,
    "krylov_rtol": 1e-5,  # linear solver relative tolerance
    "newton_atol": 1e-35,
    "newton_rtol": 1e-5,  # newton method relative tolerance
    "dofs_per_node": 15000,  # calc_boundary_values degrees of freedom per mpi node for pde solving
    "max_mpi_nodes": int(os.cpu_count() / 2),  # max nodes for fenics solver
    "cells_per_worker": 50,
    "max_pool_size": 4,  # max number of worker to extract boundary conditions at runtime
    "min_char_length": 1,  # mesh resolution
    "max_char_length": 5,  # mesh resolution
    "unit_length_exponent": -6  # for concentration conversion

}


class MyTest(unittest.TestCase):

    def get_artifact_path(self):

        artifact_path = os.path.expanduser(
            "/tmp/spatial_sim_unit_test/box_grid/{ts}_{cls}/".format(ts=time.time(), cls=self.__class__.__name__))
        return artifact_path

    def setUp(self) -> None:

        path = self.get_artifact_path()
        if os.path.exists(path):
            rmtree(path)
        else:
            os.makedirs(path)
        self.path = path

    def tearDown(self):
        pass
        # rmtree(os.path.abspath(os.path.join(artifact_path, "..")))


class PostProcessingTest(MyTest):

    def set_up_statemanager(self):
        scenario = scenario_setup(deepcopy(cytokines), deepcopy(cell_types_dict), [], deepcopy(geometry),
                                  deepcopy(numeric))

        stMan = StateManager(self.path)
        stMan.scenario = scenario
        return stMan

    def get_scan_container(self, scenario):
        pool = scenario.parameter_pool
        t_q = pool.get_template("q")
        q = ScannableParameter(t_q(10), lambda x, v: v * x)
        sec = scenario.get_entity_type_by_name("sec")
        scan_space = np.linspace(5, 15, 4)
        sec_q_def = ScanDefintion(q, "IL-2", scan_space, ScanType.ENTITY, field_quantity="il2",
                                  entity_type=sec)

        scan_container = ScanContainer()
        scan_container.add_single_parameter_scan([sec_q_def], scan_name="sec_q")

        return scan_container

    def test_static(self):
        stMan = self.set_up_statemanager()
        stMan.scan_tree.compress_xml_log_file = True

        stMan.T = [0, 1]
        stMan.clear_log_files()
        stMan.run(model_names=["pde_model", "ode_model"], number_of_replicats=3)

        post_processor = PostProcessor(self.path)
        post_processor.unit_length_exponent = -6
        post_processor.run_post_process(2, kde=False)

        plotter = Plotter(self.path)
        plotter.subplots(1, 1)
        plotter.global_steady_state_plot("Concentration", plotter.scan_index_key, style=plotter.model_name_key,
                                         legend="brief", ci="sd", xlog=False)
        plotter.make_legend()
        plotter.savefig(os.path.join(self.path, "images/") + "plot.png")

    def test_timeseries(self):
        stMan = self.set_up_statemanager()
        stMan.scan_tree.compress_xml_log_file = True
        t_unit = 3600
        stMan.T = np.arange(0, t_unit * 4, t_unit)
        stMan.clear_log_files()
        stMan.run(model_names=["pde_model", "ode_model"], number_of_replicats=3)

        post_processor = PostProcessor(self.path)
        post_processor.unit_length_exponent = -6
        post_processor.run_post_process(2, kde=True)

        plotter = Plotter(self.path)
        plotter.subplots(1, 1)
        plotter.global_time_series_plot("Concentration", style=plotter.model_name_key, legend="brief", ci="sd")
        plotter.make_legend()
        plotter.savefig(os.path.join(self.path, "images/") + "plot.png")

    def test_scan(self):
        stMan = self.set_up_statemanager()
        stMan.scan_container = self.get_scan_container(stMan.scenario)
        stMan.scan_tree.compress_xml_log_file = True

        stMan.T = [0, 1]
        stMan.clear_log_files()
        stMan.run(model_names=["pde_model", "ode_model"], number_of_replicats=3)

        post_processor = PostProcessor(self.path)
        post_processor.unit_length_exponent = -6
        post_processor.run_post_process(2, kde=False)

        plotter = Plotter(self.path)
        plotter.subplots(1, 1)
        plotter.global_steady_state_plot("Concentration", plotter.scan_index_key, style=plotter.model_name_key,
                                         legend="brief", ci="sd", xlog=False)
        plotter.make_legend()
        plotter.savefig(os.path.join(self.path, "images/") + "plot.png")
