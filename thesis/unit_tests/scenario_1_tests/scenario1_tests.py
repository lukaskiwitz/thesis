import os
import time
import unittest
from shutil import rmtree
from typing import List

import numpy as np

from thesis.main.MyPlotter import Plotter
from thesis.main.ParameterSet import ScannableParameter
from thesis.main.PostProcess import PostProcessor
from thesis.main.ScanContainer import ScanDefintion, ScanType, ScanContainer
from thesis.main.StateManager import StateManager
from thesis.unit_tests.test_scenario1 import setup as scenario_setup


class MyTest(unittest.TestCase):

    def get_artifact_path(self):

        artifact_path = os.path.expanduser(
            "/tmp/spatial_sim_unit_test/scenario1/{ts}_{cls}/".format(ts=time.time(), cls=self.__class__.__name__))
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


class StaticSimContainerTest(MyTest):

    def static_pde_test(self):

        scenario = scenario_setup()
        sc = scenario.get_sim_container(None, model_index=0)
        sc.path = self.path
        sc.top_path = self.path

        def post_step(sc: 'SimContainer', time_index: int, replicat_index: int, t: float, T: List[float]) -> None:
            for f in sc.global_problems:
                f.get_result_element(replicat_index, time_index, 0,
                                     sc.markers, sc.marker_lookup,
                                     sc.top_path,
                                     os.path.join(sc.path, "replicat_{}".format(replicat_index)))

        sc.post_step = post_step
        sc.initialize()
        sc.run([0, 1])

    def static_ode_test(self):

        scenario = scenario_setup()
        sc = scenario.get_sim_container(None, model_index=1)
        sc.path = self.path
        sc.top_path = self.path

        def post_step(sc: 'SimContainer', time_index: int, replicat_index: int, t: float, T: List[float]) -> None:
            for f in sc.global_problems:
                f.get_result_element(replicat_index, time_index, 0,
                                     sc.markers, sc.marker_lookup,
                                     sc.top_path,
                                     os.path.join(sc.path, "replicat_{}".format(replicat_index)))

        sc.post_step = post_step
        sc.initialize()
        sc.run([0, 1])


class TimeSeriesSimContainerTest(MyTest):

    def pde_test(self):

        scenario = scenario_setup()
        sc = scenario.get_sim_container(None, model_index=0)
        sc.path = self.path
        sc.top_path = self.path

        def post_step(sc: 'SimContainer', time_index: int, replicat_index: int, t: float, T: List[float]) -> None:
            for f in sc.global_problems:
                f.get_result_element(replicat_index, time_index, 0,
                                     sc.markers, sc.marker_lookup,
                                     sc.top_path,
                                     os.path.join(sc.path, "replicat_{}".format(replicat_index)))

        sc.post_step = post_step
        sc.initialize()
        tu = 3600
        T = np.arange(0, tu * 10, tu)
        sc.run(T)

    def ode_test(self):

        scenario = scenario_setup()
        sc = scenario.get_sim_container(None, model_index=1)
        sc.path = self.path
        sc.top_path = self.path

        def post_step(sc: 'SimContainer', time_index: int, replicat_index: int, t: float, T: List[float]) -> None:
            for f in sc.global_problems:
                f.get_result_element(replicat_index, time_index, 0,
                                     sc.markers, sc.marker_lookup,
                                     sc.top_path,
                                     os.path.join(sc.path, "replicat_{}".format(replicat_index)))

        sc.post_step = post_step
        sc.initialize()
        tu = 3600
        T = np.arange(0, tu * 10, tu)
        sc.run(T)


class StateManagerTest(MyTest):

    def set_up_statemanager(self):
        scenario = scenario_setup()
        stMan = StateManager(self.path)
        stMan.scenario = scenario
        return stMan

    def get_scan_container(self, scenario):
        pool = scenario.parameter_pool
        t_q = pool.get_template("q")
        q = ScannableParameter(t_q(10), lambda x, v: v * x)
        sec = scenario.get_entity_type_by_name("sec")
        scan_space = np.linspace(0.5, 1.5, 10)
        sec_q_def = ScanDefintion(q, "IL-2", scan_space, ScanType.ENTITY, field_quantity="il2",
                                  entity_type=sec)

        scan_container = ScanContainer()
        scan_container.add_single_parameter_scan([sec_q_def], scan_name="sec_q")

        return scan_container

    def timeseries_pde_test(self):
        stMan = self.set_up_statemanager()
        stMan.scan_tree.compress_log_file = True
        t_unit = 3600
        stMan.T = np.arange(0, t_unit * 10, t_unit)
        stMan.clear_log_files()
        stMan.run(model_names=["pde_model"])

    def timeseries_ode_test(self):
        stMan = self.set_up_statemanager()
        stMan.scan_tree.compress_log_file = True
        t_unit = 3600
        stMan.T = np.arange(0, t_unit * 10, t_unit)
        stMan.clear_log_files()
        stMan.run(model_names=["ode_model"])

    def scan_pde_test(self):
        stMan = self.set_up_statemanager()
        stMan.scan_container = self.get_scan_container(stMan.scenario)
        stMan.scan_tree.compress_log_file = True

        stMan.T = [0, 1]
        stMan.clear_log_files()
        stMan.run(model_names=["pde_model"])

    def scan_ode_test(self):
        stMan = self.set_up_statemanager()
        stMan.scan_container = self.get_scan_container(stMan.scenario)
        stMan.scan_tree.compress_log_file = True

        stMan.T = [0, 1]
        stMan.clear_log_files()
        stMan.run(model_names=["ode_model"])


class PostProcessingTest(MyTest):

    def set_up_statemanager(self):
        scenario = scenario_setup()
        stMan = StateManager(self.path)
        stMan.scenario = scenario
        return stMan

    def get_scan_container(self, scenario):
        pool = scenario.parameter_pool
        t_q = pool.get_template("q")
        q = ScannableParameter(t_q(10), lambda x, v: v * x)
        sec = scenario.get_entity_type_by_name("sec")
        scan_space = np.linspace(5, 15, 10)
        sec_q_def = ScanDefintion(q, "cytokine", scan_space, ScanType.ENTITY, field_quantity="cyt",
                                  entity_type=sec)

        scan_container = ScanContainer()
        scan_container.add_single_parameter_scan([sec_q_def], scan_name="sec_q")

        return scan_container

    def timeseries_ode_test(self):
        stMan = self.set_up_statemanager()
        stMan.scan_tree.compress_log_file = True
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

    def scan_ode_test(self):
        stMan = self.set_up_statemanager()
        stMan.scan_container = self.get_scan_container(stMan.scenario)
        stMan.scan_tree.compress_log_file = True

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
