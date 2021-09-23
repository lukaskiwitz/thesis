import os
import time
import unittest
from shutil import rmtree
from typing import List

import numpy as np

from thesis.main.ParameterSet import ScannablePhysicalParameter
from thesis.main.ScanContainer import ScanDefintion, ScanType, ScanContainer
from thesis.main.StateManager import StateManager
from thesis.unit_tests.test_scenario1 import setup as scenario_setup

artifact_path = os.path.expanduser("/tmp/spatial_sim_unit_test/unit_tests_{ts}/".format(ts=time.time()))


class MyTest(unittest.TestCase):

    def setUp(self) -> None:
        path = os.path.abspath(os.path.join(artifact_path, ".."))
        if os.path.exists(path):
            rmtree(path)
        else:
            os.makedirs(path)

    def tearDown(self):
        rmtree(os.path.abspath(os.path.join(artifact_path, "..")))


class StaticSimContainerTest(MyTest):

    def static_pde_test(self):

        scenario = scenario_setup()
        sc = scenario.get_sim_container(None, model_index=0)
        sc.path = artifact_path

        def post_step(sc: 'SimContainer', time_index: int, t: float, T: List[float]) -> None:
            for f in sc.global_problems:
                f.get_result_element(time_index, 0, sc.path, [], {})

        sc.post_step = post_step
        sc.initialize()
        sc.run([0, 1])

    def static_ode_test(self):

        scenario = scenario_setup()
        sc = scenario.get_sim_container(None, model_index=1)
        sc.path = artifact_path

        def post_step(sc: 'SimContainer', time_index: int, t: float, T: List[float]) -> None:
            for f in sc.global_problems:
                f.get_result_element(time_index, 0, sc.path, [], {})

        sc.post_step = post_step
        sc.initialize()
        sc.run([0, 1])


class TimeSeriesSimContainerTest(MyTest):

    def pde_test(self):

        scenario = scenario_setup()
        sc = scenario.get_sim_container(None, model_index=0)
        sc.path = artifact_path

        def post_step(sc: 'SimContainer', time_index: int, t: float, T: List[float]) -> None:
            for f in sc.global_problems:
                f.get_result_element(time_index, 0, sc.path, [], {})

        sc.post_step = post_step
        sc.initialize()
        tu = 3600
        T = np.arange(0, tu * 10, tu)
        sc.run(T)

    def ode_test(self):

        scenario = scenario_setup()
        sc = scenario.get_sim_container(None, model_index=1)
        sc.path = artifact_path

        def post_step(sc: 'SimContainer', time_index: int, t: float, T: List[float]) -> None:
            for f in sc.global_problems:
                f.get_result_element(time_index, 0, sc.path, [], {})

        sc.post_step = post_step
        sc.initialize()
        tu = 3600
        T = np.arange(0, tu * 10, tu)
        sc.run(T)


class StateManagerTest(MyTest):

    def set_up_statemanager(self):
        scenario = scenario_setup()
        stMan = StateManager(artifact_path)
        stMan.scenario = scenario
        return stMan

    def get_scan_container(self, scenario):
        pool = scenario.parameter_pool
        t_q = pool.get_template("q")
        q = ScannablePhysicalParameter(t_q(10), lambda x, v: v * x)
        sec = scenario.get_entity_type_by_name("sec")
        scan_space = np.linspace(5, 15, 10)
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
        stMan.run(model_index=0)

    def timeseries_ode_test(self):
        stMan = self.set_up_statemanager()
        stMan.scan_tree.compress_log_file = True
        t_unit = 3600
        stMan.T = np.arange(0, t_unit * 10, t_unit)
        stMan.clear_log_files()
        stMan.run(model_index=1)

    def scan_pde_test(self):
        stMan = self.set_up_statemanager()
        stMan.scan_container = self.get_scan_container(stMan.scenario)
        stMan.scan_tree.compress_log_file = True

        stMan.T = [0, 1]
        stMan.clear_log_files()
        stMan.run(model_index=0)

    def scan_ode_test(self):
        stMan = self.set_up_statemanager()
        stMan.scan_container = self.get_scan_container(stMan.scenario)
        stMan.scan_tree.compress_log_file = True

        stMan.T = [0, 1]
        stMan.clear_log_files()
        stMan.run(model_index=1)
