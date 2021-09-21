from thesis.unit_tests.test_scenario1 import setup as scenario_setup
import unittest
from typing import List
import os
import time
from shutil import rmtree

artifact_path = os.path.expanduser("~/tmp/spatial_sim_unit_test/unit_tests_{ts}/".format(ts=time.time()))

class TestSimContainer(unittest.TestCase):

    def scenario1_static_test(self):

        scenario = scenario_setup()
        sc = scenario.get_sim_container(None, model_index=0)
        sc.path = artifact_path


        def post_step(sc: 'SimContainer', time_index: int, t: float, T: List[float]) -> None:

            for f in sc.fields:
                f.get_result_element(time_index, 0, sc.path, [],{})

        sc.post_step = post_step
        sc.initialize()
        sc.run([0,1])

    def setUp(self) -> None:
        path = os.path.abspath(os.path.join(artifact_path,".."))
        if os.path.exists(path):
            rmtree(path)
    # def tearDown(self):
    #     rmtree(os.path.abspath(os.path.join(artifact_path,"..")))




