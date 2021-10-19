import logging
import os

from parameters import cytokines, cell_types_dict, geometry, numeric, boundary
from thesis.main.ParameterSet import ScannableParameter
from thesis.main.PostProcess import PostProcessor, ParaviewRender
from thesis.main.ScanContainer import ScanContainer, ScanDefintion, ScanType
from thesis.main.StateManager import StateManager
from thesis.scenarios.box_grid import setup, assign_fractions

"""Define paths for result folder structure"""
solution_path = "/extra/kiwitz/boundary_pieces_example/test_1/"
os.makedirs(solution_path, exist_ok=True)
logging.basicConfig(
    filename=os.path.join(solution_path, "sim.log"),
    level=logging.INFO,
    filemode="w",
    format='%(levelname)s::%(asctime)s %(message)s',
    datefmt='%I:%M:%S')

scenario = setup(
    cytokines,
    cell_types_dict,
    boundary,
    geometry,
    numeric)

"""scan over secretion rate for sec-cells"""
pool = scenario.parameter_pool

"""scans over secretion rate on named boundary piece"""
t_q = pool.get_template("q")
q = ScannableParameter(t_q(1), lambda x, v: v)
sec_q_def = ScanDefintion(q, "IL-2", [1, 5, 10], ScanType.BOUNDARY, field_quantity="il2",
                          boundary_pieces_name="top_boundary")

scan_container = ScanContainer()
scan_container.add_single_parameter_scan([sec_q_def], scan_name="top_il2_q")

state_manager = StateManager(solution_path)
state_manager.scan_container = scan_container
state_manager.scenario = scenario
state_manager.T = [0, 1]


def pre_scan(state_manager, scan_index):
    assign_fractions(state_manager.sim_container, 0)


state_manager.pre_scan = pre_scan
state_manager.run()

pp = PostProcessor(solution_path)
pp.unit_length_exponent = -6

pp.computations.append(ParaviewRender)
pp.run_post_process(4)
