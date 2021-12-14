import logging
import os

from parameters import cytokines, cell_types_dict, geometry, numeric, boundary
from thesis.scenarios.box_grid import setup, assign_fractions

"""Define paths for result folder structure"""
pde_solution_path = "/extra/kiwitz/simcontainer_example/test_1/pde_model"
ode_solution_path = "/extra/kiwitz/simcontainer_example/test_1/ode_model"

os.makedirs("/extra/kiwitz/simcontainer_example/test_1/", exist_ok=True)


scenario = setup(
    cytokines,
    cell_types_dict,
    boundary,
    geometry,
    numeric)

from thesis.main.ParameterSet import ParameterSet

"""pde model"""
sim_container = scenario.get_sim_container(ParameterSet("dummy", []), 0)
sim_container.path = pde_solution_path
sim_container.top_path = pde_solution_path

assign_fractions(sim_container, 0)
sim_container.initialize()
sim_container.run([0, 1])

for field in sim_container.global_problems:
    field.save_result_to_file(0, sim_container.path)

"""ode model"""
ode_solution_path = "/extra/kiwitz/simcontainer_example/test_1/ode_model"
sim_container = scenario.get_sim_container(ParameterSet("dummy", []), 1)
sim_container.path = ode_solution_path
sim_container.top_path = ode_solution_path

assign_fractions(sim_container, 0)
sim_container.initialize()
sim_container.run([0, 1])

for field in sim_container.global_problems:
    field.save_result_to_file(0, sim_container.path)
