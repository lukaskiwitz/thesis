import logging
import os

import numpy as np

import thesis.main.StateManager as StateManager
from parameters import cytokines, cell_types_dict, geometry, numeric, path, ext_cache
from thesis.main.ParameterSet import MiscParameter, ScannableParameter, PhysicalParameter
from thesis.main.ScanContainer import ScanContainer, ScanDefintion, ScanType
from thesis.scenarios.box_grid import setup

os.makedirs(path, exist_ok=True)
logging.basicConfig(
    filename=os.path.join(path, "sim.log"),
    level=logging.INFO,
    filemode="w",
    format='%(levelname)s::%(asctime)s %(message)s',
    datefmt='%I:%M:%S')


def update_state(sc, replicat_index):
    """
    Grid clustering is implemented using cellBehaviourUtilities.grid_clustering.
    Note that cell fractions are controlled by the clustering function.
    Arbitrary apc positions can be specified.

    """

    from thesis.cellBehaviourUtilities.grid_clustering import make_clusters
    np.random.seed(replicat_index)  # fix numpy seed for reproducibility

    """extract grid/cell positions from entity list"""
    cell_grid_positions = []
    for i, e in enumerate(sc.entity_list):
        cell_grid_positions.append(e.center)
        e.p.add_parameter_with_collection(MiscParameter("id", int(i)))

    cell_grid_positions = np.array(cell_grid_positions)

    """specify apc (clustering centers) position"""
    apcs = np.array([[80, 80, 80]])
    # apcs = get_apc_positions(cell_grid_positions, no_apcs=5)

    """extract parameters for clustering for simulation parameter set(s)"""
    fractions_dict = sc.p.get_collection("fractions").get_as_dictionary()
    del fractions_dict["naive"]  # does not cluster
    fractions = [sc.p.get_physical_parameter(k, "fractions").get_in_sim_unit() for k in fractions_dict.keys()]
    cluster_strengths = [sc.get_entity_type_by_name(k).p.get_physical_parameter("strength", "clustering") for k in
                         fractions_dict.keys()]
    cluster_strengths = [i.get_in_sim_unit() if i is not None else 1 for i in cluster_strengths]

    """generates clustered cell type assignments"""
    cell_type = make_clusters(cell_grid_positions, apcs, fractions, cluster_strengths)

    """assigns cell types to simulation entities"""
    for i, e in enumerate(sc.entity_list):
        if cell_type[i] == 0:
            continue
        e.change_type = list(fractions_dict.keys())[cell_type[i] - 1]


scan_container = ScanContainer()

scenario = setup(cytokines, cell_types_dict, [], geometry, numeric, path)
pool = scenario.parameter_pool

"""Retrieves entity types from sim container"""
naive = scenario.get_entity_type_by_name("naive")
abs = scenario.get_entity_type_by_name("abs")
sec = scenario.get_entity_type_by_name("sec")

s = 10
scan_space = np.linspace(1, 0, s)

cs = ScannableParameter(PhysicalParameter("strength", 0.5), lambda x, v: v)

sec_cs_def = ScanDefintion(cs, "clustering", scan_space, ScanType.ENTITY, entity_type=sec)
abs_cs_def = ScanDefintion(cs, "clustering", scan_space, ScanType.ENTITY, entity_type=abs)
scan_container.add_single_parameter_scan([sec_cs_def, abs_cs_def])

state_manager = StateManager.StateManager(path)
state_manager.scenario = scenario
state_manager.scan_container = scan_container
state_manager.dt = 1
state_manager.T = [0, 1]


def pre_replicat(sc, time_index, replicat_index, t, T):
    update_state(sc, replicat_index)


state_manager.pre_replicat = pre_replicat

"""Runs the ParameterScan"""
state_manager.run(ext_cache=ext_cache, model_names=["pde_model"], number_of_replicats=1)
