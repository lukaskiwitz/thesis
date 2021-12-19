import logging

import numpy as np

from thesis.main.ParameterSet import MiscParameter
from thesis.main.my_debug import message

module_logger = logging.getLogger(__name__)


def update_state(sc, replicat_index, parameter_pool, apcs=np.array([[140, 140, 140]])):
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
    # apcs = get_apc_positions(cell_grid_positions, no_apcs=5)

    """extract parameters for clustering for simulation parameter set(s)"""
    fractions_dict = sc.p.get_collection("fractions").get_as_dictionary()
    del fractions_dict["Tnaive"]  # does not cluster
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




    sc.apply_type_changes(replicat_index)

    type_numbers = [i + 1 for i,_ in enumerate(fractions_dict.keys())]
    type_names = [x for _, x in enumerate(fractions_dict.keys())]
    type_dict = {}
    for n,name in enumerate(type_names):
        type_dict[name] = type_numbers[n]

    Tsec_draws = np.where(cell_type == type_dict["Tsec"])[0]
    Th_draws = np.where(cell_type == type_dict["Th"])[0]
    Treg_draws = np.where(cell_type == type_dict["Treg"])[0]

    from thesis.scripts.paper_models.utilities.states import updateState as updateStatePatrick
    updateStatePatrick.parameter_pool = parameter_pool
    updateStatePatrick.set_R_lognorm_parameters(updateStatePatrick, sc,
                                                Tsec_draws=Tsec_draws,
                                                Th_draws=Th_draws,
                                                Treg_draws=Treg_draws,
                                                q_il2_sum=None
                                                )

    message("Number of secreting cells: " + str(len(Tsec_draws)), module_logger)
    message("Number of Ths: " + str(len(Th_draws)), module_logger)
    message("Number of Tregs: " + str(len(Treg_draws)), module_logger)
