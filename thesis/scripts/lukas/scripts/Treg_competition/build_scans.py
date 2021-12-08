import logging
import os
from copy import deepcopy

import numpy as np
from matplotlib import pyplot as plt

from thesis.main.MyPlotter import Plotter
from thesis.main.ParameterSet import MiscParameter, PhysicalParameter, ParameterSet, \
    ParameterCollection, PhysicalParameterTemplate, MiscParameterTemplate
from thesis.main.ScanContainer import ScanContainer, ScanSample
from thesis.main.StateManager import StateManager


def scan_1(scenario, q_scale, cs_strength=1, pSTAT_ec50=1, treg_k=860, treg_half_time=1):
    scenario = deepcopy(scenario)
    q_scale = deepcopy(q_scale)

    pool = scenario.parameter_pool
    naive = scenario.get_entity_type_by_name("naive")
    sec = scenario.get_entity_type_by_name("sec")
    treg = scenario.get_entity_type_by_name("treg")
    cs_t = PhysicalParameterTemplate(PhysicalParameter("strength", 0.5))

    pSTAT_k_t = MiscParameterTemplate(MiscParameter("pSTAT_k", 860))
    pSTAT_ec50_t = MiscParameterTemplate(MiscParameter("pSTAT_N", 1))
    half_time_t = MiscParameterTemplate(MiscParameter("pos_half_time", 1))

    scan_samples = []

    scan_name = "cs_{cs}_ec50_{N}_k_{k}_ht_{ht}".format(cs=cs_strength, N=pSTAT_ec50, k=treg_k, ht=treg_half_time)

    for q in q_scale:
        sample = ScanSample(
            [
                ParameterCollection("scan", [
                    MiscParameterTemplate(MiscParameter("value", 0, is_global=True))(q)
                ])
            ],
            [
                treg.get_updated(
                    ParameterSet("dummy", [
                        ParameterCollection("clustering", [cs_t(cs_strength)]),
                        ParameterCollection("misc",
                                            [pSTAT_k_t(treg_k), pSTAT_ec50_t(pSTAT_ec50), half_time_t(treg_half_time)])
                    ])),
                naive.get_updated(
                    ParameterSet("dummy", [
                        ParameterCollection("misc", [pSTAT_k_t(860), pSTAT_ec50_t(pSTAT_ec50)])
                    ])),
                sec.get_updated(
                    ParameterSet("dummy", [
                        ParameterCollection("misc", [pSTAT_ec50_t(pSTAT_ec50)]),
                        ParameterCollection("clustering", [cs_t(cs_strength)]),
                        ParameterCollection("IL-2", [pool.get_template("q")(q)], field_quantity="il2"),
                    ]
                                 )
                )

            ],
            {},
            scan_name=scan_name
        )
        scan_samples.append(sample)
    return scan_samples, scan_name


def scan_2(scenario, f_scale, cs_strength=1, pSTAT_ec50=1, treg_k=860, treg_half_time=1):
    scenario = deepcopy(scenario)
    f_scale = deepcopy(f_scale)

    pool = scenario.parameter_pool
    naive = scenario.get_entity_type_by_name("naive")
    sec = scenario.get_entity_type_by_name("sec")
    treg = scenario.get_entity_type_by_name("treg")
    cs_t = PhysicalParameterTemplate(PhysicalParameter("strength", 0.5))

    pSTAT_k_t = MiscParameterTemplate(MiscParameter("pSTAT_k", 860))
    pSTAT_ec50_t = MiscParameterTemplate(MiscParameter("pSTAT_N", 1))
    half_time_t = MiscParameterTemplate(MiscParameter("pos_half_time", 1))

    scan_samples = []

    scan_name = "cs_{cs}_ec50_{N}_k_{k}_ht_{ht}".format(cs=cs_strength, N=pSTAT_ec50, k=treg_k, ht=treg_half_time)

    for fsec in f_scale:
        sample = ScanSample(
            [
                ParameterCollection("scan", [
                    MiscParameterTemplate(MiscParameter("value", 0, is_global=True))(fsec)
                ]),
                ParameterCollection("fractions", [
                    PhysicalParameterTemplate(PhysicalParameter("sec", 0, is_global=True))(fsec)
                ])
            ],
            [
                treg.get_updated(
                    ParameterSet("dummy", [
                        ParameterCollection("clustering", [cs_t(cs_strength)]),
                        ParameterCollection("misc",
                                            [pSTAT_k_t(treg_k), pSTAT_ec50_t(pSTAT_ec50), half_time_t(treg_half_time)])
                    ])),
                naive.get_updated(
                    ParameterSet("dummy", [
                        ParameterCollection("misc", [pSTAT_k_t(860), pSTAT_ec50_t(pSTAT_ec50)])
                    ])),
                sec.get_updated(
                    ParameterSet("dummy", [
                        ParameterCollection("misc", [pSTAT_ec50_t(pSTAT_ec50)]),
                        ParameterCollection("clustering", [cs_t(cs_strength)]),
                    ]
                                 )
                )

            ],
            {},
            scan_name=scan_name
        )
        scan_samples.append(sample)
    return scan_samples, scan_name


def update_state(sc, replicat_index):
    from thesis.cellBehaviourUtilities.grid_clustering import make_clusters

    np.random.seed(replicat_index)
    cell_grid_positions = []
    for i, e in enumerate(sc.entity_list):
        cell_grid_positions.append(e.center)
        e.p.add_parameter_with_collection(MiscParameter("id", int(i)))

    cell_grid_positions = np.array(cell_grid_positions)
    apcs = np.array([[55.0, 55.0, 55.0]])
    # apcs = get_apc_positions(cell_grid_positions, no_apcs=5)

    fractions_dict = sc.p.get_collection("fractions").get_as_dictionary()
    del fractions_dict["naive"]

    fractions = [sc.p.get_physical_parameter(k, "fractions").get_in_sim_unit() for k in fractions_dict.keys()]

    cluster_strengths = [sc.get_entity_type_by_name(k).p.get_physical_parameter("strength", "clustering") for k in
                         fractions_dict.keys()]
    cluster_strengths = [i.get_in_sim_unit() if i is not None else 1 for i in cluster_strengths]

    cell_type = make_clusters(cell_grid_positions, apcs, fractions, cluster_strengths)

    for i, e in enumerate(sc.entity_list):
        if cell_type[i] == 0: continue
        e.change_type = list(fractions_dict.keys())[cell_type[i] - 1]

    # distribute_receptors(sc.entity_list,replicat_index=replicat_index, type_name="naive")


def compute_samples_points(scan_samples, scenario, path, time_range, ext_cache=""):
    scenario = deepcopy(scenario)
    scan_samples = deepcopy(scan_samples)
    path = deepcopy(path)
    time_range = deepcopy(time_range)
    ext_cache = deepcopy(ext_cache)

    os.makedirs(path, exist_ok=True)

    logging.basicConfig(
        filename=os.path.join(path, "sim.log"),
        level=logging.INFO,
        filemode="w",
        format='%(levelname)s::%(asctime)s %(message)s',
        datefmt='%I:%M:%S')

    scan_container = ScanContainer()
    for s in scan_samples:
        scan_container.add_sample(s)

    path = path + "/"

    scenario.marker_lookup = {"naive": 1, "sec": 2, "treg": 3}
    scenario.markers.append("IL-2_R")
    from thesis.cellBehaviourUtilities.cell_solver import kineticSolver

    scenario.internal_solvers = [kineticSolver]

    state_manager = StateManager(path)
    state_manager.scenario = scenario
    state_manager.scan_container = scan_container
    state_manager.compress_log_file = True

    state_manager.T = time_range

    def pre_replicat(sc, time_index, replicat_index, t, T):
        update_state(state_manager.sim_container, replicat_index)

    state_manager.pre_replicat = pre_replicat

    """Runs the ParameterScan"""
    state_manager.run(ext_cache=ext_cache, model_names=["pde_model"], number_of_replicats=1)


def post_process(path, time_indices, n_procs=2):
    path = path + "/"
    from thesis.main.PostProcess import PostProcessor
    pp = PostProcessor(path)
    pp.unit_length_exponent = -6

    pp.image_settings = {
        "cell_colors": "Dark2",
        "cell_color_key": "type_name",
        "round_legend_labels": 4,
        "legend_title": "",
        "dpi": 350,
    }
    pp.run_post_process(n_procs, time_indices=time_indices)


def get_new_samples_points(path_list, grad_max=10):
    path_list = deepcopy(path_list)
    #
    # with open('./path_list.txt', 'w') as f:
    #     print(path_list, file=f)

    path = path_list[-1]
    IMGPATH = os.path.join(path, "images/")
    os.makedirs(IMGPATH, exist_ok=True)

    plotter = Plotter(path_list)
    plotter.cell_df[plotter.time_key] = plotter.cell_df[plotter.time_key] / 3600
    plotter.global_df[plotter.time_key] = plotter.global_df[plotter.time_key] / 3600
    plotter.t_max = plotter.global_df[plotter.time_key].max()
    plotter.label_replacement.update({

        "D": "Diffusion coefficient",
        "sec_amax": "$a_{max}$",
        "sec_q": "IL-2 / (cell * s)",
        "f_sec": "% secretors",
        "f_abs": "% consumers",
        "abs_R": "IL-2R / cell",
        "Kc": "IL-2R EC50",
        "kd": "IL-2 decay",
        plotter.scan_index_key: "q (molec/(cell*s))",
        "Concentration": "Concentration (nM)",
        "run": "total",
        "run:scan_sample:SimContainer:run": "time series",
        "run:scan_sample:SimContainer:run:step": "timestep",
        "run:scan_sample:update_sim_container": "update sim_container",
        "run:write_element_tree": "write xml file",
        "run:scan_sample": "scan sample"

    })
    activation_threshold = 3e3
    os.makedirs(IMGPATH, exist_ok=True)
    a = 8.3
    b = np.sqrt(2) * a

    df = plotter.cell_df
    time_key = "time"
    t_max = plotter.t_max
    scan_index_key = "scan_value"
    time_index_key = "time_index"
    scan_name_key = "scan_name_scan_name"

    total_N = \
        df.groupby([scan_index_key, time_index_key, scan_name_key, "path_name"]).count().reset_index()[
            "IL-2_R"].unique()[0]

    cf = lambda df: df.loc[(df.type_name.isin(["naive"])) & (df[time_key] == t_max)]

    df["activated"] = df["IL-2_R"] > activation_threshold
    df = cf(df)

    def cf2(df):
        return len(df.loc[df.activated == True])

    df_total = df.groupby(
        [plotter.scan_index_key, plotter.time_index_key, plotter.time_key, plotter.scan_name_key]).count().reset_index()
    df = df.groupby([scan_index_key, time_index_key, scan_name_key, "path_name"]).apply(cf2).reset_index()

    df["N"] = df[0] / df_total["IL-2_R"]

    plt.figure()

    x = np.array(df[scan_index_key])
    y = np.array(df["N"])

    plt.plot(x, y, "-x")
    y = np.abs(np.diff(y))
    plt.step(x[:-1], y)

    p = np.argwhere((y > grad_max) & (y > 1 / total_N))[:, 0]

    new_points = np.array([(x[i] + x[i + 1]) / 2 for i in p])

    new_scan_scale = np.sort(np.concatenate([x, new_points]))

    plt.scatter(x, [0] * len(x), color="gray")
    plt.scatter(new_points, [0.05] * len(new_points), color="red")

    plt.plot(plt.gca().get_xlim(), [grad_max, grad_max])
    plt.tight_layout()
    plt.savefig(IMGPATH + "subdiv_plot.pdf")

    # plt.show()
    # plt.close(plt.gcf())
    return new_points
