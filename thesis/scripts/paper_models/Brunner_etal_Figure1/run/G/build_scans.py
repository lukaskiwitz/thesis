import glob
import numbers
import os
from copy import deepcopy
from enum import Enum
from typing import List

import numpy as np
from matplotlib import pyplot as plt

from thesis.main.MyPlotter import Plotter
from thesis.main.ParameterSet import MiscParameter, PhysicalParameter, ParameterSet, \
    ParameterCollection, PhysicalParameterTemplate, MiscParameterTemplate
from thesis.main.ScanContainer import ScanContainer, ScanSample
from thesis.main.StateManager import StateManager
from thesis.scripts.paper_models.utilities.states import updateState


def chunks(lst, n):
    """Yield successive n-sized chunks from lst.
    https://stackoverflow.com/questions/312443/how-do-you-split-a-list-into-evenly-sized-chunks
    """
    for i in range(0, len(lst), n):
        yield lst[i:i + n]


def get_existing_path_list(path, start=0):
    if start is not None:
        assert start >= 0
    import re
    path_list = []
    glob_list = glob.glob(os.path.join(path, "refine*"))
    refine_paths = [re.search(r"refine_(?P<n>\d)+_\d+$", rp) for rp in glob_list]
    previous_iterations = np.unique([int(rp.groupdict()["n"]) for rp in refine_paths])
    if start is not None:
        previous_iterations = previous_iterations[np.where(previous_iterations < start)]

    for i in previous_iterations:
        path_list.append([rp.string for rp in refine_paths if int(rp.groupdict()["n"]) == i])

    if len(path_list) == 0:
        return path_list, 0
    if start is not None and start > previous_iterations.max() + 1:
        print("Largest found refine iteration is {pi}".format(pi=previous_iterations.max()))

    return path_list, previous_iterations.max() + 1


class ScanAxis(Enum):
    ftreg = "ftreg"
    fsec = "fsec"
    fth = "fth"
    treg_k = "treg_k"
    treg_half_time = "treg_half_time"
    Tsec_q = "Tsec_q"
    pSTAT_ec50 = "pSTAT_ec50"
    sigma = "sigma"
    treg_cs_strength = "treg_cs_strength"
    Tsec_cs_strength = "Tsec_cs_strength"
    gamma = "gamma"
    v = "v"
    R_per_cell = "R_per_cell"


    def __str__(self):
        return self.value


def get_conditioned_scan_sample(pool, Th, Tsec, treg, parameter_dict, scan_value_dict):
    cs_t = PhysicalParameterTemplate(PhysicalParameter("strength", 0.5))
    pSTAT_k_t = MiscParameterTemplate(MiscParameter("pSTAT_k", 860))
    pSTAT_ec50_t = MiscParameterTemplate(MiscParameter("pSTAT_N", 1))
    half_time_t = MiscParameterTemplate(MiscParameter("pos_half_time", 1))
    sigma_t = MiscParameterTemplate(MiscParameter("sigma", 1))
    gamma_t = MiscParameterTemplate(MiscParameter("gamma", 1))
    v_t = MiscParameterTemplate(MiscParameter("v", 1))
    R_per_cell_t = MiscParameterTemplate(MiscParameter("R_per_cell", 1))

    value_dict = deepcopy(parameter_dict)
    value_dict.update(scan_value_dict)

    fractions_collection = ParameterCollection("fractions", [])
    my_collection = ParameterCollection("my_c", [])

    if "ftreg" in value_dict.keys():
        ftreg = value_dict["ftreg"]
        fractions_collection.set_physical_parameter(
            PhysicalParameterTemplate(PhysicalParameter("Treg", 0, is_global=True))(ftreg)
        )

    if "fsec" in value_dict.keys():
        fsec = value_dict["fsec"]
        fractions_collection.set_physical_parameter(
            PhysicalParameterTemplate(PhysicalParameter("TTsec", 0, is_global=True))(fsec)
        )

    if "fth" in value_dict.keys():
        fth = value_dict["fth"]
        fractions_collection.set_physical_parameter(
            PhysicalParameterTemplate(PhysicalParameter("Th", 1, is_global=True))(fth)
        )

    treg_parameter_set = ParameterSet("dummy", [ParameterCollection("clustering", []), ParameterCollection("misc", [])])
    if "treg_k" in value_dict.keys():
        treg_k = value_dict["treg_k"]
        treg_parameter_set.get_collection("misc").set_misc_parameter(pSTAT_k_t(treg_k))

    if "treg_half_time" in value_dict.keys():
        treg_half_time = value_dict["treg_half_time"]
        treg_parameter_set.get_collection("misc").set_misc_parameter(half_time_t(treg_half_time))

    Th_parameter_set = ParameterSet("dummy",
                                       [ParameterCollection("clustering", []), ParameterCollection("misc", [])])
    Tsec_parameter_set = ParameterSet("dummy", [ParameterCollection("clustering", []), ParameterCollection("misc", [])])

    if "Tsec_q" in value_dict.keys():
        Tsec_q = value_dict["Tsec_q"]
        Tsec_parameter_set.update(ParameterCollection("IL-2", [pool.get_template("q")(Tsec_q)], field_quantity="il2"))

    if "pSTAT_ec50" in value_dict.keys():
        pSTAT_ec50 = value_dict["pSTAT_ec50"]
        Th_parameter_set.get_collection("misc").set_misc_parameter(pSTAT_ec50_t(pSTAT_ec50))
        treg_parameter_set.get_collection("misc").set_misc_parameter(pSTAT_ec50_t(pSTAT_ec50))

    if "sigma" in value_dict.keys():
        sigma = value_dict["sigma"]
        Th_parameter_set.get_collection("misc").set_misc_parameter(sigma_t(sigma))
        treg_parameter_set.get_collection("misc").set_misc_parameter(sigma_t(sigma))

    if "treg_cs_strength" in value_dict.keys():
        cs_strength = value_dict["treg_cs_strength"]
        treg_parameter_set.get_collection("clustering").set_physical_parameter(cs_t(cs_strength))
        # Tsec_parameter_set.get_collection("clustering").set_physical_parameter(cs_t(cs_strength))
    if "Tsec_cs_strength" in value_dict.keys():
        cs_strength = value_dict["Tsec_cs_strength"]
        Tsec_parameter_set.get_collection("clustering").set_physical_parameter(cs_t(cs_strength))
    if "gamma" in value_dict.keys():
        gamma = value_dict["gamma"]
        Th_parameter_set.get_collection("misc").set_misc_parameter(gamma_t(gamma))
        Tsec_parameter_set.get_collection("misc").set_misc_parameter(gamma_t(gamma))

    if "v" in value_dict.keys():
        v = value_dict["v"]
        Th_parameter_set.get_collection("misc").set_misc_parameter(v_t(v))
        my_collection.set_physical_parameter(
            PhysicalParameterTemplate(PhysicalParameter("v", 0, is_global=True))(v)
        )
    if "R_per_cell" in value_dict.keys():
        R_per_cell = value_dict["R_per_cell"]
        Th_parameter_set.get_collection("misc").set_misc_parameter(R_per_cell_t(R_per_cell))
        my_collection.set_physical_parameter(
            PhysicalParameterTemplate(PhysicalParameter("R_per_cell", 0, is_global=True))(R_per_cell)
        )

    scan_name = get_scan_name(scan_value_dict, **parameter_dict)
    sample = ScanSample(
        [fractions_collection, my_collection],
        [
            treg.get_updated(treg_parameter_set),
            Th.get_updated(Th_parameter_set),
            Tsec.get_updated(Tsec_parameter_set)
        ],
        {},
        scan_name=scan_name
    )
    return sample


def get_scan_name(scan_value_dict, **kwargs):
    scan_name = ""
    for k, v in scan_value_dict.items():
        scan_name = scan_name + "{k}_scan_".format(k=str(k), v=str(v))
    scan_name = scan_name[0:-1]
    for k, v in kwargs.items():
        scan_name = scan_name + "{k}_{v}_".format(k=str(k), v=str(v))

    return scan_name[0:-1]


def build_scan_function(scenario, scan_value_name, scan_scale, **kwargs):
    scenario = deepcopy(scenario)
    scan_scale = deepcopy(scan_scale)
    pool = scenario.parameter_pool
    Th = scenario.get_entity_type_by_name("Th")
    Tsec = scenario.get_entity_type_by_name("Tsec")
    treg = scenario.get_entity_type_by_name("Treg")

    scan_samples = []
    for scan_value in scan_scale:
        sample = get_conditioned_scan_sample(pool, Th, Tsec, treg, kwargs, {scan_value_name: scan_value})

        sample.p.update(
            ParameterSet("update",
                         [
                             ParameterCollection("scan", [
                                 MiscParameterTemplate(MiscParameter("value", 0, is_global=True))(scan_value)]),
                         ]), overwrite=True
        )

        scan_samples.append(sample)
    if len(scan_samples) == 0:
        scan_name = "empty"
    else:
        scan_name = scan_samples[0].p.get_misc_parameter("scan_name", "scan_name").get_in_sim_unit()
    return scan_samples, scan_name

def get_draws(sc, offset_ids, excluded_ids, fraction, at_least_one = False):
    possible_draws = np.setdiff1d(range(len(sc.entity_list)), offset_ids)
    free_draws = np.setdiff1d(possible_draws, excluded_ids)

    if fraction == 0:
        adj_fraction = 0
    else:
        adj_fraction = 1 / (len(free_draws) / (fraction * len(sc.entity_list)))
    if np.abs(1.0 - adj_fraction) < 0.009:
        adj_fraction = 1

    amount_of_draws = int(round(len(free_draws) * adj_fraction))
    if at_least_one == True:
        if amount_of_draws == 0:
            amount_of_draws = 1
            print("set to 1")
    try:
        draws = np.random.choice(free_draws, amount_of_draws, replace=False)
    except ValueError:
        draws = np.random.choice(free_draws, amount_of_draws - 1, replace=False)
    return draws

def update_state(sc, pool, replicat_index):
    offset_ids = []
    varying_Tsec_fraction = False
    try:
        Tsec_fraction = sc.entity_list[0].p.get_physical_parameter("Tsec_fraction", "IL-2").get_in_sim_unit()
        varying_Tsec_fraction = True
    except:
        Tsec_fraction = sc.entity_list[0].p.get_as_dictionary()["fractions_Tsec"]

    try:
        Treg_fraction = sc.entity_list[0].p.get_as_dictionary()["fractions_Treg"]
        Treg_fraction = sc.entity_list[0].p.get_physical_parameter("Treg_fraction", "IL-2").get_in_sim_unit()
    except AttributeError:
        Treg_fraction = sc.entity_list[0].p.get_as_dictionary()["fractions_Treg"]
    except KeyError:
        Treg_fraction = 0

    Tsec_draws = get_draws(sc, offset_ids=offset_ids, excluded_ids=[], fraction=Tsec_fraction,
                                at_least_one=True)

    if varying_Tsec_fraction == True or Treg_fraction != 0:
        Th_fraction = 1 - Tsec_fraction - Treg_fraction
    elif Treg_fraction == 0:
        Th_fraction = sc.entity_list[0].p.get_as_dictionary()["fractions_Th"]

    Treg_draws = get_draws(sc, offset_ids=offset_ids, excluded_ids=list(Tsec_draws), fraction=Treg_fraction)

    Th_draws = get_draws(sc, offset_ids=offset_ids, excluded_ids=list(Tsec_draws) + list(Treg_draws),
                              fraction=Th_fraction)

    no_of_Th = 0
    no_of_Treg = 0

    for i, e in enumerate(sc.entity_list):
        e.change_type = "Th"
        if i in Tsec_draws:
            e.change_type = "Tsec"
        if i in Th_draws:
            e.change_type = "Th"
            no_of_Th += 1
        if i in offset_ids:
            e.change_type = "Th"
            Th_draws = np.append(Th_draws, i)
            no_of_Th += 1
        if i in Treg_draws:
            e.change_type = "Treg"
            no_of_Treg += 1
        e.p.add_parameter_with_collection(MiscParameter("id", int(i)))
    sc.apply_type_changes(replicat_index)

    print("Number of secreting cells: " + str(len(Tsec_draws)))
    print("Number of Ths: " + str(no_of_Th))
    print("Number of Tregs: " + str(no_of_Treg))
    if len(Tsec_draws) > 0:
        try:
            global_q = sc.entity_list[Tsec_draws[0]].p.get_physical_parameter("global_q", "IL-2").get_in_sim_unit()
        except:
            global_q = sc.entity_list[Tsec_draws[0]].p.get_as_dictionary()["IL-2_global_q"]

        q = sc.entity_list[Tsec_draws[0]].p.get_physical_parameter("q", "IL-2").get_in_post_unit()
        if global_q == True:
            q_il2_sum = len(sc.entity_list) * 0.1 * q * N_A ** -1 * 1e9
            print("Cells q: " + str(q_il2_sum / len(Tsec_draws) / (N_A ** -1 * 1e9)))
        else:
            for draw in Tsec_draws:
                try:
                    half_time = sc.entity_list[draw].p.get_misc_parameter("sec_start", "misc").get_in_post_unit()
                    beta = half_time / np.log(2)
                    sc.entity_list[draw].p.get_misc_parameter("sec_start", "misc").set_in_post_unit(
                        float(np.random.exponential(beta, 1)))
                    sc.entity_list[draw].p.get_misc_parameter("tmp_q", "misc").set_in_post_unit(
                        sc.entity_list[draw].p.get_physical_parameter("q", "IL-2").get_in_post_unit())
                    if beta != 0:
                        sc.entity_list[draw].p.get_physical_parameter("q", "IL-2").set_in_post_unit(0)
                except AttributeError:
                    # print("sec_start or tmp_q not set, skipping")
                    pass
            print("Cells q: " + str(q))
            q_il2_sum = None
    else:
        q_il2_sum = None

    if len(Th_draws) == 0:
        Th_draws = np.setdiff1d(range(len(sc.entity_list)), Tsec_draws)
        assert len(Th_draws) == no_of_Th
        print("Th draws was zero, set to diff. between all cells and Tsec draws")

    my_c = sc.p.get_collection("my_c")
    v = my_c.get_as_dictionary()["v"]
    R_per_cell = my_c.get_as_dictionary()["R_per_cell"]

    n_entities = sc.get_number_of_entities()
    Tsec_n = n_entities["Tsec"] if "Tsec" in n_entities.keys() else 0
    Th_n = n_entities["Th"] if "Th" in n_entities.keys() else 0

    n = Tsec_n + Th_n
    total_R = R_per_cell * n

    Tsec_R = (v * total_R / Tsec_n) if Tsec_n > 0 else 0
    if Tsec_R == 0:
        Tsec_R = 1
    Th_R = (total_R - Tsec_R * Tsec_n) / Th_n if Th_n > 0 else 0

    Th = sc.get_entity_type_by_name("Th")
    Tsec = sc.get_entity_type_by_name("Tsec")
    R_t = pool.get_template("R")
    sc.add_entity_type(Th.get_updated(ParameterCollection("IL-2", [R_t(Th_R)], field_quantity="il2")))
    sc.add_entity_type(Tsec.get_updated(ParameterCollection("IL-2", [R_t(Tsec_R)], field_quantity="il2")))
    sc.reapply_entity_types(replicat_index)


def compute_samples_points(scan_samples, scenario, path, time_range, ext_cache=""):
    scenario = deepcopy(scenario)
    scan_samples = deepcopy(scan_samples)
    path = deepcopy(path)
    time_range = deepcopy(time_range)
    ext_cache = deepcopy(ext_cache)


    scan_container = ScanContainer()
    for s in scan_samples:
        scan_container.add_sample(s)

    path = path + "/"

    scenario.marker_lookup = {"Th": 1, "Tsec": 2, "Treg": 3}
    scenario.markers += ["IL-2_R", "IL-2_pSTAT5", "type_name", "IL-2_surf_c"]
    from thesis.cellBehaviourUtilities.cell_solver import kineticSolver

    scenario.internal_solvers = [kineticSolver]

    state_manager = StateManager(path)
    state_manager.scenario = scenario
    state_manager.scan_container = scan_container
    state_manager.compress_log_file = True

    state_manager.T = time_range
    kineticSolver.T = time_range

    # uS = updateState(0, 0, None, scenario.parameter_pool, [], [], [], offset=0)

    def pre_replicat(sc, time_index, replicat_index, t, T):
        update_state(state_manager.sim_container, scenario.parameter_pool, replicat_index)
        # uS.step(state_manager.sim_container)

    state_manager.pre_replicat = pre_replicat

    """Runs the ParameterScan"""
    state_manager.run(ext_cache=ext_cache, model_names=["pde_model"], number_of_replicats=30)


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
    pp.computations = []
    pp.run_post_process(n_procs, time_indices=time_indices)


def get_new_samples_points(path_list, refine_step, grad_max=10):
    path_list = deepcopy(path_list)
    path = path_list[-1]

    plotter = Plotter(path_list)
    plotter.cell_df[plotter.time_key] = plotter.cell_df[plotter.time_key] / 3600
    plotter.global_df[plotter.time_key] = plotter.global_df[plotter.time_key] / 3600
    plotter.t_max = plotter.global_df[plotter.time_key].max()
    df = plotter.cell_df

    def objective_function(df, fold_change_threshold):

        cf = lambda df: df.loc[(df.type_name.isin(["Th"])) & (
            (df[plotter.time_key] == plotter.t_max)|
            (df[plotter.time_index_key] == 0)
        )]

        df = cf(df)

        df_t = df.loc[df[plotter.time_key] == plotter.t_max]

        df_0 = df.loc[df[plotter.time_index_key] == 0]
        df_0 = df_0.set_index(df_t.index)

        fold_change = df_t["IL-2_R"]/df_0["IL-2_R"]

        df = df.loc[ (df[plotter.time_key] == plotter.t_max)]
        df["activated"] = (fold_change <= 1 / (0.5*df_t["misc_gamma"])) | (fold_change >= 0.5*df_t["misc_gamma"])

        def cf2(df):
            return len(df.loc[df.activated == True])

        df_total = df.groupby(
            [plotter.scan_index_key, plotter.time_index_key, plotter.time_key,
             plotter.scan_name_key]).count().reset_index()
        df = df.groupby([plotter.scan_index_key, plotter.time_index_key, plotter.scan_name_key, "path_name"]).apply(
            cf2).reset_index()
        df["N"] = df[0] / df_total["IL-2_R"]
        x = np.array(df[plotter.scan_index_key])
        y = np.array(df["N"])

        return x, y

    def minimum_step(df):
        total_N = get_total_number_of_cells(df)
        return 1 / total_N

    def get_total_number_of_cells(df):
        return df.groupby(
            [plotter.scan_index_key, plotter.time_index_key, plotter.scan_name_key, "path_name"]).count().reset_index()[
            "IL-2_R"].unique()[0]

    x, y = objective_function(df, 2)
    plt.figure()
    plt.plot(x, y, "-x")
    y = np.abs(np.diff(y))
    plt.step(x[:-1], y)

    min_step = minimum_step(df)
    p = np.argwhere((y > grad_max) & (y > min_step))[:, 0]
    new_points = np.array([(x[i] + x[i + 1]) / 2 for i in p])

    plt.scatter(x, [0] * len(x), color="gray")
    plt.scatter(new_points, [0.05] * len(new_points), color="red")

    plt.plot(plt.gca().get_xlim(), [grad_max, grad_max])
    plt.tight_layout()

    IMGPATH = os.path.join(path, "../images_{i}/".format(i=refine_step))
    os.makedirs(IMGPATH, exist_ok=True)

    plt.savefig(IMGPATH + "subdiv_plot.pdf")
    return new_points


def extract_scan_axis_from_dictionary(scan_kwargs):
    scan_kwargs = deepcopy(scan_kwargs)

    scan_value_name = None
    scan_scale = None

    for k, v in scan_kwargs.items():
        if isinstance(v, List):
            scan_scale = np.linspace(v[0], v[1], int(v[2]))
            scan_value_name = k
            break

    assert scan_value_name is not None and scan_scale is not None
    del scan_kwargs[scan_value_name]
    for k, v in scan_kwargs.items():
        assert isinstance(v, numbers.Number)

    return str(scan_value_name), scan_scale, {str(k): v for k, v in scan_kwargs.items()}
