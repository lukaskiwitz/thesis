import logging
import random
from copy import deepcopy
from typing import List, Dict

import fenics as fcs
import numpy as np
from scipy.constants import N_A

from thesis.main.BC import OuterIntegral
from thesis.main.EntityType import CellType
from thesis.main.MyDomainTemplate import MyBoundingBoxTemplate
from thesis.main.MyEntityLocator import MyCellGridLocator, MyRandomCellLocator, MyBridsonCellLocator
from thesis.main.MyFieldTemplate import MyCytokineTemplate, MyMeanCytokineTemplate
from thesis.main.MyGlobalModel import MyPDEModel, MyODEModel
from thesis.main.MyInteractionTemplate import FieldInteractionType, MyFieldInteractionTemplate
from thesis.main.MyParameterPool import MyParameterPool
from thesis.main.MyScenario import MyScenario
from thesis.main.ParameterSet import ParameterSet, PhysicalParameterTemplate, PhysicalParameter, ParameterCollection, \
    MiscParameter
from thesis.main.PostProcessUtil import get_concentration_conversion as get_cc
from thesis.main.bcFunctions import cellBC

module_logger = logging.getLogger(__name__)

ule = -6

templates = {
    "R": PhysicalParameterTemplate(PhysicalParameter("R", 0, to_sim=N_A ** -1 * 1e9)),
    "k_on": PhysicalParameterTemplate(PhysicalParameter("k_on", 111.6, to_sim=1e15 / 60 ** 2, is_global=True)),
    "k_off": PhysicalParameterTemplate(PhysicalParameter("k_off", 0.83, to_sim=1 / 60 ** 2, is_global=True)),
    "k_endo": PhysicalParameterTemplate(PhysicalParameter("k_endo", 1.1e-3, to_sim=1, is_global=True)),
    "q": PhysicalParameterTemplate(PhysicalParameter("q", 0, to_sim=N_A ** -1 * 1e9)),
    "alpha": PhysicalParameterTemplate(PhysicalParameter("alpha", 0, to_sim=N_A ** -1 * 1e9)),
    "D": PhysicalParameterTemplate(PhysicalParameter("D", 10, to_sim=1, is_global=True)),
    "kd": PhysicalParameterTemplate(PhysicalParameter("kd", 0.1, to_sim=1 / (60 ** 2), is_global=True)),
    "threshold": PhysicalParameterTemplate(PhysicalParameter("ths", 0.1, to_sim=1 / get_cc(ule))),
    "Kc": PhysicalParameterTemplate(PhysicalParameter("Kc", 0.01, to_sim=1 / get_cc(ule))),
    "bw": PhysicalParameterTemplate(PhysicalParameter("bw", 10, to_sim=10 ** (6 + ule))),
    "cluster_strength": PhysicalParameterTemplate(PhysicalParameter("strength", 0, to_sim=1)),
    "rho": PhysicalParameterTemplate(PhysicalParameter("rho", 0, to_sim=10 ** (-6 - ule))),
    "amax": PhysicalParameterTemplate(PhysicalParameter("amax", 0, to_sim=N_A ** -1 * 1e9)),
    "mc": PhysicalParameterTemplate(PhysicalParameter("mc", 0, to_sim=1)),

    "gamma": PhysicalParameterTemplate(PhysicalParameter("gamma", 0.1, to_sim=1, is_global=True)),
    "R_start": PhysicalParameterTemplate(PhysicalParameter("R_start", 20000, to_sim=N_A ** -1 * 1e9)),
    "pSTAT5": PhysicalParameterTemplate(PhysicalParameter("pSTAT5", 0, to_sim=1)),
    "EC50": PhysicalParameterTemplate(PhysicalParameter("EC50", 0, to_sim=1 / get_cc(ule))),
    "global_q": PhysicalParameterTemplate(PhysicalParameter("global_q", True, to_sim=1, is_global=True)),
    "KD": PhysicalParameterTemplate(PhysicalParameter("KD", 7.437e-3, to_sim=1 / get_cc(ule))),  # post = nM
    "nu": PhysicalParameterTemplate(PhysicalParameter("nu", 1, to_sim=1 / 60 ** 2, is_global=True)),  # receptors/s
    "eta": PhysicalParameterTemplate(PhysicalParameter("eta", 1, to_sim=1 / 60 ** 2, is_global=True)),
}


def setup(cytokines, cell_types, boundary, geometry_dict, numeric, custom_pool=None) -> MyScenario:
    """args need to be copied for thread saftey"""
    cytokines = deepcopy(cytokines)
    cell_types = deepcopy(cell_types)
    boundary = deepcopy(boundary)
    geometry_dict = deepcopy(geometry_dict)
    numeric = deepcopy(numeric)
    custom_pool = deepcopy(custom_pool)

    parameter_pool = get_standard_pool()

    if isinstance(custom_pool, MyParameterPool):
        parameter_pool.join(custom_pool, overwrite=False)

    numeric = ParameterCollection("numeric", [MiscParameter(k, v, is_global=True) for k, v in numeric.items()])
    cell_types, fractions = _make_cell_types(cell_types, cytokines, parameter_pool)

    for ct in cell_types:
        ct.p.add_collection(ParameterCollection("rho", [parameter_pool.get_template("rho")(5)]))

    domain_template = MyBoundingBoxTemplate()

    pde_model = MyPDEModel("pde_model")
    ode_model = MyODEModel("ode_model")

    pde_model.domain_template = domain_template

    for c in cytokines:
        cytokine_template = MyCytokineTemplate()
        cytokine_template.name = c["name"]
        cytokine_template.field_quantity = c["field_quantity"]
        cytokine_template.collection = parameter_pool.get_as_collection(c, name=c["name"],
                                                                        field_quantity=c["field_quantity"])

        pde_model.add_field_template(cytokine_template)

        mean_cytokine_template = MyMeanCytokineTemplate()

        mean_cytokine_template.name = c["name"]
        mean_cytokine_template.field_quantity = c["field_quantity"]
        mean_cytokine_template.collection = parameter_pool.get_as_collection(c, name=c["name"],
                                                                             field_quantity=c["field_quantity"])
        ode_model.add_field_template(mean_cytokine_template)

    na = geometry_dict["norm_area"]
    geometry = ParameterCollection("geometry", [MiscParameter(k, v) for k, v in geometry_dict.items()])
    geometry.set_parameter(PhysicalParameter("norm_area", na, to_sim=1), overwrite=True)

    domain_parameter_set = ParameterSet("domain", [])
    domain_parameter_set.add_collection(geometry)

    domain_template.bc_list = _make_domain_bc(cytokines, boundary, numeric, domain_parameter_set, parameter_pool)

    randomize = geometry.get_misc_parameter("randomize")
    if randomize is not None and randomize.get_in_post_unit():
        if randomize.get_in_post_unit() == "random_walk":
            locator = MyRandomCellLocator()
        elif randomize.get_in_post_unit() == "bridson":
            locator = MyBridsonCellLocator()

    else:
        locator = MyCellGridLocator()

    scenario = MyScenario(parameter_pool)

    scenario.global_models = [pde_model, ode_model]

    scenario.entity_types = cell_types
    scenario.entity_locators = [locator]

    scenario.global_parameters.update(geometry)
    scenario.global_parameters.update(numeric)
    scenario.global_parameters.update(fractions)

    return scenario


def get_standard_pool():
    parameter_pool = MyParameterPool()
    for i, t in templates.items():
        parameter_pool.add_template(t)
    return parameter_pool


def _make_cell_types(cell_types, cytokines, parameter_pool) -> (List[CellType], ParameterCollection):
    fractions = ParameterCollection("fractions", [], is_global=True)
    cell_types_list = []

    for ct in cell_types:

        fractions.set_parameter(PhysicalParameter(ct["name"], ct["fraction"], is_global=True))
        cell_p_set = ParameterSet(ct["name"], [])
        interactions = []

        for c in cytokines:
            if c["field_quantity"] in ct.keys():

                ct_dict = deepcopy(c)
                ct_dict.update(ct[c["field_quantity"]])
                del ct_dict["name"]
                del ct_dict["field_quantity"]
                p = []

                for k, v in ct_dict.items():
                    if parameter_pool.get_template(k) is not None:
                        p.append(parameter_pool.get_template(k)(v))
                    else:
                        p.append(MiscParameter(k, v))
                interactions.append(MyFieldInteractionTemplate(c["field_quantity"], FieldInteractionType.INTEGRAL))
                del ct[c["field_quantity"]]

                collection = ParameterCollection(c["name"], p)
                collection.field_quantity = c["field_quantity"]
                cell_p_set.add_collection(collection)

        for k, v in ct.items():
            if isinstance(v, Dict):
                plist = []
                for kk, vv in v.items():
                    if parameter_pool.get_template(kk) is not None:
                        plist.append(parameter_pool.get_template(kk)(vv))
                    else:
                        plist.append(MiscParameter(kk, vv))
                col = ParameterCollection(k, plist)
                cell_p_set.add_collection(col)

            else:
                if parameter_pool.get_template(k) is not None:
                    cell_p_set.add_parameter_with_collection(parameter_pool.get_template(k)(v))
                else:
                    cell_p_set.add_parameter_with_collection(MiscParameter(k, v))

        cell_type = CellType(cell_p_set, ct["name"], ct["internal_solver"])
        cell_type.interactions = interactions
        cell_types_list.append(cell_type)

    return cell_types_list, fractions


def _make_domain_bc(cytokines, boundary, numeric, domain_parameter_set, parameter_pool):
    domainBC = []

    for piece in boundary:
        for key, piece_line in piece.items():
            if key in [c["field_quantity"] for c in cytokines]:
                expr = piece["expr"]
                outer_integral = OuterIntegral(cellBC, expr, p=deepcopy(domain_parameter_set),
                                               field_quantity=key, name=piece["name"])

                i = [c["field_quantity"] for c in cytokines].index(key)
                parameters = []
                line = {c["field_quantity"]: c for c in deepcopy(cytokines)}[key]
                del line["name"]
                del line["field_quantity"]

                line.update(piece_line)
                for name, value in line.items():
                    if parameter_pool.get_template(name) is not None:
                        parameters.append(parameter_pool.get_template(name)(value))
                    else:
                        parameters.append(MiscParameter(name, value))

                s = ParameterSet("update", [ParameterCollection(cytokines[i]["name"], parameters, field_quantity=key)])
                outer_integral.p.update(s, overwrite=True)
                domainBC.append(outer_integral)
    if len(domainBC) == 0:
        c = cytokines[0]

        def dummy_bc(u, p, fq, area=1):
            return fcs.Constant(0)

        outer_integral = OuterIntegral(dummy_bc, "true", p=domain_parameter_set, field_quantity=c["field_quantity"],
                                       name="box")
        domainBC.append(outer_integral)

    return domainBC


def assign_fractions(sc, t):
    """sets cell types according to the values given in fractions.
        The pseudo random seed depends on t, so that cell placement is repeatable. """

    ran = random.Random()
    # ran.seed(t)

    for i, e in enumerate(sc.entity_list):

        fractions = sc.p.get_collection("fractions")
        e.change_type = fractions.parameters[0].name

        draw = ran.random()
        s = 0
        for f in fractions.parameters[1:]:
            s = s + f.get_in_sim_unit()
            if draw < s:
                e.change_type = f.name
                break
        e.p.add_parameter_with_collection(MiscParameter("id", int(i)))


def distribute_receptors(entity_list, replicat_index, type_name, var=1):
    R = np.unique(
        [e.p.get_physical_parameter("R", "IL-2").get_in_post_unit() for e in entity_list if e.type_name == type_name])
    # R_start = np.unique([e.p.get_physical_parameter("R_start", "R_start").get_in_post_unit() for e in entity_list if
    #                      e.type_name == type_name])

    if len(R) == 0: return None
    assert len(R) == 1

    E = R[0]
    if E == 0:
        return None
    # np.random.seed(replicat_index)

    tmp_sigma = np.sqrt(np.log((var * E) ** 2 / E ** 2 + 1))
    mean = np.log(E) - 1 / 2 * tmp_sigma ** 2
    for e in entity_list:
        if e.type_name == type_name:
            R_draw = np.random.lognormal(mean, tmp_sigma)
            e.p.get_physical_parameter("R", "IL-2").set_in_post_unit(R_draw)
            e.p.get_physical_parameter("R_start", "R_start").set_in_post_unit(R_draw)
            e.p.get_misc_parameter("R_start_pos", "misc").set_in_post_unit(R_draw)
            e.p.get_misc_parameter("R_start_neg", "misc").set_in_post_unit(R_draw)
