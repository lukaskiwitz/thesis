import random
from copy import deepcopy
from typing import List, Dict

import numpy as np
from scipy.constants import N_A

import thesis.main.BC as bc
import thesis.main.Entity as Entity
import thesis.main.FieldProblem as fp
import thesis.main.MySolver
import thesis.main.SimContainer as SC
from thesis.main.EntityType import CellType
from thesis.main.ParameterSet import ParameterSet, ParameterCollection, PhysicalParameter, MiscParameter, \
    PhysicalParameterTemplate
from thesis.main.PostProcess import get_concentration_conversion as get_cc
from thesis.main.bcFunctions import cellBC
from thesis.main.my_debug import message
import fenics as fcs
"""Sets up parameter templates. This are callable object, which return a full copy of themselves 
with a new value (set in post units). This is so that conversion information has to be specified only one."""


def makeCellListGrid(global_parameter, cytokines, xLine, yLine, zLine):
    cellList = []
    ran = random.Random()
    ran.seed(1)

    for x in xLine:
        for y in yLine:
            for z in zLine:
                # np.random.seed(1)
                max_dist = 10
                r = 0  # np.random.normal(5,2)
                if r > max_dist:
                    r = 0.01 * max_dist
                else:
                    r *= 0.01

                a = np.random.uniform(0, 2 * np.pi)

                # b = np.random.uniform(0,2*np.pi)

                pos = [
                    x + r * np.cos(a),
                    y + r * np.sin(a),
                    z
                ]
                # z += pos_ran[2]
                cell_bcs = []

                for c in cytokines:
                    cell_bcs.append(bc.Integral(cellBC, field_quantity=c["field_quantity"]))

                cell = Entity.Cell(pos,
                                   global_parameter.get_physical_parameter("rho", "rho").get_in_sim_unit(), cell_bcs)
                cell.name = "cell"
                cellList.append(cell)
    return cellList


def get_cuboid_domain(x, y, z, margin):
    p1 = [
        x[0] - margin,
        y[0] - margin,
        z[0] - margin
    ]
    p2 = [
        x[-1] + margin,
        y[-1] + margin,
        z[-1] + margin
    ]
    return p1, p2


def setup(cytokines, cell_types, boundary, geometry_dict, numeric, path, ext_cache="", cell_list = None):
    message("---------------------------------------------------------------")
    message("Setup")
    message("---------------------------------------------------------------")

    templates = get_parameter_templates(numeric["unit_length_exponent"])

    def make_grid(gd: Dict, key: str):
        if key in gd.keys():
            return np.round(np.arange(gd["margin"], gd[key], gd["distance"]), 2)
        else:
            return [0]

    x = make_grid(geometry_dict, "x_grid")
    y = make_grid(geometry_dict, "y_grid")
    z = make_grid(geometry_dict, "z_grid")

    margin = geometry_dict["margin"]

    global_parameter = make_global_parameters(cytokines, geometry_dict, numeric, templates)
    na = geometry_dict["norm_area"]
    geometry = ParameterCollection("geometry", [PhysicalParameter("norm_area", na, to_sim=1)])
    domain_parameter_set = deepcopy(global_parameter)
    domain_parameter_set.add_collection(geometry)
    domain_parameter_set.name = "domain"

    cell_type_list, fractions = make_cell_types(cell_types, cytokines, templates)
    global_parameter.add_collection(fractions)


    domain_bc = make_domain_bc(cytokines, boundary, numeric, domain_parameter_set, margin, x)

    p1, p2 = get_cuboid_domain(x,y,z, margin)


    domain = thesis.main.Entity.DomainCube(
        p1,
        p2, domain_bc, periodic = False)



    sc = SC.SimContainer(global_parameter)

    """FieldProblems"""
    for c in cytokines:

        solver = thesis.main.MySolver.MyDiffusionSolver()
        solver.timeout = 24 * 60**2

        fieldProblem = fp.FieldProblem()
        fieldProblem.field_name = c["name"]
        fieldProblem.field_quantity = c["field_quantity"]

        fieldProblem.set_solver(solver)

        if not ext_cache == "":
            fieldProblem.ext_cache = ext_cache

        fieldProblem.set_outer_domain(domain)
        sc.add_problem(fieldProblem)


    """top level path"""
    sc.path = path

    """adds cell to simulation"""

    if cell_list is None:
        for i in makeCellListGrid(global_parameter, cytokines, x, y, z):
            sc.add_entity(i)
    else:
        for i in cell_list:
            cell_bcs = []
            for c in cytokines:
                cell_bcs.append(
                    bc.Integral(cellBC, field_quantity=c["field_quantity"])
                )
            p = [i["x"],i["y"],i["z"]]

            cell = Entity.Cell(p,i["radius"], cell_bcs)
            cell.name = i["name"]
            sc.add_entity(cell)

    """adds entity types"""
    for ct in cell_type_list:
        sc.add_entity_type(ct)

    message("initializing sim container")
    """sets external path for subdomain markers"""
    if not ext_cache == "":
        sc.initialize(load_subdomain=True, file_name="mesh")
    else:
        sc.initialize()
    message("initialization complete")

    from thesis.main.ScanContainer import ScanSample
    default = deepcopy(ScanSample(global_parameter.collections, cell_type_list, {}))
    sc.default_sample = default

    return sc


def make_domain_bc(cytokines, boundary, numeric, domain_parameter_set, margin, x):
    domainBC = []
    templates = get_parameter_templates(numeric["unit_length_exponent"])

    for piece in boundary:
        for key, line in piece.items():
            if key in [c["field_quantity"] for c in cytokines]:
                outer_integral = bc.OuterIntegral(
                    cellBC, piece["expr"].format(d=x[0] - margin), field_quantity=key,#todo only works in x direction!
                    p = deepcopy(domain_parameter_set), name = piece["name"]
                )

                i = [c["field_quantity"] for c in cytokines].index(key)
                parameters = []
                for name, value in line.items():
                    if name in templates:
                        parameters.append(templates[name](value))
                    else:
                        parameters.append(MiscParameter(name, value))

                s = ParameterSet("update", [ParameterCollection(cytokines[i]["name"], parameters, field_quantity=key)])
                outer_integral.p.update(s, override=True)

                domainBC.append(outer_integral)
    if len(domainBC) == 0:
        c = cytokines[0]

        def dummy_bc(u,p,fq,area = 1):
            return fcs.Constant(0)

        outer_integral = bc.OuterIntegral(
                    dummy_bc, "true", field_quantity=c["field_quantity"],
                    p = deepcopy(domain_parameter_set), name = "box"
                )
        domainBC.append(outer_integral)


    return domainBC


def make_cell_types(cell_types, cytokines, templates) -> (List[CellType], ParameterCollection):
    fractions = ParameterCollection("fractions", [], is_global=True)
    clustering = ParameterCollection("clustering", [])
    motility = ParameterCollection("motility", [])

    cell_types_list = []

    t_threshold = templates["threshold"]
    t_Kc = templates["Kc"]
    t_amax = templates["amax"]

    t_bw = templates["bw"]
    t_cluster_strength = templates["cluster_strength"]
    t_q = templates["q"]
    t_R = templates["R"]
    t_mc = templates["mc"]

    for ct in cell_types:

        fractions.set_parameter(PhysicalParameter(ct["name"], ct["fraction"], is_global=True))

        cell_p_set = ParameterSet(ct["name"], [])
        misc = ParameterCollection("misc", [])
        cell_p_set.add_collection(clustering)
        cell_p_set.add_collection(motility)
        cell_p_set.add_collection(misc)

        if "bw" in ct.keys():
            bw = t_bw(ct["bw"])
            clustering.set_parameter(bw)
        if "cluster_strength" in ct.keys():
            s = t_cluster_strength(ct["cluster_strength"])
            clustering.set_parameter(s)
        if "mc" in ct.keys():
            mc = t_mc(ct["mc"])
            motility.set_parameter(mc)

        if "misc" in ct.keys():
            for k, v in ct["misc"].items():
                misc.set_parameter(MiscParameter(k, v))

        for c in cytokines:
            if c["field_quantity"] in ct.keys():
                ct_dict = ct[c["field_quantity"]]
                p = []

                if "R" in ct_dict.keys():
                    p.append(t_R(ct_dict["R"]))
                    del ct_dict["R"]
                if "q" in ct_dict.keys():
                    p.append(t_q(ct_dict["q"]))
                    del ct_dict["q"]
                if "ths" in ct_dict.keys():
                    p.append(t_threshold(ct_dict["ths"]))
                    del ct_dict["ths"]
                if "Kc" in ct_dict.keys():
                    p.append(t_Kc(ct_dict["Kc"]))
                    del ct_dict["Kc"]
                if "amax" in ct_dict.keys():
                    p.append(t_amax(ct_dict["amax"]))
                    del ct_dict["amax"]
                if "k_off" in ct_dict.keys():
                    p.append(t_amax(ct_dict["k_off"]))
                    del ct_dict["k_off"]
                if "k_endo" in ct_dict.keys():
                    p.append(t_amax(ct_dict["k_endo"]))
                    del ct_dict["k_endo"]

                for k, v in ct_dict.items():
                    p.append(MiscParameter(k, v))

                collection = ParameterCollection(c["name"], p)
                collection.field_quantity = c["field_quantity"]
                cell_p_set.add_collection(collection)

        cell_type = CellType(cell_p_set, ct["name"], ct["internal_solver"])
        cell_types_list.append(cell_type)
    return cell_types_list, fractions


def make_global_parameters(cytokines: List, geometry: Dict, numeric: Dict, templates) -> ParameterSet:
    t_k_on = templates["k_on"]
    t_k_off = templates["k_off"]
    t_k_endo = templates["k_endo"]
    t_D = templates["D"]
    t_kd = templates["kd"]
    t_Kc = templates["Kc"]
    t_rho = templates["rho"]
    t_q = templates["q"]
    t_R = templates["R"]
    t_amax = templates["amax"]

    global_parameter = ParameterSet("global", [])
    global_parameter.add_collection(ParameterCollection("rho", [t_rho(geometry["rho"])]))

    numeric_c = ParameterCollection("numeric", [], is_global=True)
    global_parameter.add_collection(numeric_c)

    for c in cytokines:

        collection = ParameterCollection(c["name"], [])
        collection.field_quantity = c["field_quantity"]
        p = [t_R(0), t_q(0), t_amax(0), t_Kc(0.1)]

        if "k_on" in c.keys():
            p.append(t_k_on(c["k_on"]))
        else:
            p.append(t_k_on(111.6))
        if "k_off" in c.keys():

            p.append(t_k_off(c["k_off"]))
        else:
            p.append(t_k_off(0.83))
        if "k_endo" in c.keys():
            p.append(t_k_endo(c["k_endo"]))
        else:
            p.append(t_k_endo(1.1e-3, in_sim=True))
        if "D" in c.keys():
            p.append(t_D(c["D"]))
        else:
            p.append(t_D(10))
        if "kd" in c.keys():
            p.append(t_kd(c["kd"]))
        else:
            p.append(t_kd(0.1))
        for p in p:
            collection.set_parameter(p)
        global_parameter.add_collection(collection)

    for k, v in numeric.items():
        numeric_c.set_misc_parameter(MiscParameter(k, v))

    return global_parameter


def get_parameter_templates(ule):
    templates = {
        "R": PhysicalParameterTemplate(PhysicalParameter("R", 0, to_sim=N_A ** -1 * 1e9)),
        "k_on": PhysicalParameterTemplate(PhysicalParameter("k_on", 111.6, to_sim=1e15 / 60 ** 2, is_global=True)),
        "k_off": PhysicalParameterTemplate(PhysicalParameter("k_off", 0.83, to_sim=1 / 60 ** 2, is_global=True)),
        "k_endo": PhysicalParameterTemplate(PhysicalParameter("k_endo", 1.1e-3, to_sim=1, is_global=True)),
        "q": PhysicalParameterTemplate(PhysicalParameter("q", 0, to_sim=N_A ** -1 * 1e9)),
        "D": PhysicalParameterTemplate(PhysicalParameter("D", 10, to_sim=1, is_global=True)),
        "kd": PhysicalParameterTemplate(PhysicalParameter("kd", 0.1, to_sim=1 / (60 ** 2), is_global=True)),
        "threshold": PhysicalParameterTemplate(PhysicalParameter("ths", 0.1, to_sim=1 / get_cc(ule))),
        "Kc": PhysicalParameterTemplate(PhysicalParameter("Kc", 0.01, to_sim=1 / get_cc(ule))),
        "bw": PhysicalParameterTemplate(PhysicalParameter("bw", 10, to_sim=10 ** (6 + ule))),
        "cluster_strength": PhysicalParameterTemplate(PhysicalParameter("strength", 10, to_sim=1)),
        "rho": PhysicalParameterTemplate(PhysicalParameter("rho", 0, to_sim=10 ** (-6 - ule))),
        "amax": PhysicalParameterTemplate(PhysicalParameter("amax", 0, to_sim=N_A ** -1 * 1e9)),
        "mc":PhysicalParameterTemplate(PhysicalParameter("mc",0,to_sim = 1))
    }

    return templates