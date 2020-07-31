import random
from copy import deepcopy
from typing import List, Dict

import numpy as np
from scipy.constants import N_A

import BC as bc
import Entity
import FieldProblem as fp
import MySolver
import SimContainer as SC
from EntityType import CellType
from ParameterSet import ParameterSet, ParameterCollection, PhysicalParameter, MiscParameter, PhysicalParameterTemplate
from bcFunctions import cellBC
from my_debug import message

"""Sets up parameter templates. This are callable object, which return a full copy of themselves 
with a new value (set in post units). This is so that conversion information has to be specified only one."""
t_R = PhysicalParameterTemplate(PhysicalParameter("R", 20000, to_sim=N_A ** -1 * 1e9))
t_k_on = PhysicalParameterTemplate(PhysicalParameter("k_on", 540, to_sim=1e15 / 60 ** 2, is_global=True))
t_q = PhysicalParameterTemplate(PhysicalParameter("q", 0, to_sim=N_A ** -1 * 1e9))
t_D = PhysicalParameterTemplate(PhysicalParameter("D", 10, to_sim=1, is_global=True))
t_kd = PhysicalParameterTemplate(PhysicalParameter("kd", 0.0, to_sim=1 / (60 * 2), is_global=True))

t_fraction = PhysicalParameterTemplate(PhysicalParameter("fraction", 0, to_sim=1, is_global=True))
t_sigma = PhysicalParameterTemplate(PhysicalParameter("sigma", 10000, to_sim=N_A ** -1 * 1e9, is_global=True))
t_gamma = PhysicalParameterTemplate(PhysicalParameter("gamma", 0, to_sim=1, is_global=True))
t_k_factor = PhysicalParameterTemplate(PhysicalParameter("k_factor", 2, to_sim=1, is_global=True))

t_R_start = PhysicalParameterTemplate(PhysicalParameter("R_start", 0, to_sim=N_A ** -1 * 1e9))

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


def setup(cytokines, cell_types, geometry_dict, numeric, path, ext_cache=""):
    message("setup")

    def make_grid(gd: Dict, key: str):
        if key in gd.keys():
            return np.round(np.arange(gd["margin"], gd[key], gd["distance"]), 2)
        else:
            return [0]

    x = make_grid(geometry_dict, "x_grid")
    y = make_grid(geometry_dict, "y_grid")
    z = make_grid(geometry_dict, "z_grid")

    margin = geometry_dict["margin"]

    global_parameter = make_global_parameters(cytokines, geometry_dict, numeric)
    na = geometry_dict["norm_area"]
    geometry = ParameterCollection("geometry", [PhysicalParameter("norm_area", na, to_sim=1)])
    domain_parameter_set = deepcopy(global_parameter)
    domain_parameter_set.add_collection(geometry)
    domain_parameter_set.name = "domain"

    cell_type_list, fractions = make_cell_types(cell_types, cytokines, numeric)
    global_parameter.add_collection(fractions)

    domain_bc = make_domain_bc(cytokines, domain_parameter_set, margin, x)

    p1, p2 = get_cuboid_domain(x,y,z, margin)

    domain = Entity.DomainCube(
        p1,
        p2, domain_bc)

    sc = SC.SimContainer(global_parameter)

    """FieldProblems"""
    for c in cytokines:

        solver = MySolver.MyLinearSoler()

        fieldProblem = fp.FieldProblem()
        fieldProblem.field_name = c["name"]
        fieldProblem.field_quantity = c["field_quantity"]

        fieldProblem.set_solver(solver)

        if not ext_cache == "":
            fieldProblem.ext_cache = ext_cache

        fieldProblem.set_outer_domain(domain)
        sc.add_field(fieldProblem)


    """top level path"""
    sc.path = path

    """adds cell to simulation"""
    for i in makeCellListGrid(global_parameter, cytokines, x, y, z):
        sc.add_entity(i)

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
    return sc


def make_domain_bc(cytokines, domain_parameter_set, margin, x):
    domainBC = []
    for c in cytokines:
        left_boundary = bc.OuterIntegral(
            cellBC, "near(x[0],{d})".format(d=x[0] - margin), field_quantity=c["field_quantity"],
            p=deepcopy(domain_parameter_set), name="left_boundary"
        )
        box = bc.OuterIntegral(
            cellBC, "!near(x[0],{d})".format(d=x[0] - margin), field_quantity=c["field_quantity"],
            p=deepcopy(domain_parameter_set), name="box")

        domainBC.append(left_boundary)
        domainBC.append(box)

    return domainBC


def make_cell_types(cell_types, cytokines, numeric) -> (List[CellType], ParameterCollection):
    fractions = ParameterCollection("fractions", [])
    cell_types_list = []

    from PostProcess import get_concentration_conversion as get_cc
    t_threshold = PhysicalParameterTemplate(
        PhysicalParameter("ths", 0.1, to_sim=1 / get_cc(numeric["unit_length_exponent"])))

    for ct in cell_types:

        fractions.set_parameter(PhysicalParameter(ct["name"], ct["fraction"], is_global=True))
        cell_p_set = ParameterSet(ct["name"], [])
        for c in cytokines:
            if c["field_quantity"] in ct.keys():
                ct_dict = ct[c["field_quantity"]]
                p = []

                if "R" in ct_dict.keys():
                    p.append(t_R(ct_dict["R"]))
                if "q" in ct_dict.keys():
                    p.append(t_q(ct_dict["q"]))
                if "ths" in ct_dict.keys():
                    p.append(t_threshold(ct_dict["ths"]))

                collection = ParameterCollection(c["name"], p)
                collection.field_quantity = c["field_quantity"]
                cell_p_set.add_collection(collection)

        cell_type = CellType(cell_p_set, ct["name"], ct["internal_solver"])
        cell_types_list.append(cell_type)
    return cell_types_list, fractions


def make_global_parameters(cytokines: List, geometry: Dict, numeric: Dict) -> ParameterSet:


    global_parameter = ParameterSet("global", [])

    ule = numeric["unit_length_exponent"]

    rho = ParameterCollection("rho", [PhysicalParameter("rho", geometry["rho"], to_sim=10**(-6-ule))])
    global_parameter.add_collection(rho)

    numeric_c = ParameterCollection("numeric", [])
    global_parameter.add_collection(numeric_c)

    for c in cytokines:

        collection = ParameterCollection(c["name"], [])
        collection.field_quantity = c["field_quantity"]
        p = [t_R(0), t_q(0)]

        if "k_on" in c.keys():
            p.append(t_k_on(c["k_on"]))
        else:
            p.append(t_k_on(111.6))

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

    for k,v in numeric.items():
        numeric_c.set_misc_parameter(MiscParameter(k,v))

    return global_parameter

