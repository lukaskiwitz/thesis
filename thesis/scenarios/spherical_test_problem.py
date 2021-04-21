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


def setup(cytokine, boundary, geometry_dict, numeric, path, ext_cache=""):

    message("---------------------------------------------------------------")
    message("Setup")
    message("---------------------------------------------------------------")

    templates = get_parameter_templates(numeric["unit_length_exponent"])

    global_parameter = make_global_parameters(cytokine, geometry_dict, numeric, templates)

    na = geometry_dict["norm_area"]
    geometry = ParameterCollection("geometry", [PhysicalParameter("norm_area", na, to_sim=1)])

    domain_parameter_set = deepcopy(global_parameter)
    domain_parameter_set.add_collection(geometry)
    domain_parameter_set.name = "domain"


    def q(u, p, fq, area = 1):

        if p["scenario"] == "high_density":

            R = p["R"]
            k_on = p["k_on"]
            D = p["D"]
            N = p["N"]


            return -k_on*u*N*R/(area*D)

        elif p["scenario"] == "low_density":
            
            return 0
        else:
            raise  Exception

    parameters = []

    for name, value in boundary[cytokine["field_quantity"]].items():

        if name in templates:
            parameters.append(templates[name](value))
        else:
            parameters.append(MiscParameter(name, value))


    s = ParameterSet("update",[ParameterCollection(cytokine["name"], parameters, field_quantity=cytokine["field_quantity"])])
    outer_integral = bc.OuterIntegral(q, "true", field_quantity=cytokine["field_quantity"])
    outer_integral.p.update(s, override=True)
    outer_integral.p.add_collection(geometry)
    domain_bc = [
        outer_integral
    ]

    domain = thesis.main.Entity.DomainSphere(
        [0,0,0],
        geometry_dict["radius"], domain_bc)

    sc = SC.SimContainer(global_parameter)

    """FieldProblems"""


    solver = thesis.main.MySolver.MyDiffusionSolver()
    solver.timeout = 60 ** 2

    fieldProblem = fp.FieldProblem()
    fieldProblem.field_name = cytokine["name"]
    fieldProblem.field_quantity = cytokine["field_quantity"]

    fieldProblem.set_solver(solver)

    if not ext_cache == "":
        fieldProblem.ext_cache = ext_cache

    fieldProblem.set_outer_domain(domain)
    sc.add_field(fieldProblem)

    """top level path"""
    sc.path = path

    """adds cell to simulation"""

    cell_bcs = [bc.Integral(cellBC, field_quantity=cytokine["field_quantity"])]

    cell = Entity.Cell([0,0,0], geometry_dict["rho"], cell_bcs)
    cell.change_type = "cell"
    cell.p.add_parameter_with_collection(MiscParameter("id", 1))
    sc.add_entity(cell)

    """adds entity types"""
    cell_p_set = ParameterSet("cell", [])
    t_q = templates["q"]
    t_R = templates["R"]


    p = [t_R(1e2),t_q(10)]


    collection = ParameterCollection(cytokine["name"], p)
    collection.field_quantity = cytokine["field_quantity"]
    cell_p_set.add_collection(collection)

    cell_type = CellType(cell_p_set, "cell", "")
    sc.add_entity_type(cell_type)


    message("initializing sim container")
    """sets external path for subdomain markers"""
    if not ext_cache == "":
        sc.initialize(load_subdomain=True, file_name="mesh")
    else:
        sc.initialize()
    message("initialization complete")

    from thesis.main.ScanContainer import ScanSample
    default = deepcopy(ScanSample(global_parameter.collections, [cell_type], {}))
    sc.default_sample = default

    return sc


def make_global_parameters(cytokine: Dict, geometry: Dict, numeric: Dict, templates) -> ParameterSet:

    t_k_on = templates["k_on"]
    t_D = templates["D"]
    t_kd = templates["kd"]
    t_rho = templates["rho"]
    t_q = templates["q"]
    t_R = templates["R"]

    global_parameter = ParameterSet("global", [])
    global_parameter.add_collection(ParameterCollection("rho", [t_rho(geometry["rho"])]))

    numeric_c = ParameterCollection("numeric", [], is_global=True)
    global_parameter.add_collection(numeric_c)


    p = [t_R(0), t_q(0)]
    p.append(t_k_on(cytokine["k_on"]))
    p.append(t_D(cytokine["D"]))
    p.append(t_kd(cytokine["kd"]))

    collection = ParameterCollection(cytokine["name"], p)
    collection.field_quantity = cytokine["field_quantity"]
    global_parameter.add_collection(collection)
    for k, v in numeric.items():
        numeric_c.set_misc_parameter(MiscParameter(k, v))

    return global_parameter


def get_parameter_templates(ule):
    templates = {
        "R": PhysicalParameterTemplate(PhysicalParameter("R", 0, to_sim=N_A ** -1 * 1e9)),
        "k_on": PhysicalParameterTemplate(PhysicalParameter("k_on", 111.6, to_sim=1e15 / 60 ** 2, is_global=True)),
        "q": PhysicalParameterTemplate(PhysicalParameter("q", 0, to_sim=N_A ** -1 * 1e9)),
        "D": PhysicalParameterTemplate(PhysicalParameter("D", 10, to_sim=1, is_global=True)),
        "kd": PhysicalParameterTemplate(PhysicalParameter("kd", 0.1, to_sim=1 / (60 ** 2), is_global=True)),
        "threshold": PhysicalParameterTemplate(PhysicalParameter("ths", 0.1, to_sim=1 / get_cc(ule))),
        "Kc": PhysicalParameterTemplate(PhysicalParameter("Kc", 0.01, to_sim=1 / get_cc(ule))),
        "bw": PhysicalParameterTemplate(PhysicalParameter("bw", 10, to_sim=10 ** (6 + ule))),
        "cluster_strength": PhysicalParameterTemplate(PhysicalParameter("strength", 10, to_sim=1)),
        "rho": PhysicalParameterTemplate(PhysicalParameter("rho", 0, to_sim=10 ** (-6 - ule))),
        "amax": PhysicalParameterTemplate(PhysicalParameter("amax", 0, to_sim=N_A ** -1 * 1e9)),
        "mc": PhysicalParameterTemplate(PhysicalParameter("mc", 0, to_sim=1))
    }

    return templates
