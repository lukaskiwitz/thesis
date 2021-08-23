from thesis.main.MyFieldTemplate import MyCytokineTemplate
from thesis.main.MyDomainTemplate import MyBoxDomainTemplate, MyBoundingBoxTemplate
from thesis.main.MyGlobalModel import MyPDEModel
from thesis.main.MyEntityLocator import MyEntityLocator, MyCellGridLocator
from thesis.main.EntityType import EntityType, CellType
from thesis.main.InternalSolver import InternalSolver
from thesis.main.MyScenario import MyScenario
from thesis.main.MyParameterPool import MyParameterPool
from thesis.main.SimContainer import SimContainer
from thesis.main.MySolver import MyDiffusionSolver
from thesis.main.ParameterSet import ParameterSet, PhysicalParameterTemplate,PhysicalParameter, ParameterCollection, MiscParameter
from thesis.main.MyInteraction import MyFieldInteractionTemplate
import matplotlib.pyplot as plt
import numpy as np
from scipy.constants import N_A
from thesis.main.PostProcessUtil import get_concentration_conversion as get_cc
from typing import List, Dict
from copy import deepcopy
import os

ule = -6

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

def setup(cytokines, cell_types, boundary, geometry_dict, numeric, custom_pool = None):


    parameter_pool = MyParameterPool()

    for i, t in templates.items():
        parameter_pool.add_template(t)

    if isinstance(custom_pool,MyParameterPool):
        parameter_pool.join(custom_pool, override=False)

    numeric = ParameterCollection("numeric", [MiscParameter(k, v, is_global=True) for k, v in numeric.items()])
    cell_types, fractions = make_cell_types(cell_types, cytokines, parameter_pool)

    for ct in cell_types:
        ct.p.add_collection(ParameterCollection("rho", [parameter_pool.get_template("rho")(5)]))

    # x = make_grid(geometry_dict, "x_grid")
    # y = make_grid(geometry_dict, "y_grid")
    # z = make_grid(geometry_dict, "z_grid")

    # p1,p2 = get_cuboid_domain(x,y,z, geometry_dict["margin"])

    # domain_template = MyBoxDomainTemplate()
    domain_template = MyBoundingBoxTemplate()


    pde_model = MyPDEModel("pde_model")
    pde_model.domain_template = domain_template
    for c in cytokines:
        cytokine_template = MyCytokineTemplate()
        cytokine_template.name = c["name"]
        cytokine_template.field_quantity = c["field_quantity"]
        pde_model.add_field_template(cytokine_template)


    na = geometry_dict["norm_area"]
    geometry = ParameterCollection("geometry", [MiscParameter(k, v) for k,v in geometry_dict.items()])
    geometry.set_parameter(PhysicalParameter("norm_area", na, to_sim=1),override = True)

    domain_parameter_set = ParameterSet("domain",[])
    domain_parameter_set.add_collection(geometry)

    domain_template.bc_list = make_domain_bc(cytokines, boundary, numeric, domain_parameter_set, parameter_pool)
    locator = MyCellGridLocator()
    # internal_solver = InternalSolver()

    scenario = MyScenario(parameter_pool)

    scenario.global_models = [pde_model]
    # scenario.internal_solvers = [internal_solver]
    scenario.entity_types = cell_types
    scenario.entity_locators = [locator]

    scenario.global_parameters.update(geometry)
    scenario.global_parameters.update(numeric)
    scenario.global_parameters.update(fractions)

    return scenario

def make_cell_types(cell_types, cytokines, parameter_pool) -> (List[CellType], ParameterCollection):

    fractions = ParameterCollection("fractions", [], is_global=True)
    cell_types_list = []

    for ct in cell_types:

        fractions.set_parameter(PhysicalParameter(ct["name"], ct["fraction"], is_global=True))
        cell_p_set = ParameterSet(ct["name"], [])
        for c in cytokines:
            if c["field_quantity"] in ct.keys():
                ct_dict = ct[c["field_quantity"]]
                p = []

                for k, v in ct_dict.items():
                    if parameter_pool.get_template(k) is not None:
                        p.append(parameter_pool.get_template(k)(v))
                    else:
                        p.append(MiscParameter(k, v))
                del ct[c["field_quantity"]]

                collection = ParameterCollection(c["name"], p)
                collection.field_quantity = c["field_quantity"]
                cell_p_set.add_collection(collection)


        for k,v in ct.items():
            if isinstance(v,Dict):
                plist = []
                for kk,vv in v.items():
                    if parameter_pool.get_template(kk) is not None:
                        plist.append(parameter_pool.get_template(kk)(vv))
                    else:
                        plist.append(MiscParameter(kk, vv))
                col = ParameterCollection(k,plist)
                cell_p_set.add_collection(col)

            else:
                if parameter_pool.get_template(k) is not None:
                     cell_p_set.add_parameter_with_collection(parameter_pool.get_template(k)(v))
                else:
                    cell_p_set.add_parameter_with_collection(MiscParameter(k, v))


        from thesis.main.MyInteraction import FieldInteractionType

        cell_type = CellType(cell_p_set, ct["name"], ct["internal_solver"])
        cell_type.interactions.append(MyFieldInteractionTemplate("il2", FieldInteractionType.INTEGRAL))
        cell_types_list.append(cell_type)

    return cell_types_list, fractions

def make_domain_bc(cytokines, boundary, numeric, domain_parameter_set, parameter_pool):
    domainBC = []
    from thesis.main.BC import OuterIntegral
    from thesis.main.bcFunctions import cellBC
    import fenics as fcs

    for piece in boundary:
        for key, line in piece.items():
            if key in [c["field_quantity"] for c in cytokines]:
                outer_integral = OuterIntegral(cellBC, piece["expr"], p = deepcopy(domain_parameter_set), field_quantity=key, name = piece["name"])

                i = [c["field_quantity"] for c in cytokines].index(key)
                parameters = []
                for name, value in line.items():
                    if parameter_pool.get_template(name) is not None:
                        parameters.append(parameter_pool.get_template(name)(value))
                    else:
                        parameters.append(MiscParameter(name, value))

                s = ParameterSet("update", [ParameterCollection(cytokines[i]["name"], parameters, field_quantity=key)])
                outer_integral.p.update(s, override=True)
                domainBC.append(outer_integral)
    if len(domainBC) == 0:
        c = cytokines[0]

        def dummy_bc(u,p,fq,area = 1):
            return fcs.Constant(0)

        outer_integral = OuterIntegral(dummy_bc, "true", p = domain_parameter_set, field_quantity=c["field_quantity"], name = "box")
        domainBC.append(outer_integral)


    return domainBC


