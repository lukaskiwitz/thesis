from typing import Dict

from scipy.constants import N_A

from thesis.main.EntityType import CellType
from thesis.main.MyDomainTemplate import MySphereDomainTemplate
from thesis.main.MyEntityLocator import MyCellListLocator
from thesis.main.MyFieldTemplate import MyCytokineTemplate
from thesis.main.MyGlobalModel import MyPDEModel
from thesis.main.MyInteractionTemplate import FieldInteractionType, MyFieldInteractionTemplate
from thesis.main.MyParameterPool import MyParameterPool
from thesis.main.MyScenario import MyScenario
from thesis.main.ParameterSet import ParameterSet, PhysicalParameterTemplate, PhysicalParameter, ParameterCollection, \
    MiscParameter
from thesis.main.PostProcessUtil import get_concentration_conversion as get_cc

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
    "mc": PhysicalParameterTemplate(PhysicalParameter("mc", 0, to_sim=1)),

    "gamma": PhysicalParameterTemplate(PhysicalParameter("gamma", 0.1, to_sim=1, is_global=True)),
    "R_start": PhysicalParameterTemplate(PhysicalParameter("R_start", 20000, to_sim=N_A ** -1 * 1e9)),
    "pSTAT5": PhysicalParameterTemplate(PhysicalParameter("pSTAT5", 0, to_sim=1)),
    "EC50": PhysicalParameterTemplate(PhysicalParameter("EC50", 0, to_sim=1 / get_cc(ule))),
    "global_q": PhysicalParameterTemplate(PhysicalParameter("global_q", True, to_sim=1, is_global=True)),
    "KD": PhysicalParameterTemplate(PhysicalParameter("KD", 7.437e-12, to_sim=1 / get_cc(ule))),  # post = nM
    "nu": PhysicalParameterTemplate(PhysicalParameter("nu", 1, to_sim=1 / 60 ** 2, is_global=True)),  # receptors/s
    "eta": PhysicalParameterTemplate(PhysicalParameter("eta", 1, to_sim=1 / 60 ** 2, is_global=True)),
}


def get_standard_pool():
    parameter_pool = MyParameterPool()
    for i, t in templates.items():
        parameter_pool.add_template(t)
    return parameter_pool


def setup(cytokine, boundary, geometry_dict, numeric, custom_pool=None) -> MyScenario:
    parameter_pool = get_standard_pool()
    if isinstance(custom_pool, MyParameterPool):
        parameter_pool.join(custom_pool, overwrite=False)

    numeric = ParameterCollection("numeric", [MiscParameter(k, v, is_global=True) for k, v in numeric.items()])
    domain_template = MySphereDomainTemplate([0, 0, 0], geometry_dict["radius"])
    pde_model = MyPDEModel("pde_model")
    pde_model.domain_template = domain_template

    cytokine_template = MyCytokineTemplate()
    cytokine_template.name = cytokine["name"]
    cytokine_template.field_quantity = cytokine["field_quantity"]
    pde_model.add_field_template(cytokine_template)

    # global_parameter = make_global_parameters(cytokine, geometry_dict, numeric, templates)

    na = geometry_dict["norm_area"]
    geometry = ParameterCollection("geometry", [MiscParameter(k, v) for k, v in geometry_dict.items()])
    geometry.set_parameter(PhysicalParameter("norm_area", na, to_sim=1), overwrite=True)
    domain_parameter_set = ParameterSet("domain", [])
    domain_parameter_set.add_collection(geometry)

    def q(u, p, fq, area=1):

        if p["scenario"] == "high_density":

            R = p["R"]
            k_on = p["k_on"]
            D = p["D"]
            N = p["N"]

            return -k_on * u * N * R / (area * D)

        elif p["scenario"] == "low_density":

            return 0
        else:
            raise Exception

    parameters = []

    for name, value in boundary[cytokine["field_quantity"]].items():

        if name in templates:
            parameters.append(parameter_pool.get_template(name)(value))
        else:
            parameters.append(MiscParameter(name, value))

    s = ParameterSet("update",
                     [ParameterCollection(cytokine["name"], parameters, field_quantity=cytokine["field_quantity"])])
    from thesis.main.BC import OuterIntegral
    outer_integral = OuterIntegral(q, "true", field_quantity=cytokine["field_quantity"])
    outer_integral.p.update(s, overwrite=True)
    outer_integral.p.add_collection(geometry)
    domain_template.bc_list = [
        outer_integral
    ]

    scenario = MyScenario(parameter_pool)
    scenario.global_models = [pde_model]

    """adds entity types"""
    cell_p_set = ParameterSet("cell", [])
    cell_p_set.add_collection(parameter_pool.get_as_collection({"R": 1e2, "q": 10}, name=cytokine["name"],
                                                               field_quantity=cytokine["field_quantity"]))
    cell_p_set.add_collection(parameter_pool.get_as_collection({"rho": 5}, name="rho"))

    cell_type = CellType(cell_p_set, "cell", "")
    cell_type.interactions = [MyFieldInteractionTemplate(cytokine["field_quantity"], FieldInteractionType.INTEGRAL)]

    scenario.entity_types.append(cell_type)

    locator = MyCellListLocator([[0, 0, 0]], [cell_type])
    scenario.entity_locators = [locator]

    scenario.global_parameters.update(geometry)
    scenario.global_parameters.update(numeric)

    return scenario


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