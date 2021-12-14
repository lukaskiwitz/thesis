from scipy.constants import N_A

from thesis.main.EntityType import CellType
from thesis.main.MyDomainTemplate import MyBoundingBoxTemplate
from thesis.main.MyEntityLocator import MyCellListLocator
from thesis.main.MyFieldTemplate import MyCytokineTemplate, MyMeanCytokineTemplate
from thesis.main.MyGlobalModel import MyPDEModel, MyODEModel
from thesis.main.MyInteractionTemplate import FieldInteractionType
from thesis.main.MyInteractionTemplate import MyFieldInteractionTemplate
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
    "mc": PhysicalParameterTemplate(PhysicalParameter("mc", 0, to_sim=1))
}

numeric = {
    "linear_solver": "gmres",
    "preconditioner": "amg",
    "linear": False,
    "krylov_atol": 1e-35,
    "krylov_rtol": 1e-5,
    "newton_atol": 1e-35,
    "newton_rtol": 1e-5,
    "dofs_per_node": 15000,
    "max_mpi_nodes": 2,
    "cells_per_worker": 50,
    "max_pool_size": 1,
    "min_char_length": 0.1,
    "max_char_length": 5,
    "unit_length_exponent": -6
}

geometry = {
    "distance": 20,
    "margin": 20,
    "rho": 5
}
def setup() -> MyScenario:

    parameter_pool = MyParameterPool()

    for i, t in templates.items():
        parameter_pool.add_template(t)

    pde_model = MyPDEModel("pde_model")
    pde_model.domain_template = MyBoundingBoxTemplate()
    cytokine_template = MyCytokineTemplate()
    cytokine_template.name = "cytokine"
    cytokine_template.field_quantity = "cyt"
    pde_model.add_field_template(cytokine_template)

    ode_model = MyODEModel("ode_model")
    mean_cytokine_template = MyMeanCytokineTemplate()
    mean_cytokine_template.name = "cytokine"
    mean_cytokine_template.field_quantity = "cyt"
    ode_model.add_field_template(mean_cytokine_template)

    scenario = MyScenario(parameter_pool)

    scenario.global_models = [pde_model, ode_model]

    p = ParameterSet("global", [
        parameter_pool.get_as_collection({"D": None, "kd": None}),
        ParameterCollection("geometry", [MiscParameter(k, v) for k, v in geometry.items()]),
        ParameterCollection("numeric", [MiscParameter(k, v) for k, v in numeric.items()])
    ])
    scenario.global_parameters.update(p)


    cell_type1 = CellType(ParameterSet("cell_dummy",[
        parameter_pool.get_as_collection({"D":None,"R":1e4,"q":0,"k_on":None,"k_off":None,"k_endo":None,"Kc":None},name="cytokine",field_quantity="cyt"),
        parameter_pool.get_as_collection({"rho": geometry["rho"]}, name="rho")
    ]),"abs",None)
    cell_type1.interactions = [MyFieldInteractionTemplate("cyt", FieldInteractionType.INTEGRAL)]
    cell_type1.p.get_collection("cytokine").set_misc_parameter(MiscParameter("bc_type", "patrick_saturation"))

    cell_type2 = CellType(ParameterSet("cell_dummy", [
        parameter_pool.get_as_collection(
            {"D": None, "R": 1e2, "q": 1, "k_on": None, "k_off": None, "k_endo": None, "Kc": None}, name="cytokine",
            field_quantity="cyt"),
        parameter_pool.get_as_collection({"rho": geometry["rho"]}, name="rho")
    ]), "sec", None)
    cell_type2.interactions = [MyFieldInteractionTemplate("cyt", FieldInteractionType.INTEGRAL)]
    cell_type2.p.get_collection("cytokine").set_misc_parameter(MiscParameter("bc_type","patrick_saturation"))
    scenario.entity_types = [cell_type1, cell_type2]
    locator = MyCellListLocator([[-geometry["distance"] / 2, 0, 0], [geometry["distance"] / 2, 0, 0]],
                                [cell_type1, cell_type2])
    scenario.entity_locators = [locator]

    return scenario


