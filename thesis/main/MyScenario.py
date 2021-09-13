from copy import deepcopy

from thesis.main.MyParameterPool import MyParameterPool
from thesis.main.ParameterSet import ParameterSet
from thesis.main.ScanContainer import ScanSample
from thesis.main.SimContainer import SimContainer
from typing import Union

class MyScenario:

    def __init__(self, parameter_pool):

        assert isinstance(parameter_pool, MyParameterPool)

        self.global_models = []
        self.internal_solvers = []
        self.entity_types = []
        self.entity_locators = []
        self.parameter_pool = parameter_pool
        self.global_parameters: ParameterSet = ParameterSet("scenario_dummy", [])

    def get_sim_container(self, p:Union[ParameterSet,None]):

        parameter_set = deepcopy(self.global_parameters)
        if p is not None:
            parameter_set.update(p, override=True)

        for locator in self.entity_locators:
            cell_list = locator.get_entity_list(self.entity_types[0], parameter_set)

        for model in self.global_models:
            parameter_set.update(model.build_parameter_set(self.parameter_pool))

        sc = SimContainer(parameter_set)

        for model in self.global_models:
            for p in model.get_problem_list(parameter_set):
                sc.add_problem(p)

        for c in cell_list:
            sc.add_entity(c)

        for entity_type in self.entity_types:
            sc.add_entity_type(entity_type)

        default = deepcopy(ScanSample(parameter_set.collections, self.entity_types, {}))
        sc.default_sample = default

        return sc

    def get_entity_type_by_name(self, name: str):
        for e in self.entity_types:
            if e.name == name:
                return e

        return None

    def serialize_to_xml(self):
        raise NotImplementedError

    def deserialize_from_xml(self):
        raise NotImplementedError
