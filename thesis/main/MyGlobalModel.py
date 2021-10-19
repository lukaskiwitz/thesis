from abc import ABC, abstractmethod
from typing import List

from thesis.main.FieldProblem import GlobalProblem, MeanFieldProblem, FieldProblem
from thesis.main.MyDomainTemplate import MyDomainTemplate
from thesis.main.MyFieldTemplate import MyCytokineTemplate, MyMeanCytokineTemplate, MyFieldTemplate
from thesis.main.ParameterSet import ParameterSet


class MyGlobalModel(ABC):

    def __init__(self, name: str):
        self.name: str = name
        self.field_templates: List[MyFieldTemplate] = []

    def build_parameter_set(self, parameter_pool):
        p = ParameterSet(self.name, [])

        for template in self.field_templates:
            p.add_collection(template.build_parameter_collection(parameter_pool))

        return p

    @abstractmethod
    def add_field_template(self, field_template: MyFieldTemplate): pass

    @abstractmethod
    def get_problem_list(self, p: ParameterSet) -> List[GlobalProblem]: pass


class MyODEModel(MyGlobalModel):

    def __init__(self, name: str):
        self.name: str = name
        self.field_templates: List[MyMeanCytokineTemplate] = []

    def add_field_template(self, field_template: MyMeanCytokineTemplate):
        assert isinstance(field_template, MyMeanCytokineTemplate)
        self.field_templates.append(field_template)

    def get_problem_list(self, p: ParameterSet) -> List[MeanFieldProblem]:
        problem_list: List[MeanFieldProblem] = []

        for templates in self.field_templates:
            field_problem = templates.get_problem()
            problem_list.append(field_problem)

        return problem_list


class MyPDEModel(MyGlobalModel):

    def __init__(self, name: str):
        super().__init__(name)
        self.field_templates: List[MyCytokineTemplate] = []
        self.domain_template: MyDomainTemplate = None

    def add_field_template(self, field_template: MyCytokineTemplate):
        assert isinstance(field_template, MyFieldTemplate)
        self.field_templates.append(field_template)

    def get_problem_list(self, p: ParameterSet) -> List[FieldProblem]:
        problem_list: List[FieldProblem] = []

        for templates in self.field_templates:
            field_problem = templates.get_problem()
            problem_list.append(field_problem)
            field_problem.domain_template = self.domain_template

        return problem_list
