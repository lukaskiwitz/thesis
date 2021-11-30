from abc import ABC, abstractmethod

from thesis.main.FieldProblem import FieldProblem, MeanFieldProblem, GlobalProblem
from thesis.main.MyParameterPool import MyParameterPool
from thesis.main.MySolver import MyDiffusionSolver, MyMeanFieldSolver
from thesis.main.ParameterSet import ParameterCollection


class MyFieldTemplate(ABC):

    def __init__(self):
        self.name: str = None
        self.field_quantity: str = None
        self.collection: ParameterCollection = ParameterCollection("FieldTemplateDummy", [])

    @abstractmethod
    def get_problem(self) -> GlobalProblem: pass

    @abstractmethod
    def build_parameter_collection(self, parameter_pool: MyParameterPool) -> ParameterCollection: pass


class MyMeanCytokineTemplate(MyFieldTemplate):

    def get_problem(self) -> MeanFieldProblem:
        mean_field_problem = MeanFieldProblem()
        mean_field_problem.field_name = self.name
        mean_field_problem.field_quantity = self.field_quantity
        mean_field_problem.solver = MyMeanFieldSolver()

        return mean_field_problem

    def build_parameter_collection(self, parameter_pool: MyParameterPool) -> ParameterCollection:
        collection = ParameterCollection(self.name, [])
        collection.field_quantity = self.field_quantity

        collection.set_parameter(parameter_pool.get_template("k_on")())
        collection.set_parameter(parameter_pool.get_template("Kc")())
        collection.set_parameter(parameter_pool.get_template("k_off")())
        collection.set_parameter(parameter_pool.get_template("k_endo")())
        collection.set_parameter(parameter_pool.get_template("D")())
        collection.set_parameter(parameter_pool.get_template("kd")())

        collection.update(self.collection, overwrite=True)
        return collection


class MyCytokineTemplate(MyFieldTemplate):

    def __init__(self):
        super().__init__()
        self.ext_cache: str = ""

    def get_problem(self) -> FieldProblem:
        fieldProblem = FieldProblem()
        fieldProblem.field_name = self.name
        fieldProblem.field_quantity = self.field_quantity
        fieldProblem.ext_cache = self.ext_cache
        fieldProblem.solver = MyDiffusionSolver()
        fieldProblem.timeout = 24 * 60 ** 2

        return fieldProblem

    def build_parameter_collection(self, parameter_pool: MyParameterPool) -> ParameterCollection:
        collection = ParameterCollection(self.name, [])
        collection.field_quantity = self.field_quantity

        collection.set_parameter(parameter_pool.get_template("k_on")())
        collection.set_parameter(parameter_pool.get_template("k_off")())
        collection.set_parameter(parameter_pool.get_template("k_endo")())
        collection.set_parameter(parameter_pool.get_template("D")())
        collection.set_parameter(parameter_pool.get_template("kd")())

        collection.update(self.collection, overwrite=True)
        return collection
