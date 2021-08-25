from thesis.main.FieldProblem import FieldProblem
from thesis.main.ParameterSet import ParameterCollection
from thesis.main.MySolver import MyDiffusionSolver
class MyFieldTemplate:

    def __init__(self):
        self.name = None
        self.field_quantity = None
        self.p = None
        self.ext_cache = ""

    def get_problem(self):
        pass
    def build_parameter_collection(self):
        pass


class MyCytokineTemplate(MyFieldTemplate):

    def get_problem(self):

        fieldProblem = FieldProblem()
        fieldProblem.field_name = self.name
        fieldProblem.field_quantity = self.field_quantity
        fieldProblem.ext_cache = self.ext_cache
        fieldProblem.solver = MyDiffusionSolver()
        fieldProblem.timeout = 24 * 60 ** 2

        return fieldProblem

    def build_parameter_collection(self, parameter_pool):

        collection = ParameterCollection(self.name, [])
        collection.field_quantity = self.field_quantity


        collection.set_parameter(parameter_pool.get_template("k_on")())
        collection.set_parameter(parameter_pool.get_template("k_off")())
        collection.set_parameter(parameter_pool.get_template("k_endo")())
        collection.set_parameter(parameter_pool.get_template("D")())
        collection.set_parameter(parameter_pool.get_template("kd")())

        return collection


