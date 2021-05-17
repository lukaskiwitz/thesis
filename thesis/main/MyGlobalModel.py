from thesis.main.MyFieldTemplate import MyFieldTemplate, MyCytokineTemplate
from thesis.main.ParameterSet import ParameterSet

class MyGlobalModel:

    def __init__(self,name: str):

        self.name = name


class MyODEModel(MyGlobalModel):
    pass

class MyPDEModel(MyGlobalModel):

    def __init__(self, name: str):

        self.name = name
        self.field_templates = []
        self.domain_template = None

    def add_field_template(self,field_template: MyFieldTemplate):

        assert isinstance(field_template,MyFieldTemplate)
        self.field_templates.append(field_template)

    def get_problem_list(self, p):

        problem_list = []

        for templates in self.field_templates:
            field_problem = templates.get_problem()
            problem_list.append(field_problem)
            field_problem.domain_template = self.domain_template

        return problem_list

    def build_parameter_set(self, parameter_pool):

        p = ParameterSet(self.name,[])

        for template in self.field_templates:
            p.add_collection(template.build_parameter_collection(parameter_pool))

        return p


