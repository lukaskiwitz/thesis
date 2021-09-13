from thesis.main.ParameterSet import ParameterTemplate,ParameterCollection
from typing import List, Mapping

class MyParameterPool:

    def __init__(self):
        self.parameter_templates: List[ParameterTemplate] = []

    def join(self, pool: 'MyParameterPool', override=True):

        own_names = [i.name for i in self.parameter_templates]

        for t in pool.parameter_templates:
            if (t.name in own_names and override):
                for own_t in self.parameter_templates:
                    if t.name == own_t.name:
                        del own_t
                        break
                self.parameter_templates.append(t)

            elif t.name not in own_names:
                self.parameter_templates.append(t)

    def get_template(self, name: str):
        for t in self.parameter_templates:
            if t.name == name:
                return t
        return None

    def add_template(self, template: ParameterTemplate):
        assert isinstance(template, ParameterTemplate)
        self.parameter_templates.append(template)

    def get_as_collection(self, parameter_name_dict: Mapping[str,float], name = "dummy",field_quantity=""):

        c = ParameterCollection(name,[],field_quantity=field_quantity)
        for k,v in parameter_name_dict.items():
            t = self.get_template(k)
            if t is not None:
                    if v is None:
                        c.set_parameter(t(in_sim = False),override=True)
                    else:
                        c.set_parameter(t(value = v, in_sim = False))
        return c

