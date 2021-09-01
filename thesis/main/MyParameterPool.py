from thesis.main.ParameterSet import ParameterTemplate


class MyParameterPool:

    def __init__(self):
        self.parameter_templates = []

    def join(self, pool, override=True):

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
